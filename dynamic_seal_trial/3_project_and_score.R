suppressMessages({library(mizer); library(therMizer); library(dplyr); library(tidyr); library(parallel)})
source("Ess_Fishmip_data/Useful_functions.R")
quiet <- function(x) suppressWarnings(suppressMessages(x)); area <- 146360555008; sc <- Sys.getenv("SC"); forage <- c("Atlantic Herring","Sand Lance","Capelin")
C <- readRDS(file.path(sc, "seal_calibrated.rds")); vm <- C$vm; ex <- C$ex; u217 <- C$u217; tmat <- C$tmat; target <- C$target_seals
fish <- setdiff(as.character(C$params@species_params$species), "Grey Seal")
dg <- readRDS("/Users/juliab6/Desktop/Ez June 2026/ESS_Restructured/Ess_data/seal_abundance_digitised.RDS")
ef <- read.csv("Ess_Fishmip_data/Effort_spin_1861_2010.csv", check.names=FALSE); emat <- as.matrix(ef[, setdiff(names(ef),"year")]); rownames(emat) <- as.character(ef$year)
seal_k <- function(q) sum(q@initial_n["Grey Seal",]*q@dw)/1000
fish_eaten_kt <- function(q) { idx <- (length(q@w_full)-length(q@w)+1):length(q@w_full); pr <- getPredRate(q)["Grey Seal", idx]; sum(sapply(fish, function(s) sum(q@interaction["Grey Seal", s]*pr*q@initial_n[s,]*q@w*q@dw)))/1e9 }
ref <- as.character(1978:1981); ny <- 2000; lab <- as.character(seq_len(ny))
const_inputs <- function(q) { r <- quiet(upgradeTherParams(q, ocean_temp=matrix(rep(colMeans(tmat[ref,]), each=ny), nrow=ny, dimnames=list(time=lab, realm=colnames(tmat))),
   n_pp_array=matrix(rep(colMeans(u217[ref,]), each=ny), nrow=ny, dimnames=list(time=lab, w=names(q@initial_n_pp))), vertical_migration_array=vm, exposure_array=ex, aerobic_effect=TRUE, metabolism_effect=TRUE))
   initialNResource(r) <- q@initial_n_pp; resource_dynamics(r) <- "resource_constant"; r }
solve <- function(q) quiet(steady(q, t_max=300, progress_bar=FALSE, preserve="erepro"))
# ---- variant B: seals eat fish and outside food only, scaled to the old gear's removal ----
make_fish_eating <- function(q, target_kt = 22.5) {
  q@species_params$interaction_resource[q@species_params$species == "Grey Seal"] <- 0
  for (i in 1:9) { now <- quiet(fish_eaten_kt(q)); if (!is.finite(now) || now <= 0) stop("fish eaten is ", now, " at pass ", i); step <- min(2.5, max(0.4, target_kt/now)); sv <- q@search_vol; sv["Grey Seal",] <- sv["Grey Seal",] * step; q <- setSearchVolume(q, search_vol=sv); if (anyNA(q@initial_n)) stop("NA in abundances before steady, pass ", i)
    q <- solve(q); rm <- q@species_params$R_max; rm[q@species_params$species=="Grey Seal"] <- rm[q@species_params$species=="Grey Seal"]*target/seal_k(q); q@species_params$R_max <- rm
    q@initial_n["Grey Seal",] <- q@initial_n["Grey Seal",]*target/seal_k(q); q <- quiet(matchBiomasses(q, species=fish))
    cat(sprintf("  fish-eating seal, pass %d: fish eaten %.1f kt/yr | seals %.1f thousand | fish %.3f g/m2 | seal feeding level %.2f\n", i, quiet(fish_eaten_kt(q)), seal_k(q), sum(getBiomass(q)[fish])/area, mean(getFeedingLevel(q)["Grey Seal", q@w>=16000]))) }
  q <- solve(q); q }
# ---- two-phase projection: few seals through the spin-up, then free to grow from 1961 ----
run <- function(q) {
  gsp <- q@given_species_params; gsp$temp_min <- q@species_params$temp_min; gsp$temp_max <- q@species_params$temp_max; q@given_species_params <- gsp; stopifnot(!anyNA(q@species_params$temp_min))
  qd <- quiet(upgradeTherParams(q, ocean_temp=tmat, n_pp_array=u217, vertical_migration_array=vm, exposure_array=ex, aerobic_effect=TRUE, metabolism_effect=TRUE)); stopifnot(!anyNA(scaled_temp_effect(qd, 1861)), isTRUE(all.equal(qd@search_vol, q@search_vol)))
  E <- emat[, names(q@initial_effort)]; s60 <- dg$scotian_mean[dg$year==1960]; is_seal <- q@species_params$species=="Grey Seal"
  a <- qd; a@species_params$R_max[is_seal] <- q@species_params$R_max[is_seal]*s60/target; a@initial_n["Grey Seal",] <- q@initial_n["Grey Seal",]*s60/seal_k(q)
  s1 <- quiet(project(a, effort=E[as.character(1861:1961),], progress_bar=FALSE))
  b <- qd; b@species_params$R_max[is_seal] <- 8100; n1 <- s1@n[dim(s1@n)[1],,]; b@initial_n[] <- n1; b@initial_n_pp[] <- s1@n_pp[dim(s1@n_pp)[1],]
  s2 <- quiet(project(b, effort=E[as.character(1961:2010),], progress_bar=FALSE))
  YG <- getYieldGear(s2); list(B=getBiomass(s2), Y=apply(YG, c(1,3), sum), seals=getN(s2)[, "Grey Seal"]/1000, B1=getBiomass(s1), sim=s2) }
qA <- C$params; cat("variant A (as Paper 2): fish eaten", round(quiet(fish_eaten_kt(qA)),1), "kt/yr\n")
qB <- readRDS(file.path(sc, "seal_qB.rds")); spB <- qB@species_params; cat("variant B fish with erepro > 1:", paste(sprintf("%s %.2f", fish[spB[fish,"erepro"]>1], spB[fish,"erepro"][spB[fish,"erepro"]>1]), collapse=", "), "\n")
rB <- getBiomass(qB)[fish]/spB[fish,"biomass_observed"]; cat(sprintf("variant B calibration: %d/19 fish within 2x, median ratio %.3f, max fish erepro %.2f\n", sum(rB>0.5 & rB<2), median(rB), max(spB[fish,"erepro"])))
res <- mclapply(list(A=qA, B=qB), run, mc.cores=2)
saveRDS(list(res=lapply(res, function(r) r[c("B","Y","seals","B1")]), qA=qA, qB=qB), file.path(sc, "seal_runs.rds")); saveRDS(res$B$sim, file.path(sc, "seal_simB.rds")); saveRDS(res$A$sim, file.path(sc, "seal_simA.rds"))
sg <- readRDS("Ess_Fishmip_data/sims_GFDL_recalibrated.RDS"); Bg <- getBiomass(sg); YGg <- getYieldGear(sg); Yg <- apply(YGg[, setdiff(dimnames(YGg)[[2]],"Seals"), , drop=FALSE], c(1,3), sum)
obsB <- readRDS("Ess_Fishmip_data/biomass_ts.RDS") |> transmute(Year=as.integer(Year), Species=as.character(Species), obs=Biomass)
obsC <- readRDS("Ess_Fishmip_data/new_NAFO_4VW_Catches.RDS") |> pivot_longer(-Year, names_to="Species", values_to="obs") |> filter(!is.na(obs), obs>0)
long <- function(m) as.data.frame(as.table(m)) |> setNames(c("Year","Species","mod")) |> mutate(Year=as.integer(as.character(Year)), Species=as.character(Species))
ps <- function(mod, obs, years) inner_join(mod, obs, by=c("Year","Species")) |> filter(Year %in% years, mod>0, obs>0, !Species %in% forage) |> group_by(Species) |> filter(n()>=4) |> summarise(bias=median(log10(mod/obs)), r=suppressWarnings(cor(log(mod),log(obs))), rmse=sqrt(mean((log10(mod)-log10(obs))^2)), .groups="drop")
f <- function(d) sprintf("x%.2f %2d/%2d r%+.2f rmse %.2f", 10^median(d$bias), sum(abs(d$bias)<log10(2)), nrow(d), median(d$r,na.rm=TRUE), median(d$rmse))
M <- list("seal gear (current model)"=list(Bg,Yg), "dynamic seal A: as Paper 2"=list(res$A$B[,fish],res$A$Y[,fish]), "dynamic seal B: eats fish only"=list(res$B$B[,fish],res$B$Y[,fish]))
cat(sprintf("\n%-32s | %-29s | %-29s | %-29s | %-29s\n", "16 spp without multipliers", "biomass 1970-2004", "biomass 2005-2010 held out", "catch 1970-2004", "catch 2005-2010 held out"))
for (m in names(M)) cat(sprintf("%-32s | %s | %s | %s | %s\n", m, f(ps(long(M[[m]][[1]]), obsB, 1970:2004)), f(ps(long(M[[m]][[1]]), obsB, 2005:2010)), f(ps(long(M[[m]][[2]]), obsC, 1970:2004)), f(ps(long(M[[m]][[2]]), obsC, 2005:2010))))
rel <- function(B, s) mean(B[as.character(1994:2004), s])/mean(B[as.character(1978:1981), s]); ob <- function(s) with(obsB[obsB$Species==s,], mean(obs[Year %in% 1994:2004])/mean(obs[Year %in% 1978:1981]))
sps <- c("Atlantic Cod","Haddock","American Plaice","Witch Flounder","Winter Skate","Atlantic Herring","Sand Lance")
cat(sprintf("\n%-32s | %s | total 16spp 1961 / 1980 / 2010\n", "1994-2004 rel. to 1978-81", paste(sprintf("%9s", c("cod","haddock","plaice","witch","w.skate","herring","sandlance")), collapse=" ")))
for (m in names(M)) { B <- M[[m]][[1]]; t16 <- function(y) sum(B[as.character(y), setdiff(colnames(B), forage)])/area; cat(sprintf("%-32s | %s | %.2f / %.2f / %.2f\n", m, paste(sprintf("%9.2f", sapply(sps, function(s) rel(B,s))), collapse=" "), t16(1961), t16(1980), t16(2010))) }
cat(sprintf("%-32s | %s | survey 16spp 2005-10 mean 3.96\n", "observed", paste(sprintf("%9.2f", sapply(sps, ob)), collapse=" ")))
cat("\nseal numbers, thousand:  year  data | A | B\n"); for (y in c(1961,1970,1980,1990,2000,2010)) cat(sprintf("  %d  %5.0f | %6.0f | %6.0f\n", y, dg$scotian_mean[dg$year==y], res$A$seals[as.character(y)], res$B$seals[as.character(y)]))
