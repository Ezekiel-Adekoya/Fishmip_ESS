suppressMessages({library(mizer); library(therMizer); library(dplyr); library(tidyr)})
source("Ess_Fishmip_data/Useful_functions.R")
quiet <- function(x) suppressWarnings(suppressMessages(x)); area <- 146360555008; sc <- Sys.getenv("SC"); forage <- c("Atlantic Herring","Sand Lance","Capelin")
C <- readRDS(file.path(sc, "seal_calibrated.rds")); vm <- C$vm; ex <- C$ex; u217 <- C$u217; tmat <- C$tmat; target <- C$target_seals
q <- C$params; sp_all <- as.character(q@species_params$species); fish <- setdiff(sp_all, "Grey Seal"); is_seal <- sp_all == "Grey Seal"
dg <- readRDS("/Users/juliab6/Desktop/Ez June 2026/ESS_Restructured/Ess_data/seal_abundance_digitised.RDS"); counts <- setNames(dg$scotian_mean, dg$year)
ef <- read.csv("Ess_Fishmip_data/Effort_spin_1861_2010.csv", check.names=FALSE); emat <- as.matrix(ef[, setdiff(names(ef),"year")]); rownames(emat) <- as.character(ef$year)
seal_k <- function(p) sum(p@initial_n["Grey Seal",]*p@dw)/1000
fish_eaten_by_species <- function(p, n=p@initial_n, n_pp=p@initial_n_pp, t=1) { idx <- (length(p@w_full)-length(p@w)+1):length(p@w_full); pr <- getPredRate(p, n=n, n_pp=n_pp, t=t)["Grey Seal", idx]; sapply(fish, function(s) sum(p@interaction["Grey Seal", s]*pr*n[s,]*p@w*p@dw))/1e9 }
solve <- function(p) quiet(steady(p, t_max=300, progress_bar=FALSE, preserve="erepro"))
# keep the seal thermal limits where therMizer re-reads them
g <- q@given_species_params; g$temp_min <- q@species_params$temp_min; g$temp_max <- q@species_params$temp_max
# ---- 1. diet: no plankton resource; 2. prey size: preferred prey about 200 g ----
q@species_params$interaction_resource[is_seal] <- 0; g$interaction_resource[g$species=="Grey Seal"] <- 0
seal_w <- 225000; q@species_params$beta[is_seal] <- seal_w/200; q@species_params$sigma[is_seal] <- 1.5; g$beta[g$species=="Grey Seal"] <- seal_w/200; g$sigma[g$species=="Grey Seal"] <- 1.5; q@given_species_params <- g
old_kernel <- getPredKernel(q)["Grey Seal", which.min(abs(q@w - seal_w)), ]; q <- quiet(setPredKernel(q)); new_kernel <- getPredKernel(q)["Grey Seal", which.min(abs(q@w - seal_w)), ]
half <- function(k) range(q@w_full[k > 0.5*max(k)]); cat(sprintf("prey a 225 kg seal prefers (kernel above half its peak): before %.0f-%.0f g, now %.0f-%.0f g\n", half(old_kernel)[1], half(old_kernel)[2], half(new_kernel)[1], half(new_kernel)[2]))
cat("seal search volume still the compensated one after resetting the kernel:", isTRUE(all.equal(q@search_vol["Grey Seal",], C$params@search_vol["Grey Seal",])), "\n")
# ---- 3. scale feeding to the old gear's 1980 removal, re-matching the fish each pass ----
target_kt <- 22.5
for (i in 1:10) { now <- sum(quiet(fish_eaten_by_species(q))); if (!is.finite(now) || now <= 0) stop("fish eaten is ", now, " at pass ", i)
  step <- min(2.5, max(0.4, target_kt/now)); sv <- q@search_vol; sv["Grey Seal",] <- sv["Grey Seal",]*step; q <- setSearchVolume(q, search_vol=sv); q <- solve(q)
  q@species_params$R_max[is_seal] <- q@species_params$R_max[is_seal]*target/seal_k(q); q@initial_n["Grey Seal",] <- q@initial_n["Grey Seal",]*target/seal_k(q); q <- quiet(matchBiomasses(q, species=fish))
  cat(sprintf("  pass %2d: fish eaten %5.1f kt/yr | fish %.3f g/m2 | seal feeding level %.2f\n", i, sum(quiet(fish_eaten_by_species(q))), sum(getBiomass(q)[fish])/area, mean(getFeedingLevel(q)["Grey Seal", q@w>=16000]))) }
q <- solve(q); sp <- q@species_params; r <- getBiomass(q)[fish]/sp[fish,"biomass_observed"]; yg <- colSums(getYieldGear(q))[fish]; yo <- sp[fish,"yield_observed"]; ok <- !is.na(yo) & yo>0
eaten <- quiet(fish_eaten_by_species(q)); cat(sprintf("\ncalibrated: %d/19 fish within 2x (median %.3f) | yield %d/%d within 2x | fish with erepro > 1: %s\n", sum(r>0.5&r<2), median(r), sum((yg/yo)[ok]>0.5&(yg/yo)[ok]<2), sum(ok), paste(sprintf("%s %.2f", fish[sp[fish,"erepro"]>1], sp[fish,"erepro"][sp[fish,"erepro"]>1]), collapse=", ")))
cat(sprintf("seals take %.1f kt of fish a year at %.1f thousand = %.2f t per seal; top prey: %s\n", sum(eaten), seal_k(q), sum(eaten)*1e3/(seal_k(q)*1e3), paste(sprintf("%s %.1f", names(sort(eaten, decreasing=TRUE))[1:6], sort(eaten, decreasing=TRUE)[1:6]), collapse=", ")))
saveRDS(q, file.path(sc, "seal_fixed_params.rds"))
# ---- projection: spin-up with few seals, then seal numbers set from the counts every year ----
qd <- quiet(upgradeTherParams(q, ocean_temp=tmat, n_pp_array=u217, vertical_migration_array=vm, exposure_array=ex, aerobic_effect=TRUE, metabolism_effect=TRUE)); stopifnot(!anyNA(scaled_temp_effect(qd, 1861)), isTRUE(all.equal(qd@search_vol, q@search_vol)))
E <- emat[, names(q@initial_effort)]
a <- qd; a@species_params$R_max[is_seal] <- q@species_params$R_max[is_seal]*counts["1960"]/target; a@initial_n["Grey Seal",] <- q@initial_n["Grey Seal",]*counts["1960"]/seal_k(q)
s1 <- quiet(project(a, effort=E[as.character(1861:1961),], progress_bar=FALSE)); n <- s1@n[dim(s1@n)[1],,]; npp <- s1@n_pp[dim(s1@n_pp)[1],]
years <- 1961:2010; B <- matrix(NA, length(years), length(sp_all), dimnames=list(years, sp_all)); Y <- B[, fish]; eat <- B[, fish]; sealN <- setNames(rep(NA, length(years)), years); juvM <- matrix(NA, length(years), 3, dimnames=list(years, c("Atlantic Cod","Haddock","American Plaice")))
idx <- (length(q@w_full)-length(q@w)+1):length(q@w_full)
for (y in years) { n["Grey Seal",] <- n["Grey Seal",] * counts[as.character(y)] / (sum(n["Grey Seal",]*q@dw)/1000)
  b <- qd; b@initial_n[] <- n; b@initial_n_pp[] <- npp; Ey <- rbind(E[as.character(y),], E[as.character(y),]); rownames(Ey) <- c(y, y+1)
  s <- quiet(project(b, effort=Ey, progress_bar=FALSE)); B[as.character(y),] <- getBiomass(s)[1,]; yg <- getYieldGear(s); Y[as.character(y),] <- colSums(yg[1,,])[fish]; sealN[as.character(y)] <- sum(n["Grey Seal",]*q@dw)/1000
  eat[as.character(y),] <- quiet(fish_eaten_by_species(b, n=n, n_pp=npp, t=y)); pr <- getPredRate(b, n=n, n_pp=npp, t=y)["Grey Seal", idx]
  for (s2 in colnames(juvM)) { win <- q@w >= 90 & q@w <= 480; wt <- n[s2,win]*q@w[win]*q@dw[win]; juvM[as.character(y), s2] <- sum(q@interaction["Grey Seal", s2]*pr[win]*wt)/sum(wt) }
  n <- s@n[2,,]; npp <- s@n_pp[2,] }
saveRDS(list(B=B, Y=Y, eat=eat, sealN=sealN, juvM=juvM), file.path(sc, "seal_fixed_run.rds"))
cat("\nyear | seals (thousand) | fish eaten kt | t per seal | seal mortality on 90-480 g cod / haddock / plaice, per year | total fish 16 spp g/m2\n")
for (y in c(1961,1970,1980,1990,1995,2000,2005,2010)) { yc <- as.character(y); cat(sprintf("%d | %5.0f | %6.1f | %.2f | %.3f / %.3f / %.3f | %.2f\n", y, sealN[yc], sum(eat[yc,]), sum(eat[yc,])/sealN[yc], juvM[yc,1], juvM[yc,2], juvM[yc,3], sum(B[yc, setdiff(fish, forage)])/area)) }
obsB <- readRDS("Ess_Fishmip_data/biomass_ts.RDS") |> transmute(Year=as.integer(Year), Species=as.character(Species), obs=Biomass)
obsC <- readRDS("Ess_Fishmip_data/new_NAFO_4VW_Catches.RDS") |> pivot_longer(-Year, names_to="Species", values_to="obs") |> filter(!is.na(obs), obs>0)
long <- function(m) as.data.frame(as.table(m)) |> setNames(c("Year","Species","mod")) |> mutate(Year=as.integer(as.character(Year)), Species=as.character(Species))
ps <- function(mod, obs, years) inner_join(mod, obs, by=c("Year","Species")) |> filter(Year %in% years, mod>0, obs>0, !Species %in% forage) |> group_by(Species) |> filter(n()>=4) |> summarise(bias=median(log10(mod/obs)), r=suppressWarnings(cor(log(mod),log(obs))), rmse=sqrt(mean((log10(mod)-log10(obs))^2)), .groups="drop")
f <- function(d) sprintf("x%.2f %2d/%2d r%+.2f rmse %.2f", 10^median(d$bias), sum(abs(d$bias)<log10(2)), nrow(d), median(d$r,na.rm=TRUE), median(d$rmse))
cat(sprintf("\nrevised dynamic seal | biomass 1970-2004 %s | biomass 2005-2010 %s | catch 1970-2004 %s | catch 2005-2010 %s\n", f(ps(long(B[,fish]), obsB, 1970:2004)), f(ps(long(B[,fish]), obsB, 2005:2010)), f(ps(long(Y), obsC, 1970:2004)), f(ps(long(Y), obsC, 2005:2010))))
rel <- function(s) mean(B[as.character(1994:2004), s])/mean(B[as.character(1978:1981), s]); ob <- function(s) with(obsB[obsB$Species==s,], mean(obs[Year %in% 1994:2004])/mean(obs[Year %in% 1978:1981]))
for (s in c("Atlantic Cod","Haddock","American Plaice","Witch Flounder","Winter Skate","Atlantic Herring","Sand Lance")) cat(sprintf("  %-18s 1994-2004 rel. to 1978-81: model %.2f | observed %.2f\n", s, rel(s), ob(s)))
