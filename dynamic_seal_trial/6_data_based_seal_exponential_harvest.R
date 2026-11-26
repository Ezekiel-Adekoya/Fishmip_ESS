suppressMessages({library(mizer); library(therMizer); library(dplyr); library(tidyr)})
source("Ess_Fishmip_data/Useful_functions.R")
quiet <- function(x) suppressWarnings(suppressMessages(x)); area <- 146360555008; sc <- Sys.getenv("SC"); forage <- c("Atlantic Herring","Sand Lance","Capelin")
C <- readRDS(file.path(sc, "seal_calibrated.rds")); vm <- C$vm; ex <- C$ex; u217 <- C$u217; tmat <- C$tmat; target <- C$target_seals
D <- readRDS(file.path(sc, "seal_demography.rds")); dem <- D$with_harvest
q <- C$params; sp_all <- as.character(q@species_params$species); fish <- setdiff(sp_all, "Grey Seal"); is_seal <- sp_all == "Grey Seal"
dg <- readRDS("/Users/juliab6/Desktop/Ez June 2026/ESS_Restructured/Ess_data/seal_abundance_digitised.RDS"); counts <- setNames(dg$scotian_mean, dg$year)
ef <- read.csv("Ess_Fishmip_data/Effort_spin_1861_2010.csv", check.names=FALSE); emat <- as.matrix(ef[, setdiff(names(ef),"year")]); rownames(emat) <- as.character(ef$year)
idx <- (length(q@w_full)-length(q@w)+1):length(q@w_full); big <- q@w >= 16000
seal_k <- function(n) sum(n["Grey Seal",]*q@dw)/1000
eaten <- function(p, n=p@initial_n, n_pp=p@initial_n_pp, t=1) { pr <- getPredRate(p, n=n, n_pp=n_pp, t=t)["Grey Seal", idx]; sapply(fish, function(s) sum(p@interaction["Grey Seal", s]*pr*n[s,]*p@w*p@dw))/1e9 }
solve <- function(p) quiet(steady(p, t_max=300, progress_bar=FALSE, preserve="erepro"))
g <- q@given_species_params; g$temp_min <- q@species_params$temp_min; g$temp_max <- q@species_params$temp_max
q@species_params$interaction_resource[is_seal] <- 0; g$interaction_resource[g$species=="Grey Seal"] <- 0
q@species_params$beta[is_seal] <- 225000/20; q@species_params$sigma[is_seal] <- 1.5; g$beta[g$species=="Grey Seal"] <- 225000/20; g$sigma[g$species=="Grey Seal"] <- 1.5; q@given_species_params <- g; q <- quiet(setPredKernel(q))
# ---- consumption target from the folder: 12 x 170 kg x 85% fish, of which 19.1% is not sand lance ----
per_seal_other_fish <- 12*170*0.85/1000 * (1 - 0.809); target_kt <- per_seal_other_fish * target
cat(sprintf("target: %.2f t per seal of fish other than sand lance = %.1f kt at %.0f thousand seals\n", per_seal_other_fish, target_kt, target))
other <- function(e) sum(e[names(e) != "Sand Lance"])
for (i in 1:30) { now <- other(quiet(eaten(q))); if (!is.finite(now) || now <= 0) stop("fish eaten is ", now, " at pass ", i); if (abs(now/target_kt - 1) < 0.03 && i > 1) break
  step <- min(if (now < 0.2*target_kt) 2.5 else 1.25, max(0.7, (target_kt/now)^0.6)); sv <- q@search_vol; sv["Grey Seal",] <- sv["Grey Seal",]*step; q <- setSearchVolume(q, search_vol=sv); q <- solve(q)
  q@species_params$R_max[is_seal] <- q@species_params$R_max[is_seal]*target/seal_k(q@initial_n); q@initial_n["Grey Seal",] <- q@initial_n["Grey Seal",]*target/seal_k(q@initial_n); q <- quiet(matchBiomasses(q, species=fish))
  cat(sprintf("  pass %2d: other fish eaten %5.2f kt/yr | fish %.3f g/m2 | seal feeding level %.2f\n", i, other(quiet(eaten(q))), sum(getBiomass(q)[fish])/area, mean(getFeedingLevel(q)["Grey Seal", big]))) }
# outside food set so the seal feeds at the level its demography was fitted at
for (i in 1:4) { enc <- getEncounter(q)["Grey Seal", ]; ee <- ext_encounter(q); want <- D$feeding_level/(1-D$feeding_level) * q@species_params$h[is_seal] * q@w^q@species_params$n[is_seal]; ee["Grey Seal", big] <- pmax(0, ee["Grey Seal", big] + (want - enc)[big]); ext_encounter(q) <- ee }
q <- solve(q); q@initial_n["Grey Seal",] <- q@initial_n["Grey Seal",]*target/seal_k(q@initial_n); q <- quiet(matchBiomasses(q, species=fish)); q <- solve(q)
sp <- q@species_params; r <- getBiomass(q)[fish]/sp[fish,"biomass_observed"]; yg <- colSums(getYieldGear(q))[fish]; yo <- sp[fish,"yield_observed"]; ok <- !is.na(yo) & yo>0; e0 <- quiet(eaten(q))
cat(sprintf("\ncalibrated: %d/19 fish within 2x (median %.3f) | yield %d/%d within 2x | fish with erepro > 1: %s | seal feeding level %.2f\n", sum(r>0.5&r<2), median(r), sum((yg/yo)[ok]>0.5&(yg/yo)[ok]<2), sum(ok), paste(fish[sp[fish,"erepro"]>1], collapse=", "), mean(getFeedingLevel(q)["Grey Seal", big])))
cat(sprintf("seals take %.2f kt of fish a year at %.1f thousand (%.2f t per seal); top prey: %s\n", sum(e0), seal_k(q@initial_n), sum(e0)/seal_k(q@initial_n), paste(sprintf("%s %.2f", names(sort(e0, decreasing=TRUE))[1:6], sort(e0, decreasing=TRUE)[1:6]), collapse=", ")))
saveRDS(q, file.path(sc, "seal_fixed3_params.rds"))
# ---- projection: spin-up with few seals, then exponential growth with a harvest in 1967-1990 ----
qd <- quiet(upgradeTherParams(q, ocean_temp=tmat, n_pp_array=u217, vertical_migration_array=vm, exposure_array=ex, aerobic_effect=TRUE, metabolism_effect=TRUE)); stopifnot(!anyNA(scaled_temp_effect(qd, 1861)), isTRUE(all.equal(qd@search_vol, q@search_vol)))
E <- emat[, names(q@initial_effort)]
a <- qd; a@species_params$R_max[is_seal] <- q@species_params$R_max[is_seal]*counts["1960"]/target; a@initial_n["Grey Seal",] <- q@initial_n["Grey Seal",]*counts["1960"]/seal_k(q@initial_n)
s1 <- quiet(project(a, effort=E[as.character(1861:1961),], progress_bar=FALSE)); n <- s1@n[dim(s1@n)[1],,]; npp <- s1@n_pp[dim(s1@n_pp)[1],]
n["Grey Seal",] <- n["Grey Seal",]*counts["1961"]/seal_k(n)
base <- qd; base@species_params$erepro[is_seal] <- dem["erepro"]; base@species_params$R_max[is_seal] <- dem["R_max"]
phases <- list(list(1961, 1967, 0), list(1967, 1991, dem["harvest_rate"]), list(1991, 2010, 0)); sims <- list()
for (ph in phases) { b <- base; b@mu_b["Grey Seal", ] <- base@mu_b["Grey Seal", ] + ph[[3]]; b@initial_n[] <- n; b@initial_n_pp[] <- npp
  s <- quiet(project(b, effort=E[as.character(ph[[1]]:ph[[2]]),], progress_bar=FALSE)); sims[[length(sims)+1]] <- s; n <- s@n[dim(s@n)[1],,]; npp <- s@n_pp[dim(s@n_pp)[1],] }
take <- function(fun) { parts <- lapply(seq_along(sims), function(i) { m <- fun(sims[[i]]); if (i < length(sims)) m[-nrow(m), , drop=FALSE] else m }); do.call(rbind, parts) }
B <- take(getBiomass); Y <- take(function(s) apply(getYieldGear(s), c(1,3), sum))[, fish]; N <- take(function(s) sapply(dimnames(s@n)[[1]], function(tt) s@n[tt,,], simplify="array") |> aperm(c(3,1,2)) |> (\(a) { m <- matrix(apply(a, 1, function(x) sum(x["Grey Seal",]*q@dw)/1000), ncol=1); rownames(m) <- dimnames(a)[[1]]; m })())
sealN <- setNames(N[,1], rownames(N)); years <- as.integer(rownames(B))
state <- function(y) { for (s in sims) { tt <- dimnames(s@n)[[1]]; if (as.character(y) %in% tt) return(list(n=s@n[as.character(y),,], npp=s@n_pp[as.character(y),])) } }
eat <- t(sapply(years, function(y) { st <- state(y); quiet(eaten(base, n=st$n, n_pp=st$npp, t=y)) })); rownames(eat) <- years
juv <- t(sapply(years, function(y) { st <- state(y); pr <- getPredRate(base, n=st$n, n_pp=st$npp, t=y)["Grey Seal", idx]; sapply(c("Atlantic Cod","Haddock","American Plaice"), function(s2) { win <- q@w >= 90 & q@w <= 480; wt <- st$n[s2,win]*q@w[win]*q@dw[win]; sum(q@interaction["Grey Seal", s2]*pr[win]*wt)/sum(wt) }) })); rownames(juv) <- years
era <- 1990:1995; n_era <- Reduce(`+`, lapply(era, function(y) state(y)$n))/length(era); pr_era <- Reduce(`+`, lapply(era, function(y) { st <- state(y); getPredRate(base, n=st$n, n_pp=st$npp, t=y)["Grey Seal", idx] }))/length(era)
prey_w <- sapply(fish, function(s) { m <- q@interaction["Grey Seal", s]*pr_era; sum(m*n_era[s,]*q@w*q@dw)/sum(m*n_era[s,]*q@dw) })
saveRDS(list(B=B, Y=Y, eat=eat, sealN=sealN, juv=juv, prey_w=prey_w, dem=dem), file.path(sc, "seal_fixed3_run.rds"))
cat("\nyear | seals model / count (thousand) | fish eaten kt | t per seal | seal mortality on 90-480 g cod / haddock / plaice | total fish 16 spp g/m2\n")
for (y in c(1961,1967,1970,1980,1990,1995,2000,2005,2010)) { yc <- as.character(y); cat(sprintf("%d | %5.0f / %3.0f | %6.1f | %.2f | %.3f / %.3f / %.3f | %.2f\n", y, sealN[yc], counts[yc], sum(eat[yc,]), sum(eat[yc,])/sealN[yc], juv[yc,1], juv[yc,2], juv[yc,3], sum(B[yc, setdiff(fish, forage)])/area)) }
obsB <- readRDS("Ess_Fishmip_data/biomass_ts.RDS") |> transmute(Year=as.integer(Year), Species=as.character(Species), obs=Biomass)
obsC <- readRDS("Ess_Fishmip_data/new_NAFO_4VW_Catches.RDS") |> pivot_longer(-Year, names_to="Species", values_to="obs") |> filter(!is.na(obs), obs>0)
long <- function(m) as.data.frame(as.table(m)) |> setNames(c("Year","Species","mod")) |> mutate(Year=as.integer(as.character(Year)), Species=as.character(Species))
ps <- function(mod, obs, yrs) inner_join(mod, obs, by=c("Year","Species")) |> filter(Year %in% yrs, mod>0, obs>0, !Species %in% forage) |> group_by(Species) |> filter(n()>=4) |> summarise(bias=median(log10(mod/obs)), r=suppressWarnings(cor(log(mod),log(obs))), rmse=sqrt(mean((log10(mod)-log10(obs))^2)), .groups="drop")
f <- function(d) sprintf("x%.2f %2d/%2d r%+.2f rmse %.2f", 10^median(d$bias), sum(abs(d$bias)<log10(2)), nrow(d), median(d$r,na.rm=TRUE), median(d$rmse))
cat(sprintf("\ndata-based dynamic seal | biomass 1970-2004 %s | biomass 2005-2010 %s | catch 1970-2004 %s | catch 2005-2010 %s\n", f(ps(long(B[,fish]), obsB, 1970:2004)), f(ps(long(B[,fish]), obsB, 2005:2010)), f(ps(long(Y), obsC, 1970:2004)), f(ps(long(Y), obsC, 2005:2010))))
rel <- function(s) mean(B[as.character(1994:2004), s])/mean(B[as.character(1978:1981), s]); ob <- function(s) with(obsB[obsB$Species==s,], mean(obs[Year %in% 1994:2004])/mean(obs[Year %in% 1978:1981]))
for (s in c("Atlantic Cod","Haddock","American Plaice","Witch Flounder","Winter Skate","Atlantic Herring","Sand Lance")) cat(sprintf("  %-18s 1994-2004 rel. to 1978-81: model %.2f | observed %.2f\n", s, rel(s), ob(s)))
