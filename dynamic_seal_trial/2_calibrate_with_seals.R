suppressMessages({library(mizer); library(therMizer)})
source("Ess_Fishmip_data/Useful_functions.R")
quiet <- function(x) suppressWarnings(suppressMessages(x)); area <- 146360555008; sc <- Sys.getenv("SC")
S <- readRDS(file.path(sc, "seal_structure.rds")); p <- S$p; L <- S$L; fish <- S$fish; n0 <- S$n0; N0 <- S$N0
p0 <- readRDS(Sys.getenv("BASE_PARAMS", "Ess_Fishmip_data/params_GFDL_steady.RDS"))
dg <- readRDS("/Users/juliab6/Desktop/Ez June 2026/ESS_Restructured/Ess_data/seal_abundance_digitised.RDS")
cat("digitised seals (thousand):", paste(sprintf("%d %.0f", c(1960,1970,1978,1981,1990,2000,2010), dg$scotian_mean[match(c(1960,1970,1978,1981,1990,2000,2010), dg$year)]), collapse=" | "), "\n")
target_seals <- mean(dg$scotian_mean[dg$year %in% 1978:1981]); cat("calibration target, 1978-81 mean:", round(target_seals,1), "thousand\n")
seal_thousand <- function(q) sum(q@initial_n["Grey Seal", ] * q@dw)/1000
# ---- state: fish and resource exactly as calibrated; seal spectrum shape from Paper 2 ----
nn <- p@initial_n; nn[] <- 0; nn[fish, 1:n0] <- p0@initial_n; stopifnot(isTRUE(all.equal(unname(L@w), unname(p@w))))
nn["Grey Seal", ] <- L@initial_n["Grey Seal", ] * target_seals / (sum(L@initial_n["Grey Seal", ] * L@dw)/1000); p@initial_n <- nn
npp <- p@initial_n_pp; npp[] <- 0; npp[1:N0] <- p0@initial_n_pp; p@initial_n_pp <- npp; resource_dynamics(p) <- "resource_constant"
ee <- ext_encounter(p); ee["Grey Seal", ] <- ext_encounter(L)["Grey Seal", ]; ext_encounter(p) <- ee
# ---- therMizer inputs padded for the seal and the six new size bins ----
u <- as.matrix(read.table("Ess_Fishmip_data/GFDL_resource_spectra_unanchored.dat", header=TRUE, check.names=FALSE))
u217 <- cbind(u, matrix(-100, nrow(u), length(p@w_full)-N0)); dimnames(u217) <- list(time=rownames(u), w=names(p@initial_n_pp))
tm <- read.csv("Ess_Fishmip_data/Temp_by_depth_spin.csv"); rownames(tm) <- as.character(tm$year); tmat <- as.matrix(tm[, setdiff(names(tm),"year"), drop=FALSE])
e <- new.env(); load("Ess_Fishmip_data/GFDL_ESS_vertical_migration_array.Rdata", envir=e); vm0 <- e$vertical_migration_array; ex0 <- as.matrix(read.csv("Ess_Fishmip_data/GFDL_ESS_exposure_array.csv", row.names=1))
sp_all <- as.character(p@species_params$species)
vm <- array(0, dim=c(2, length(sp_all), length(p@w)), dimnames=list(realm=dimnames(vm0)[[1]], sp=sp_all, w=signif(p@w,3)))
vm[, fish, 1:n0] <- vm0; for (j in (n0+1):length(p@w)) vm[, fish, j] <- vm0[, , n0]; vm[, "Grey Seal", ] <- 0.5
ex <- cbind(ex0, `Grey Seal`=1); colnames(ex) <- sp_all
p@species_params$temp_min <- c(p0@species_params$temp_min, -50); p@species_params$temp_max <- c(p0@species_params$temp_max, 80)
ref <- as.character(1978:1981); ny <- 2000; lab <- as.character(seq_len(ny))
const_inputs <- function(q) quiet(upgradeTherParams(q, ocean_temp=matrix(rep(colMeans(tmat[ref,]), each=ny), nrow=ny, dimnames=list(time=lab, realm=colnames(tmat))),
   n_pp_array=matrix(rep(colMeans(u217[ref,]), each=ny), nrow=ny, dimnames=list(time=lab, w=names(q@initial_n_pp))), vertical_migration_array=vm, exposure_array=ex, aerobic_effect=TRUE, metabolism_effect=TRUE))
# ---- seals are endotherms: cancel therMizer's scaling at the reference temperature ----
pt <- const_inputs(p); a <- scaled_temp_effect(pt, 1)["Grey Seal", 1]; m <- therMizer:::thermizer_metabolism_effect(pt, 1)["Grey Seal", 1]
cat(sprintf("therMizer scaling on the seal at the reference temperature: search volume x%.4f, metabolism x%.5f\n", a, m))
for (T in c(3, 5, 7.5)) { tt <- matrix(rep(c(T+1.3, T-1.3), each=ny), nrow=ny, dimnames=list(time=lab, realm=colnames(tmat))); q <- quiet(upgradeTherParams(p, ocean_temp=tt, n_pp_array=other_params(pt)$n_pp_array, vertical_migration_array=vm, exposure_array=ex, aerobic_effect=TRUE, metabolism_effect=TRUE))
  cat(sprintf("  mean temperature %.1f C: search volume factor %.4f (%+.1f%% vs reference), metabolism factor %.5f (%+.0f%%)\n", T, scaled_temp_effect(q,1)["Grey Seal",1], 100*(scaled_temp_effect(q,1)["Grey Seal",1]/a-1), therMizer:::thermizer_metabolism_effect(q,1)["Grey Seal",1], 100*(therMizer:::thermizer_metabolism_effect(q,1)["Grey Seal",1]/m-1))) }
sv <- p@search_vol; sv["Grey Seal", ] <- sv["Grey Seal", ]/a; p <- setSearchVolume(p, search_vol=sv); mb <- p@metab; mb["Grey Seal", ] <- mb["Grey Seal", ]/m; p <- setMetabolicRate(p, metab=mb)
ee <- ext_encounter(p); ee["Grey Seal", ] <- ext_encounter(L)["Grey Seal", ]; ext_encounter(p) <- ee
pt <- const_inputs(p); initialNResource(pt) <- { x <- 10^colMeans(u217[ref,])/pt@dw_full; x[pt@w_full >= pt@resource_params$w_pp_cutoff] <- 0; x }; resource_dynamics(pt) <- "resource_constant"
cat("after compensation, seal search volume kept:", isTRUE(all.equal(pt@search_vol["Grey Seal",], sv["Grey Seal",])), "| resource equals the gear model's:", isTRUE(all.equal(unname(pt@initial_n_pp[1:N0]), unname(p0@initial_n_pp))), "\n")
# ---- calibrate: seals held near the 1978-81 number, fish re-matched to their observed biomass ----
solve <- function(q) quiet(steady(q, t_max=300, progress_bar=FALSE, preserve="erepro"))
rep_line <- function(lbl, q) { b <- getBiomass(q); cat(sprintf("%-30s fish %.3f g/m2 | seals %.1f thousand | max fish erepro %.3g\n", lbl, sum(b[fish])/area, seal_thousand(q), max(q@species_params[fish,"erepro"]))); invisible(q) }
rep_line("start", pt)
q <- pt
for (i in 1:5) { q <- solve(q); rep_line(sprintf("steady %d", i), q)
  rm <- q@species_params$R_max; rm[sp_all=="Grey Seal"] <- rm[sp_all=="Grey Seal"] * target_seals/seal_thousand(q); q@species_params$R_max <- rm; q@initial_n["Grey Seal", ] <- q@initial_n["Grey Seal", ] * target_seals/seal_thousand(q)
  q <- quiet(matchBiomasses(q, species=fish)) }
q <- solve(q); rep_line("final", q)
sp <- q@species_params; ob <- sp[fish,"biomass_observed"]; mo <- getBiomass(q)[fish]; r <- mo/ob
yg <- getYieldGear(q); yo <- sp[fish,"yield_observed"]; ym <- colSums(yg)[fish]; ok <- !is.na(yo) & yo > 0
cat(sprintf("biomass: %d/%d fish within 2x, median ratio %.3f | yield: %d/%d within 2x, median ratio %.3f | fish erepro > 1: %s\n", sum(r>0.5 & r<2), length(r), median(r), sum((ym/yo)[ok]>0.5 & (ym/yo)[ok]<2), sum(ok), median((ym/yo)[ok]), paste(fish[sp[fish,"erepro"]>1], collapse=", ")))
cat("seal R_max at the 1978-81 level:", signif(sp["Grey Seal","R_max"],4), "| seal erepro", signif(sp["Grey Seal","erepro"],3), "| seal feeding level range", paste(round(range(getFeedingLevel(q)["Grey Seal", q@w>=16000]),2), collapse="-"), "\n")
pm <- getPredMort(q); cat("consumption by seals of model fish, kt/yr:", round(sum(sapply(fish, function(s) sum(q@interaction["Grey Seal",s] * getPredRate(q)["Grey Seal", (length(q@w_full)-length(q@w)+1):length(q@w_full)] * q@initial_n[s,] * q@w * q@dw)))/1e9, 1), "\n")
saveRDS(list(params=q, vm=vm, ex=ex, u217=u217, tmat=tmat, target_seals=target_seals), file.path(sc, "seal_calibrated.rds"))
