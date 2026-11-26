suppressMessages(library(mizer)); quiet <- function(x) suppressWarnings(suppressMessages(x)); sc <- Sys.getenv("SC")
# The full model to take the seal from, and the feeding level to fit its demography at.
# Defaults are the first trial's; script 8 writes seal_on_bowen_diet.rds with both.
seal_source <- Sys.getenv("SEAL_SOURCE", "seal_calibrated.rds")
C <- readRDS(file.path(sc, seal_source)); full <- C$params
dg <- readRDS("/Users/juliab6/Desktop/Ez June 2026/ESS_Restructured/Ess_data/seal_abundance_digitised.RDS"); counts <- setNames(dg$scotian_mean, dg$year)[as.character(1961:2010)]
# ---- a seal-only model on the same size grid, fed at a fixed feeding level ----
row <- as.data.frame(full@species_params["Grey Seal", ]); keep <- c("species","w_min","w_mat","w_max","beta","sigma","h","ks","k","z0","alpha","n","p","w_mat25","w_repro_max","m","erepro","R_max")
row <- row[, intersect(keep, names(row))]; row$interaction_resource <- 0
gp <- data.frame(species="Grey Seal", gear="Seal_Harvest", sel_func="knife_edge", knife_edge_size=16000, catchability=1)
dx <- log10(full@w[2]/full@w[1]); first_bin <- full@w[which(full@w <= 16000 * (1 + 1e-9)) |> max()]; row$w_min <- first_bin
p <- quiet(newMultispeciesParams(row, gear_params=transform(gp, knife_edge_size=first_bin), initial_effort=c(Seal_Harvest=0), min_w=first_bin, max_w=max(full@w), no_w=round(log10(max(full@w)/first_bin)/dx) + 1, n=3/4))
cat("seal-only grid:", length(p@w), "bins from", round(min(p@w)), "to", round(max(p@w)), "g, same spacing as the full model:", isTRUE(all.equal(log10(p@w[2]/p@w[1]), dx)), "\n")
feeding_level <- if (is.null(C$seal_feeding_level)) 0.40 else C$seal_feeding_level; ee <- ext_encounter(p); ee[1, ] <- feeding_level/(1-feeding_level) * row$h * p@w^row$n; ext_encounter(p) <- ee
cat("seal-only model: feeding level", round(mean(getFeedingLevel(p)[1, p@w >= 16000]), 2), "| background mortality", round(row$z0, 3), "per year | matures at", row$w_mat/1000, "kg\n")
numbers <- function(n) sum(n * p@dw)/1000
set_demography <- function(erepro, R_max) { q <- p; q@species_params$erepro <- erepro; q@species_params$R_max <- R_max; q }
# growth rate with no ceiling and no harvest, as a function of erepro
growth_rate <- function(erepro) { q <- set_demography(erepro, Inf); q@initial_n[1, ] <- q@initial_n[1, ] * 1e-6; s <- quiet(project(q, t_max=80, effort=0, progress_bar=FALSE)); N <- getN(s)[, 1]; unname(log(N["80"]/N["50"])/30) }
cat("unharvested growth rate per year at erepro 1:", round(growth_rate(1), 3), "| at 0.3:", round(growth_rate(0.3), 3), "| at 0.1:", round(growth_rate(0.1), 3), "| at 0.03:", round(growth_rate(0.03), 3), "\n")
# stable size structure for a growing population, used as the 1961 starting shape
shape <- { q <- set_demography(0.3, Inf); q@initial_n[1, ] <- q@initial_n[1, ]*1e-6; s <- quiet(project(q, t_max=60, effort=0, progress_bar=FALSE)); s@n[dim(s@n)[1], 1, ] }
years <- 1961:2010; harvest_years <- 1967:1990
simulate <- function(erepro, harvest_rate, R_max) { q <- set_demography(erepro, R_max); q@initial_n[1, ] <- shape * counts["1961"]/numbers(shape)
  eff <- matrix(ifelse(years %in% harvest_years, harvest_rate, 0), ncol=1, dimnames=list(time=years, gear="Seal_Harvest")); s <- quiet(project(q, effort=eff, progress_bar=FALSE)); setNames(getN(s)[, 1]/1000, years) }
sse <- function(th) { N <- tryCatch(simulate(exp(th[1]), max(th[2], 0), exp(th[3])), error=function(e) NULL); if (is.null(N) || any(!is.finite(N)) || any(N <= 0)) return(1e6); sum((log(N) - log(counts))^2) + ifelse(th[2] < 0, 1e3*th[2]^2, 0) + ifelse(th[1] > 0, 1e3*th[1]^2, 0) }
fit_h  <- optim(c(log(0.2), 0.03, log(5e3)), sse, control=list(maxit=400, reltol=1e-6))
fit_0  <- optim(c(log(0.2), log(5e3)), function(t2) sse(c(t2[1], 0, t2[2])), control=list(maxit=400, reltol=1e-6))
est <- c(erepro=exp(fit_h$par[1]), harvest_rate=max(fit_h$par[2],0), R_max=exp(fit_h$par[3])); est0 <- c(erepro=exp(fit_0$par[1]), R_max=exp(fit_0$par[2]))
Nh <- simulate(est["erepro"], est["harvest_rate"], est["R_max"]); N0 <- simulate(est0["erepro"], 0, est0["R_max"])
cat(sprintf("\nwith harvest 1967-1990:  erepro %.3f | harvest rate %.3f per year | R_max %.0f | unharvested growth rate %.3f per year | log RMSE %.3f\n", est["erepro"], est["harvest_rate"], est["R_max"], growth_rate(est["erepro"]), sqrt(fit_h$value/length(years))))
cat(sprintf("no harvest:              erepro %.3f | R_max %.0f | growth rate %.3f per year | log RMSE %.3f\n", est0["erepro"], est0["R_max"], growth_rate(est0["erepro"]), sqrt(fit_0$value/length(years))))
cat("\nyear   count   with harvest   no harvest   (thousand)\n"); for (y in c(1961,1967,1970,1975,1980,1985,1990,1995,2000,2005,2010)) cat(sprintf("%d  %6.0f  %10.0f  %10.0f\n", y, counts[as.character(y)], Nh[as.character(y)], N0[as.character(y)]))
taken <- est["harvest_rate"] * Nh[as.character(harvest_years)] * 1000; cat(sprintf("\nimplied seals removed per year: %.0f in 1967, %.0f in 1980, %.0f in 1990; total 1967-1990 about %.0f\n", taken["1967"], taken["1980"], taken["1990"], sum(taken)))
saveRDS(list(with_harvest=est, no_harvest=est0, Nh=Nh, N0=N0, shape=shape, feeding_level=feeding_level), file.path(sc, "seal_demography.rds"))
