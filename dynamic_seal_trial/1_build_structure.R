suppressMessages({library(mizer); library(therMizer)})
source("Ess_Fishmip_data/Useful_functions.R")
quiet <- function(x) suppressWarnings(suppressMessages(x))
p0 <- readRDS(Sys.getenv("BASE_PARAMS", "Ess_Fishmip_data/params_GFDL_steady.RDS"))           # FishMIP model, seals as a gear
L <- readRDS("/Users/juliab6/Desktop/Ez June 2026/ESS_Paper2_climate/output/paper2_base_LOCKED.rds")
g <- L@given_species_params; g$w_repro_max[is.na(g$w_repro_max)] <- g$w_max[is.na(g$w_repro_max)]; L@given_species_params <- g; L <- quiet(validParams(L))
fish <- as.character(p0@species_params$species)
seal_row <- L@species_params["Grey Seal", ]
cat("Paper 2 seal parameters:\n"); print(seal_row[, c("w_min","w_mat","w_max","beta","sigma","h","gamma","ks","k","z0","R_max","erepro","alpha","n","p")])
cat("seal prey preference (interaction row):\n"); print(round(L@interaction["Grey Seal", ], 2))
ee <- ext_encounter(L)["Grey Seal", ]; hw <- seal_row$h * L@w^seal_row$n
cat("external food as a fraction of maximum intake (sizes >= 16 kg):", round(unique(round((ee/hw)[L@w >= 16000], 3)), 3), "\n")
cat("exposure array dims:", paste(dim(as.matrix(read.csv("Ess_Fishmip_data/GFDL_ESS_exposure_array.csv", row.names=1))), collapse="x"), "\n")
e <- new.env(); load("Ess_Fishmip_data/GFDL_ESS_vertical_migration_array.Rdata", envir=e); vm <- e$vertical_migration_array; cat("vertical migration dims:", paste(dim(vm), collapse="x"), "| dimnames:", paste(names(dimnames(vm)), collapse=","), "\n")
ex <- as.matrix(read.csv("Ess_Fishmip_data/GFDL_ESS_exposure_array.csv", row.names=1)); print(round(ex[, 1:4], 2)); cat("vm for cod, realm sums over sizes range:", paste(round(range(colSums(vm[, 1, ])), 3), collapse="-"), "\n")
# ---- species table: fish exactly as calibrated (background mortality carried in z0), plus the seal ----
s0 <- p0@species_params; s0$z0 <- ext_mort(p0)[, 1]
cols <- union(names(s0), names(seal_row)); for (cn in setdiff(cols, names(s0))) s0[[cn]] <- NA; for (cn in setdiff(cols, names(seal_row))) seal_row[[cn]] <- NA
sp2 <- rbind(s0[, cols], seal_row[, cols]); sp2$species <- as.character(sp2$species)
inter2 <- rbind(cbind(p0@interaction, `Grey Seal` = 0), `Grey Seal` = c(L@interaction["Grey Seal", fish], 0)); dimnames(inter2) <- list(sp2$species, sp2$species)
gp <- as.data.frame(gear_params(p0)); gp <- gp[gp$gear != "Seals", ]; eff <- p0@initial_effort[names(p0@initial_effort) != "Seals"]
dx <- log10(max(p0@w)/min(p0@w))/(length(p0@w)-1); k <- ceiling(log10(300000/max(p0@w))/dx); no_w <- length(p0@w)+k; max_w <- min(p0@w)*10^(dx*(no_w-1))
rp <- p0@resource_params
p <- quiet(newMultispeciesParams(sp2, interaction=inter2, gear_params=gp, initial_effort=eff, min_w=min(p0@w), max_w=max_w, no_w=no_w, kappa=rp$kappa, lambda=rp$lambda, w_pp_cutoff=rp$w_pp_cutoff, n=3/4))
N0 <- length(p0@w_full); n0 <- length(p0@w)
cat(sprintf("\ngrid: %d -> %d bins (max %.0f g); original bins preserved to %.1e; w_full %d -> %d\n", n0, length(p@w), max(p@w), max(abs(p@w[1:n0]-p0@w)), N0, length(p@w_full)))
# ---- do the rebuilt fish have the same rates as the calibrated ones? ----
chk <- function(nm, a, b) cat(sprintf("  %-22s max relative difference %.2e\n", nm, max(abs(a-b)/pmax(abs(b), 1e-300))))
chk("search volume", p@search_vol[fish, 1:n0], p0@search_vol); chk("max intake", p@intake_max[fish, 1:n0], p0@intake_max); chk("metabolism", p@metab[fish, 1:n0], p0@metab)
chk("background mortality", p@mu_b[fish, 1:n0], p0@mu_b); chk("maturity x repro", p@psi[fish, 1:n0], p0@psi); chk("fishing selectivity", p@selectivity[names(eff), fish, 1:n0], p0@selectivity[names(eff), , ])
chk("catchability", p@catchability[names(eff), fish], p0@catchability[names(eff), ]); chk("erepro", p@species_params[fish,"erepro"], p0@species_params$erepro); chk("R_max", p@species_params[fish,"R_max"], p0@species_params$R_max)
saveRDS(list(p=p, L=L, fish=fish, n0=n0, N0=N0), file.path(Sys.getenv("SC"), "seal_structure.rds"))
