# Measured grey seal diets at Sable Island, as percent of prey weight
#
# Two periods, both from scats collected on Sable Island:
#
#   1991-1993  Bowen & Harrison (1994), Mar. Ecol. Prog. Ser. 112: 1-11,
#              Table 2, samples from July 1991 to January 1993. The figures
#              are as copied into make_S6_seal_diet.R in the ESS_Restructured
#              project; the paper itself was not to hand.
#   1993-1998  Bowen & Harrison (2007), NAMMCO Sci. Publ. 6: 123-134,
#              Table 2, samples from May 1993 to January 1998.
#
# Only species that are in the model are listed. NA means the table gives no
# figure for that species in that period. "White Hake" in 1993-1998 is the
# paper's "Hake".
seal_diet_percent <- rbind(
  "1991-1993" = c("Sand Lance" = 80.9, "Atlantic Cod" = 11.0, "Silver Hake" = 2.6,
                  "Yellowtail Flounder" = 1.1, "American Plaice" = 1.4,
                  "Longhorn Sculpin" = NA, "Capelin" = 0.4, "Redfish Unseparated" = 1.3,
                  "Haddock" = 0.05, "Pollock" = 0.2, "Atlantic Herring" = 0.2,
                  "White Hake" = 0.2, "Witch Flounder" = 0.4),
  "1993-1998" = c("Sand Lance" = 61.7, "Atlantic Cod" = 13.3, "Silver Hake" = 7.6,
                  "Yellowtail Flounder" = 4.4, "American Plaice" = 4.3,
                  "Longhorn Sculpin" = 2.9, "Capelin" = 2.0, "Redfish Unseparated" = 1.0,
                  "Haddock" = 0.8, "Pollock" = 0.7, "Atlantic Herring" = 0.5,
                  "White Hake" = 0.4, "Witch Flounder" = 0.2))
seal_diet_years <- list("1991-1993" = 1991:1993, "1993-1998" = 1993:1998)

# Shares of the fish eaten, among the model species in the tables.
seal_diet_share_all_listed <- seal_diet_percent / rowSums(seal_diet_percent, na.rm = TRUE)

# The seal is calibrated species by species only for prey that make up at
# least 1% of the diet in one of the periods and have a stock in the model
# large enough to supply it. The rest share one preference, so the seal takes
# them in proportion to what is there. Capelin is left out although it was 2%
# of the 1993-1998 diet: the model holds 0.2 kt of it (the raw trawl figure),
# and asking for 2% forced the seal's highest preference onto it and upset the
# whole calibration. Herring, haddock, pollock, white hake and witch are each
# under 1% in both periods.
calibrated_prey <- c("Sand Lance", "Atlantic Cod", "Silver Hake", "Yellowtail Flounder",
                     "American Plaice", "Longhorn Sculpin", "Redfish Unseparated")
seal_diet_share <- seal_diet_share_all_listed[, calibrated_prey]
share_of_unnamed_species_together <- 1 - mean(rowSums(seal_diet_share, na.rm = TRUE))

# The diet to start the calibration from: the average of the two periods.
starting_seal_diet_share <- colMeans(seal_diet_share, na.rm = TRUE)
share_of_unnamed_species_together <- 1 - sum(starting_seal_diet_share)
