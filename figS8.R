################################################################################
#################################### fig s8 ####################################
################################################################################
library(dplyr)

spmat.df <- read.csv("skyline_spmat.csv", check.names = F)
taxo.df <- read.csv("skyline_taxo.csv")
# shared.df <- read.csv("skyline_sharedmat.csv")
meta <- c("site_no", "type", "sp_richness")
# roof.covariates.df <- read.csv("skyline_covmat_greenroof.csv")
# ground.covariates.df <- read.csv("skyline_covmat_groundlevel.csv")

spmat.df <- spmat.df %>%
  filter(type == "Roof") %>%
  select(-any_of(meta)) %>%
  select(which(colSums(.) > 0))

spmat.df[spmat.df > 1] <- 1
spmat.long <- data.frame(taxa = colnames(spmat.df),
                         no_sites = colSums(spmat.df))
spmat.long$mobility <- taxo.df[taxo.df$taxa %in% spmat.long$taxa, ]$mobility
spmat.long$mobility <- factor(spmat.long$mobility, levels = c("flying", "epigeic", "edaphic"),
                              labels = c("Flying", "Epigeic", "Edaphic"))

ggplot(spmat.long, aes(x = no_sites, fill = mobility)) +
  geom_density(alpha = 0.3, linetype = "solid") +
  geom_histogram(aes(y = stat(density)), alpha = 0.5, bins = 54, color = "black", position = "identity") +
  scale_fill_manual(name = "",
                    values = c("Flying" = "#2AB7CA", "Epigeic" = "#FE4A49", "Edaphic" = "#D4B483")) +
  labs(x = "Number of sites", y = "Species density") +
  theme_pubr(base_size = 14) +
  facet_wrap(~ mobility) +
  theme(strip.background = element_blank(), 
        strip.text = element_blank())