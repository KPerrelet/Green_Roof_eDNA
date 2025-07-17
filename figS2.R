################################################################################
#################################### fig s2 ####################################
################################################################################
library(vegan)
library(ggplot2)
library(ggpubr)

set.seed(1) # for reproducibility 

# Import data
cntrl.df <- read.csv("skyline_controlmat.csv")

# Perform Non-metric Multidimensional Scaling (NMDS) on the dataset, excluding "year" and "site_no" columns
cntrl.nmds <- metaMDS(cntrl.df[, -which(colnames(cntrl.df) %in% c("year", "site_no"))],
                      distance = "jaccard",
                      k = 2,
                      maxit = 1000,
                      try = 1000,
                      trymax = 1000,
                      wascores = TRUE)

# Extract NMDS site scores (coordinates for each sample in NMDS space)
cntrl.scrs <- as.data.frame(scores(cntrl.nmds, display = "sites"))
cntrl.scrs$year <- as.factor(cntrl.df$year)
cntrl.scrs$site_no <- as.factor(cntrl.df$site_no)

# Perform PERMANOVA (Permutational Multivariate Analysis of Variance) to test the effect of 'year' on community composition
permanova_type <- adonis2(cntrl.df[, -which(colnames(cntrl.df) %in% c("year", "site_no"))] ~ cntrl.df$year, method = "jaccard")

plot1 <- ggplot(cntrl.scrs, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = year, shape = site_no), size = 3) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, fill = year),
               type = "norm", level = 0.9, geom = "polygon",
               linetype = "solid", alpha = 0.3) +
  annotate(geom = "text", label = paste0("stress = ", signif(cntrl.nmds$stress, 4)), x = Inf, y = Inf, hjust = 1, vjust = 1, size = 5) + 
  scale_color_brewer(name = "Control", palette = "Set1") +
  scale_fill_brewer(name = "Control", palette = "Set1") +
  scale_shape_manual(values = c(15, 17, 19, 21)) + 
  theme_pubr(base_size = 18) +
  annotate("text", 
           label = paste("Permutation test: p-value =", round(permanova_type$`Pr(>F)`[1], digits = 3)), size = 5, 
           x = Inf, y = Inf, hjust = 1, vjust = 3) + 
  guides(shape = guide_legend(title = "Site number"),
         color = guide_legend(title = "Year"), 
         fill = guide_legend(title = "Year"))

# Perform PERMANOVA grouped by 'site_no' to test if samples cluster by site number
permanova_type <- adonis2(cntrl.df[, -which(colnames(cntrl.df) %in% c("year", "site_no"))] ~ cntrl.df$site_no, method = "jaccard")

plot2 <- ggplot(cntrl.scrs, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(color = site_no, shape = site_no), size = 3) +
  stat_ellipse(aes(x = NMDS1, y = NMDS2, fill = site_no),
               type = "norm", level = 0.9, geom = "polygon",
               linetype = "solid", alpha = 0.3) + 
  annotate(geom = "text", label = paste0("stress = ", signif(cntrl.nmds$stress, 4)), x = Inf, y = Inf, hjust = 1, vjust = 1, size = 5) + 
  scale_color_brewer(name = "Control", palette = "Dark2") +
  scale_fill_brewer(name = "Control", palette = "Dark2") +
  scale_shape_manual(values = c(15, 17, 19, 21)) + 
  theme_pubr(base_size = 18) +
  annotate("text", 
           label = paste("Permutation test: p-value =", round(permanova_type$`Pr(>F)`[1], digits = 3)), size = 5, 
           x = Inf, y = Inf, hjust = 1, vjust = 3) + 
  guides(shape = guide_legend(title = "Site number"), 
         fill = guide_legend(title = "Site number"), 
         color = guide_legend(title = "Site number"))

ggarrange(plot1, plot2, common.legend = F)

# Calculate species richness for each sample and add it as a new column to the dataframe
cntrl.df$sp_richness <- specnumber(cntrl.df[, -which(colnames(cntrl.df) %in% c("year", "site_no"))])

# Summarize species richness by 'site_no' and 'year', calculating the mean and standard deviation
summary.df <- cntrl.df %>%
  group_by(site_no, year) %>%
  dplyr::summarise(mean_sp_richness = mean(sp_richness, na.rm = TRUE),
                   sd_sp_richness = sd(sp_richness, na.rm = TRUE),
                   .groups = "drop")

plot3 <- ggplot(summary.df, aes(x = factor(site_no), y = mean_sp_richness, fill = factor(year))) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  geom_errorbar(aes(ymin = mean_sp_richness - sd_sp_richness, ymax = mean_sp_richness + sd_sp_richness),
                position = position_dodge(width = 0.9), width = 0.25) +
  scale_fill_brewer(palette = "Set1") +  
  scale_x_discrete(labels=c("Control A" = "A", "Control B" = "B", "Control C" = "C", "Control D" = "D")) +
  labs(x = "Control", y = "Species richness", fill = "Year") +
  theme_pubr(base_size = 18) 

ggarrange(plot1, plot2, plot3, widths = c(1, 1, 0.7), nrow = 1, common.legend = T, labels = c("A", "B", "C"))
