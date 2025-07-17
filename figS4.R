################################################################################
#################################### fig s3 ####################################
################################################################################
library(lme4)
library(lmerTest)
library(ggplot2)
library(ggforce)
library(ggpubr)
library(vegan)

# Import data
samplmat.df <- read.csv("skyline_spmat_per_sample.csv", check.names = F)
meta <- c("rowid", "site_no", "sample", "replicate", "type", "sp_richness")
samplmat.df$type <- factor(samplmat.df$type, levels = c("Roof", "Ground"))  # Ensure correct order

# Fit a negative binomial mixed model to predict species richness, with random intercepts for site_no
glmm <- glmer.nb(sp_richness ~ type + (1 | site_no), samplmat.df)
summary(glmm)

# Extract the p-value for the effect of type "Ground" from the model summary
p_value <- round(coef(summary(glmm))[, "Pr(>|z|)"]["typeGround"], 3)
p_value <- ifelse(p_value == 0, "<0.001", p_value)

alpha.plot <- ggplot(samplmat.df, aes(x = type, y = sp_richness, color = type)) +
  geom_violin(trim = T, alpha = 0, size = 1) +
  geom_boxplot(width = 0.2, lwd = 0.8) +
  labs(x = "", y = "Species richness") +
  geom_jitter(shape = 16, position = position_jitter(0.2), aes(color = type), alpha = 0.4) +
  scale_color_manual(values = c("forestgreen", "sienna4"), labels = c("Green roof", "Ground-level")) +
  scale_fill_manual(values = c("forestgreen", "sienna4"), labels = c("Green roof", "Ground-level")) +  
  scale_x_discrete(labels=c("Roof" = "Green\nroof", "Ground" = "Ground-\nlevel")) +
  annotate("text", x = Inf, y = Inf,
           label = paste("p-value = ", p_value),
           hjust = 1, vjust = 1, size = 7) +
  theme_pubr(base_size = 18) +
  guides(fill = "none", color = "none")

# Perform Principal Coordinates Analysis (PCoA) using Jaccard distance for community composition
sample.pcoa <- cmdscale(vegdist(samplmat.df[, -which(colnames(samplmat.df) %in% meta)], method = "jaccard", binary = T), eig = T)

# Calculate the percentage variance explained by the first two principal coordinates
sample.pcoa1 <- sample.pcoa$eig[1]/sum(sample.pcoa$eig) * 100
sample.pcoa2 <- sample.pcoa$eig[2]/sum(sample.pcoa$eig) * 100
sample.pcoa <- as.data.frame(sample.pcoa$points)
sample.pcoa$type <- samplmat.df$type

sample.pcoa$site_no <- as.factor(paste0(samplmat.df$site_no, samplmat.df$type))

beta.plot <- ggplot(sample.pcoa, aes(x=V1, y=V2)) +
  geom_point(aes(color = type), size = 3) +
  geom_mark_hull(aes(fill = site_no), alpha = 0, lwd = 0.7, expand = 0, radius = 0, alpha = 0.1) +
  stat_ellipse(aes(group = type, fill = type),
               type = "norm", level = 0.95, geom = "polygon",
               linetype = "solid", alpha = 0.3)  +
  scale_color_manual(values = c("forestgreen", "sienna4"), labels = c("Green roof", "Ground-level")) +
  scale_fill_manual(values = c(rep(c("forestgreen", "sienna4"), 53))) +
  labs(x = paste0("PCoA1 (", round(sample.pcoa1, 2), "%)"), y = paste0("PCoA2 (", round(sample.pcoa2, 2), "%)")) +
  theme_pubr(base_size = 18) +
  theme(text = element_text(size = 18)) +
  guides(fill = "none")

# Perform PERMANOVA (permutational multivariate analysis of variance) to assess differences between 'type' (Roof vs. Ground)
adonis2(samplmat.df[, -which(colnames(samplmat.df) %in% meta)] ~ samplmat.df$type, method = "jaccard")

# Calculate the distance to centroid for each group ('Roof' vs 'Ground') using the betadisper function
sample.dist <- betadisper(vegdist(samplmat.df[, -which(colnames(samplmat.df) %in% meta)], method = "jaccard", binary = T), 
                          samplmat.df$type, 
                          type = "centroid") # distance to centroid in each group (soil / GR)

# Prepare the data for visualization
sample.dist <- data.frame(sample.dist$distances, sample.dist$group,"Within type")
colnames(sample.dist) <- c("distance", "type", "scale")
sample.dist$site_no <- samplmat.df$site_no
sample.dist <- sample.dist %>%
  dplyr::group_by(type) %>%
  dplyr::mutate(mean = mean(distance), sd = sd(distance)) # Calculate mean and standard deviation for distances

# Fit a linear mixed model (LMM) to assess differences in distances to centroids by 'type'
lmm <- lmerTest::lmer(distance ~ type + (1 | site_no), data = sample.dist)
summary(lmm)

# Extract the p-value for the 'typeGround' effect
p_value <- round(coef(summary(lmm))[, "Pr(>|t|)"]["typeGround"], 3)
p_value <- ifelse(p_value == 0, "<0.001", p_value)

dist.plot <- ggplot(sample.dist, aes(x = type, y = distance, color = type)) +
  geom_violin(trim = T, alpha = 0, size = 1) +
  labs(x = "", y = "Distance to type centroid") +
  geom_boxplot(width = 0.2, lwd = 0.8) +
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.3) +
  scale_color_manual(values = c("forestgreen", "sienna4")) +
  scale_x_discrete(labels=c("Roof" = "Green\nroof", "Ground" = "Ground-\nlevel")) +
  theme_pubr(base_size = 18) +
  annotate("text", x = Inf, y = Inf,
           label = paste("p-value = ", p_value),
           hjust = 1, vjust = 1, size = 7) + 
  guides(fill = "none", color = "none")

# Calculate distance to centroid for 'Roof' and 'Ground' groups separately
roof.dist <- betadisper(vegdist(samplmat.df[samplmat.df$type == "Roof", -which(colnames(samplmat.df) %in% meta)], method = "jaccard", binary = T),
                        samplmat.df[samplmat.df$type == "Roof", ]$site_no, 
                        type = "centroid")
ground.dist <- betadisper(vegdist(samplmat.df[samplmat.df$type == "Ground", -which(colnames(samplmat.df) %in% meta)], method = "jaccard", binary = T),
                          samplmat.df[samplmat.df$type == "Ground", ]$site_no, 
                          type = "centroid")

# Prepare distance data for 'Roof' and 'Ground' groups
roof.dist <- data.frame(distance = roof.dist$distances,
                        type = "Roof",
                        scale = "Within sites",
                        site_no = samplmat.df[samplmat.df$type == "Roof", ]$site_no) %>%
  dplyr::group_by(type) %>%
  dplyr::mutate(mean = mean(distance), sd = sd(distance))

ground.dist <- data.frame(distance = ground.dist$distances,
                          type = "Ground",
                          scale = "Within sites",
                          site_no = samplmat.df[samplmat.df$type == "Ground", ]$site_no) %>%
  dplyr::group_by(type) %>%
  dplyr::mutate(mean = mean(distance), sd = sd(distance))

# Combine the distance data for both 'Roof' and 'Ground' groups
all.dist <- rbind(roof.dist, ground.dist) 
all.dist$type <- factor(all.dist$type, levels = c("Roof", "Ground"))  # Ensure correct order

# Fit a linear mixed model to assess differences in distances to site centroids by 'type'
lmm <- lmer(distance ~ type + (1 | site_no), data = all.dist)
summary(lmm)

# Extract the p-value for the 'typeGround' effect in the new model
p_value <- round(coef(summary(lmm))[, "Pr(>|t|)"]["typeGround"], 3)
p_value <- ifelse(p_value == 0, "<0.001", p_value)

dist.plot2 <- ggplot(all.dist, aes(x = type, y = distance, color = type)) +
  geom_violin(trim = T, alpha = 0, size = 1) +
  labs(x = "", y = "Distance to site centroid") +
  geom_boxplot(width = 0.2, lwd = 0.8) +
  geom_jitter(shape = 16, position = position_jitter(0.2), alpha = 0.3) +
  scale_color_manual(values = c("forestgreen", "sienna4")) +
  scale_x_discrete(labels=c("Roof" = "Green\nroof", "Ground" = "Ground-\nlevel")) +
  theme_pubr(base_size = 18) +
  annotate("text", x = Inf, y = Inf,
           label = paste("p-value = ", p_value),
           hjust = 1, vjust = 1, size = 7) +
  guides(fill = "none", color = "none")

ggarrange(beta.plot, 
          ggarrange(alpha.plot, dist.plot, dist.plot2, common.legend = T, nrow = 1, labels = c("B", "C", "D"), font.label = list(size = 20)),
          nrow = 2,
          labels = c("A", ""),
          font.label = list(size = 20), 
          heights = c(1, 0.8),
          common.legend = T)
