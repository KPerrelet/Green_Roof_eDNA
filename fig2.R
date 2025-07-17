################################################################################
##################################### fig 2 ####################################
################################################################################
library(vegan)
library(tidyr)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(ggsignif)

# Function to compute pairwise Jaccard distances and format as dataframe
make_dist_df <- function(df) {
  df.dist <- df %>%
    dplyr::select(-any_of(meta)) %>%  # Remove metadata columns
    vegdist(method = "jaccard", binary = T) # Compute Jaccard distance
  
  df.dist.df <- data.frame(as.vector(df.dist),
                           which(lower.tri(df.dist),
                                 arr.ind = TRUE))
  colnames(df.dist.df) <- c("distance", "site2", "site1")
  df.dist.df$type1 <- df$type[df.dist.df$site1]
  df.dist.df$type2 <- df$type[df.dist.df$site2]
  df.dist.df$type1to2 <- paste(df.dist.df$type1, df.dist.df$type2, sep = "-")
  
  return(df.dist.df)
}

# Import data
spmat.df <- read.csv("skyline_spmat.csv", check.names = F)
meta <- c("site_no", "type", "sp_richness")
taxo.df <- read.csv("skyline_taxo.csv")

# Subset species data based on mobility modes (i.e., flying, epigeic, edaphic)
flying.df <- spmat.df[, which(colnames(spmat.df) %in% taxo.df[taxo.df$mobility == "flying", ]$taxa)]
flying.df$sp_richness <- specnumber(flying.df)
flying.df$type <- spmat.df$type

epigeic.df <- spmat.df[, which(colnames(spmat.df) %in% taxo.df[taxo.df$mobility == "epigeic", ]$taxa)]
epigeic.df$sp_richness <- specnumber(epigeic.df)
epigeic.df$type <- spmat.df$type

edaphic.df <- spmat.df[, which(colnames(spmat.df) %in% taxo.df[taxo.df$mobility == "edaphic", ]$taxa)]
edaphic.df$sp_richness <- specnumber(edaphic.df)
edaphic.df$type <- spmat.df$type

# Combine richness data into a single dataframe for visualization
summary.df <- rbind(data.frame(richness = spmat.df$sp_richness, type = spmat.df$type, mobility = "All"), 
                     data.frame(richness = flying.df$sp_richness, type = flying.df$type, mobility = "Flying"), 
                     data.frame(richness = epigeic.df$sp_richness, type = epigeic.df$type, mobility = "Epigeic"), 
                     data.frame(richness = edaphic.df$sp_richness, type = edaphic.df$type, mobility = "Edaphic"))

# Compute Wilcoxon test p-values per mobility group and apply BH correction
p_values <- summary.df %>%
  dplyr::group_by(mobility) %>%
  dplyr::summarise(p = wilcox.test(richness ~ type, data = cur_data())$p.value)

p_values$p_adj <- p.adjust(p_values$p, method = "BH")
p_values$signif_label <- ifelse(p_values$p_adj < 0.1, ifelse(p_values$p_adj < 0.05, ifelse(p_values$p_adj < 0.01, ifelse(p_values$p_adj < 0.001, "***","**"),"*"),"."), "")

summary.df$type <- factor(summary.df$type, levels = c("Roof", "Ground"))  # Ensure correct order
alpha.plot <- ggplot(summary.df, aes(x = type, y = richness, color = type)) +
  geom_boxplot(width = 0.5, lwd = 0.8) +
  ylab("Species richness") +
  scale_color_manual(values = c("forestgreen", "sienna4"), name = "", labels = c("Green roof", "Ground-level")) +
  theme_pubr(base_size = 18) +
  expand_limits(y = max(summary.df$richness) * 1.2) + 
  stat_signif( # The results do not change with or without Benjamini-Hochberg correction so we can use a non-corrected wilcoxon test for visualization purposes
    test = "wilcox.test",
    comparisons = list(c("Roof", "Ground")),
    map_signif_level = TRUE,
    margin_top = -0.1,
    color = "black",
    size = 1,
    textsize = 6
  ) +
  facet_wrap(~factor(mobility, c("All", "Flying", "Epigeic", "Edaphic")), nrow = 1) +
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(face = "bold", size = 16))

# Principal Coordinates Analysis (PCoA) based on Jaccard distance for different mobility groups
all.pcoa <- cmdscale(vegdist(spmat.df[, -which(colnames(spmat.df) %in% meta)], method = "jaccard", binary = T), eig = T)
all.pcoa1 <- all.pcoa$eig[1]/sum(all.pcoa$eig) * 100
all.pcoa2 <- all.pcoa$eig[2]/sum(all.pcoa$eig) * 100
all.pcoa <- as.data.frame(all.pcoa$points)
all.pcoa$type <- spmat.df$type
all.pcoa$type <- factor(all.pcoa$type, levels = c("Roof", "Ground"))  # Ensure correct order

flying.pcoa <- cmdscale(vegdist(flying.df[, -which(colnames(flying.df) %in% meta)], method = "jaccard", binary = T), eig = T)
flying.pcoa1 <- flying.pcoa$eig[1]/sum(flying.pcoa$eig) * 100
flying.pcoa2 <- flying.pcoa$eig[2]/sum(flying.pcoa$eig) * 100
flying.pcoa <- as.data.frame(flying.pcoa$points)
flying.pcoa$type <- spmat.df$type
flying.pcoa$type <- factor(flying.pcoa$type, levels = c("Roof", "Ground"))  # Ensure correct order

epigeic.df$site_No <- spmat.df$site_No
epigeic.df <- epigeic.df[epigeic.df$sp_richness >= 1, ] #We remove sites with no epigeic species as we cannot calculate a Jaccard distance matrix otherwise
epigeic.pcoa <- cmdscale(vegdist(epigeic.df[, -which(colnames(epigeic.df) %in% meta)], method = "jaccard", binary = T), eig = T)
epigeic.pcoa1 <- epigeic.pcoa$eig[1]/sum(epigeic.pcoa$eig) * 100
epigeic.pcoa2 <- epigeic.pcoa$eig[2]/sum(epigeic.pcoa$eig) * 100
epigeic.pcoa <- as.data.frame(epigeic.pcoa$points)
epigeic.pcoa$type <- epigeic.df$type
epigeic.pcoa$type <- factor(epigeic.pcoa$type, levels = c("Roof", "Ground"))  # Ensure correct order


edaphic.pcoa <- cmdscale(vegdist(edaphic.df[, -which(colnames(edaphic.df) %in% meta)], method = "jaccard", binary = T), eig = T)
edaphic.pcoa1 <- edaphic.pcoa$eig[1]/sum(edaphic.pcoa$eig) * 100
edaphic.pcoa2 <- edaphic.pcoa$eig[2]/sum(edaphic.pcoa$eig) * 100
edaphic.pcoa <- as.data.frame(edaphic.pcoa$points)
edaphic.pcoa$type <- spmat.df$type
edaphic.pcoa$type <- factor(edaphic.pcoa$type, levels = c("Roof", "Ground"))  # Ensure correct order

all.pcoa.plot <- ggplot(all.pcoa, aes(x=V1, y=V2)) +    
  geom_point(aes(color = type), size = 3) +
  stat_ellipse(aes(group = type, fill = type),
               type = "norm", level = 0.95, geom = "polygon",
               linetype = "solid", alpha = 0.3) + 
  scale_color_manual(values = c("forestgreen", "sienna4"), labels = c("Green roof", "Ground-level")) +
  scale_fill_manual(values = c("forestgreen", "sienna4"), labels = c("Green roof", "Ground-level")) +
  labs(x = paste0("PCoA1 (", round(all.pcoa1, 2), "%)"), y = paste0("PCoA2 (", round(all.pcoa2, 2), "%)")) + 
  theme_pubr(base_size = 18)

flying.pcoa.plot <- ggplot(flying.pcoa, aes(x=V1, y=V2)) +    
  geom_point(aes(color = type), size = 3) +
  stat_ellipse(aes(group = type, fill = type),
               type = "norm", level = 0.95, geom = "polygon",
               linetype = "solid", alpha = 0.3) + 
  scale_color_manual(values = c("forestgreen", "sienna4"), labels = c("Green roof", "Ground-level")) +
  scale_fill_manual(values = c("forestgreen", "sienna4"), labels = c("Green roof", "Ground-level")) +
  labs(x = paste0("PCoA1 (", round(flying.pcoa1, 2), "%)"), y = paste0("PCoA2 (", round(flying.pcoa2, 2), "%)")) + 
  theme_pubr(base_size = 18)

epigeic.pcoa.plot <- ggplot(epigeic.pcoa, aes(x=V1, y=V2)) +    
  geom_point(aes(color = type), size = 3) +
  stat_ellipse(aes(group = type, fill = type),
               type = "norm", level = 0.95, geom = "polygon",
               linetype = "solid", alpha = 0.3) + 
  scale_color_manual(values = c("forestgreen", "sienna4"), labels = c("Green roof", "Ground-level")) +
  scale_fill_manual(values = c("forestgreen", "sienna4"), labels = c("Green roof", "Ground-level")) +
  labs(x = paste0("PCoA1 (", round(epigeic.pcoa1, 2), "%)"), y = paste0("PCoA2 (", round(epigeic.pcoa2, 2), "%)")) + 
  theme_pubr(base_size = 18)

edaphic.pcoa.plot <- ggplot(edaphic.pcoa, aes(x=V1, y=V2)) +    
  geom_point(aes(color = type), size = 3) +
  stat_ellipse(aes(group = type, fill = type),
               type = "norm", level = 0.95, geom = "polygon",
               linetype = "solid", alpha = 0.3) + 
  scale_color_manual(values = c("forestgreen", "sienna4"), labels = c("Green roof", "Ground-level")) +
  scale_fill_manual(values = c("forestgreen", "sienna4"), labels = c("Green roof", "Ground-level")) +
  labs(x = paste0("PCoA1 (", round(edaphic.pcoa1, 2), "%)"), y = paste0("PCoA2 (", round(edaphic.pcoa2, 2), "%)")) + 
  theme_pubr(base_size = 18)

# Generate pairwise distance matrices for different mobility groups
all.dist.df <- make_dist_df(spmat.df)
flying.dist.df <- make_dist_df(flying.df)
epigeic.dist.df <- make_dist_df(epigeic.df)
edaphic.dist.df <- make_dist_df(edaphic.df)

# Combine distances into a single dataframe for visualization
summary.dist.df <- rbind(data.frame(all.dist.df, mobility = "All"), 
                         data.frame(flying.dist.df, mobility = "Flying"), 
                         data.frame(epigeic.dist.df, mobility = "Epigeic"), 
                         data.frame(edaphic.dist.df, mobility = "Edaphic"))

# Compute Wilcoxon test p-values per mobility group and apply BH correction
p_values <- summary.dist.df %>%
  dplyr::filter(type1to2 %in% c("Roof-Roof", "Ground-Ground")) %>%
  dplyr::group_by(mobility) %>%
  dplyr::summarise(p = wilcox.test(distance ~ type1to2, data = cur_data())$p.value)

p_values$p_adj <- p.adjust(p_values$p, method = "BH")
p_values$signif_label <- ifelse(p_values$p_adj < 0.1, ifelse(p_values$p_adj < 0.05, ifelse(p_values$p_adj < 0.01, ifelse(p_values$p_adj < 0.001, "***","**"),"*"),"."), "")

beta.plot <- ggplot(summary.dist.df[summary.dist.df$type1to2 != "Roof-Ground", ], aes(x = type1to2, y = distance, color = type1to2)) +
  labs(y = "Dissimilarity", x = "Pairwise comparison") + 
  geom_boxplot(width = 0.5, lwd = 0.8) +
  scale_color_manual(values = c("forestgreen", "sienna4"), 
                     name = "", labels = c("Roof-roof", "Ground-ground")) +
  theme_pubr(base_size = 18) +
  expand_limits(y = max(edaphic.dist.df$distance) * 1.3) +  
  stat_signif( # Again, the results do not change with the correction so we can use uncorrected wilcoxon tests for visualization
    test = "wilcox.test",
    comparisons = list(
      c("Roof-Roof", "Ground-Ground")
    ),
    map_signif_level = TRUE,
    margin_top = -0.1,    
    color = "black",
    step_increase = 0.1,
    size = 1,
    textsize = 6
  ) + 
  facet_wrap(~factor(mobility, c("All", "Flying", "Epigeic", "Edaphic")), nrow = 1) +
  theme(axis.text.x = element_blank(), 
        axis.ticks = element_blank(), 
        strip.background = element_blank(), 
        strip.text = element_text(face = "bold", size = 16),
        legend.position = "bottom")

ggarrange(NULL, 
          ggarrange(all.pcoa.plot, flying.pcoa.plot, epigeic.pcoa.plot, edaphic.pcoa.plot, 
                    legend = "none",
                    labels = c("All species", "Flying species", "Epigeic species", "Edaphic species"), 
                    label.y = 1.05, 
                    font.label = list(size = 16, face = "bold"), 
                    hjust = -0.2, vjust = 1.5), 
          ggarrange(alpha.plot, beta.plot, align = "hv", common.legend = T, legend = "bottom", labels = c("B", "C"), 
                    label.y = 1.05, 
                    font.label = list(size = 20)),
          labels = c("", "A", ""),
          # label.y = 1.3,
          vjust = -2,
          font.label = list(size = 20), 
          legend = "bottom", 
          align = "hv",
          ncol = 1, nrow = 3, 
          # heights = c(0.1, 1, 0.6)
          heights = c(0.1, 1, 0.5)
)
