################################################################################
################################ fig 3 + fig s4 ################################
################################################################################
library(dplyr)
library(grid)
library(VennDiagram)
library(stringr)
library(ggnewscale)
library(ggplot2)
library(ggpubr)
library(geomtextpath)

# Function to compute species richness at the class and order level
compute_class_order_counts <- function(site, taxo.df) {
  
  # Initialize empty vectors for class-level counts
  class_counts <- setNames(rep(0, length(unique(taxo.df$class))), unique(taxo.df$class))
  class_fracs <- class_counts
  
  # Compute class-level richness and read fractions
  for (clas in unique(taxo.df$class)) {
    class_counts[clas]  <- nrow(taxo.df[taxo.df$taxa %in% colnames(site) & taxo.df$class == clas, ])
    class_fracs[clas] <- sum(site[, which(colnames(site) %in% taxo.df[taxo.df$class == clas, ]$taxa)])
  }
  
  # Initialize empty vectors for order-level counts (only for Insecta)
  order_counts <- setNames(rep(0, length(unique(taxo.df[taxo.df$class == "Insecta", ]$order))), unique(taxo.df[taxo.df$class == "Insecta", ]$order))
  order_fracs <- order_counts
  
  # Compute order-level richness and read fractions
  for (ordr in unique(taxo.df[taxo.df$class == "Insecta", ]$order)) {
    order_counts[ordr]  <- nrow(taxo.df[taxo.df$taxa %in% colnames(site) & taxo.df$order == ordr, ])
    order_fracs[ordr] <- sum(site[, which(colnames(site) %in% taxo.df[taxo.df$order == ordr, ]$taxa)])
  }
  
  # Combine class and order-level counts, excluding 'Insecta' from the class-level summary
  class_counts <- append(class_counts[setdiff(names(class_counts), "Insecta")], order_counts)
  class_fracs <- append(class_fracs[setdiff(names(class_fracs), "Insecta")], order_fracs)
  return(list(class_counts, class_fracs))
}

# Function to compute species richness based on mobility traits
compute_mobility_counts <- function(site, taxo.df) {
  mobility_counts <- setNames(rep(0, length(unique(taxo.df$mobility))), unique(taxo.df$mobility))
  mobility_fracs <- mobility_counts
  
  # Compute richness and read fractions for each mobility mode
  for (mobi in unique(taxo.df$mobility)) {
    mobility_counts[mobi]  <- nrow(taxo.df[taxo.df$taxa %in% colnames(site) & taxo.df$mobility == mobi, ])
    mobility_fracs[mobi] <- sum(site[, which(colnames(site) %in% taxo.df[taxo.df$mobility == mobi, ]$taxa)])
  }
  return(list(mobility_counts, mobility_fracs))
}

# Function to aggregate taxonomic composition at each site
composition_agg_row <- function(df) {
  res <- data.frame()
  res2 <- data.frame()
  
  # Iterate over each site in the dataset
  for (i in 1:nrow(df)) {
    site <- df[i, ]
    site_no <- site$site_no
    
    # Remove metadata columns and empty species columns
    site <- site %>%
      dplyr::select(-any_of(meta)) %>%
      dplyr::select_if(colSums(.) > 0)
    
    # Compute class/order and mobility-based richness
    class_order <- compute_class_order_counts(site, taxo.df)
    class_counts <- class_order[[1]]
    class_fracs <- class_order[[2]]
    
    mobility <- compute_mobility_counts(site, taxo.df)
    mobility_counts <- mobility[[1]]
    mobility_fracs <- mobility[[2]]
    
    # Create data frames for class/order-level and mobility-level summaries
    temp <- data.frame(site_no = site_no, 
                       taxa = names(class_counts), 
                       richness = unlist(class_counts), 
                       reads = unlist(class_fracs))
    temp2 <- data.frame(site_no = site_no, 
                        taxa_mobility = names(mobility_counts), 
                        richness_mobility = unlist(mobility_counts), 
                        reads_mobility = unlist(mobility_fracs))
    res <- rbind(res, temp)
    res2 <- rbind(res2, temp2)
  }
  
  # Regroup certain taxa and add a postion (ymin, ymax), as well as a label for the doughnut plot
  mobility_counts <- setNames(rep(0, 3 * length(unique(res$taxa))), paste(rep(unique(res$taxa), each = 3), unique(taxo.df$mobility), sep = "_"))
  mobility_fracs <- mobility_counts
  for (taxa in unique(res$taxa)) {
    for (mobi in unique(res2$taxa_mobility)) {
      mobility_counts[paste(taxa, mobi, sep = "_")] <- nrow(taxo.df[(taxo.df$class == taxa | taxo.df$order == taxa) & taxo.df$mobility == mobi & taxo.df$taxa %in% colnames(df), ])
    }
  }
  final_mobility_df <- data.frame(taxa = str_split_fixed(names(mobility_counts), "_", 2)[, 1], mobility = str_split_fixed(names(mobility_counts), "_", 2)[, 2], count = unlist(mobility_counts))
  
  res3 <- res %>%
    dplyr::group_by(taxa) %>%
    dplyr::summarise(total_count = nrow(taxo.df[taxo.df$taxa %in% colnames(df) & (taxo.df$class == taxa | taxo.df$order == taxa), ])) %>%
    ungroup() %>%
    mutate(taxa = case_when(
      taxa %in% c("Ostracoda", "Branchiopoda", "Malacostraca", "Hexanauplia") ~ "Crustacea", 
      taxa %in% c("Chilopoda", "Diplopoda", "Symphyla") ~ "Myriapoda",
      taxa %in% c("Hymenoptera", "Blattodea", "Mecoptera", "Neuroptera", "Orthoptera", "Psocoptera",
                  "Thysanoptera", "Megaloptera", "Odonata", "Plecoptera", "Ephemeroptera", "Trichoptera") ~ "Others",
      TRUE ~ taxa
    )) %>%
    dplyr::group_by(taxa) %>%
    dplyr::summarise(total_count = sum(total_count, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    dplyr::arrange(match(taxa, sort(final_mobility_df$taxa))) %>%
    dplyr::mutate(
      fraction = total_count / sum(total_count),
      ymax = cumsum(fraction),
      ymin = c(0, head(ymax, n = -1)), 
      label_y = (ymin + ymax) / 2,
      label_x = if_else(dplyr::row_number() %% 2 == 0, 2.3, 2)
    )
  
  final_mobility_df <- final_mobility_df %>%
    dplyr::filter(count > 0) %>%
    mutate(taxa = case_when(
      taxa %in% c("Ostracoda", "Branchiopoda", "Malacostraca", "Hexanauplia") ~ "Crustacea",
      taxa %in% c("Chilopoda", "Diplopoda", "Symphyla") ~ "Myriapoda",
      taxa %in% c("Hymenoptera", "Blattodea", "Mecoptera", "Neuroptera", "Orthoptera", "Psocoptera",
                  "Thysanoptera", "Megaloptera", "Odonata", "Plecoptera", "Ephemeroptera", "Trichoptera") ~ "Others",
      TRUE ~ taxa
    )) %>%
    arrange(match(taxa, c("Arachnida", "Clitellata", "Collembola", "Crustacea", 
                          "Eutardigrada", "Myriapoda", "Coleoptera", "Diptera", 
                          "Hemiptera", "Lepidoptera", "Others"))) %>%
    dplyr::group_by(taxa, mobility) %>%
    dplyr::summarise(count = sum(count, na.rm = TRUE)) %>%
    dplyr::ungroup() %>%
    left_join(res3, by = "taxa") %>%
    dplyr::group_by(taxa) %>%
    dplyr::mutate(
      mob_proportion = count / total_count) %>%
    dplyr::ungroup()
  
  final_mobility_df <- final_mobility_df[order(final_mobility_df$ymin), ]
  
  final_mobility_df$ymin_role <- 0
  final_mobility_df$ymax_role <- 0
  for (i in 1:nrow(final_mobility_df)) {
    if (final_mobility_df[i, ]$count != 0 ) {
      final_mobility_df[i, ]$ymax_role <- final_mobility_df[i, ]$ymin_role + ((final_mobility_df[i, ]$ymax - final_mobility_df[i, ]$ymin) * (final_mobility_df[i, ]$count / final_mobility_df[i, ]$total_count))
      final_mobility_df[i + 1, ]$ymin_role <- final_mobility_df[i, ]$ymax_role
    } else { 
      final_mobility_df[i, ]$ymax_role <- final_mobility_df[i, ]$ymin_role
      final_mobility_df[i + 1, ]$ymin_role <- final_mobility_df[i, ]$ymax_role
    }
  }
  
  res <- res %>%
    mutate(taxa = case_when(
      taxa %in% c("Ostracoda", "Branchiopoda", "Malacostraca", "Hexanauplia") ~ "Crustacea",
      taxa %in% c("Chilopoda", "Diplopoda", "Symphyla") ~ "Myriapoda",
      taxa %in% c("Hymenoptera", "Blattodea", "Mecoptera", "Neuroptera", "Orthoptera", "Psocoptera",
                  "Thysanoptera", "Megaloptera", "Odonata", "Plecoptera", "Ephemeroptera", "Trichoptera") ~ "Others",
      TRUE ~ taxa
    )) %>%
    group_by(taxa, site_no) %>%
    dplyr::summarise(richness = sum(richness)) %>%
    arrange(match(taxa, c("Arachnida", "Clitellata", "Collembola", "Crustacea", 
                          "Myriapoda", "Coleoptera", "Diptera", "Hemiptera", 
                          "Lepidoptera", "Others", "Psocoptera")))
  
  return(list(res, res2, final_mobility_df, res3))
}


# Import data
spmat.df <- read.csv("skyline_spmat.csv", check.names = FALSE)
taxo.df <- read.csv("skyline_taxo.csv")
shared.df <- read.csv("skyline_sharedmat.csv", check.names = FALSE)
meta <- c("site_no", "type", "sp_richness")

# Filter data for green roof sites and remove metadata columns
roof.df <- spmat.df %>%
  dplyr::filter(type == "Roof") %>%
  dplyr::select(-any_of(meta)) %>%
  select_if(colSums(.) != 0)
roof.df$site_no <-  spmat.df[spmat.df$type == "Roof", ]$site_no

# Same for ground-level data
ground.df <- spmat.df %>%
  dplyr::filter(type == "Ground") %>%
  dplyr::select(-any_of(meta)) %>%
  select_if(colSums(.) != 0)
ground.df$site_no <-  spmat.df[spmat.df$type == "Ground", ]$site_no

shared.df <- shared.df %>%
  dplyr::select(-any_of(c("jaccard", "turnover", "nestedness"))) 

summary.df <- data.frame(
  site_no = factor(),
  richness_ground = numeric(),
  richness_roof = numeric(),
  richness_shared = numeric()
)

# Compute species richness for ground, roof, and shared species at each site
for (site in unique(spmat.df$site_no)) { 
  site_ij <- spmat.df[spmat.df$site_no == site, ]
  
  # Extract and clean data for ground sites
  soil_ij <- site_ij[site_ij$type == "Ground", -which(colnames(site_ij) %in% meta)]
  soil_ij <- soil_ij[, -which(colSums(soil_ij) == 0)]
  
  # Extract and clean data for roof sites
  sky_ij <- site_ij[site_ij$type == "Roof",  -which(colnames(site_ij) %in% meta)]
  sky_ij <- sky_ij[, -which(colSums(sky_ij) == 0)]
  
  # Compute richness values
  richness_shared <- length(which(colnames(soil_ij) %in% colnames(sky_ij)))
  richness_ground <- length(colnames(soil_ij)[-which(colnames(soil_ij) %in% colnames(sky_ij))])
  richness_roof <- length(colnames(sky_ij)[-which(colnames(sky_ij) %in% colnames(soil_ij))])
  
  summary.df <- rbind(summary.df, data.frame(
    site_no = unique(site_ij$site_no),
    richness_ground = richness_ground,
    richness_roof = richness_roof,
    richness_shared = richness_shared))
}

mean_ground_unique <- round(mean(summary.df$richness_ground), 0)
mean_roof_unique <- round(mean(summary.df$richness_roof), 0)
mean_shared <- round(mean(summary.df$richness_shared), 0)
sd_ground_unique <- sd(summary.df$richness_ground)
sd_roof_unique <- sd(summary.df$richness_roof)
sd_shared <- sd(summary.df$richness_shared)

grid.newpage()

### To save in .svg format to open in Inkscape
# svg("VennDiagram.svg", width = 8, height = 8)

draw.pairwise.venn(area1 = mean_ground_unique + mean_shared ,
                   area2 = mean_roof_unique + mean_shared ,
                   cross.area = mean_shared + sd_shared, # + standard deviation
                   fill = c("sienna4", "forestgreen"),
                   alpha = 0.5, #le plus transparent
                   lty = "blank",
                   label.col = "transparent")

draw.pairwise.venn(area1 = mean_ground_unique + mean_shared ,
                   area2 = mean_roof_unique + mean_shared ,
                   cross.area = mean_shared - sd_shared,# - standard deviation
                   fill = c("sienna4", "forestgreen"),
                   alpha = 0.2,
                   lty = "blank",
                   label.col = "transparent")

draw.pairwise.venn(area1 = mean_ground_unique + mean_shared,
                   area2 = mean_roof_unique + mean_shared,
                   cross.area = mean_shared,
                   col = c("sienna4", "forestgreen"),
                   category = c("Ground-level", "Green Roof"),
                   cat.fontfamily = rep("sans", 2),
                   cat.cex = 1.25,
                   alpha = 0,
                   lwd = 2,
                   fontfamily = rep("sans", 3),
                   cex = 1.2 )

# dev.off()

# Compute the composition of major taxonomic groups and mobility modes at each site
# This is repeated for green roof and ground-level sites, as well as for species that are shared. 
roof.compo <- composition_agg_row(roof.df)[[1]]
roof.mobility <- composition_agg_row(roof.df)[[3]] 
roof.compo.tot <- composition_agg_row(roof.df)[[4]] 

ground.compo <- composition_agg_row(ground.df)[[1]]
ground.mobility <- composition_agg_row(ground.df)[[3]] 
ground.compo.tot <- composition_agg_row(ground.df)[[4]] 

shared.compo <- composition_agg_row(shared.df)[[1]]
shared.mobility <- composition_agg_row(shared.df)[[3]] 
shared.compo.tot <- composition_agg_row(shared.df)[[4]] 

# Create doughnut plots
plot1 <- ggplot() +
  geom_textpath(data = roof.compo.tot[roof.compo.tot$total_count > 0, ],
                aes(x = label_x, y = label_y, label = taxa), 
                angle = 90,
                size = 4, 
                color = "black", 
                spacing = 5, 
                # hjust = 1,  # Adjust horizontal alignment to center more
                vjust = 0.9,  # Adjust vertical alignment to center better
                text_only = TRUE) +
  geom_rect(data = roof.mobility,
            aes(xmin = 1.7, xmax = 1.9, ymin = ymin_role, ymax = ymax_role, fill = mobility),
  ) +
  scale_fill_manual(values =   c("#000000", "#2AB7CA", "#FE4A49", "#D4B483"), breaks = c("Aquatic", "flying", "edaphic", "epigeic"), labels = c("Aquatic", "Flying", "Epigeic", "Edaphic"), name = "Mobility") +
  new_scale_fill() + # This lets us add another external layer
  geom_rect(data = roof.compo.tot,
            aes(xmin = 0, xmax = 1.65, ymin = ymin, ymax = ymax, fill = taxa)) +
  geom_text(data = roof.compo.tot %>%
              dplyr::filter(total_count > 15), # Only larger slices
            aes(x = 1.3,
                y = (ymin + ymax) / 2, # Center the label in each slice
                label = total_count), # Display only the count, not taxa names
            size = 4,
            color = "black") +
  scale_fill_brewer(palette = "Set3") +
  guides(fill = "none") +
  theme_void() +
  coord_polar(theta = "y")

plot2 <- ggplot() +
  geom_textpath(data = ground.compo.tot[ground.compo.tot$total_count > 0, ],
                aes(x = label_x, y = label_y, label = taxa), 
                angle = 90,
                size = 4, 
                color = "black", 
                spacing = 5, 
                # hjust = 1,  # Adjust horizontal alignment to center more
                vjust =0.9,  # Adjust vertical alignment to center better
                text_only = TRUE) +
  geom_rect(data = ground.mobility,
            aes(xmin = 1.7, xmax = 1.9, ymin = ymin_role, ymax = ymax_role, fill = mobility),
  ) +
  scale_fill_manual(values =   c("#000000", "#2AB7CA", "#FE4A49", "#D4B483"), breaks = c("Aquatic", "flying", "edaphic", "epigeic"), labels = c("Aquatic", "Flying", "Epigeic", "Edaphic"), name = "Mobility") +
  new_scale_fill() +
  geom_rect(data = ground.mobility,
            aes(xmin = 0, xmax = 1.65, ymin = ymin, ymax = ymax, fill = taxa)) +
  geom_text(data = ground.compo.tot %>%
              dplyr::filter(total_count > 15), # Only larger slices
            aes(x = 1.3, 
                y = (ymin + ymax) / 2, # Center the label in each slice
                label = total_count), # Display only the count, not taxa names
            size = 4, 
            color = "black") +
  scale_fill_brewer(palette = "Set3") +
  guides(fill = "none") +
  theme_void() +
  coord_polar(theta = "y")

plot3 <- ggplot() +
  geom_textpath(data = shared.compo.tot[shared.compo.tot$total_count > 0, ],
                aes(x = label_x, y = label_y, label = taxa), 
                angle = 90,
                size = 4, 
                color = "black", 
                spacing = 5, 
                # hjust = 1,  # Adjust horizontal alignment to center more
                vjust = 0.9,  # Adjust vertical alignment to center better
                text_only = TRUE) +
  geom_rect(data = shared.mobility,
            aes(xmin = 1.7, xmax = 1.9, ymin = ymin_role, ymax = ymax_role, fill = mobility),
  ) +
  scale_fill_manual(values =   c("#000000", "#2AB7CA", "#FE4A49", "#D4B483"), breaks = c("Aquatic", "Winged", "Ground-dwelling", "Soil-living"), labels = c("Aquatic", "Flying", "Epigeic", "Edaphic")) +
  new_scale_fill() +
  geom_rect(data = shared.compo.tot,
            aes(xmin = 0, xmax = 1.65, ymin = ymin, ymax = ymax, fill = taxa)) +
  geom_text(data = shared.compo.tot %>%
              dplyr::filter(total_count > 15), # Only larger slices
            aes(x = 1.3, 
                y = (ymin + ymax) / 2, # Center the label in each slice
                label = total_count), # Display only the count, not Taxa names
            size = 4, 
            color = "black") +
  scale_fill_brewer(palette = "Set3") +
  guides(fill = "none") +
  theme_void() +
  coord_polar(theta = "y")


ggarrange(plot2, plot3, plot1, nrow = 1, 
          labels = c("Ground-level", "Shared", "Green roof"),
          common.legend = T)

# Fig S4
roof.compo$site_no <- as.factor(rep(seq(1, 52), 11))
ground.compo$site_no <- as.factor(rep(seq(1, 52), 11))
shared.compo$site_no <- as.factor(rep(seq(1, 52), 11))

plot1 <- ggplot(roof.compo, aes(x = site_no, y = richness, fill = taxa)) +
  geom_bar(position = "fill", stat = "identity", color = "black", lwd = 0.2) +
  labs(x = "", y = "Fraction of species") +
  scale_fill_brewer(palette = "Set3") +
  theme_pubclean(base_size = 16) +
  theme(axis.text.x = element_blank())

plot2 <- ggplot(ground.compo, aes(x = as.factor(site_no), y = richness, fill = taxa)) +
  geom_bar(position = "fill", stat = "identity", color = "black", lwd = 0.2) +
  labs(x = "", y = "Fraction of species") +
  scale_fill_brewer(palette = "Set3") +
  theme_pubclean(base_size = 16) +
  theme(axis.text.x = element_blank())

plot3 <- ggplot(shared.compo, aes(x = as.factor(site_no), y = richness, fill = taxa)) +
  geom_bar(position = "fill", stat = "identity", color = "black", lwd = 0.2) +
  labs(x = "Site", y = "Fraction of species") +
  scale_fill_brewer(palette = "Set3") +
  theme_pubclean(base_size = 16)

ggarrange(plot1, plot2, plot3, nrow = 3, common.legend = T,
          labels = c("Green roof", "Ground-level", "Shared species"),
          vjust = 0.5)
