################################################################################
#################################### fig s5 ####################################
################################################################################
library(iNEXT)
library(ggplot2)
library(ggpubr)

# Function to generate lists of taxa by different categories (e.g., class, order, mobility)
generate_taxa_list <- function(df, taxo.df) {
  taxa.list <- list()
  
  # Add overall 'Total' taxa category
  taxa <- list(t(as.data.frame(df[, which(colnames(df) %in% taxo.df$taxa)]))) 
  names(taxa) <- "Total"
  taxa.list <- append(taxa.list, taxa)
  
  # Add taxa for each class and order (if applicable)
  for (clas in unique(taxo.df$class)) {
    if (clas == "Insecta") {
      for (ordr in unique(taxo.df[taxo.df$class == clas, ]$order)) {
        taxa <- df[, which(colnames(df) %in% taxo.df[taxo.df$order == ordr, ]$taxa)]
        taxa <- list(t(as.data.frame(taxa))) 
        names(taxa) <- ordr
        taxa.list <- append(taxa.list, taxa)
      }
    } else {
      taxa <- df[, which(colnames(df) %in% taxo.df[taxo.df$class == clas, ]$taxa)]
      taxa <- list(t(as.data.frame(taxa))) 
      names(taxa) <- clas
      taxa.list <- append(taxa.list, taxa)
    }
  } 
  
  # Add taxa by mobility (Epigeic, Edaphic, Flying)
  mobility.taxa.list <- list()
  for (mobi in unique(taxo.df$mobility)) {
    if (mobi == "Aquatic") {
      next # Skip "Aquatic" taxa
    }
    mobility.taxa <- df[, which(colnames(df) %in% taxo.df[taxo.df$mobility == mobi, ]$taxa)]
    mobility.taxa <- list(t(as.data.frame(mobility.taxa))) 
    names(mobility.taxa) <- mobi
    mobility.taxa.list <- append(mobility.taxa.list, mobility.taxa)
    
  }
  return(list(taxa.list, mobility.taxa.list))
}

# Import data
spmat.df <- read.csv("skyline_spmat.csv", check.names = F)
meta <- c("site_no", "type", "sp_richness")
taxo.df <- read.csv("skyline_taxo.csv")

# Split the data into Roof and Ground types by filtering based on the 'type' column
roof.df <- spmat.df[spmat.df$type == "Roof", -which(colnames(spmat.df) %in% meta)]
ground.df <- spmat.df[spmat.df$type == "Ground", -which(colnames(spmat.df) %in% meta)]

# Combine certain taxa into broader groups to simplify analysis
taxo.df[taxo.df$class %in% c("Symphyla", "Diplopoda", "Chilopoda"), ]$class <-"Myriapoda"
taxo.df[taxo.df$class %in% c("Ostracoda", "Branchiopoda", "Malacostraca", "Hexanauplia"), ]$class <- "Crustacea"
taxo.df[taxo.df$order %in% c("Hymenoptera", "Blattodea", "Mecoptera", "Neuroptera", "Orthoptera", "Thysanoptera", 
                             "Megaloptera", "Odonata", "Plecoptera", "Ephemeroptera", "Trichoptera", "Psocoptera"), ]$order <- "Others"

# Generate taxa lists for Roof and Ground types
roof.taxa.list <- generate_taxa_list(spmat.df[spmat.df$type == "Roof", ], taxo.df)[[1]]

# Perform rarefaction analysis for Roof taxa using iNEXT
roof.taxa.inext <- iNEXT(roof.taxa.list, q = 0, datatype = "incidence_raw")
roof.taxa.df <- fortify(roof.taxa.inext, type = 1) 
roof.taxa.df$type <- "Roof"

# Repeat for ground-level sites, excluding the very rare "Eutardigrada" taxa for avoid computational errors (as they are so rare)
ground.taxa.list <- generate_taxa_list(spmat.df[spmat.df$type == "Ground", ], taxo.df)[[1]]
ground.taxa.list <- ground.taxa.list[names(ground.taxa.list) %in% "Eutardigrada" == F] 
ground.taxa.inext <- iNEXT(ground.taxa.list, q = 0, datatype = "incidence_raw")
ground.taxa.df <- fortify(ground.taxa.inext, type = 1) 
ground.taxa.df$type <- "Ground"

# Combine both Roof and Ground taxa data frames
all.taxa.df <- rbind(roof.taxa.df, ground.taxa.df)

# Filter to only include observed values for plotting
all.taxa.point <- all.taxa.df[which(all.taxa.df$Method=="Observed"),]

# Repeat the same process for mobility categories (Epigeic, Edaphic, Flying)
roof.mobi.list <- generate_taxa_list(spmat.df[spmat.df$type == "Roof", ], taxo.df)[[2]]
roof.mobi.inext <- iNEXT(roof.mobi.list, q = 0, datatype = "incidence_raw")
roof.mobi.df <- fortify(roof.mobi.inext, type = 1) 
roof.mobi.df$type <- "Roof"

ground.mobi.list <- generate_taxa_list(spmat.df[spmat.df$type == "Ground", ], taxo.df)[[2]]
ground.mobi.inext <- iNEXT(ground.mobi.list, q = 0, datatype = "incidence_raw")
ground.mobi.df <- fortify(ground.mobi.inext, type = 1) 
ground.mobi.df$type <- "Ground"

all.mobi.df <- rbind(roof.mobi.df, ground.mobi.df)
all.mobi.point <- all.mobi.df[which(all.mobi.df$Method=="Observed"),]

plot1 <- ggplot() +
  coord_cartesian(ylim = c(0, max(log(all.taxa.point$y) + 0.1))) + 
  geom_point(data = all.taxa.point, aes(x = x, y = log(y), shape = type, color = Assemblage), size = 4) +
  geom_line(data = roof.taxa.df[roof.taxa.df$Method != "Observed", ], aes(x = x, y = log(y), linetype = Method, color = Assemblage), lwd = 1) +
  geom_line(data = ground.taxa.df[ground.taxa.df$Method != "Observed", ], aes(x = x, y = log(y), linetype = Method, color = Assemblage), lwd = 1) +
  geom_point(data = roof.taxa.df[roof.taxa.df$Method != "Observed", ] %>% dplyr::group_by(Assemblage) %>% dplyr::slice_tail(n = 1),
             aes(x = x + 3, y = log(y), color = Assemblage), shape = 2, size = 4) +
  geom_point(data = ground.taxa.df[ground.taxa.df$Method != "Observed", ] %>% dplyr::group_by(Assemblage) %>% dplyr::slice_tail(n = 1),
             aes(x = x + 3, y = log(y), color = Assemblage), shape = 1, size = 4) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Number of sites", y = "log(Species richness)") +
  theme_pubr(base_size = 18) +
  scale_shape_manual(labels = c("Green roof", "Ground-level"), values = c(17, 16)) +
  scale_color_manual(name = "Taxonomy", values = {  
    colors <- RColorBrewer::brewer.pal(n = length(unique(all.taxa.point$Assemblage)), "Set3")
    colors[length(colors)] <- "black"
    colors
  }) + 
  theme(legend.position="right", legend.box="vertical", legend.text = element_text(size = 18)) 

plot2 <- ggplot() +
  geom_point(data = all.mobi.point, aes(x = x, y = log(y), shape = type, color = Assemblage), size = 4) +
  geom_line(data = roof.mobi.df[roof.mobi.df$Method != "Observed", ], aes(x = x, y = log(y), linetype = Method, color = Assemblage), lwd = 1) +
  geom_line(data = ground.mobi.df[ground.mobi.df$Method != "Observed", ], aes(x = x, y = log(y), linetype = Method, color = Assemblage), lwd = 1) +
  geom_point(data = roof.mobi.df[roof.mobi.df$Method != "Observed", ] %>% dplyr::group_by(Assemblage) %>% dplyr::slice_tail(n = 1),
             aes(x = x + 3, y = log(y), color = Assemblage), shape = 2, size = 4) +
  geom_point(data = ground.mobi.df[ground.mobi.df$Method != "Observed", ] %>% dplyr::group_by(Assemblage) %>% dplyr::slice_tail(n = 1),
             aes(x = x + 3, y = log(y), color = Assemblage), shape = 1, size = 4) +
  scale_linetype_manual(values = c("dashed", "solid")) +
  labs(x = "Number of sites", y = "log(Species richness)") +
  theme_pubr(base_size = 18) +
  scale_shape_manual(labels = c("Green roof", "Ground-level"), values = c(17, 16)) +
  scale_color_manual(name = "Mobility", labels = c("Epigeic", "Edaphic", "Flying"), values = c("#FE4A49", "#D4B483", "#2AB7CA")) +
  theme(legend.position="right", legend.box="vertical", legend.text = element_text(size = 18)) +
  guides(linetype = "none", 
         shape = "none")

ggarrange(plot1, plot2, common.legend = F, legend = "bottom", align = "hv")
