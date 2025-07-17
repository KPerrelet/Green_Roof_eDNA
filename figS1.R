################################################################################
#################################### fig s1 ####################################
################################################################################
library(ggcorrplot)
library(dplyr)
library(reshape2)

# Function to generate asterisks based on p-values
labs.function = function(x){
  case_when(x >= 0.05 ~ "",
            x < 0.05 & x >= 0.01 ~ "*",
            x < 0.01 & x >= 0.001 ~ "**",
            x < 0.001 ~ "***")
}

# Import data
roof.covariates.df <- read.csv("skyline_covmat_greenroof.csv")
roof.covariates.df$pv <- ifelse(roof.covariates.df$pv == 0, FALSE, TRUE)
roof.covariates.df$substrate_type <- as.factor(roof.covariates.df$substrate_type)

# Define the list of covariates (columns) to include in the correlation matrix
covariates <- c("vegetation", "height", "depth", "age_roof", "area",
                "green_frac_50", "green_frac_500", "distance_gr")

# Scale the selected covariates to have zero mean and unit variance (standardization)
roof.covariates.df[, covariates] <- scale(roof.covariates.df[, covariates], center = T, scale = T) # scale to zero

# Compute the correlation matrix for the scaled covariates, rounding to 2 decimal places
corr_mat <- round(cor(scale(roof.covariates.df[, covariates]), use = "na.or.complete"), 2)
colnames(corr_mat) <- c("Vegetation", 
                        "Height",
                        "Depth",
                        "Roof age",
                        "Area",
                        "% green (50 m)", 
                        "% green (500 m)", 
                        "Distance to roof")
rownames(corr_mat) <- c("Vegetation", 
                        "Height",
                        "Depth",
                        "Roof age",
                        "Area",
                        "% green (50 m)", 
                        "% green (500 m)", 
                        "Distance to roof")

# Generate a p-value matrix from the correlation matrix (using ggcorrplot)
p.df <- as.data.frame(ggcorrplot::cor_pmat(corr_mat))

# Apply the function to create a matrix of asterisks based on p-values
p.labs = p.df  %>%                      
  mutate_all(labs.function)

p.labs$Var1 = as.factor(rownames(p.labs))
p.labs = reshape2::melt(p.labs, id.vars = "Var1", variable.name = "Var2", value.name = "lab")

# Initial ggcorrplot
cor.plot <- ggcorrplot(corr_mat, method = "circle", type = "upper", insig = "blank", show.diag = T, lab = F, lab_size = 5) +
  labs(x = "", y = "") +
  theme_void(base_size = 16) + 
  theme(strip.background = element_blank(),
        panel.spacing = unit(1, "lines"), 
        legend.text = element_text(size = 16), 
        axis.text.y = element_text(size = 16, hjust = 1),
        axis.title = element_text(size = 18)) +
  scale_fill_gradient2(low = "#b2182b", 
                       high = "#2166ac", 
                       name = "Correlation") 

# Subsetting asteriks matrix to only those rows within ggcorrplot data
p.labs$in.df = ifelse(is.na(match(paste0(p.labs$Var1, p.labs$Var2), 
                                  paste0(cor.plot[["data"]]$Var1, cor.plot[["data"]]$Var2))),
                      "No", "Yes")

p.labs = select(filter(p.labs, in.df == "Yes"), -in.df)

cor.plot + 
  geom_text(aes(x = p.labs$Var1, 
                y = p.labs$Var2), 
            label = p.labs$lab, 
            nudge_y = -0.08,
            size = 6)
