#####Ambient Lake Data Mixed Effect Models incorporating land type#####
#R version 4.2.2 
## 23 September 2026

#Updated to scale after LCBD 

###Clear working environment
rm(list=ls())


##Reading in Data 
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) #set working directory
getwd()
WQ.zoo <- read.csv("WQ.taxa25.csv", header = TRUE)

#install.packages("remotes")
#library(remotes)
#install_version("ggeffects", "1.4.0")


library(tidyverse) #everyday data Wrangling
library(lme4) #used to run lmer() function, perform 
library(GGally)#used for ggpairs correlation plots
library(car)#used for QQ plot generation
library(ggeffects)#used for ggpredict
library(brms)#baysian mixed models
library(bayestestR)#used for summary stats
library(vegan)
library(hillR)#calculating hill numbers
library(adespatial)


#Isolating WQ from taxa
WQ_unclean <- WQ.zoo[, c(1,51:61,63:75)]

#Prepping the data 
WQ.DoY <- WQ_unclean %>%
  separate_wider_delim(lake.date, delim = "_", names=c("lakeID", "year", "DoY"),cols_remove = FALSE)

WQ.DoY$DoY <- as.numeric(WQ.DoY$DoY)
WQ.DoY$year <- as.numeric(WQ.DoY$year)

summary(WQ.zoo)

#Calculating z scores of each variable 
WQ_zscores <- WQ.DoY %>%  
  mutate(pH_z = (pH - mean(pH, na.rm = TRUE))/sd(pH, na.rm = TRUE))%>%  
  mutate(DO_z = (DO_mg.l - mean(DO_mg.l, na.rm = TRUE))/sd(DO_mg.l, na.rm = TRUE))%>%  
  mutate(Secchi_z = (Secchi_m - mean(Secchi_m, na.rm = TRUE))/sd(Secchi_m, na.rm = TRUE))%>%  
  mutate(BottomDepth_z = (Bottom.Depth_m - mean(Bottom.Depth_m, na.rm = TRUE))/sd(Bottom.Depth_m, na.rm = TRUE))%>%  
  mutate(WaterTemp_z = (Temp.water.degC - mean(Temp.water.degC, na.rm = TRUE))/sd(Temp.water.degC, na.rm = TRUE))%>%  
  mutate(FSS_z = (Fixed.Suspended.Solids_mg.l - mean(Fixed.Suspended.Solids_mg.l, na.rm = TRUE))/sd(Fixed.Suspended.Solids_mg.l, na.rm = TRUE))%>%  
  mutate(Turb_z = (Turbidity_NTU - mean(Turbidity_NTU, na.rm = TRUE))/sd(Turbidity_NTU, na.rm = TRUE))%>%  
  mutate(Conduct_z = (Specific_Conductance_uS.cm - mean(Specific_Conductance_uS.cm, na.rm = TRUE))/sd(Specific_Conductance_uS.cm, na.rm = TRUE))%>%  
  mutate(TSS_z = (Total_Suspended_Solids_mg.l - mean(Total_Suspended_Solids_mg.l, na.rm = TRUE))/sd(Total_Suspended_Solids_mg.l, na.rm = TRUE))%>%  
  mutate(Phosphorus_z = (Total_Phorphorus_Mixed_Forms_mg.l.as.P. - mean(Total_Phorphorus_Mixed_Forms_mg.l.as.P., na.rm = TRUE))/sd(Total_Phorphorus_Mixed_Forms_mg.l.as.P., na.rm = TRUE))%>%  
  mutate(Orthophosphate_z = (Orthophosphate_.as_P._mg.l - mean(Orthophosphate_.as_P._mg.l, na.rm = TRUE))/sd(Orthophosphate_.as_P._mg.l, na.rm = TRUE))%>%  
  mutate(NitrateNitrite_z = (Nitrate.Nitrite_.as_N._mg.l - mean(Nitrate.Nitrite_.as_N._mg.l, na.rm = TRUE))/sd(Nitrate.Nitrite_.as_N._mg.l, na.rm = TRUE))%>%  
  mutate(DOSaturation_z = (Dissolved_Oxygen_Saturation - mean(Dissolved_Oxygen_Saturation, na.rm = TRUE))/sd(Dissolved_Oxygen_Saturation, na.rm = TRUE))%>%  
  mutate(Chlorophylla_z = (Chlorophyll_a - mean(Chlorophyll_a, na.rm = TRUE))/sd(Chlorophyll_a, na.rm = TRUE))%>%  
  mutate(Carbon_z = (Carbon - mean(Carbon, na.rm = TRUE))/sd(Carbon, na.rm = TRUE))%>%  
  mutate(year_z = (year - mean(year, na.rm = TRUE))/sd(year, na.rm = TRUE))

WQ_zscores$pH_z2 <- scale(WQ_zscores$pH, center = TRUE, scale = TRUE)
summary(WQ_zscores)


#Histogram of DOY
hist(WQ_zscores$DoY, breaks = 60)


WQ <- WQ_zscores %>%
  mutate(Season = case_when(
    between(DoY, 134, 183) ~ "early",
    between(DoY, 187, 236) ~ "middle",
    between(DoY, 237, 278) ~"late"))

summary(WQ)


#####Creating N:P ratios based on Grace's original Code#####

# Summing columns Kjeldahl_nitrogen_asN_mg.l and NitrateNitrite_asN_mg.l, replacing NA values in Kjeldahl_nitrogen_asN_mg.l with values from Nutrient.nitrogen_mg.l
WQ$Total_Nitrogen <- ifelse(is.na(WQ$Kjeldahl_nitrogen.as_N.), 
                                   WQ$Nitrate.Nitrite_.as_N._mg.l + WQ$Nutrient.nitrogen_mg.l, WQ$Kjeldahl_nitrogen.as_N. + WQ$Nitrate.Nitrite_.as_N._mg.l)

#Create a column with all phosphorous. If Total phosphorus isnt available.

# Create a new column incorporating all values from Total.Phosphorus..mixed.forms..as.P._mg.l
WQ$Total_Phosphorus <- WQ$Total_Phorphorus_Mixed_Forms_mg.l.as.P.

WQ$N.P_Ratio <- WQ$Total_Nitrogen / WQ$Total_Phosphorus
print(WQ$N.P_Ratio)


#####Diversity Metrics#####
taxa <- WQ.zoo[, c(1:49)]


#Prepping the data 
taxa.DoY <- taxa %>%
  separate_wider_delim(lake.date, delim = "_", names=c("lakeID", "year", "DoY"),cols_remove = FALSE)

table(taxa.DoY$year)




#Replacing NAs with zeros
taxa.fix <- taxa %>%
  replace(is.na(.),0)

#Hill number calculation 
hill.rich <- hill_taxa(taxa.fix[,-1], q = 0, MARGIN = 1, base = exp(1))
hill.Shannon <- hill_taxa(taxa.fix[,-1], q = 1, MARGIN = 1, base = exp(1))
hill.InverseSimpson <- hill_taxa(taxa.fix[,-1], q = 2, MARGIN = 1, base = exp(1))

#Hill number beta calculations
  #Removing rows that don't contain any taxa
taxa.fix.beta <- taxa.fix[rowSums(taxa.fix[,-1]) !=0, ] %>%
  separate_wider_delim(lake.date, delim = "_", names=c("lakeID", "year", "DoY"),cols_remove = FALSE) %>%
  filter(!lakeID == "213")%>% #These lakes were removed because they were only sampled once
  filter(!lakeID == "239")%>%      #and cannot calculate beta by lake with only one sample
  filter(!lakeID == "243")%>%
  filter(!lakeID == "250")%>%
  filter(!lakeID == "263")%>%
  filter(!lakeID == "267")

#####Alpha Diversity
#adding LD
lake.date <- taxa$lake.date #taking lake date from taxa df 
LD.df <- as.data.frame(lake.date) #making lake date a df
LD.df$hill.rich <- hill.rich #adding rich to the df
LD.df$hill.Shannon <- hill.Shannon #adding rich to the df
LD.df$hill.InverseSimpson <- hill.InverseSimpson #adding rich to the df




#How to preserve the lake data as the row names when running through the HillR code  
taxa.fix.beta <- as.data.frame(taxa.fix.beta)
rownames(taxa.fix.beta) <- taxa.fix.beta$lake.date

#creating list of unique lakes
lakes <-unique(taxa.fix.beta$lakeID)


#####Beta Function#####
#####By lakes#####
betafunction.q0 <- function(lake_id) {
  # Subset data for the given lake ID
  subset_data <- taxa.fix.beta[taxa.fix.beta$lakeID == lake_id, -(1:4)]
  
  # Remove columns that contain only zeros
  subset_data <- subset_data[, colSums(subset_data != 0) > 0]
  
  # Call the function with the subsetted data and run pairwise
  pairwise_result <- hill_taxa_parti_pairwise(subset_data, 
                                              q = 0, 
                                              rel_then_pool = TRUE, 
                                              pairs = c("unique"),
                                              .progress = TRUE)
  
  # Convert pairwise result to a symmetric similarity matrix
  wide_pair <- pairwise_result %>%
    select(site1, site2, local_similarity) %>%
    pivot_wider(
      names_from = site2,
      values_from = local_similarity,
      values_fill = NA  # Replace NAs with 0
    )
  
  # Convert to numeric symmetric matrix
  wide_matrix <- as.matrix(wide_pair[, -1])  # Exclude row names column
  rownames(wide_matrix) <- wide_pair[[1]]    # Set row names
  
  # Add site names as rows and columns (ensure all sites are included)
  all_sites <- unique(c(pairwise_result$site1, pairwise_result$site2))
  
  wide_pair <- wide_pair %>%
    complete(site1 = all_sites) %>%  # Add any missing rows
    arrange(site1)
  
  # Convert to matrix
  wide_matrix <- as.matrix(wide_pair[, -1])
  rownames(wide_matrix) <- wide_pair$site1  # Set row names from site1
  
  # Add any missing columns to ensure symmetry
  missing_cols <- setdiff(all_sites, colnames(wide_matrix))
  wide_matrix <- cbind(wide_matrix, matrix(0, nrow = nrow(wide_matrix), ncol = length(missing_cols)))
  colnames(wide_matrix)[(ncol(wide_matrix) - length(missing_cols) + 1):ncol(wide_matrix)] <- missing_cols
  
  # Reorder rows and columns to match
  wide_matrix <- wide_matrix[all_sites, all_sites]
  
  # Replace NAs with 0 and ensure symmetry
  wide_matrix[is.na(wide_matrix)] <- 0
  wide_matrix <- (wide_matrix + t(wide_matrix)) / 2  # Ensure matrix is symmetric
  
  # Add the diagonal trace of 0
  diag(wide_matrix) <- 0
  
  # Run LCBD calculation
  lcbd_result <- LCBD.comp(wide_matrix, sqrt.D = FALSE)
  
  # Combine LCBD results with site identifiers
  lcbd_output <- data.frame(
    site = rownames(wide_matrix),
    LCBD = lcbd_result$LCBD
  )
  
  # Scaling lcbd
  lcbd_output$LCBD <- (lcbd_output$LCBD - min(lcbd_output$LCBD)) / 
    (max(lcbd_output$LCBD) - min(lcbd_output$LCBD))
  
  # Return results as a list
  return(list(
    pairwise_result = pairwise_result,
    lcbd_result = lcbd_output
  ))
}

# Run the function for multiple lakes
lakes.q0 <- lapply(lakes, betafunction.q0)

#####lakes Q=1
betafunction.q1 <- function(lake_id) {
  # Subset data for the given lake ID
  subset_data <- taxa.fix.beta[taxa.fix.beta$lakeID == lake_id, -(1:4)]
  
  # Remove columns that contain only zeros
  subset_data <- subset_data[, colSums(subset_data != 0) > 0]
  
  # Call the function with the subsetted data and run pairwise
  pairwise_result <- hill_taxa_parti_pairwise(subset_data, 
                                              q = 1, 
                                              rel_then_pool = TRUE, 
                                              pairs = c("unique"),
                                              .progress = TRUE)
  
  # Convert pairwise result to a symmetric similarity matrix
  wide_pair <- pairwise_result %>%
    select(site1, site2, local_similarity) %>%
    pivot_wider(
      names_from = site2,
      values_from = local_similarity,
      values_fill = NA  # Replace NAs with 0
    )
  
  # Convert to numeric symmetric matrix
  wide_matrix <- as.matrix(wide_pair[, -1])  # Exclude row names column
  rownames(wide_matrix) <- wide_pair[[1]]    # Set row names
  
  # Add site names as rows and columns (ensure all sites are included)
  all_sites <- unique(c(pairwise_result$site1, pairwise_result$site2))
  
  wide_pair <- wide_pair %>%
    complete(site1 = all_sites) %>%  # Add any missing rows
    arrange(site1)
  
  # Convert to matrix
  wide_matrix <- as.matrix(wide_pair[, -1])
  rownames(wide_matrix) <- wide_pair$site1  # Set row names from site1
  
  # Add any missing columns to ensure symmetry
  missing_cols <- setdiff(all_sites, colnames(wide_matrix))
  wide_matrix <- cbind(wide_matrix, matrix(0, nrow = nrow(wide_matrix), ncol = length(missing_cols)))
  colnames(wide_matrix)[(ncol(wide_matrix) - length(missing_cols) + 1):ncol(wide_matrix)] <- missing_cols
  
  # Reorder rows and columns to match
  wide_matrix <- wide_matrix[all_sites, all_sites]
  
  # Replace NAs with 0 and ensure symmetry
  wide_matrix[is.na(wide_matrix)] <- 0
  wide_matrix <- (wide_matrix + t(wide_matrix)) / 2  # Ensure matrix is symmetric
  
  # Add the diagonal trace of 0
  diag(wide_matrix) <- 0
  
  # Run LCBD calculation
  lcbd_result <- LCBD.comp(wide_matrix, sqrt.D = FALSE)
  
  # Combine LCBD results with site identifiers
  lcbd_output <- data.frame(
    site = rownames(wide_matrix),
    LCBD = lcbd_result$LCBD
  )
  
  # Scaling lcbd
  lcbd_output$LCBD <- (lcbd_output$LCBD - min(lcbd_output$LCBD)) / 
    (max(lcbd_output$LCBD) - min(lcbd_output$LCBD))
  
  # Return results as a list
  return(list(
    pairwise_result = pairwise_result,
    lcbd_result = lcbd_output
  ))
}

# Run the function for multiple lakes
lakes.q1 <- lapply(lakes, betafunction.q1)

#####lakes Q=2
betafunction.q2 <- function(lake_id) {
  # Subset data for the given lake ID
  subset_data <- taxa.fix.beta[taxa.fix.beta$lakeID == lake_id, -(1:4)]
  
  # Remove columns that contain only zeros
  subset_data <- subset_data[, colSums(subset_data != 0) > 0]
  
  # Call the function with the subsetted data and run pairwise
  pairwise_result <- hill_taxa_parti_pairwise(subset_data, 
                                              q = 2, 
                                              rel_then_pool = TRUE, 
                                              pairs = c("unique"),
                                              .progress = TRUE)
  
  # Convert pairwise result to a symmetric similarity matrix
  wide_pair <- pairwise_result %>%
    select(site1, site2, local_similarity) %>%
    pivot_wider(
      names_from = site2,
      values_from = local_similarity,
      values_fill = NA  # Replace NAs with 0
    )
  
  # Convert to numeric symmetric matrix
  wide_matrix <- as.matrix(wide_pair[, -1])  # Exclude row names column
  rownames(wide_matrix) <- wide_pair[[1]]    # Set row names
  
  # Add site names as rows and columns (ensure all sites are included)
  all_sites <- unique(c(pairwise_result$site1, pairwise_result$site2))
  
  wide_pair <- wide_pair %>%
    complete(site1 = all_sites) %>%  # Add any missing rows
    arrange(site1)
  
  # Convert to matrix
  wide_matrix <- as.matrix(wide_pair[, -1])
  rownames(wide_matrix) <- wide_pair$site1  # Set row names from site1
  
  # Add any missing columns to ensure symmetry
  missing_cols <- setdiff(all_sites, colnames(wide_matrix))
  wide_matrix <- cbind(wide_matrix, matrix(0, nrow = nrow(wide_matrix), ncol = length(missing_cols)))
  colnames(wide_matrix)[(ncol(wide_matrix) - length(missing_cols) + 1):ncol(wide_matrix)] <- missing_cols
  
  # Reorder rows and columns to match
  wide_matrix <- wide_matrix[all_sites, all_sites]
  
  # Replace NAs with 0 and ensure symmetry
  wide_matrix[is.na(wide_matrix)] <- 0
  wide_matrix <- (wide_matrix + t(wide_matrix)) / 2  # Ensure matrix is symmetric
  
  # Add the diagonal trace of 0
  diag(wide_matrix) <- 0
  
  # Run LCBD calculation
  lcbd_result <- LCBD.comp(wide_matrix, sqrt.D = FALSE)
  
  # Combine LCBD results with site identifiers
  lcbd_output <- data.frame(
    site = rownames(wide_matrix),
    LCBD = lcbd_result$LCBD
  )
  
  # Scaling lcbd 
  lcbd_output$LCBD <- (lcbd_output$LCBD - min(lcbd_output$LCBD)) / 
    (max(lcbd_output$LCBD) - min(lcbd_output$LCBD))
  
  # Return results as a list
  return(list(
    pairwise_result = pairwise_result,
    lcbd_result = lcbd_output
  ))
}

# Run the function for multiple lakes
lakes.q2 <- lapply(lakes, betafunction.q2)


saveRDS(lakes.q0, file = "lakes.q0.rds")
saveRDS(lakes.q1, file = "lakes.q1.rds")
saveRDS(lakes.q2, file = "lakes.q2.rds")

lakes.q0 <- readRDS(file = "lakes.q0.rds")
lakes.q1 <- readRDS(file = "lakes.q1.rds")
lakes.q2 <- readRDS(file = "lakes.q2.rds")


# Save the list object
save(lakes.q0, file = "lakes.q0.RData")

# Load the list object later
load("lakes.q0.RData")

#####By year#####
#creating list of unique lakes
year <-unique(taxa.fix.beta$year)

#####betafunction.q0
Year.betafunction.q0 <- function(year_id) {
  # Subset data for the given lake ID
  subset_data <- taxa.fix.beta[taxa.fix.beta$year == year_id, -(1:4)]
  
  # Remove columns that contain only zeros
  subset_data <- subset_data[, colSums(subset_data != 0) > 0]
  
  # Call the function with the subsetted data and run pairwise
  pairwise_result <- hill_taxa_parti_pairwise(subset_data, 
                                              q = 0, 
                                              rel_then_pool = TRUE, 
                                              pairs = c("unique"),
                                              .progress = TRUE)
  
  # Convert pairwise result to a symmetric similarity matrix
  wide_pair <- pairwise_result %>%
    select(site1, site2, local_similarity) %>%
    pivot_wider(
      names_from = site2,
      values_from = local_similarity,
      values_fill = NA  # Replace NAs with 0
    )
  
  # Convert to numeric symmetric matrix
  wide_matrix <- as.matrix(wide_pair[, -1])  # Exclude row names column
  rownames(wide_matrix) <- wide_pair[[1]]    # Set row names
  
  # Add site names as rows and columns (ensure all sites are included)
  all_sites <- unique(c(pairwise_result$site1, pairwise_result$site2))
  
  wide_pair <- wide_pair %>%
    complete(site1 = all_sites) %>%  # Add any missing rows
    arrange(site1)
  
  # Convert to matrix
  wide_matrix <- as.matrix(wide_pair[, -1])
  rownames(wide_matrix) <- wide_pair$site1  # Set row names from site1
  
  # Add any missing columns to ensure symmetry
  missing_cols <- setdiff(all_sites, colnames(wide_matrix))
  wide_matrix <- cbind(wide_matrix, matrix(0, nrow = nrow(wide_matrix), ncol = length(missing_cols)))
  colnames(wide_matrix)[(ncol(wide_matrix) - length(missing_cols) + 1):ncol(wide_matrix)] <- missing_cols
  
  # Reorder rows and columns to match
  wide_matrix <- wide_matrix[all_sites, all_sites]
  
  # Replace NAs with 0 and ensure symmetry
  wide_matrix[is.na(wide_matrix)] <- 0
  wide_matrix <- (wide_matrix + t(wide_matrix)) / 2  # Ensure matrix is symmetric
  
  # Add the diagonal trace of 0
  diag(wide_matrix) <- 0
  
  # Run LCBD calculation
  lcbd_result <- LCBD.comp(wide_matrix, sqrt.D = FALSE)
  
  # Combine LCBD results with site identifiers
  lcbd_output <- data.frame(
    site = rownames(wide_matrix),
    LCBD = lcbd_result$LCBD
  )
  
  # Scaling lcbd 
  lcbd_output$LCBD <- (lcbd_output$LCBD - min(lcbd_output$LCBD)) / 
    (max(lcbd_output$LCBD) - min(lcbd_output$LCBD))
  
  # Return results as a list
  return(list(
    pairwise_result = pairwise_result,
    lcbd_result = lcbd_output
  ))
}

# Run the function for multiple years
year.q0 <- lapply(year, Year.betafunction.q0)


#####betafunction.q1
Year.betafunction.q1 <- function(year_id) {
  # Subset data for the given lake ID
  subset_data <- taxa.fix.beta[taxa.fix.beta$year == year_id, -(1:4)]
  
  # Remove columns that contain only zeros
  subset_data <- subset_data[, colSums(subset_data != 0) > 0]
  
  # Call the function with the subsetted data and run pairwise
  pairwise_result <- hill_taxa_parti_pairwise(subset_data, 
                                              q = 1, 
                                              rel_then_pool = TRUE, 
                                              pairs = c("unique"),
                                              .progress = TRUE)
  
  print(pairwise_result)
  
  # Convert pairwise result to a symmetric similarity matrix
  wide_pair <- pairwise_result %>%
    select(site1, site2, local_similarity) %>%
    pivot_wider(
      names_from = site2,
      values_from = local_similarity,
      values_fill = NA  # Replace NAs with 0
    )
  
  # Convert to numeric symmetric matrix
  wide_matrix <- as.matrix(wide_pair[, -1])  # Exclude row names column
  rownames(wide_matrix) <- wide_pair[[1]]    # Set row names
  
  # Add site names as rows and columns (ensure all sites are included)
  all_sites <- unique(c(pairwise_result$site1, pairwise_result$site2))
  
  wide_pair <- wide_pair %>%
    complete(site1 = all_sites) %>%  # Add any missing rows
    arrange(site1)
  
  # Convert to matrix
  wide_matrix <- as.matrix(wide_pair[, -1])
  rownames(wide_matrix) <- wide_pair$site1  # Set row names from site1
  
  # Add any missing columns to ensure symmetry
  missing_cols <- setdiff(all_sites, colnames(wide_matrix))
  wide_matrix <- cbind(wide_matrix, matrix(0, nrow = nrow(wide_matrix), ncol = length(missing_cols)))
  colnames(wide_matrix)[(ncol(wide_matrix) - length(missing_cols) + 1):ncol(wide_matrix)] <- missing_cols
  
  # Reorder rows and columns to match
  wide_matrix <- wide_matrix[all_sites, all_sites]
  
  # Replace NAs with 0 and ensure symmetry
  wide_matrix[is.na(wide_matrix)] <- 0
  wide_matrix <- (wide_matrix + t(wide_matrix)) / 2  # Ensure matrix is symmetric
  
  # Add the diagonal trace of 0
  diag(wide_matrix) <- 0
  
  # Run LCBD calculation
  lcbd_result <- LCBD.comp(wide_matrix, sqrt.D = FALSE)
  
  # Combine LCBD results with site identifiers
  lcbd_output <- data.frame(
    site = rownames(wide_matrix),
    LCBD = lcbd_result$LCBD
  )
  
  # Scaling lcbd 
  lcbd_output$LCBD <- (lcbd_output$LCBD - min(lcbd_output$LCBD)) / 
    (max(lcbd_output$LCBD) - min(lcbd_output$LCBD))
  
  # Return results as a list
  return(list(
    pairwise_result = pairwise_result,
    lcbd_result = lcbd_output
  ))
}

# Run the function for multiple years
year.q1 <- lapply(year, Year.betafunction.q1)

#####betafunction.q2
Year.betafunction.q2 <- function(year_id) {
  # Subset data for the given lake ID
  subset_data <- taxa.fix.beta[taxa.fix.beta$year == year_id, -(1:4)]
  
  # Remove columns that contain only zeros
  subset_data <- subset_data[, colSums(subset_data != 0) > 0]
  
  # Call the function with the subsetted data and run pairwise
  pairwise_result <- hill_taxa_parti_pairwise(subset_data, 
                                              q = 2, 
                                              rel_then_pool = TRUE, 
                                              pairs = c("unique"),
                                              .progress = TRUE)
  
  # Convert pairwise result to a symmetric similarity matrix
  wide_pair <- pairwise_result %>%
    select(site1, site2, local_similarity) %>%
    pivot_wider(
      names_from = site2,
      values_from = local_similarity,
      values_fill = NA  # Replace NAs with 0
    )
  
  # Convert to numeric symmetric matrix
  wide_matrix <- as.matrix(wide_pair[, -1])  # Exclude row names column
  rownames(wide_matrix) <- wide_pair[[1]]    # Set row names
  
  # Add site names as rows and columns (ensure all sites are included)
  all_sites <- unique(c(pairwise_result$site1, pairwise_result$site2))
  
  wide_pair <- wide_pair %>%
    complete(site1 = all_sites) %>%  # Add any missing rows
    arrange(site1)
  
  # Convert to matrix
  wide_matrix <- as.matrix(wide_pair[, -1])
  rownames(wide_matrix) <- wide_pair$site1  # Set row names from site1
  
  # Add any missing columns to ensure symmetry
  missing_cols <- setdiff(all_sites, colnames(wide_matrix))
  wide_matrix <- cbind(wide_matrix, matrix(0, nrow = nrow(wide_matrix), ncol = length(missing_cols)))
  colnames(wide_matrix)[(ncol(wide_matrix) - length(missing_cols) + 1):ncol(wide_matrix)] <- missing_cols
  
  # Reorder rows and columns to match
  wide_matrix <- wide_matrix[all_sites, all_sites]
  
  # Replace NAs with 0 and ensure symmetry
  wide_matrix[is.na(wide_matrix)] <- 0
  wide_matrix <- (wide_matrix + t(wide_matrix)) / 2  # Ensure matrix is symmetric
  
  # Add the diagonal trace of 0
  diag(wide_matrix) <- 0
  
  # Run LCBD calculation
  lcbd_result <- LCBD.comp(wide_matrix, sqrt.D = FALSE)
  
  # Combine LCBD results with site identifiers
  lcbd_output <- data.frame(
    site = rownames(wide_matrix),
    LCBD = lcbd_result$LCBD
  )
  
  # Scaling lcbd 
  lcbd_output$LCBD <- (lcbd_output$LCBD - min(lcbd_output$LCBD)) / 
    (max(lcbd_output$LCBD) - min(lcbd_output$LCBD))
  
  # Return results as a list
  return(list(
    pairwise_result = pairwise_result,
    lcbd_result = lcbd_output
  ))
}

# Run the function for multiple years
year.q2 <- lapply(year, Year.betafunction.q2)


saveRDS(year.q0, file = "year.q0.rds")
saveRDS(year.q1, file = "year.q1.rds")
saveRDS(year.q2, file = "year.q2.rds")

year.q0 <- readRDS(file = "year.q0.rds")
year.q1 <- readRDS(file = "year.q1.rds")
year.q2 <- readRDS(file = "year.q2.rds")

######extracting lcbd from large lists#####
#extracting lcbd by lakes
lcbd_results <- lapply(lakes.q0, function(x) x$lcbd_result)
lcbd.L.q0 <- do.call(rbind, lcbd_results)
lcbd.L.q0 <- lcbd.L.q0 %>%
  rename(lake.LCBD.q0 = LCBD)%>%
  rename(lake.date = site)

lcbd_results <- lapply(lakes.q1, function(x) x$lcbd_result)
lcbd.L.q1 <- do.call(rbind, lcbd_results)
lcbd.L.q1 <- lcbd.L.q1 %>%
  rename(lake.LCBD.q1 = LCBD)%>%
  rename(lake.date = site)

lcbd_results <- lapply(lakes.q2, function(x) x$lcbd_result)
lcbd.L.q2 <- do.call(rbind, lcbd_results)
lcbd.L.q2 <- lcbd.L.q2 %>%
  rename(lake.LCBD.q2 = LCBD)%>%
  rename(lake.date = site)

#extracting lcbd by year
lcbd_results <- lapply(year.q0, function(x) x$lcbd_result)
lcbd.Y.q0 <- do.call(rbind, lcbd_results)
lcbd.Y.q0 <- lcbd.Y.q0 %>%
  rename(year.LCBD.q0 = LCBD)%>%
  rename(lake.date = site)

lcbd_results <- lapply(year.q1, function(x) x$lcbd_result)
lcbd.Y.q1 <- do.call(rbind, lcbd_results)
lcbd.Y.q1 <- lcbd.Y.q1 %>%
  rename(year.LCBD.q1 = LCBD)%>%
  rename(lake.date = site)

lcbd_results <- lapply(year.q2, function(x) x$lcbd_result)
lcbd.Y.q2 <- do.call(rbind, lcbd_results)
lcbd.Y.q2 <- lcbd.Y.q2 %>%
  rename(year.LCBD.q2 = LCBD)%>%
  rename(lake.date = site)



#joining dfs
LD.df <- left_join(LD.df, lcbd.L.q0, by = "lake.date")
LD.df <- left_join(LD.df, lcbd.L.q1, by = "lake.date")
LD.df <- left_join(LD.df, lcbd.L.q2, by = "lake.date")
LD.df <- left_join(LD.df, lcbd.Y.q0, by = "lake.date")
LD.df <- left_join(LD.df, lcbd.Y.q1, by = "lake.date")
LD.df <- left_join(LD.df, lcbd.Y.q2, by = "lake.date")

#converting lcbd values less than zero to 0
#lcbd.columns <- c("lake.LCBD.q0","lake.LCBD.q1", "lake.LCBD.q2","year.LCBD.q0", "year.LCBD.q1", "year.LCBD.q2")
#LD.df[lcbd.columns] <- lapply(LD.df[lcbd.columns], function(x) ifelse(x < 0, 0, x))


#Adding WQ diversity metrics to WQ df
WQ.div <- left_join(WQ, LD.df, by = "lake.date")

#####LOAD data saved from directory####
saveRDS(WQ.div, file = "WQ.div.rds")

WQ.div <- readRDS(file = "WQ.div.rds")


summary(WQ.div)


##Looking at Richness
hist(WQ.div$hill.rich)
qqPlot(WQ.div$hill.rich)
#appears normal

#Looking at hill shannon
hist(WQ.div$hill.Shannon)
qqPlot(WQ.div$hill.Shannon)

#Looking at hill inverse simpson
###4 lake dates did not have any taxa and so inverse 
WQ.div[WQ.div$lake.date == "12_2018_176" & WQ.div$hill.InverseSimpson == Inf, "hill.InverseSimpson"] <- 0
WQ.div[WQ.div$lake.date == "16_2018_226" & WQ.div$hill.InverseSimpson == Inf, "hill.InverseSimpson"] <- 0
WQ.div[WQ.div$lake.date == "75_2018_262" & WQ.div$hill.InverseSimpson == Inf, "hill.InverseSimpson"] <- 0
WQ.div[WQ.div$lake.date == "140_2019_140" & WQ.div$hill.InverseSimpson == Inf, "hill.InverseSimpson"] <- 0


hist(WQ.div$hill.InverseSimpson)
qqPlot(WQ.div$hill.InverseSimpson)

hist(WQ.div$lake.LCBD.q0)
qqPlot(WQ.div$lake.LCBD.q0)
hist(WQ.div$lake.LCBD.q1)
qqPlot(WQ.div$lake.LCBD.q1)
hist(WQ.div$lake.LCBD.q2)
qqPlot(WQ.div$lake.LCBD.q2)
hist(WQ.div$year.LCBD.q0)
qqPlot(WQ.div$year.LCBD.q0)
hist(WQ.div$year.LCBD.q1)
qqPlot(WQ.div$year.LCBD.q1)
hist(WQ.div$year.LCBD.q2)
qqPlot(WQ.div$year.LCBD.q2)



#Creating new column of early sample dates with a new column "lake.year"
WQ.div.early <- WQ.div %>%
  filter(Season == "early")%>%
  mutate(lake.year = paste(lakeID, year, sep = "_"))

#looking at duplicated lake.years
  #and then choosing the earlier DoY samples preferentially
duplicates <- duplicated(WQ.div.early$lake.year)
WQ.div.early <- WQ.div.early[!(duplicates & WQ.div.early$DoY !=min(WQ.div.early$DoY[duplicates])),]

#the distribution of samples per year
table(WQ.div.early$year)

#looking at new richness distribution and it is the same
qqPlot(WQ.div.early$hill.rich)
hist(WQ.div.early$hill.rich) 
#plot of richness over time
plot(data = WQ.div.early, hill.rich ~ year)

#looking at new shannon distribution and it is the same
qqPlot(WQ.div.early$hill.Shannon)
hist(WQ.div.early$hill.Shannon) 
#plot of richness over time
plot(data = WQ.div.early, hill.Shannon ~ year)

hist(WQ.div.early$lake.LCBD.q0)
qqPlot(WQ.div.early$lake.LCBD.q0)
hist(WQ.div.early$lake.LCBD.q1)
qqPlot(WQ.div$lake.LCBD.q1)
hist(WQ.div.early$lake.LCBD.q2)
qqPlot(WQ.div.early$lake.LCBD.q2)
hist(WQ.div.early$year.LCBD.q0)
qqPlot(WQ.div.early$year.LCBD.q0)
hist(WQ.div.early$year.LCBD.q1)
qqPlot(WQ.div.early$year.LCBD.q1)
hist(WQ.div.early$year.LCBD.q2)
qqPlot(WQ.div.early$year.LCBD.q2)



#####TEMP TIME MODEL#####
######Richness####


##Determine the total number of NAs within Water Temp Column
na_count <- sum(is.na(WQ.div.early$WaterTemp_z))
##Removing NAs
Temp.div.early <- WQ.div.early[complete.cases(WQ.div.early$WaterTemp_z), ]  
##There are no NAs in Rich
na_count <- sum(is.na(Temp.div.early$hill.rich))

#Looking at how many times each lake was sampled
lakeID.df <- as.data.frame(table(Temp.div.early$lakeID))
table(lakeID.df$Freq)
lakeID.df <- lakeID.df %>%
  rename(lakeID = Var1)

mean(Temp.div.early$Temp.water.degC)



#Filtering by lakeID to remove lakes which were not sampled as frequently
Temp.div.early <- Temp.div.early %>%
  left_join(lakeID.df, by = "lakeID") %>%
  filter(Freq >= 9) #Freq >= is selecting how many observations each lake needed to be included  

#Calculating Z
Temp.div.early$WaterTemp_z <- scale(Temp.div.early$Temp.water.degC, center = TRUE, scale = TRUE)

hist(Temp.div.early$hill.rich)

##building model
time.model <- brm(formula = hill.rich ~ WaterTemp_z + (1 + WaterTemp_z | lakeID),
                  data = Temp.div.early,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .95))


#Model Checks
pairs(time.model)
summary(time.model)
pp_check(time.model, ndraws = 100)

time_model_ppcheck_Alpha_q0 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/temp_time_Alpha_q0.png",
  plot = time_model_ppcheck_Alpha_q0,
  width = 3,
  height = 2,
  units = "in"
)

diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)



pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

summary(pred.a)
summary(pred.b)

saveRDS(time.model, file = "time_model_temp_rich.rds")
saveRDS(pred.a, file = "pred.A_timemodel_temp_rich.rds")
saveRDS(pred.b, file = "pred.B_timemodel_temp_rich.rds")

time.model <- readRDS(file = "time_model_temp_rich.rds")
pred.a <- readRDS(file = "pred.A_timemodel_temp_rich.rds")
pred.b <- readRDS(file = "pred.B_timemodel_temp_rich.rds")



#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x2 <- pred.a$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)
pred.b$x2 <- pred.b$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)

shared_xlim <- range(Temp.div.early$Temp.water.degC, na.rm = TRUE)
fixed_margins <- theme(plot.margin = unit(c(0.5, 1, 0.5, 1.5), "cm"))

pdf(file = "lakeID.RandomEffect_WaterTemp_Rich.pdf")
ggpred_plot <- (ggplot() + 
    geom_ribbon(data = pred.a, 
                aes(ymin = conf.low, ymax = conf.high, x = x2), 
                alpha = 0.2) +  # error band
    geom_point(data = Temp.div.early,                      # adding the raw data (scaled values)
                aes(x = Temp.water.degC, y = hill.rich, colour = lakeID),
               alpha = 0.1,
               shape = 16) +
    geom_line(data = pred.b, aes(x2, predicted, colour = group),
              alpha = 0.6) +      
    geom_line(data = pred.a, aes(x2, predicted), linewidth = 1.5, linetype = "longdash") +
    labs(y = "α (q = 0)") + 
    theme_bw()) +
  theme(
    axis.title.x = element_blank(),  # Removes X-axis label
    axis.text.x = element_blank(),   # Removes X-axis numbers
    axis.ticks.x = element_blank(),  # Removes X-axis tick marks
    axis.title.y = element_text(size = 14),  # Y-axis label size
    axis.text = element_text(size = 12)  # Axis tick labels size 
    ) +
  guides(color = "none") +
  fixed_margins + 
  coord_cartesian(xlim = shared_xlim)

print(ggpred_plot)
dev.off()
ggpred_plot
TempPlot.A.q0 <- ggpred_plot 

######Hill.Shannon####
##building model
hist(Temp.div.early$hill.Shannon)


time.model <- brm(formula = hill.Shannon ~ WaterTemp_z + (1 + WaterTemp_z | lakeID),
                  data = Temp.div.early,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .8))
#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 5)

pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

time_model_ppcheck_Alpha_q1 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/temp_time_Alpha_q1.png",
  plot = time_model_ppcheck_Alpha_q1,
  width = 3,
  height = 2,
  units = "in")


pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_temp_shan")
saveRDS(pred.a, file = "pred.A_timemodel_temp_shan")
saveRDS(pred.b, file = "pred.B_timemodel_temp_shan")

time.model <- readRDS(file = "time_model_temp_shan")
pred.a <- readRDS(file = "pred.A_timemodel_temp_shan")
pred.b <- readRDS(file = "pred.B_timemodel_temp_shan")

summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)


pdf(file = "lakeID.RandomEffect_WaterTemp_Hill.Shannon.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Temp.div.early,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = hill.Shannon, colour = lakeID),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Water Temperature (C°)", y = "α (q = 1)") + 
                  theme_bw()) +
  theme(
    axis.title.x = element_blank(),  # Removes X-axis label
    axis.text.x = element_blank(),   # Removes X-axis numbers
    axis.ticks.x = element_blank(),  # Removes X-axis tick marks
    axis.title.y = element_text(size = 14),  # Y-axis label size
    axis.text = element_text(size = 12)  # Axis tick labels size 
  ) +
  guides(color = "none")+
  fixed_margins + 
  coord_cartesian(xlim = shared_xlim)
ggpred_plot
dev.off()
ggpred_plot
TempPlot.A.q1 <- ggpred_plot

######Hill.Inverse Simpson#####
##building model
hist(Temp.div.early$hill.InverseSimpson)


time.model <- brm(formula = hill.InverseSimpson ~ WaterTemp_z + (1 + WaterTemp_z | lakeID),
                  data = Temp.div.early,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .8))
#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 4)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

time_model_ppcheck_Alpha_q2 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/temp_time_Alpha_q2.png",
  plot = time_model_ppcheck_Alpha_q2,
  width = 3,
  height = 2,
  units = "in")


pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

#saving R objects
saveRDS(time.model, file = "time_model_temp_simp")
saveRDS(pred.a, file = "pred.A_timemodel_temp_simp")
saveRDS(pred.b, file = "pred.B_timemodel_temp_simp")

time.model <- readRDS(file = "time_model_temp_simp")
pred.a <- readRDS(file = "pred.A_timemodel_temp_simp")
pred.b <- readRDS(file = "pred.B_timemodel_temp_simp")


summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)




pdf(file = "lakeID.RandomEffect_WaterTemp_Hill.InverseSimpson.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Temp.div.early,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = hill.InverseSimpson, colour = lakeID),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Water Temperature (C°)", y = "α (q = 2)") + 
                  theme_bw()) +
  theme(
    axis.title.x = element_blank(),  # Removes X-axis label
    axis.text.x = element_blank(),   # Removes X-axis numbers
    axis.ticks.x = element_blank(),  # Removes X-axis tick marks
    axis.title.y = element_text(size = 14),  # Y-axis label size
    axis.text = element_text(size = 12)  # Axis tick labels size 
  ) +
  guides(color = "none")+
  fixed_margins
ggpred_plot
dev.off()
ggpred_plot
TempPlot.A.q2 <- ggpred_plot


######LCBD Lake q0#####
##building model
hist(Temp.div.early$lake.LCBD.q0)
qqPlot(Temp.div.early$lake.LCBD.q0)

summary(Temp.div.early)

#standard “Smithson & Verkuilen” (2006) transformation so that we can use zero inflated beta with zeros and 1s
n <- nrow(Temp.div.early)
Temp.div.early$lake.LCBD.q0 <- (Temp.div.early$lake.LCBD.q0 * (n - 1) + 0.5) / n

time.model <- brm(formula = lake.LCBD.q0 ~ WaterTemp_z + (1 + WaterTemp_z | lakeID),
                  data = Temp.div.early,
                    family = Beta(link = "cloglog"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .99))


#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 4)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

time_model_ppcheck_LCBD_q0 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/temp_time_LCBD_q0.png",
  plot = time_model_ppcheck_LCBD_q0,
  width = 3,
  height = 2,
  units = "in")

pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

#saving R objects
saveRDS(time.model, file = "time_model_temp_LCBD.L.q0")
saveRDS(pred.a, file = "pred.A_timemodel_temp_LCBD.L.q0")
saveRDS(pred.b, file = "pred.B_timemodel_temp_LCBD.L.q0")

time.model <- readRDS(file = "time_model_temp_LCBD.L.q0")
pred.a <- readRDS(file = "pred.A_timemodel_temp_LCBD.L.q0")
pred.b <- readRDS(file = "pred.B_timemodel_temp_LCBD.L.q0")


summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)




pdf(file = "lakeID.RandomEffect_WaterTemp_LCBD.L.q0.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Temp.div.early,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = lake.LCBD.q0, colour = lakeID),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Water Temperature (C)", y = "LCBD (q = 0)") + 
                  theme_bw()) +
  theme(
    axis.title.x = element_blank(),  # Removes X-axis label
    axis.text.x = element_blank(),   # Removes X-axis numbers
    axis.ticks.x = element_blank(),  # Removes X-axis tick marks
    axis.title.y = element_text(size = 14),  # Y-axis label size
    axis.text = element_text(size = 12)  # Axis tick labels size 
                  ) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot
TempPlot.B.q0 <- ggpred_plot




######LCBD Lake q1#####
##building model
hist(Temp.div.early$lake.LCBD.q1)
qqPlot(Temp.div.early$lake.LCBD.q1)
summary(Temp.div.early)

#standard “Smithson & Verkuilen” (2006) transformation so that we can use zero inflated beta with zeros and 1s
n <- nrow(Temp.div.early)
Temp.div.early$lake.LCBD.q1 <- (Temp.div.early$lake.LCBD.q1 * (n - 1) + 0.5) / n

time.model <- brm(formula = lake.LCBD.q1 ~ WaterTemp_z + (1 + WaterTemp_z | lakeID),
                  data = Temp.div.early,
                  family = Beta(link = "cloglog"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .9))
#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 4)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

time_model_ppcheck_LCBD_q1 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/temp_time_LCBD_q1.png",
  plot = time_model_ppcheck_LCBD_q1,
  width = 3,
  height = 2,
  units = "in")

pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

#saving R objects
saveRDS(time.model, file = "time_model_temp_LCBD.L.q1")
saveRDS(pred.a, file = "pred.A_timemodel_temp_LCBD.L.q1")
saveRDS(pred.b, file = "pred.B_timemodel_temp_LCBD.L.q1")

time.model <- readRDS(file = "time_model_temp_LCBD.L.q1")
pred.a <- readRDS(file = "pred.A_timemodel_temp_LCBD.L.q1")
pred.b <- readRDS(file = "pred.B_timemodel_temp_LCBD.L.q1")


summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)




pdf(file = "lakeID.RandomEffect_WaterTemp_LCBD.L.q1.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Temp.div.early,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = lake.LCBD.q1, colour = lakeID),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Water Temperature (C)", y = "LCBD (q = 1)") + 
                  theme_bw()) +
  theme(
    axis.title.x = element_blank(),  # Removes X-axis label
    axis.text.x = element_blank(),   # Removes X-axis numbers
    axis.ticks.x = element_blank(),  # Removes X-axis tick marks
    axis.title.y = element_text(size = 14),  # Y-axis label size
    axis.text = element_text(size = 12)  # Axis tick labels size 
  ) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot
TempPlot.B.q1 <- ggpred_plot

######LCBD Lake q2#####
##building model
hist(Temp.div.early$lake.LCBD.q2)
qqPlot(Temp.div.early$lake.LCBD.q2)

#standard “Smithson & Verkuilen” (2006) transformation so that we can use beta with zeros and 1s
n <- nrow(Temp.div.early)
Temp.div.early$lake.LCBD.q2 <- (Temp.div.early$lake.LCBD.q2 * (n - 1) + 0.5) / n

summary(Temp.div.early)

time.model <- brm(formula = lake.LCBD.q2 ~ WaterTemp_z + (1 + WaterTemp_z | lakeID),
                  data = Temp.div.early,
                  family = Beta(link = "cloglog"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .8))
#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 4)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

time_model_ppcheck_LCBD_q2 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/temp_time_LCBD_q2.png",
  plot = time_model_ppcheck_LCBD_q2,
  width = 3,
  height = 2,
  units = "in")


pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

#saving R objects
saveRDS(time.model, file = "time_model_temp_LCBD.L.q2")
saveRDS(pred.a, file = "pred.A_timemodel_temp_LCBD.L.q2")
saveRDS(pred.b, file = "pred.B_timemodel_temp_LCBD.L.q2")

time.model <- readRDS(file = "time_model_temp_LCBD.L.q2")
pred.a <- readRDS(file = "pred.A_timemodel_temp_LCBD.L.q2")
pred.b <- readRDS(file = "pred.B_timemodel_temp_LCBD.L.q2")


summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)




pdf(file = "lakeID.RandomEffect_WaterTemp_LCBD.L.q2.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Temp.div.early,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = lake.LCBD.q2, colour = lakeID),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Water Temperature (C°)", y = "LCBD (q = 2)") + 
                  theme_bw()) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),  # Y-axis label size
    axis.text = element_text(size = 12)  # Axis tick labels size 
  ) +
  guides(color = "none") 
ggpred_plot
dev.off()
ggpred_plot
TempPlot.B.q2 <- ggpred_plot

#####TEMP SPACE MODEL#####
######Richness##########
##building model
Temp.div.early$year <- as.factor(Temp.div.early$year)
summary(Temp.div.early)

space.model <- brm(formula = hill.rich ~ WaterTemp_z + (1 + WaterTemp_z | year),
                  data = Temp.div.early,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .8))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence

#model checks
pairs(space.model)
summary(space.model)
print(space.model, digits = 4)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_Alpha_q0 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/temp_space_Alpha_q0.png",
  plot = space_model_ppcheck_Alpha_q0,
  width = 3,
  height = 2,
  units = "in")

pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_rich")
saveRDS(pred.a, file = "pred.a_spacemodel_temp_rich")
saveRDS(pred.b, file = "pred.b_spacemodel_temp_rich")

space.model <- readRDS(file = "space_model_temp_rich")
pred.a <- readRDS(file = "pred.a_spacemodel_temp_rich")
pred.b <- readRDS(file = "pred.b_spacemodel_temp_rich")

summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_Rich.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Temp.div.early,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = hill.rich, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Water Temperature (C)", y = "Richness (as Hill numbers)") + 
                  theme_bw()) +
  theme(
      axis.title.x = element_blank(),  # Removes X-axis label
      axis.text.x = element_blank(),   # Removes X-axis numbers
      axis.ticks.x = element_blank(),  # Removes X-axis tick marks
      axis.title.y = element_blank(),  # Y-axis label size
      axis.text = element_blank()  # Axis tick labels size 
                  ) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot
TempPlot.A.q0.S <- ggpred_plot


#####TEMP SPACE MODEL
#######Hill.Shannon####
hist(Temp.div.early$hill.Shannon)
qqPlot(Temp.div.early$hill.Shannon)


##building model
space.model <- brm(formula = hill.Shannon ~ WaterTemp_z + (1 + WaterTemp_z | year),
                  data = Temp.div.early,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .8))
#model checks
summary(space.model)
print(space.model, digits = 4)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_Alpha_q1 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/temp_space_Alpha_q1.png",
  plot = space_model_ppcheck_Alpha_q1,
  width = 3,
  height = 2,
  units = "in")

pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_shan.rds")
saveRDS(pred.a, file = "pred.A_spacemodel_temp_shan.rds")
saveRDS(pred.b, file = "pred.B_spacemodel_temp_shan.rds")

space.model <- readRDS(file = "space_model_temp_shan.rds")
pred.a <- readRDS(file = "pred.A_spacemodel_temp_shan.rds")
pred.b <- readRDS(file = "pred.B_spacemodel_temp_shan.rds")

summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_Shannon.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Temp.div.early,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = hill.Shannon, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Water Temperature (C)", y = "Richness (as Hill numbers)") + 
                  theme_bw()) +
  theme(
    axis.title.x = element_blank(),  # Removes X-axis label
    axis.text.x = element_blank(),   # Removes X-axis numbers
    axis.ticks.x = element_blank(),  # Removes X-axis tick marks
    axis.title.y = element_blank(),  # Y-axis label size
    axis.text = element_blank()  # Axis tick labels size 
  ) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot
TempPlot.A.q1.S <- ggpred_plot

#####TEMP SPACE MODEL
########Hill.Inverse Simpson####

space.model <- brm(formula = hill.InverseSimpson ~ WaterTemp_z + (1 + WaterTemp_z | year),
                   data = Temp.div.early,
                   family = Gamma(link = "log"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .8))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence
#model checks
pairs(space.model)
summary(space.model)
print(space.model, digit = 4)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_Alpha_q2 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/temp_space_Alpha_q2.png",
  plot = space_model_ppcheck_Alpha_q2,
  width = 3,
  height = 2,
  units = "in")

pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_simp")
saveRDS(pred.a, file = "pred.A_spacemodel_temp_simp")
saveRDS(pred.b, file = "pred.B_spacemodel_temp_simp")

space.model <- readRDS(file = "space_model_temp_simp")
pred.a <- readRDS(file = "pred.A_spacemodel_temp_simp")
pred.b <- readRDS(file = "pred.B_spacemodel_temp_simp")

summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_Inverse Simpson.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Temp.div.early,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = hill.InverseSimpson, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Water Temperature (C)", y = "Inverse Simpson (as Hill numbers)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()

ggpred_plot
TempPlot.A.q2.S <- ggpred_plot


######LCBD Year q0##########
##building model
Temp.div.early$year <- as.factor(Temp.div.early$year)
summary(Temp.div.early)

hist(Temp.div.early$year.LCBD.q0)
qqPlot(Temp.div.early$year.LCBD.q0)

#standard “Smithson & Verkuilen” (2006) transformation so that we can use zero inflated beta with zeros and 1s
n <- nrow(Temp.div.early)
Temp.div.early$year.LCBD.q0 <- (Temp.div.early$year.LCBD.q0 * (n - 1) + 0.5) / n
Temp.div.early$year.LCBD.q1 <- (Temp.div.early$year.LCBD.q1 * (n - 1) + 0.5) / n
Temp.div.early$year.LCBD.q2 <- (Temp.div.early$year.LCBD.q2 * (n - 1) + 0.5) / n

space.model <- brm(formula = year.LCBD.q0 ~ WaterTemp_z + (1 + WaterTemp_z | year),
                   data = Temp.div.early,
                   family = Beta(link = "cloglog"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = 0.8, max_treedepth = 12))
#model checks
pairs(space.model)
summary(space.model)
print(space.model, digits = 7)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_LCBD_q0 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/temp_space_LCBD_q0.png",
  plot = space_model_ppcheck_LCBD_q0,
  width = 3,
  height = 2,
  units = "in")

pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_year.LCBD.q0")
saveRDS(pred.a, file = "pred.a_spacemodel_temp_year.LCBD.q0")
saveRDS(pred.b, file = "pred.b_spacemodel_temp_year.LCBD.q0")

space.model <- readRDS(file = "space_model_temp_year.LCBD.q0")
pred.a <- readRDS(file = "pred.a_spacemodel_temp_year.LCBD.q0")
pred.b <- readRDS(file = "pred.b_spacemodel_temp_year.LCBD.q0")

summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_year.LCBD.q0.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Temp.div.early,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = year.LCBD.q0, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Water Temperature (C)", y = "LCBD (q = 0)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot
TempPlot.B.q0.S <- ggpred_plot


######LCBD Year q1##########
##building model
Temp.div.early$year <- as.factor(Temp.div.early$year)
summary(Temp.div.early)

hist(Temp.div.early$year.LCBD.q1)
qqPlot(Temp.div.early$year.LCBD.q1)

space.model <- brm(formula = year.LCBD.q1 ~ WaterTemp_z + (1 + WaterTemp_z | year),
                   data = Temp.div.early,
                   family = Beta(link = "cloglog"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .8, max_treedepth = 12))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence

#model checks
pairs(space.model)
summary(space.model)
print(space.model, digits = 6)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_LCBD_q1 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/temp_space_LCBD_q1.png",
  plot = space_model_ppcheck_LCBD_q1,
  width = 3,
  height = 2,
  units = "in")

pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_year.LCBD.q1")
saveRDS(pred.a, file = "pred.a_spacemodel_temp_year.LCBD.q1")
saveRDS(pred.b, file = "pred.b_spacemodel_temp_year.LCBD.q1")

space.model <- readRDS(file = "space_model_temp_year.LCBD.q1")
pred.a <- readRDS(file = "pred.a_spacemodel_temp_year.LCBD.q1")
pred.b <- readRDS(file = "pred.b_spacemodel_temp_year.LCBD.q1")

summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_year.LCBD.q1.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Temp.div.early,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = year.LCBD.q1, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Water Temperature (C)", y = "LCBD (q = 1)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot
TempPlot.B.q1.S <- ggpred_plot

######LCBD Year q2##########
##building model
Temp.div.early$year <- as.factor(Temp.div.early$year)
summary(Temp.div.early)

hist(Temp.div.early$year.LCBD.q2)
qqPlot(Temp.div.early$year.LCBD.q2)

space.model <- brm(formula = year.LCBD.q2 ~ WaterTemp_z + (1 + WaterTemp_z | year),
                   data = Temp.div.early,
                   family = Beta(link = "cloglog"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .85, max_treedepth = 12))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence

#model checks
pairs(space.model)
summary(space.model)
print(space.model, digits = 6)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_LCBD_q2 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/temp_space_LCBD_q2.png",
  plot = space_model_ppcheck_LCBD_q2,
  width = 3,
  height = 2,
  units = "in")

pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_year.LCBD.q2")
saveRDS(pred.a, file = "pred.a_spacemodel_temp_year.LCBD.q2")
saveRDS(pred.b, file = "pred.b_spacemodel_temp_year.LCBD.q2")

space.model <- readRDS(file = "space_model_temp_year.LCBD.q2")
pred.a <- readRDS(file = "pred.a_spacemodel_temp_year.LCBD.q2")
pred.b <- readRDS(file = "pred.b_spacemodel_temp_year.LCBD.q2")

summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.early$Temp.water.degC) + mean(Temp.div.early$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_year.LCBD.q2.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Temp.div.early,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = year.LCBD.q2, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Water Temperature (C)", y = "LCBD (q = 2)") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot

dev.off()
ggpred_plot
TempPlot.B.q2.S <- ggpred_plot

#####Nitrogen TIME MODEL#####
######Richness#####



##TN from grace's code using TN instead of nitrate-nitrite
Nit.div.early <- WQ.div.early %>%
  group_by(lake.date) %>%
  mutate(
    Total_Nitrogen = if_else(is.na(Nutrient.nitrogen_mg.l), Kjeldahl_nitrogen.as_N. + Nitrate.Nitrite_.as_N._mg.l, Nutrient.nitrogen_mg.l)
  ) %>%
  ungroup()  
summary(WQ.div.early)

#Calculating Z
Nit.div.early$TN_z <- scale(Nit.div.early$Total_Nitrogen, center = TRUE, scale = TRUE)
summary(Nit.div.early)

##Determine the total number of NAs within Water Nitrate Column
na_count <- sum(is.na(Nit.div.early$TN_z))




##Removing NAs
Nit.div.early <- Nit.div.early[complete.cases(Nit.div.early$TN_z), ]  
##There are no NAs in Rich
na_count <- sum(is.na(Nit.div.early$hill.rich))

#Looking at how many times each lake was sampled
lakeID.df <- as.data.frame(table(Nit.div.early$lakeID))
table(lakeID.df$Freq)
lakeID.df <- lakeID.df %>%
  rename(lakeID = Var1)

#Filtering by lakeID to remove lakes which were not sampled as frequently
Nit.div.early <- Nit.div.early %>%
  left_join(lakeID.df, by = "lakeID") %>%
  filter(Freq >= 9) #Freq >= is selecting how many observations each lake needed to be included  



##building model
time.model <- brm(formula = hill.rich ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.early,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .95))


getwd()
#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 4)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

time_model_ppcheck_Alpha_q0 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TN_time_Alpha_q0.png",
  plot = time_model_ppcheck_Alpha_q0,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(time.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("TN_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_TN_rich")
saveRDS(pred.a, file = "pred.A_timemodel_TN_rich")
saveRDS(pred.b, file = "pred.B_timemodel_TN_rich")

time.model <- readRDS(file = "time_model_TN_rich")
pred.a <- readRDS(file = "pred.A_timemodel_TN_rich")
pred.b <- readRDS(file = "pred.B_timemodel_TN_rich")


summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitriogen_Rich.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Nit.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.rich, colour = lakeID),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "α (q = 0)") + 
                  theme_bw() +
  guides(color = "none"))
ggpred_plot
dev.off()
ggpred_plot

TNPlot.A.q0.T <- ggpred_plot



######Hill.Shannon#####
##building model
hist(Nit.div.early$hill.Shannon)


time.model <- brm(formula = hill.Shannon ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.early,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .9))
#Model Checks
pairs(time.model)
summary(time.model)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

time_model_ppcheck_Alpha_q1 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TN_time_Alpha_q1.png",
  plot = time_model_ppcheck_Alpha_q1,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(time.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("TN_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_TN_shan")
saveRDS(pred.a, file = "pred.A_timemodel_TN_shan")
saveRDS(pred.b, file = "pred.B_timemodel_TN_shan")

time.model <- readRDS(file = "time_model_TN_shan")
pred.a <- readRDS(file = "pred.A_timemodel_TN_shan")
pred.b <- readRDS(file = "pred.B_timemodel_TN_shan")


summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitrogen_Hill.Shannon.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Nit.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.Shannon, colour = lakeID),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha =0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "α (q = 1)") + 
                  theme_bw() +
  guides(color = "none"))
ggpred_plot
dev.off()
ggpred_plot

TNPlot.A.q1.T <- ggpred_plot

######Hill.Inverse Simpson####
## building model
hist(Nit.div.early$hill.InverseSimpson)


time.model <- brm(formula = hill.InverseSimpson ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.early,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .97))
#Model Checks
pairs(time.model)
summary(time.model)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

time_model_ppcheck_Alpha_q2 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TN_time_Alpha_q2.png",
  plot = time_model_ppcheck_Alpha_q2,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(time.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("TN_z [all]", "lakeID"), type = "re", ppd = FALSE)

summary(pred.a)
summary(pred.b)

saveRDS(time.model, file = "time_model_TN_simp")
saveRDS(pred.a, file = "pred.A_timemodel_TN_simp")
saveRDS(pred.b, file = "pred.B_timemodel_TN_simp")

time.model <- readRDS(file = "time_model_TN_simp")
pred.a <- readRDS(file = "pred.A_timemodel_TN_simp")
pred.b <- readRDS(file = "pred.B_timemodel_TN_simp")


#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitrogen_Hill.InverseSimpson.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Nit.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.InverseSimpson, colour = lakeID),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  
                  labs(x = "Total Nitrogen (mg/l)", y = "α (q = 2)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot

TNPlot.A.q2.T <- ggpred_plot
######LCBD Lake q0######

#standard “Smithson & Verkuilen” (2006) transformation so that we can use zero inflated beta with zeros and 1s
n <- nrow(Nit.div.early)
Nit.div.early$lake.LCBD.q0 <- (Nit.div.early$lake.LCBD.q0 * (n - 1) + 0.5) / n
Nit.div.early$lake.LCBD.q1 <- (Nit.div.early$lake.LCBD.q1 * (n - 1) + 0.5) / n
Nit.div.early$lake.LCBD.q2 <- (Nit.div.early$lake.LCBD.q2 * (n - 1) + 0.5) / n

##building model
time.model <- brm(formula = lake.LCBD.q0 ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.early,
                  family = Beta(link = "cloglog"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .85))



#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 4)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

time_model_ppcheck_LCBD_q0 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TN_time_LCBD_q0.png",
  plot = time_model_ppcheck_LCBD_q0,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(time.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("TN_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_TN_lake.LCBD.q0")
saveRDS(pred.a, file = "pred.A_timemodel_TN_lake.LCBD.q0")
saveRDS(pred.b, file = "pred.B_timemodel_TN_lake.LCBD.q0")

time.model <- readRDS(file = "time_model_TN_lake.LCBD.q0")
pred.a <- readRDS(file = "pred.A_timemodel_TN_lake.LCBD.q0")
pred.b <- readRDS(file = "pred.B_timemodel_TN_lake.LCBD.q0")


summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitrogen_lake.LCBD.q0.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Nit.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = lake.LCBD.q0, colour = lakeID),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "LCBD (q = 0)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot

TNPlot.B.q0.T <- ggpred_plot

######LCBD Lake q1######

##building model
time.model <- brm(formula = lake.LCBD.q1 ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.early,
                  family = Beta(link = "cloglog"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .85))



#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 4)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

pred.a <- ggpredict(time.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("TN_z [all]", "lakeID"), type = "re", ppd = FALSE)

time_model_ppcheck_LCBD_q1 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TN_time_LCBD_q1.png",
  plot = time_model_ppcheck_LCBD_q1,
  width = 3,
  height = 2,
  units = "in"
)

saveRDS(time.model, file = "time_model_TN_lake.LCBD.q1")
saveRDS(pred.a, file = "pred.A_timemodel_TN_lake.LCBD.q1")
saveRDS(pred.b, file = "pred.B_timemodel_TN_lake.LCBD.q1")

time.model <- readRDS(file = "time_model_TN_lake.LCBD.q1")
pred.a <- readRDS(file = "pred.A_timemodel_TN_lake.LCBD.q1")
pred.b <- readRDS(file = "pred.B_timemodel_TN_lake.LCBD.q1")


summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitrogen_lake.LCBD.q1.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Nit.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = lake.LCBD.q1, colour = lakeID),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "LCBD (q = 1)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot

TNPlot.B.q1.T <- ggpred_plot

######LCBD Lake q2######

##building model
time.model <- brm(formula = lake.LCBD.q2 ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.early,
                  family = Beta(link = "cloglog"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .85))



#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 4)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

time_model_ppcheck_LCBD_q2 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TN_time_LCBD_q2.png",
  plot = time_model_ppcheck_LCBD_q2,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(time.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("TN_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_TN_lake.LCBD.q2")
saveRDS(pred.a, file = "pred.A_timemodel_TN_lake.LCBD.q2")
saveRDS(pred.b, file = "pred.B_timemodel_TN_lake.LCBD.q2")

time.model <- readRDS(file = "time_model_TN_lake.LCBD.q2")
pred.a <- readRDS(file = "pred.A_timemodel_TN_lake.LCBD.q2")
pred.b <- readRDS(file = "pred.B_timemodel_TN_lake.LCBD.q2")


summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitrogen_lake.LCBD.q2.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Nit.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = lake.LCBD.q2, colour = lakeID),
                             alpha = 0.1, 
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Nitrate Nitrite (mg/l)", y = "LCBD (q = 2)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot

TNPlot.B.q2.T <- ggpred_plot



#####Nitrogen SPACE MODEL#####
######Hill.Richness####
##building model
Nit.div.early$year <- as.factor(Nit.div.early$year)
summary(Nit.div.early)

space.model <- brm(formula = hill.rich ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.early,
                   family = Gamma(link = "log"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .8))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence

#model checks
pairs(space.model)
summary(space.model)
print(space.model, digit = 4)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_Alpha_q0 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TN_space_Alpha_q0.png",
  plot = space_model_ppcheck_Alpha_q0,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_TN_rich")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_rich")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_rich")

space.model <- readRDS(file = "space_model_TN_rich")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_rich")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_rich")


summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_Rich.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Nit.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.rich, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "α (q = 0)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot

TNPlot.A.q0.S <- ggpred_plot

#######Hill.Shannon########
hist(Nit.div.early$hill.Shannon)
qqPlot(Nit.div.early$hill.Shannon)


##building model
space.model <- brm(formula = hill.Shannon ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.early,
                   family = Gamma(link = "log"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .8))
#model checks
pairs(space.model)
print(space.model, digits = 4)
summary(space.model)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_Alpha_q1 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TN_space_Alpha_q1.png",
  plot = space_model_ppcheck_Alpha_q1,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_TN_shan")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_shan")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_shan")

space.model <- readRDS(file = "space_model_TN_shan")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_shan")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_shan")


summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_Shannon.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Nit.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.Shannon, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "α (q = 1)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot

TNPlot.A.q1.S <- ggpred_plot

#######hill.Inverse Simpson####

space.model <- brm(formula = hill.InverseSimpson ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.early,
                   family = Gamma(link = "log"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .8))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence
#model checks
pairs(space.model)
summary(space.model)

pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_Alpha_q2 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TN_space_Alpha_q2.png",
  plot = space_model_ppcheck_Alpha_q2,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]","year"), type = "re", ppd = FALSE)

summary(pred.a)
summary(pred.b)

saveRDS(space.model, file = "space_model_TN_simp")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_simp")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_simp")

space.model <- readRDS(file = "space_model_TN_simp")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_simp")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_simp")


#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_InverseSimpson.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Nit.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.InverseSimpson, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x,predicted), linewidth = 1.5, linetype = "longdash") +
                  
                  labs(x = "Total Nitrogen (mg/l)", y = "α (q = 2)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()

ggpred_plot

TNPlot.A.q2.S <- ggpred_plot

######LCBD Year q0####
##building model
Nit.div.early$year <- as.factor(Nit.div.early$year)
summary(Nit.div.early)

#standard “Smithson & Verkuilen” (2006) transformation so that we can use zero inflated beta with zeros and 1s
n <- nrow(Nit.div.early)
Nit.div.early$year.LCBD.q0 <- (Nit.div.early$year.LCBD.q0 * (n - 1) + 0.5) / n
Nit.div.early$year.LCBD.q1 <- (Nit.div.early$year.LCBD.q1 * (n - 1) + 0.5) / n
Nit.div.early$year.LCBD.q2 <- (Nit.div.early$year.LCBD.q2 * (n - 1) + 0.5) / n

space.model <- brm(formula = year.LCBD.q0 ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.early,
                   family = Beta(link = "cloglog"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .8, max_treedepth = 12))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence

#model checks
pairs(space.model)
summary(space.model)
print(space.model, digits = 7)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_LCBD_q0 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TN_space_LCBD_q0.png",
  plot = space_model_ppcheck_LCBD_q0,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_TN_year.LCBD.q0")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_year.LCBD.q0")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_year.LCBD.q0")

space.model <- readRDS(file = "space_model_TN_year.LCBD.q0")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_year.LCBD.q0")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_year.LCBD.q0")


summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_year.LCBD.q0.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Nit.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = year.LCBD.q0, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Nitrate Nitrite (mg/l)", y = "LCBD (q = 1)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot

TNPlot.B.q0.S <- ggpred_plot

######LCBD Year q1####
##building model
Nit.div.early$year <- as.factor(Nit.div.early$year)
summary(Nit.div.early)

space.model <- brm(formula = year.LCBD.q1 ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.early,
                   family = Beta(link = "cloglog"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .9, max_treedepth = 12))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence

#model checks
pairs(space.model)
print(space.model, digits = 7)
summary(space.model)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_LCBD_q1 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TN_space_LCBD_q1.png",
  plot = space_model_ppcheck_LCBD_q1,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_TN_year.LCBD.q1")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_year.LCBD.q1")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_year.LCBD.q1")

space.model <- readRDS(file = "space_model_TN_year.LCBD.q1")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_year.LCBD.q1")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_year.LCBD.q1")


summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_year.LCBD.q1.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Nit.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = year.LCBD.q1, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "LCBD (q = 1)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot

TNPlot.B.q1.S <- ggpred_plot

######LCBD Year q2####
##building model
Nit.div.early$year <- as.factor(Nit.div.early$year)
summary(Nit.div.early)

space.model <- brm(formula = year.LCBD.q2 ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.early,
                   family = Beta(link = "cloglog"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .9, max_treedepth = 12))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence

#model checks
pairs(space.model)
summary(space.model)
print(space.model, digits = 5)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_LCBD_q2 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TN_space_LCBD_q2.png",
  plot = space_model_ppcheck_LCBD_q2,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_TN_year.LCBD.q2")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_year.LCBD.q2")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_year.LCBD.q2")

space.model <- readRDS(file = "space_model_TN_year.LCBD.q2")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_year.LCBD.q2")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_year.LCBD.q2")


summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.early$Total_Nitrogen) + mean(Nit.div.early$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_year.LCBD.q2.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = Nit.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = year.LCBD.q2, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "LCBD (q = 2)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot

TNPlot.B.q2.S <- ggpred_plot

#####Phosphorus TIME MODEL#####
######Hill.Richness####

##Determine the total number of NAs within Water Nitrate Column
na_count <- sum(is.na(WQ.div.early$Phosphorus_z))
##Removing NAs
TP.div.early <- WQ.div.early[complete.cases(WQ.div.early$Phosphorus_z), ]  
##There are no NAs in Rich
na_count <- sum(is.na(TP.div.early$hill.rich))

#Looking at how many times each lake was sampled
lakeID.df <- as.data.frame(table(TP.div.early$lakeID))
table(lakeID.df$Freq)
lakeID.df <- lakeID.df %>%
  rename(lakeID = Var1)

#Filtering by lakeID to remove lakes which were not sampled as frequently
TP.div.early <- TP.div.early %>%
  left_join(lakeID.df, by = "lakeID") %>%
  filter(Freq >= 9) #Freq >= is selecting how many observations each lake needed to be included  

#Calculating Z
TP.div.early$Phosphorus_z <- scale(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P., center = TRUE, scale = TRUE)
summary(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
summary(TP.div.early$Phosphorus_z)
summary(TP.div.early$Total_Phosphorus)

##building model
time.model <- brm(formula = hill.rich ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.early,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .8))

#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 3)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

time_model_ppcheck_Alpha_q0 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TP_time_Alpha_q0.png",
  plot = time_model_ppcheck_Alpha_q0,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_Phos_rich")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_rich")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_rich")

time.model <- readRDS(file = "time_model_Phos_rich")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_rich")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_rich")

summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


summary(pred.a$predicted)
summary(pred.a$conf.low)
summary(pred.a$conf.high)

png(file = "lakeID.RandomEffect_Total_Phosphorus_Rich.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = TP.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.rich, colour = lakeID),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "α (q = 0)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot

TPPlot.A.q0.T <- ggpred_plot

######Hill.Shannon####
##building model
hist(TP.div.early$hill.Shannon)


time.model <- brm(formula = hill.Shannon ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.early,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .9))
#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digit = 4)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z[all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z[all]", "lakeID"), type = "re", ppd = FALSE)

time_model_ppcheck_Alpha_q1 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TP_time_Alpha_q1.png",
  plot = time_model_ppcheck_Alpha_q1,
  width = 3,
  height = 2,
  units = "in"
)

saveRDS(time.model, file = "time_model_Phos_shan")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_shan")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_shan")

time.model <- readRDS(file = "time_model_Phos_shan")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_shan")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_shan")

summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)



png(file = "lakeID.RandomEffect_Phosphorus_z_Hill.Shannon.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = TP.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.Shannon, colour = lakeID),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "α (q = 1)") + 
                  theme_bw()) +

  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot

TPPlot.A.q1.T <- ggpred_plot

######Hill.Inverse Simpson#####
## building model
hist(TP.div.early$hill.InverseSimpson)


time.model <- brm(formula = hill.InverseSimpson ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.early,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .85))
#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits= 4)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z [all]", "lakeID"), type = "re", ppd = FALSE)

time_model_ppcheck_Alpha_q2 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TP_time_Alpha_q2.png",
  plot = time_model_ppcheck_Alpha_q2,
  width = 3,
  height = 2,
  units = "in"
)

saveRDS(time.model, file = "time_model_Phos_simp")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_simp")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_simp")

time.model <- readRDS(file = "time_model_Phos_simp")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_simp")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_simp")

summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


png(file = "lakeID.RandomEffect_Phosphorus_z_Hill.InverseSimpson.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = TP.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.InverseSimpson, colour = lakeID),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  
                  labs(x = "Phosphorus (mg/l)", y = "α (q = 2)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot

TPPlot.A.q2.T <- ggpred_plot

######LCBD Lake q0####

n <- nrow(TP.div.early)
TP.div.early$lake.LCBD.q0 <- (TP.div.early$lake.LCBD.q0 * (n - 1) + 0.5) / n
TP.div.early$lake.LCBD.q1 <- (TP.div.early$lake.LCBD.q1 * (n - 1) + 0.5) / n
TP.div.early$lake.LCBD.q2 <- (TP.div.early$lake.LCBD.q2 * (n - 1) + 0.5) / n


##building model
time.model <- brm(formula = lake.LCBD.q0 ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.early,
                  family = Beta(link = "cloglog"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .8))

#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 7)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

time_model_ppcheck_LCBD_q0 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TP_time_LCBD_q0.png",
  plot = time_model_ppcheck_LCBD_q0,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_Phos_lake.LCBD.q0")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_lake.LCBD.q0")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_lake.LCBD.q0")

time.model <- readRDS(file = "time_model_Phos_lake.LCBD.q0")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_lake.LCBD.q0")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_lake.LCBD.q0")

summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)

png(file = "lakeID.RandomEffect_Total_Phosphorus_lake.LCBD.q0.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = TP.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = lake.LCBD.q0, colour = lakeID),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 0)") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot
dev.off()
ggpred_plot

TPPlot.B.q0.T <- ggpred_plot

######LCBD Lake q1####

##building model
time.model <- brm(formula = lake.LCBD.q1 ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.early,
                  family = Beta(link = "cloglog"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .93))

#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 7)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

time_model_ppcheck_LCBD_q1 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TP_time_LCBD_q1.png",
  plot = time_model_ppcheck_LCBD_q1,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_Phos_lake.LCBD.q1")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_lake.LCBD.q1")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_lake.LCBD.q1")

time.model <- readRDS(file = "time_model_Phos_lake.LCBD.q1")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_lake.LCBD.q1")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_lake.LCBD.q1")

summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)

png(file = "lakeID.RandomEffect_Total_Phosphorus_lake.LCBD.q1.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = TP.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = lake.LCBD.q1, colour = lakeID),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 1)") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot
dev.off()
ggpred_plot

TPPlot.B.q1.T <- ggpred_plot

######LCBD Lake q2####

##building model
time.model <- brm(formula = lake.LCBD.q2 ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.early,
                  family = Beta(link = "cloglog"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .9))

#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 7)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

time_model_ppcheck_LCBD_q2 <- pp_check(time.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TP_time_LCBD_q2.png",
  plot = time_model_ppcheck_LCBD_q2,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_Phos_lake.LCBD.q2")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_lake.LCBD.q2")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_lake.LCBD.q2")

time.model <- readRDS(file = "time_model_Phos_lake.LCBD.q2")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_lake.LCBD.q2")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_lake.LCBD.q2")

summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)

png(file = "lakeID.RandomEffect_Total_Phosphorus_lake.LCBD.q2.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = TP.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = lake.LCBD.q2, colour = lakeID),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 2)") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot
dev.off()
ggpred_plot

TPPlot.B.q2.T <- ggpred_plot

#####Phosphorus SPACE MODEL#####
######Richness######
##building model
TP.div.early$year <- as.factor(TP.div.early$year)
summary(TP.div.early)

space.model <- brm(formula = hill.rich ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.early,
                   family = Gamma(link = "log"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .85))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence

#model checks
pairs(space.model)
summary(space.model)
print(space.model, digits = 2)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_Alpha_q0 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TP_space_Alpha_q0.png",
  plot = space_model_ppcheck_Alpha_q0,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(space.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_rich")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_rich")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_rich")

space.model <- readRDS(file = "space_model_Phos_rich")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_rich")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_rich")

summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


png(file = "year.RandomEffect_Phosphorus_z_Rich.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = TP.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.rich, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "α (q = 0)") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot

dev.off()
ggpred_plot

TPPlot.A.q0.S <- ggpred_plot

######Hill.Shannon####
hist(TP.div.early$hill.Shannon)
qqPlot(TP.div.early$hill.Shannon)


##building model
space.model <- brm(formula = hill.Shannon ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.early,
                   family = Gamma(link = "log"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .8))
#model checks
pairs(space.model)
summary(space.model)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_Alpha_q1 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TP_space_Alpha_q1.png",
  plot = space_model_ppcheck_Alpha_q1,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(space.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_shan")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_shan")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_shan")

space.model <- readRDS(file = "space_model_Phos_shan")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_shan")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_shan")


summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)



summary(pred.a)

png(file = "year.RandomEffect_Phosphorus_Shannon.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = TP.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.Shannon, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "α (q = 1)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot

TPPlot.A.q1.S <- ggpred_plot


######Inverse Simpson####

space.model <- brm(formula = hill.InverseSimpson ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.early,
                   family = Gamma(link = "log"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .8))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence
#model checks
pairs(space.model)
summary(space.model)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_Alpha_q2 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TP_space_Alpha_q2.png",
  plot = space_model_ppcheck_Alpha_q2,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(space.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_simp")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_simp")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_simp")

space.model <- readRDS(file = "space_model_Phos_simp")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_simp")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_simp")

summary(pred.a)
summary(pred.b)


#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)

png(file = "year.RandomEffect_Phosphorus_InverseSimpson.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = TP.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.InverseSimpson, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  
                  labs(x = "Phosphorus (mg/l)", y = "α (q = 2)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot

TPPlot.A.q2.S <- ggpred_plot

######LCBD Year q0######
##building model
TP.div.early$year <- as.factor(TP.div.early$year)
summary(TP.div.early)

n <- nrow(TP.div.early)
TP.div.early$year.LCBD.q0 <- (TP.div.early$year.LCBD.q0 * (n - 1) + 0.5) / n
TP.div.early$year.LCBD.q1 <- (TP.div.early$year.LCBD.q1 * (n - 1) + 0.5) / n
TP.div.early$year.LCBD.q2 <- (TP.div.early$year.LCBD.q2 * (n - 1) + 0.5) / n

space.model <- brm(formula = year.LCBD.q0 ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.early,
                   family = Beta(link = "cloglog"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .9, max_treedepth = 12))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence

#model checks
pairs(space.model)
summary(space.model)
print(space.model, digits = 7)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_LCBD_q0 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TP_space_LCBD_q0.png",
  plot = space_model_ppcheck_LCBD_q2,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(space.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_year.LCBD.q0")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_year.LCBD.q0")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_year.LCBD.q0")

space.model <- readRDS(file = "space_model_Phos_year.LCBD.q0")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_year.LCBD.q0")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_year.LCBD.q0")

summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


png(file = "year.RandomEffect_Phosphorus_z_year.LCBD.q0.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = TP.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = year.LCBD.q0, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 0)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot

TPPlot.B.q0.S <- ggpred_plot

######LCBD Year q1######
##building model
TP.div.early$year <- as.factor(TP.div.early$year)
summary(TP.div.early)

space.model <- brm(formula = year.LCBD.q1 ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.early,
                   family = Beta(link = "cloglog"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .9, max_treedepth = 12))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence

#model checks
pairs(space.model)
summary(space.model)
print(space.model, digits = 7)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_LCBD_q1 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TP_space_LCBD_q1.png",
  plot = space_model_ppcheck_LCBD_q1,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(space.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_year.LCBD.q1")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_year.LCBD.q1")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_year.LCBD.q1")

space.model <- readRDS(file = "space_model_Phos_year.LCBD.q1")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_year.LCBD.q1")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_year.LCBD.q1")

summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


png(file = "year.RandomEffect_Phosphorus_z_year.LCBD.q1.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = TP.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = year.LCBD.q1, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 1)") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

TPPlot.B.q1.S <- ggpred_plot

######LCBD Year q2######
##building model
TP.div.early$year <- as.factor(TP.div.early$year)
summary(TP.div.early)

space.model <- brm(formula = year.LCBD.q2 ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.early,
                   family = Beta(link = "cloglog"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .9, max_treedepth = 12))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence

#model checks
pairs(space.model)
summary(space.model)
print(space.model, digits = 7)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

space_model_ppcheck_LCBD_q2 <- pp_check(space.model, ndraws = 100)
# Save
ggsave(
  filename = "PP-Checks/TP_space_LCBD_q2.png",
  plot = space_model_ppcheck_LCBD_q2,
  width = 3,
  height = 2,
  units = "in"
)

pred.a <- ggpredict(space.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_year.LCBD.q2")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_year.LCBD.q2")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_year.LCBD.q2")

space.model <- readRDS(file = "space_model_Phos_year.LCBD.q2")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_year.LCBD.q2")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_year.LCBD.q2")

summary(pred.a)
summary(pred.b)

#Changing pred.b so that plots are shortened to length of observation
raw_data <- attr(pred.b, "rawdata", exact = TRUE) #attr gets the specific attributes of an object
!is.null(raw_data)
ranges <- lapply(split(raw_data, raw_data$group), function(i) round(range(i$x, na.rm = TRUE),3))#the 3 is the number of decimal places to round to
for (i in names(ranges)) {
  remove <- pred.b$group == i & pred.b$x < ranges[[i]][1]
  pred.b$x[remove] <- NA
  remove <- pred.b$group == i & pred.b$x > ranges[[i]][2]
  pred.b$x[remove] <- NA
}

pred.a$x <- pred.a$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


png(file = "year.RandomEffect_Phosphorus_z_year.LCBD.q2.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.2) +  # error band
                  geom_point(data = TP.div.early,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = year.LCBD.q2, colour = year),
                             alpha = 0.1,
                             shape = 16) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),
                            alpha = 0.6) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 2)") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot


TPPlot.B.q2.S <- ggpred_plot
 
#####Comparing models####
######temperature Rich####
time.model <- readRDS(file = "time_model_temp_rich.rds")
space.model <- readRDS(file = "space_model_temp_rich")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["WaterTemp_z", 1])
time.error <- as.numeric(fixed_effects.time["WaterTemp_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["WaterTemp_z", 1])
space.error <- as.numeric(fixed_effects.space["WaterTemp_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value


######temperature Shannon####
time.model <- readRDS(file = "time_model_temp_shan")
space.model <- readRDS(file = "space_model_temp_shan.rds")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["WaterTemp_z", 1])
time.error <- as.numeric(fixed_effects.time["WaterTemp_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["WaterTemp_z", 1])
space.error <- as.numeric(fixed_effects.space["WaterTemp_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value



######temperature Inverse Simpson####
time.model <- readRDS(file = "time_model_temp_simp")
space.model <- readRDS(file = "space_model_temp_simp")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["WaterTemp_z", 1])
time.error <- as.numeric(fixed_effects.time["WaterTemp_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["WaterTemp_z", 1])
space.error <- as.numeric(fixed_effects.space["WaterTemp_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value


######temperature LCBD q = 0####
time.model <- readRDS(file = "time_model_temp_LCBD.L.q0")
space.model <- readRDS(file = "space_model_temp_year.LCBD.q0")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["WaterTemp_z", 1])
time.error <- as.numeric(fixed_effects.time["WaterTemp_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["WaterTemp_z", 1])
space.error <- as.numeric(fixed_effects.space["WaterTemp_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value


######temperature LCBD q = 1####
time.model <- readRDS(file = "time_model_temp_LCBD.L.q1")
space.model <- readRDS(file = "space_model_temp_year.LCBD.q1")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["WaterTemp_z", 1])
time.error <- as.numeric(fixed_effects.time["WaterTemp_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["WaterTemp_z", 1])
space.error <- as.numeric(fixed_effects.space["WaterTemp_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value






######temperature LCBD q = 2####
time.model <- readRDS(file = "time_model_temp_LCBD.L.q2")
space.model <- readRDS(file = "space_model_temp_year.LCBD.q2")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["WaterTemp_z", 1])
time.error <- as.numeric(fixed_effects.time["WaterTemp_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["WaterTemp_z", 1])
space.error <- as.numeric(fixed_effects.space["WaterTemp_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value


######Total Nitrogen Rich####
time.model <- readRDS(file = "time_model_TN_rich")
space.model <- readRDS(file = "space_model_TN_rich")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["TN_z", 1])
time.error <- as.numeric(fixed_effects.time["TN_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["TN_z", 1])
space.error <- as.numeric(fixed_effects.space["TN_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value


######Total Nitrogen Shannon####
time.model <- readRDS(file = "time_model_TN_shan")
space.model <- readRDS(file = "space_model_TN_shan")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["TN_z", 1])
time.error <- as.numeric(fixed_effects.time["TN_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["TN_z", 1])
space.error <- as.numeric(fixed_effects.space["TN_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value


######Total Nitrogen Inverse Simpson####
time.model <- readRDS(file = "time_model_TN_simp")
space.model <- readRDS(file = "space_model_TN_simp")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["TN_z", 1])
time.error <- as.numeric(fixed_effects.time["TN_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["TN_z", 1])
space.error <- as.numeric(fixed_effects.space["TN_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value


######Total Nitrogen LCBD q = 0####
time.model <- readRDS(file = "time_model_TN_lake.LCBD.q0")
space.model <- readRDS(file = "space_model_TN_year.LCBD.q0")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["TN_z", 1])
time.error <- as.numeric(fixed_effects.time["TN_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["TN_z", 1])
space.error <- as.numeric(fixed_effects.space["TN_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value


######Total Nitrogen LCBD q = 1####
time.model <- readRDS(file = "time_model_TN_lake.LCBD.q1")
space.model <- readRDS(file = "space_model_TN_year.LCBD.q1")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["TN_z", 1])
time.error <- as.numeric(fixed_effects.time["TN_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["TN_z", 1])
space.error <- as.numeric(fixed_effects.space["TN_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value


######Total Nitrogen LCBD q = 2####
time.model <- readRDS(file = "time_model_TN_lake.LCBD.q2")
space.model <- readRDS(file = "space_model_TN_year.LCBD.q2")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["TN_z", 1])
time.error <- as.numeric(fixed_effects.time["TN_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["TN_z", 1])
space.error <- as.numeric(fixed_effects.space["TN_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value


######Total Phosphorus Rich####
time.model <- readRDS(file = "time_model_Phos_rich")
space.model <- readRDS(file = "space_model_Phos_rich")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["Phosphorus_z", 1])
time.error <- as.numeric(fixed_effects.time["Phosphorus_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["Phosphorus_z", 1])
space.error <- as.numeric(fixed_effects.space["Phosphorus_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value


######Total Phosphorus Shannon####
time.model <- readRDS(file = "time_model_Phos_shan")
space.model <- readRDS(file = "space_model_Phos_shan")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["Phosphorus_z", 1])
time.error <- as.numeric(fixed_effects.time["Phosphorus_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["Phosphorus_z", 1])
space.error <- as.numeric(fixed_effects.space["Phosphorus_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value

######Total Phosphorus Inverse Simpson####
time.model <- readRDS(file = "time_model_Phos_simp")
space.model <- readRDS(file = "space_model_Phos_simp")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["Phosphorus_z", 1])
time.error <- as.numeric(fixed_effects.time["Phosphorus_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["Phosphorus_z", 1])
space.error <- as.numeric(fixed_effects.space["Phosphorus_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value


######Total Phosphorus LCBD q = 0 ####
time.model <- readRDS(file = "time_model_Phos_lake.LCBD.q0")
space.model <- readRDS(file = "space_model_Phos_year.LCBD.q0")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["Phosphorus_z", 1])
time.error <- as.numeric(fixed_effects.time["Phosphorus_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["Phosphorus_z", 1])
space.error <- as.numeric(fixed_effects.space["Phosphorus_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value


######Total Phosphorus LCBD q = 1 ####
time.model <- readRDS(file = "time_model_Phos_lake.LCBD.q1")
space.model <- readRDS(file = "space_model_Phos_year.LCBD.q1")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["Phosphorus_z", 1])
time.error <- as.numeric(fixed_effects.time["Phosphorus_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["Phosphorus_z", 1])
space.error <- as.numeric(fixed_effects.space["Phosphorus_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value

#
######Total Phosphorus LCBD q = 2 ####
time.model <- readRDS(file = "time_model_Phos_lake.LCBD.q2")
space.model <- readRDS(file = "space_model_Phos_year.LCBD.q2")

summary(time.model)
summary(space.model)

# Extract population-level (fixed) effects
fixed_effects.time <- fixef(time.model)

#extract population-level slope and error
time.slope <- as.numeric(fixed_effects.time["Phosphorus_z", 1])
time.error <- as.numeric(fixed_effects.time["Phosphorus_z", 2])

# Extract population-level (fixed) effects
fixed_effects.space <- fixef(space.model)

# View the population-level slope and error
space.slope <- as.numeric(fixed_effects.space["Phosphorus_z", 1])
space.error <- as.numeric(fixed_effects.space["Phosphorus_z", 2])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
z

p_value <- 2 * (1 - pnorm(abs(z)))
p_value

#########Total Phosphorus LCBD q = 2 Random############
ranef.time <- as.data.frame(ranef(time.model)$lakeID) %>%
  rownames_to_column(var = "lakeID")

# Subset data for the given lake ID
subset_data <- ranef.time[, -(2:5)]

# View the population-level slope and error
time.slope <- as.numeric(subset_data[, 2])
time.error <- as.numeric(subset_data[, 3])

z <- ((time.slope - space.slope) / (sqrt((time.error^2)+(space.error^2))))
p_value <- 2 * (1 - pnorm(abs(z)))

# Combine Z and P with site identifiers
ranef.ZPs <- data.frame(
  site = ranef.time$lakeID,
  z = z,
  p_value = p_value
)

ranef.ZPs$adjusted.p <- p.adjust(ranef.ZPs$p_value, method = "fdr")

sum(ranef.ZPs$adjusted.p < 0.05, na.rm = TRUE)


######Directionality#####
###Time Model

# Extract random effects
ranef_values <- ranef(time.model)  # Returns a list

# Check structure
str(ranef_values)

# Extract the random slopes for 'Phosphorus_z' by lakeID
random_slopes <- ranef_values$lakeID[, ,"Phosphorus_z" ]

# Convert to a data frame with credible intervals
random_slopes_summary <- data.frame(
  Estimate = random_slopes[, "Estimate"],
  CI_Lower = random_slopes[, "Q2.5"],
  CI_Upper = random_slopes[, "Q97.5"],
  Classification = ifelse(random_slopes[, "Q2.5"] > 0 & random_slopes[, "Q97.5"] > 0, "Positive",
                          ifelse(random_slopes[, "Q2.5"] < 0 & random_slopes[, "Q97.5"] < 0, "Negative",
                                 "Zero Effect"))
)

# Define all possible categories
categories <- c("Negative", "Zero Effect","Positive")

# Create the table, ensuring all categories are included
randefftable <- (table(factor(random_slopes_summary$Classification, levels = categories)))

# Print the result
print(randefftable)

randefftable_df <- as.data.frame(randefftable)
colnames(randefftable_df) <- c("Classification", "Count")


png(file = "BOX.time.RandomEffect_Phosphorus_LCBD.q2.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Phosphorus-LCBD (q = 2) Time Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot

###Space Model
# Extract random effects
ranef_values <- ranef(space.model)  # Returns a list

# Check structure
str(ranef_values)

# Extract the random slopes for 'Phosphorus_z' by year
random_slopes <- ranef_values$year[, ,"Phosphorus_z"]

# Convert to a data frame with credible intervals
random_slopes_summary <- data.frame(
  Estimate = random_slopes[, "Estimate"],
  CI_Lower = random_slopes[, "Q2.5"],
  CI_Upper = random_slopes[, "Q97.5"],
  Classification = ifelse(random_slopes[, "Q2.5"] > 0 & random_slopes[, "Q97.5"] > 0, "Positive",
                          ifelse(random_slopes[, "Q2.5"] < 0 & random_slopes[, "Q97.5"] < 0, "Negative",
                                 "Zero Effect"))
)

# Define all possible categories
categories <- c("Negative", "Zero Effect","Positive")

# Create the table, ensuring all categories are included
randefftable <- (table(factor(random_slopes_summary$Classification, levels = categories)))

# Print the result
print(randefftable)

randefftable_df <- as.data.frame(randefftable)
colnames(randefftable_df) <- c("Classification", "Count")


png(file = "BOX.space.RandomEffect_Phosphorus_LCBD.q2.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Phosphorus-LCBD (q = 2) Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot




#######Formating Graphs##########
library(patchwork)

shared_xlim <- range(Temp.div.early$Temp.water.degC, na.rm = TRUE)
TN.shared_xlim <- range(Nit.div.early$Total_Nitrogen, na.rm = TRUE)
TP.shared_xlim <- range(TP.div.early$Total_Phorphorus_Mixed_Forms_mg.l.as.P., na.rm = TRUE)

theme_no_x <- theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  plot.margin = unit(c(0.5, 1, 0.5, 1.5), "cm")  # top, right, bottom, left
)

theme_no_x_or_y <- theme(
  axis.title.x = element_blank(),
  axis.text.x = element_blank(),
  axis.ticks.x = element_blank(),
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  plot.margin = unit(c(0.5, 1, 0.5, 1.5), "cm")  # top, right, bottom, left
)

theme_with_x <- theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 12),
  axis.ticks.x = element_line(),
  plot.margin = unit(c(0.5, 1, 0.5, 1.5), "cm")
)

theme_with_x_but_no_y <- theme(
  axis.title.x = element_text(size = 14),
  axis.text.x = element_text(size = 12),
  axis.ticks.x = element_line(),  
  axis.title.y = element_blank(),
  axis.text.y = element_blank(),
  axis.ticks.y = element_blank(),
  plot.margin = unit(c(0.5, 1, 0.5, 1.5), "cm")
)

#Temperature Change plot
TempPlot.A.q0 <- TempPlot.A.q0 + coord_cartesian(xlim = shared_xlim) + theme_no_x 
TempPlot.A.q1 <- TempPlot.A.q1 + coord_cartesian(xlim = shared_xlim) + theme_no_x
TempPlot.A.q2 <- TempPlot.A.q2 + coord_cartesian(xlim = shared_xlim) + theme_no_x
TempPlot.B.q0 <- TempPlot.B.q0 + coord_cartesian(xlim = shared_xlim) + theme_no_x
TempPlot.B.q1 <- TempPlot.B.q1 + coord_cartesian(xlim = shared_xlim) + theme_no_x
TempPlot.B.q2 <- TempPlot.B.q2 + coord_cartesian(xlim = shared_xlim) + theme_with_x
TempPlot.A.q0.S <- TempPlot.A.q0.S + coord_cartesian(xlim = shared_xlim) + theme_no_x_or_y
TempPlot.A.q1.S <- TempPlot.A.q1.S + coord_cartesian(xlim = shared_xlim) + theme_no_x_or_y
TempPlot.A.q2.S <- TempPlot.A.q2.S + coord_cartesian(xlim = shared_xlim) + theme_no_x_or_y
TempPlot.B.q0.S <- TempPlot.B.q0.S + coord_cartesian(xlim = shared_xlim) + theme_no_x_or_y
TempPlot.B.q1.S <- TempPlot.B.q1.S + coord_cartesian(xlim = shared_xlim) + theme_no_x_or_y
TempPlot.B.q2.S <- TempPlot.B.q2.S + coord_cartesian(xlim = shared_xlim) + theme_with_x_but_no_y

#Land use change
TNPlot.A.q0.T <- TNPlot.A.q0.T + coord_cartesian(xlim = TN.shared_xlim) + theme_no_x 
TNPlot.A.q1.T <- TNPlot.A.q1.T + coord_cartesian(xlim = TN.shared_xlim) + theme_no_x
TNPlot.A.q2.T <- TNPlot.A.q2.T + coord_cartesian(xlim = TN.shared_xlim) + theme_no_x
TNPlot.B.q0.T <- TNPlot.B.q0.T + coord_cartesian(xlim = TN.shared_xlim) + theme_no_x
TNPlot.B.q1.T <- TNPlot.B.q1.T + coord_cartesian(xlim = TN.shared_xlim) + theme_no_x
TNPlot.B.q2.T <- TNPlot.B.q2.T + coord_cartesian(xlim = TN.shared_xlim) + theme_with_x
TNPlot.A.q0.S <- TNPlot.A.q0.S + coord_cartesian(xlim = TN.shared_xlim) + theme_no_x_or_y
TNPlot.A.q1.S <- TNPlot.A.q1.S + coord_cartesian(xlim = TN.shared_xlim) + theme_no_x_or_y
TNPlot.A.q2.S <- TNPlot.A.q2.S + coord_cartesian(xlim = TN.shared_xlim) + theme_no_x_or_y
TNPlot.B.q0.S <- TNPlot.B.q0.S + coord_cartesian(xlim = TN.shared_xlim) + theme_no_x_or_y
TNPlot.B.q1.S <- TNPlot.B.q1.S + coord_cartesian(xlim = TN.shared_xlim) + theme_no_x_or_y
TNPlot.B.q2.S <- TNPlot.B.q2.S + coord_cartesian(xlim = TN.shared_xlim) + theme_with_x_but_no_y

TPPlot.A.q0.T <- TPPlot.A.q0.T + coord_cartesian(xlim = TP.shared_xlim) + theme_no_x_or_y 
TPPlot.A.q1.T <- TPPlot.A.q1.T + coord_cartesian(xlim = TP.shared_xlim) + theme_no_x_or_y
TPPlot.A.q2.T <- TPPlot.A.q2.T + coord_cartesian(xlim = TP.shared_xlim) + theme_no_x_or_y
TPPlot.B.q0.T <- TPPlot.B.q0.T + coord_cartesian(xlim = TP.shared_xlim) + theme_no_x_or_y
TPPlot.B.q1.T <- TPPlot.B.q1.T + coord_cartesian(xlim = TP.shared_xlim) + theme_no_x_or_y
TPPlot.B.q2.T <- TPPlot.B.q2.T + coord_cartesian(xlim = TP.shared_xlim) + theme_with_x_but_no_y
TPPlot.A.q0.S <- TPPlot.A.q0.S + coord_cartesian(xlim = TP.shared_xlim) + theme_no_x_or_y
TPPlot.A.q1.S <- TPPlot.A.q1.S + coord_cartesian(xlim = TP.shared_xlim) + theme_no_x_or_y
TPPlot.A.q2.S <- TPPlot.A.q2.S + coord_cartesian(xlim = TP.shared_xlim) + theme_no_x_or_y
TPPlot.B.q0.S <- TPPlot.B.q0.S + coord_cartesian(xlim = TP.shared_xlim) + theme_no_x_or_y
TPPlot.B.q1.S <- TPPlot.B.q1.S + coord_cartesian(xlim = TP.shared_xlim) + theme_no_x_or_y
TPPlot.B.q2.S <- TPPlot.B.q2.S + coord_cartesian(xlim = TP.shared_xlim) + theme_with_x_but_no_y


combined_plot <- (
  TempPlot.A.q0 /
    TempPlot.A.q0.S /
    TempPlot.A.q1 /
    TempPlot.A.q1.S /
    TempPlot.A.q2 /
    TempPlot.A.q2.S /
    TempPlot.B.q0 /
    TempPlot.B.q0.S /
    TempPlot.B.q1 /
    TempPlot.B.q1.S /
    TempPlot.B.q2 /
    TempPlot.B.q2.S
) + plot_layout(ncol = 2)

combined_plot


Landchange_plot <- (
  TNPlot.A.q0.T / 
    TNPlot.A.q0.S /
    TPPlot.A.q0.T/
    TPPlot.A.q0.S /
    
    TNPlot.A.q1.T / 
    TNPlot.A.q1.S /
    TPPlot.A.q1.T / 
    TPPlot.A.q1.S /
    
    TNPlot.A.q2.T /
    TNPlot.A.q2.S /
    TPPlot.A.q2.T / 
    TPPlot.A.q2.S /
    
    TNPlot.B.q0.T /
    TNPlot.B.q0.S / 
    TPPlot.B.q0.T /
    TPPlot.B.q0.S /
    
    TNPlot.B.q1.T /
    TNPlot.B.q1.S /
    TPPlot.B.q1.T /
    TPPlot.B.q1.S /
    
    TNPlot.B.q2.T /
    TNPlot.B.q2.S /
    TPPlot.B.q2.T /
    TPPlot.B.q2.S
) + plot_layout(ncol = 4)

Landchange_plot

ggsave(
  filename = "Landchange_plot.pdf",   # output file name
  plot = Landchange_plot,             # your plot object
  device = "pdf",                     # file type
  width = 12,                         # width in inches
  height = 10,                        # height in inches
  units = "in",                       # specify unit type
  dpi = 300                           # resolution (useful if you convert later)
)

getwd()
