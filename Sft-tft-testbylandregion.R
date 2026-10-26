#####Code from Large Vs small lakes which is being modified for comparing land forms
#######Ambient Lake Data Mixed Effect Models incorporating lake size#####
#R version 4.2.2 
#29 August 2025

#switch to a baysian approach

###Clear working environment
rm(list=ls())


##Reading in Data 
setwd(dirname(rstudioapi::getSourceEditorContext()$path)) #set working directory
getwd()
WQ.zoo <- read.csv("WQ.taxa25.csv", header = TRUE)
landforms <- read.csv("Copy.of.lake.coords.landformregions.csv", header = TRUE)





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

summary(lakes.q2)
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
lcbd.columns <- c("lake.LCBD.q0","lake.LCBD.q1", "lake.LCBD.q2","year.LCBD.q0", "year.LCBD.q1", "year.LCBD.q2")
LD.df[lcbd.columns] <- lapply(LD.df[lcbd.columns], function(x) ifelse(x < 0, 0, x))


#Adding WQ diversity metrics to WQ df
WQ.div <- left_join(WQ, LD.df, by = "lake.date")

#####LOAD data saved from directory####
saveRDS(WQ.div, file = "WQ.div.rds")

WQ.div <- readRDS(file = "WQ.div.rds")


####Processing landform data####
landforms[landforms$Unique.Lakes == "74", "LANDFORM_R"] <- "Des Moines Lobe"
landforms$lakeID <- landforms$Unique.Lakes
landforms <- select(landforms, lakeID, LANDFORM_R) 
landforms$lakeID <- as.character(landforms$lakeID)
landforms$LANDFORM_R <- as.factor(landforms$LANDFORM_R)
landforms <- landforms %>%
  mutate(SouthernIowaDrift = if_else(LANDFORM_R == "Southern Iowa Drift Plain", 1, 0))


summary(landforms) 

WQ.div <- WQ.div %>%
  left_join(landforms, by = "lakeID")

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

#Filtering by lakeID to remove lakes which were not sampled as frequently
Temp.div.early <- Temp.div.early %>%
  left_join(lakeID.df, by = "lakeID") %>%
  filter(Freq >= 9) #Freq >= is selecting how many observations each lake needed to be included  

#Calculating Z
Temp.div.early$WaterTemp_z <- scale(Temp.div.early$Temp.water.degC, center = TRUE, scale = TRUE)

#####FILTER FOR ON or off the Southern Iowa drift plain
Temp.div.L <- Temp.div.early %>%
  filter(SouthernIowaDrift == "0")


Temp.div.S <- Temp.div.early %>%
  filter(SouthernIowaDrift == "1")

##### Not on SID LAKES #####
######Richness####


##building model
time.model <- brm(formula = hill.rich ~ WaterTemp_z + (1 + WaterTemp_z | lakeID),
                  data = Temp.div.L,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .99))


#Model Checks
pairs(time.model)
summary(time.model)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)



pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

summary(pred.a)
summary(pred.b)

saveRDS(time.model, file = "time_model_temp_rich.L")
saveRDS(pred.a, file = "pred.A_timemodel_temp_rich.L")
saveRDS(pred.b, file = "pred.B_timemodel_temp_rich.L")

time.model <- readRDS(file = "time_model_temp_rich.L")
pred.a <- readRDS(file = "pred.A_timemodel_temp_rich.L")
pred.b <- readRDS(file = "pred.B_timemodel_temp_rich.L")



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

pred.a$x2 <- pred.a$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)
pred.b$x2 <- pred.b$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)

shared_xlim <- range(Temp.div.L$Temp.water.degC, na.rm = TRUE)
fixed_margins <- theme(plot.margin = unit(c(0.5, 1, 0.5, 1.5), "cm"))

pdf(file = "lakeID.RandomEffect_WaterTemp_Rich.L.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x2), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.L,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = hill.rich, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x2, predicted, colour = group)) +      
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
                  data = Temp.div.L,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .94))
#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 5)

pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_temp_shan.L")
saveRDS(pred.a, file = "pred.A_timemodel_temp_shan.L")
saveRDS(pred.b, file = "pred.B_timemodel_temp_shan.L")

time.model <- readRDS(file = "time_model_temp_shan.L")
pred.a <- readRDS(file = "pred.A_timemodel_temp_shan.L")
pred.b <- readRDS(file = "pred.B_timemodel_temp_shan.L")

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

pred.a$x <- pred.a$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)


pdf(file = "lakeID.RandomEffect_WaterTemp_Hill.Shannon.L.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.L,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = hill.Shannon, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
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
                  data = Temp.div.L,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .95))
#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 4)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

#saving R objects
saveRDS(time.model, file = "time_model_temp_simp.L")
saveRDS(pred.a, file = "pred.A_timemodel_temp_simp.L")
saveRDS(pred.b, file = "pred.B_timemodel_temp_simp.L")

time.model <- readRDS(file = "time_model_temp_simp.L")
pred.a <- readRDS(file = "pred.A_timemodel_temp_simp.L")
pred.b <- readRDS(file = "pred.B_timemodel_temp_simp.L")


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

pred.a$x <- pred.a$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)




pdf(file = "lakeID.RandomEffect_WaterTemp_Hill.InverseSimpson.L.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.L,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = hill.InverseSimpson, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
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
hist(Temp.div.L$lake.LCBD.q0)
qqPlot(Temp.div.L$lake.LCBD.q0)

summary(Temp.div.L)

#standard “Smithson & Verkuilen” (2006) transformation so that we can use zero inflated beta with zeros and 1s
n <- nrow(Temp.div.L)
Temp.div.L$lake.LCBD.q0 <- (Temp.div.L$lake.LCBD.q0 * (n - 1) + 0.5) / n

time.model <- brm(formula = lake.LCBD.q0 ~ WaterTemp_z + (1 + WaterTemp_z | lakeID),
                  data = Temp.div.L,
                  family = Beta(link = "cloglog"),
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

pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

#saving R objects
saveRDS(time.model, file = "time_model_temp_LCBD.q0.L")
saveRDS(pred.a, file = "pred.A_timemodel_temp_LCBD.q0.L")
saveRDS(pred.b, file = "pred.B_timemodel_temp_LCBD.q0.L")

time.model <- readRDS(file = "time_model_temp_LCBD.q0.L")
pred.a <- readRDS(file = "pred.A_timemodel_temp_LCBD.q0.L")
pred.b <- readRDS(file = "pred.B_timemodel_temp_LCBD.q0.L")


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

pred.a$x <- pred.a$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)




pdf(file = "lakeID.RandomEffect_WaterTemp_LCBD.q0.L.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.L,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = lake.LCBD.q0, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
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
  guides(color = "none")+
  fixed_margins
ggpred_plot
dev.off()
ggpred_plot
TempPlot.B.q0 <- ggpred_plot




######LCBD Lake q1#####
##building model
hist(Temp.div.L$lake.LCBD.q1)
qqPlot(Temp.div.L$lake.LCBD.q1)
summary(Temp.div.L)

#standard “Smithson & Verkuilen” (2006) transformation so that we can use zero inflated beta with zeros and 1s
n <- nrow(Temp.div.L)
Temp.div.L$lake.LCBD.q1 <- (Temp.div.L$lake.LCBD.q1 * (n - 1) + 0.5) / n

time.model <- brm(formula = lake.LCBD.q1 ~ WaterTemp_z + (1 + WaterTemp_z | lakeID),
                  data = Temp.div.L,
                  family = Beta(link = "cloglog"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .8)) #Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 4)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

#saving R objects
saveRDS(time.model, file = "time_model_temp_LCBD.q1.L")
saveRDS(pred.a, file = "pred.A_timemodel_temp_LCBD.q1.L")
saveRDS(pred.b, file = "pred.B_timemodel_temp_LCBD.q1.L")

time.model <- readRDS(file = "time_model_temp_LCBD.q1.L")
pred.a <- readRDS(file = "pred.A_timemodel_temp_LCBD.q1.L")
pred.b <- readRDS(file = "pred.B_timemodel_temp_LCBD.q1.L")


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

pred.a$x <- pred.a$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)




pdf(file = "lakeID.RandomEffect_WaterTemp_LCBD.q1.L.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.L,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = lake.LCBD.q1, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
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
  guides(color = "none")+
  fixed_margins
ggpred_plot
dev.off()
ggpred_plot
TempPlot.B.q1 <- ggpred_plot

######LCBD Lake q2#####
##building model
hist(Temp.div.L$lake.LCBD.q2)
qqPlot(Temp.div.L$lake.LCBD.q2)


summary(Temp.div.L)

#standard “Smithson & Verkuilen” (2006) transformation so that we can use zero inflated beta with zeros and 1s
n <- nrow(Temp.div.L)
Temp.div.L$lake.LCBD.q2 <- (Temp.div.L$lake.LCBD.q2 * (n - 1) + 0.5) / n


time.model <- brm(formula = lake.LCBD.q2 ~ WaterTemp_z + (1 + WaterTemp_z | lakeID),
                  data = Temp.div.L,
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

pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

#saving R objects
saveRDS(time.model, file = "time_model_temp_LCBD.q2.L")
saveRDS(pred.a, file = "pred.A_timemodel_temp_LCBD.q2.L")
saveRDS(pred.b, file = "pred.B_timemodel_temp_LCBD.q2.L")

time.model <- readRDS(file = "time_model_temp_LCBD.q2.L")
pred.a <- readRDS(file = "pred.A_timemodel_temp_LCBD.q2.L")
pred.b <- readRDS(file = "pred.B_timemodel_temp_LCBD.q2.L")


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

pred.a$x <- pred.a$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)




pdf(file = "lakeID.RandomEffect_WaterTemp_LCBD.q2.L.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.L,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = lake.LCBD.q2, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Water Temperature (C°)", y = "LCBD (q = 2)") + 
                  theme_bw()) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),  # Y-axis label size
    axis.text = element_text(size = 12)  # Axis tick labels size 
  ) +
  guides(color = "none")+
  fixed_margins
ggpred_plot
dev.off()
ggpred_plot
TempPlot.B.q2 <- ggpred_plot

#####TEMP SPACE MODEL#####
######Richness##########
##building model
Temp.div.L$year <- as.factor(Temp.div.L$year)
summary(Temp.div.L)

space.model <- brm(formula = hill.rich ~ WaterTemp_z + (1 + WaterTemp_z | year),
                   data = Temp.div.L,
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


pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_rich.L")
saveRDS(pred.a, file = "pred.a_spacemodel_temp_rich.L")
saveRDS(pred.b, file = "pred.b_spacemodel_temp_rich.L")

space.model <- readRDS(file = "space_model_temp_rich.L")
pred.a <- readRDS(file = "pred.a_spacemodel_temp_rich.L")
pred.b <- readRDS(file = "pred.b_spacemodel_temp_rich.L")

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

pred.a$x <- pred.a$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_Rich.L.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.L,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = hill.rich, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
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
hist(Temp.div.L$hill.Shannon)
qqPlot(Temp.div.L$hill.Shannon)


##building model
space.model <- brm(formula = hill.Shannon ~ WaterTemp_z + (1 + WaterTemp_z | year),
                   data = Temp.div.L,
                   family = Gamma(link = "log"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .85))
#model checks
summary(space.model)
print(space.model, digits = 4)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)


pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_shan.L")
saveRDS(pred.a, file = "pred.A_spacemodel_temp_shan.L")
saveRDS(pred.b, file = "pred.B_spacemodel_temp_shan.L")

space.model <- readRDS(file = "space_model_temp_shan.L")
pred.a <- readRDS(file = "pred.A_spacemodel_temp_shan.L")
pred.b <- readRDS(file = "pred.B_spacemodel_temp_shan.L")

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

pred.a$x <- pred.a$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_Shannon.L.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.L,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = hill.Shannon, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
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
                   data = Temp.div.L,
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

pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_simp.L")
saveRDS(pred.a, file = "pred.A_spacemodel_temp_simp.L")
saveRDS(pred.b, file = "pred.B_spacemodel_temp_simp.L")

space.model <- readRDS(file = "space_model_temp_simp.L")
pred.a <- readRDS(file = "pred.A_spacemodel_temp_simp.L")
pred.b <- readRDS(file = "pred.B_spacemodel_temp_simp.L")

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

pred.a$x <- pred.a$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_Inverse Simpson.L.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.L,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = hill.InverseSimpson, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
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
Temp.div.L$year <- as.factor(Temp.div.L$year)
summary(Temp.div.L)

hist(Temp.div.L$year.LCBD.q0)
qqPlot(Temp.div.L$year.LCBD.q0)

Temp.div.L$year.LCBD.q0.t <- (Temp.div.L$year.LCBD.q0 +1)
hist(Temp.div.L$year.LCBD.q0.t)
qqPlot(Temp.div.L$year.LCBD.q0.t)

#standard “Smithson & Verkuilen” (2006) transformation so that we can use zero inflated beta with zeros and 1s
n <- nrow(Temp.div.L)
Temp.div.L$year.LCBD.q0 <- (Temp.div.L$year.LCBD.q0 * (n - 1) + 0.5) / n
Temp.div.L$year.LCBD.q1 <- (Temp.div.L$year.LCBD.q1 * (n - 1) + 0.5) / n
Temp.div.L$year.LCBD.q2 <- (Temp.div.L$year.LCBD.q2 * (n - 1) + 0.5) / n


space.model <- brm(formula = year.LCBD.q0 ~ WaterTemp_z + (1 + WaterTemp_z | year),
                   data = Temp.div.L,
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


pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_year.LCBD.q0.L")
saveRDS(pred.a, file = "pred.a_spacemodel_temp_year.LCBD.q0.L")
saveRDS(pred.b, file = "pred.b_spacemodel_temp_year.LCBD.q0.L")

space.model <- readRDS(file = "space_model_temp_year.LCBD.q0.L")
pred.a <- readRDS(file = "pred.a_spacemodel_temp_year.LCBD.q0.L")
pred.b <- readRDS(file = "pred.b_spacemodel_temp_year.LCBD.q0.L")

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

pred.a$x <- pred.a$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_year.LCBD.q0.L.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.L,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = year.LCBD.q0, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Water Temperature (C)", y = "LCBD (q = 0)") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot

dev.off()
ggpred_plot
TempPlot.B.q0.S <- ggpred_plot


######LCBD Year q1##########
##building model
Temp.div.L$year <- as.factor(Temp.div.L$year)
summary(Temp.div.L)

hist(Temp.div.L$year.LCBD.q1)
qqPlot(Temp.div.L$year.LCBD.q1)

space.model <- brm(formula = year.LCBD.q1 ~ WaterTemp_z + (1 + WaterTemp_z | year),
                   data = Temp.div.L,
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


pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_year.LCBD.q1.L")
saveRDS(pred.a, file = "pred.a_spacemodel_temp_year.LCBD.q1.L")
saveRDS(pred.b, file = "pred.b_spacemodel_temp_year.LCBD.q1.L")

space.model <- readRDS(file = "space_model_temp_year.LCBD.q1.L")
pred.a <- readRDS(file = "pred.a_spacemodel_temp_year.LCBD.q1.L")
pred.b <- readRDS(file = "pred.b_spacemodel_temp_year.LCBD.q1.L")

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

pred.a$x <- pred.a$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_year.LCBD.q1.L.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.L,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = year.LCBD.q1, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
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
Temp.div.L$year <- as.factor(Temp.div.L$year)
summary(Temp.div.L)

hist(Temp.div.L$year.LCBD.q2)
qqPlot(Temp.div.L$year.LCBD.q2)

space.model <- brm(formula = year.LCBD.q2 ~ WaterTemp_z + (1 + WaterTemp_z | year),
                   data = Temp.div.L,
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


pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_year.LCBD.q2.L")
saveRDS(pred.a, file = "pred.a_spacemodel_temp_year.LCBD.q2.L")
saveRDS(pred.b, file = "pred.b_spacemodel_temp_year.LCBD.q2.L")

space.model <- readRDS(file = "space_model_temp_year.LCBD.q2.L")
pred.a <- readRDS(file = "pred.a_spacemodel_temp_year.LCBD.q2.L")
pred.b <- readRDS(file = "pred.b_spacemodel_temp_year.LCBD.q2.L")

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

pred.a$x <- pred.a$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.L$Temp.water.degC) + mean(Temp.div.L$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_year.LCBD.q2.L.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.L,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = year.LCBD.q2, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Water Temperature (C)", y = "LCBD (q = 2)") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot

dev.off()
ggpred_plot
TempPlot.B.q2.S <- ggpred_plot


##### On SID LAKES #####
######Richness####


##building model
time.model <- brm(formula = hill.rich ~ WaterTemp_z + (1 + WaterTemp_z | lakeID),
                  data = Temp.div.S,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .99))


#Model Checks
pairs(time.model)
summary(time.model)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)



pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

summary(pred.a)
summary(pred.b)

saveRDS(time.model, file = "time_model_temp_rich.S")
saveRDS(pred.a, file = "pred.A_timemodel_temp_rich.S")
saveRDS(pred.b, file = "pred.B_timemodel_temp_rich.S")

time.model <- readRDS(file = "time_model_temp_rich.S")
pred.a <- readRDS(file = "pred.A_timemodel_temp_rich.S")
pred.b <- readRDS(file = "pred.B_timemodel_temp_rich.S")



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

pred.a$x2 <- pred.a$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)
pred.b$x2 <- pred.b$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)

shared_xlim <- range(Temp.div.S$Temp.water.degC, na.rm = TRUE)
fixed_margins <- theme(plot.margin = unit(c(0.5, 1, 0.5, 1.5), "cm"))

pdf(file = "lakeID.RandomEffect_WaterTemp_Rich.S.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x2), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.S,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = hill.rich, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x2, predicted, colour = group)) +      
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
                  data = Temp.div.S,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .94))
#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 5)

pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_temp_shan.S")
saveRDS(pred.a, file = "pred.A_timemodel_temp_shan.S")
saveRDS(pred.b, file = "pred.B_timemodel_temp_shan.S")

time.model <- readRDS(file = "time_model_temp_shan.S")
pred.a <- readRDS(file = "pred.A_timemodel_temp_shan.S")
pred.b <- readRDS(file = "pred.B_timemodel_temp_shan.S")

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

pred.a$x <- pred.a$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)


pdf(file = "lakeID.RandomEffect_WaterTemp_Hill.Shannon.S.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.S,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = hill.Shannon, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
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
                  data = Temp.div.S,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .95))
#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 4)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

#saving R objects
saveRDS(time.model, file = "time_model_temp_simp.S")
saveRDS(pred.a, file = "pred.A_timemodel_temp_simp.S")
saveRDS(pred.b, file = "pred.B_timemodel_temp_simp.S")

time.model <- readRDS(file = "time_model_temp_simp.S")
pred.a <- readRDS(file = "pred.A_timemodel_temp_simp.S")
pred.b <- readRDS(file = "pred.B_timemodel_temp_simp.S")


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

pred.a$x <- pred.a$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)




pdf(file = "lakeID.RandomEffect_WaterTemp_Hill.InverseSimpson.S.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.S,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = hill.InverseSimpson, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
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
hist(Temp.div.S$lake.LCBD.q0)
qqPlot(Temp.div.S$lake.LCBD.q0)

summary(Temp.div.S)

#standard “Smithson & Verkuilen” (2006) transformation so that we can use zero inflated beta with zeros and 1s
n <- nrow(Temp.div.S)
Temp.div.S$lake.LCBD.q0 <- (Temp.div.S$lake.LCBD.q0 * (n - 1) + 0.5) / n


time.model <- brm(formula = lake.LCBD.q0 ~ WaterTemp_z + (1 + WaterTemp_z | lakeID),
                  data = Temp.div.S,
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

pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

#saving R objects
saveRDS(time.model, file = "time_model_temp_LCBD.q0.S")
saveRDS(pred.a, file = "pred.A_timemodel_temp_LCBD.q0.S")
saveRDS(pred.b, file = "pred.B_timemodel_temp_LCBD.q0.S")

time.model <- readRDS(file = "time_model_temp_LCBD.q0.S")
pred.a <- readRDS(file = "pred.A_timemodel_temp_LCBD.q0.S")
pred.b <- readRDS(file = "pred.B_timemodel_temp_LCBD.q0.S")


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

pred.a$x <- pred.a$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)




pdf(file = "lakeID.RandomEffect_WaterTemp_LCBD.q0.S.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.S,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = lake.LCBD.q0, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
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
  guides(color = "none")+
  fixed_margins
ggpred_plot
dev.off()
ggpred_plot
TempPlot.B.q0 <- ggpred_plot




######LCBD Lake q1#####
##building model
hist(Temp.div.S$lake.LCBD.q1)
qqPlot(Temp.div.S$lake.LCBD.q1)
summary(Temp.div.S)

#standard “Smithson & Verkuilen” (2006) transformation so that we can use zero inflated beta with zeros and 1s
n <- nrow(Temp.div.S)
Temp.div.S$lake.LCBD.q1 <- (Temp.div.S$lake.LCBD.q1 * (n - 1) + 0.5) / n

time.model <- brm(formula = lake.LCBD.q1 ~ WaterTemp_z + (1 + WaterTemp_z | lakeID),
                  data = Temp.div.S,
                  family = Beta(link = "cloglog"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .8)) #Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 4)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

#saving R objects
saveRDS(time.model, file = "time_model_temp_LCBD.q1.S")
saveRDS(pred.a, file = "pred.A_timemodel_temp_LCBD.q1.S")
saveRDS(pred.b, file = "pred.B_timemodel_temp_LCBD.q1.S")

time.model <- readRDS(file = "time_model_temp_LCBD.q1.S")
pred.a <- readRDS(file = "pred.A_timemodel_temp_LCBD.q1.S")
pred.b <- readRDS(file = "pred.B_timemodel_temp_LCBD.q1.S")


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

pred.a$x <- pred.a$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)




pdf(file = "lakeID.RandomEffect_WaterTemp_LCBD.q1.S.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.S,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = lake.LCBD.q1, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
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
  guides(color = "none")+
  fixed_margins
ggpred_plot
dev.off()
ggpred_plot
TempPlot.B.q1 <- ggpred_plot

######LCBD Lake q2#####
##building model
hist(Temp.div.S$lake.LCBD.q2)
qqPlot(Temp.div.S$lake.LCBD.q2)


summary(Temp.div.S)

#standard “Smithson & Verkuilen” (2006) transformation so that we can use zero inflated beta with zeros and 1s
n <- nrow(Temp.div.S)
Temp.div.S$lake.LCBD.q2 <- (Temp.div.S$lake.LCBD.q2 * (n - 1) + 0.5) / n

time.model <- brm(formula = lake.LCBD.q2 ~ WaterTemp_z + (1 + WaterTemp_z | lakeID),
                  data = Temp.div.S,
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

pred.a <- ggpredict(time.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("WaterTemp_z [all]", "lakeID"), type = "re", ppd = FALSE)

#saving R objects
saveRDS(time.model, file = "time_model_temp_LCBD.q2.S")
saveRDS(pred.a, file = "pred.A_timemodel_temp_LCBD.q2.S")
saveRDS(pred.b, file = "pred.B_timemodel_temp_LCBD.q2.S")

time.model <- readRDS(file = "time_model_temp_LCBD.q2.S")
pred.a <- readRDS(file = "pred.A_timemodel_temp_LCBD.q2.S")
pred.b <- readRDS(file = "pred.B_timemodel_temp_LCBD.q2.S")


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

pred.a$x <- pred.a$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)




pdf(file = "lakeID.RandomEffect_WaterTemp_LCBD.q2.S.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.S,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = lake.LCBD.q2, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Water Temperature (C°)", y = "LCBD (q = 2)") + 
                  theme_bw()) +
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14),  # Y-axis label size
    axis.text = element_text(size = 12)  # Axis tick labels size 
  ) +
  guides(color = "none")+
  fixed_margins
ggpred_plot
dev.off()
ggpred_plot
TempPlot.B.q2 <- ggpred_plot

#####TEMP SPACE MODEL#####
######Richness##########
##building model
Temp.div.S$year <- as.factor(Temp.div.S$year)
summary(Temp.div.S)

space.model <- brm(formula = hill.rich ~ WaterTemp_z + (1 + WaterTemp_z | year),
                   data = Temp.div.S,
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


pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_rich.S")
saveRDS(pred.a, file = "pred.a_spacemodel_temp_rich.S")
saveRDS(pred.b, file = "pred.b_spacemodel_temp_rich.S")

space.model <- readRDS(file = "space_model_temp_rich.S")
pred.a <- readRDS(file = "pred.a_spacemodel_temp_rich.S")
pred.b <- readRDS(file = "pred.b_spacemodel_temp_rich.S")

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

pred.a$x <- pred.a$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_Rich.S.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.S,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = hill.rich, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
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
hist(Temp.div.S$hill.Shannon)
qqPlot(Temp.div.S$hill.Shannon)


##building model
space.model <- brm(formula = hill.Shannon ~ WaterTemp_z + (1 + WaterTemp_z | year),
                   data = Temp.div.S,
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


pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_shan.S")
saveRDS(pred.a, file = "pred.A_spacemodel_temp_shan.S")
saveRDS(pred.b, file = "pred.B_spacemodel_temp_shan.S")

space.model <- readRDS(file = "space_model_temp_shan.S")
pred.a <- readRDS(file = "pred.A_spacemodel_temp_shan.S")
pred.b <- readRDS(file = "pred.B_spacemodel_temp_shan.S")

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

pred.a$x <- pred.a$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_Shannon.S.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.S,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = hill.Shannon, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
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
                   data = Temp.div.S,
                   family = Gamma(link = "log"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .9)) #0.8 is default, determines how long it will "bounce", increasing can aid in convergence
#model checks
pairs(space.model)
summary(space.model)
print(space.model, digit = 4)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_simp.S")
saveRDS(pred.a, file = "pred.A_spacemodel_temp_simp.S")
saveRDS(pred.b, file = "pred.B_spacemodel_temp_simp.S")

space.model <- readRDS(file = "space_model_temp_simp.S")
pred.a <- readRDS(file = "pred.A_spacemodel_temp_simp.S")
pred.b <- readRDS(file = "pred.B_spacemodel_temp_simp.S")

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

pred.a$x <- pred.a$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_Inverse Simpson.S.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.S,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = hill.InverseSimpson, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
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
Temp.div.S$year <- as.factor(Temp.div.S$year)
summary(Temp.div.S)

hist(Temp.div.S$year.LCBD.q0)
qqPlot(Temp.div.S$year.LCBD.q0)

Temp.div.S$year.LCBD.q0.t <- (Temp.div.S$year.LCBD.q0 +1)
hist(Temp.div.S$year.LCBD.q0.t)
qqPlot(Temp.div.S$year.LCBD.q0.t)

#standard “Smithson & Verkuilen” (2006) transformation so that we can use zero inflated beta with zeros and 1s
n <- nrow(Temp.div.S)
Temp.div.S$year.LCBD.q0 <- (Temp.div.S$year.LCBD.q0 * (n - 1) + 0.5) / n
Temp.div.S$year.LCBD.q1 <- (Temp.div.S$year.LCBD.q1 * (n - 1) + 0.5) / n
Temp.div.S$year.LCBD.q2 <- (Temp.div.S$year.LCBD.q2 * (n - 1) + 0.5) / n


space.model <- brm(formula = year.LCBD.q0 ~ WaterTemp_z + (1 + WaterTemp_z | year),
                   data = Temp.div.S,
                   family = Beta(link = "cloglog"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = 0.8, max_treedepth = 12))
#model checks
pairs(space.model)
summary(space.model)
print(space.model, digits = 4)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)


pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_year.LCBD.q0.S")
saveRDS(pred.a, file = "pred.a_spacemodel_temp_year.LCBD.q0.S")
saveRDS(pred.b, file = "pred.b_spacemodel_temp_year.LCBD.q0.S")

space.model <- readRDS(file = "space_model_temp_year.LCBD.q0.S")
pred.a <- readRDS(file = "pred.a_spacemodel_temp_year.LCBD.q0.S")
pred.b <- readRDS(file = "pred.b_spacemodel_temp_year.LCBD.q0.S")

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

pred.a$x <- pred.a$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_year.LCBD.q0.S.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.S,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = year.LCBD.q0, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Water Temperature (C)", y = "LCBD (q = 0)") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot

dev.off()
ggpred_plot
TempPlot.B.q0.S <- ggpred_plot


######LCBD Year q1##########
##building model
Temp.div.S$year <- as.factor(Temp.div.S$year)
summary(Temp.div.S)

hist(Temp.div.S$year.LCBD.q1)
qqPlot(Temp.div.S$year.LCBD.q1)

space.model <- brm(formula = year.LCBD.q1 ~ WaterTemp_z + (1 + WaterTemp_z | year),
                   data = Temp.div.S,
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


pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_year.LCBD.q1.S")
saveRDS(pred.a, file = "pred.a_spacemodel_temp_year.LCBD.q1.S")
saveRDS(pred.b, file = "pred.b_spacemodel_temp_year.LCBD.q1.S")

space.model <- readRDS(file = "space_model_temp_year.LCBD.q1.S")
pred.a <- readRDS(file = "pred.a_spacemodel_temp_year.LCBD.q1.S")
pred.b <- readRDS(file = "pred.b_spacemodel_temp_year.LCBD.q1.S")

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

pred.a$x <- pred.a$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_year.LCBD.q1.S.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.S,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = year.LCBD.q1, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
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
Temp.div.S$year <- as.factor(Temp.div.S$year)
summary(Temp.div.S)

hist(Temp.div.S$year.LCBD.q2)
qqPlot(Temp.div.S$year.LCBD.q2)

space.model <- brm(formula = year.LCBD.q2 ~ WaterTemp_z + (1 + WaterTemp_z | year),
                   data = Temp.div.S,
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


pred.a <- ggpredict(space.model, terms = c("WaterTemp_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("WaterTemp_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_temp_year.LCBD.q2.S")
saveRDS(pred.a, file = "pred.a_spacemodel_temp_year.LCBD.q2.S")
saveRDS(pred.b, file = "pred.b_spacemodel_temp_year.LCBD.q2.S")

space.model <- readRDS(file = "space_model_temp_year.LCBD.q2.S")
pred.a <- readRDS(file = "pred.a_spacemodel_temp_year.LCBD.q2.S")
pred.b <- readRDS(file = "pred.b_spacemodel_temp_year.LCBD.q2.S")

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

pred.a$x <- pred.a$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)
pred.b$x <- pred.b$x * sd(Temp.div.S$Temp.water.degC) + mean(Temp.div.S$Temp.water.degC)


pdf(file = "year.RandomEffect_WaterTemp_year.LCBD.q2.S.pdf")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Temp.div.S,                      # adding the raw data (scaled values)
                             aes(x = Temp.water.degC, y = year.LCBD.q2, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Water Temperature (C)", y = "LCBD (q = 2)") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot

dev.off()
ggpred_plot
TempPlot.B.q2.S <- ggpred_plot


#####Nitrogen TIME MODEL#####
#####Not on SID#####

##TN Trial from grace's code using TN instead of nitrate-nitrite
Nit.div.L <- WQ.div.early %>%
  group_by(lake.date) %>%
  mutate(
    Total_Nitrogen = if_else(is.na(Nutrient.nitrogen_mg.l), Kjeldahl_nitrogen.as_N. + Nitrate.Nitrite_.as_N._mg.l, Nutrient.nitrogen_mg.l)
  ) %>%
  ungroup()  
summary(Nit.div.L)


##Determine the total number of NAs within Water Temp Column
na_count <- sum(is.na(WQ.div.early$Total_Nitrogen))
##Removing NAs
Nit.div.early <- WQ.div.early[complete.cases(WQ.div.early$Total_Nitrogen), ]  
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

#Calculating Z
Nit.div.early$TN_z <- scale(Nit.div.early$Total_Nitrogen, center = TRUE, scale = TRUE)

#####FILTER FOR ON or off the Southern Iowa drift plain
Nit.div.L <- Nit.div.early %>%
  filter(SouthernIowaDrift == "0")


Nit.div.S <- Nit.div.early %>%
  filter(SouthernIowaDrift == "1")




######Richness#####
##building model ### Had to adjust some parameters (treedepth for convergence)
time.model <- brm(formula = hill.rich ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.L,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3100,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .99, max_treedepth = 25))



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

saveRDS(time.model, file = "time_model_TN_rich.L")
saveRDS(pred.a, file = "pred.A_timemodel_TN_rich.L")
saveRDS(pred.b, file = "pred.B_timemodel_TN_rich.L")

time.model <- readRDS(file = "time_model_TN_rich.L")
pred.a <- readRDS(file = "pred.A_timemodel_TN_rich.L")
pred.b <- readRDS(file = "pred.B_timemodel_TN_rich.L")


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

pred.a$x <- pred.a$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitriogen_Rich.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.rich, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "Richness (as Hill numbers)", 
                       title = "Time Model of Total Nitrogen & Species Richness") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot




######Hill.Shannon#####
##building model
hist(Nit.div.L$hill.Shannon)


time.model <- brm(formula = hill.Shannon ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.L,
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
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

pred.a <- ggpredict(time.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("TN_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_TN_shan.L")
saveRDS(pred.a, file = "pred.A_timemodel_TN_shan.L")
saveRDS(pred.b, file = "pred.B_timemodel_TN_shan.L")

time.model <- readRDS(file = "time_model_TN_shan.L")
pred.a <- readRDS(file = "pred.A_timemodel_TN_shan.L")
pred.b <- readRDS(file = "pred.B_timemodel_TN_shan.L")


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

pred.a$x <- pred.a$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitrogen_Hill.Shannon.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.Shannon, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "Shannon (as Hill numbers)", 
                       title = "Time Model of Total Nitrogen & Shannon") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot


######Hill.Inverse Simpson####
## building model
hist(Nit.div.L$hill.InverseSimpson)


time.model <- brm(formula = hill.InverseSimpson ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.L,
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

pred.a <- ggpredict(time.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("TN_z [all]", "lakeID"), type = "re", ppd = FALSE)

summary(pred.a)
summary(pred.b)

saveRDS(time.model, file = "time_model_TN_simp.L")
saveRDS(pred.a, file = "pred.A_timemodel_TN_simp.L")
saveRDS(pred.b, file = "pred.B_timemodel_TN_simp.L")

time.model <- readRDS(file = "time_model_TN_simp.L")
pred.a <- readRDS(file = "pred.A_timemodel_TN_simp.L")
pred.b <- readRDS(file = "pred.B_timemodel_TN_simp.L")


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

pred.a$x <- pred.a$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitrogen_Hill.InverseSimpson.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.InverseSimpson, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  
                  labs(x = "Total Nitrogen (mg/l)", y = "Inverse Simpson (as Hill numbers)", 
                       title = "Time Model of Total Nitrogen & Inverse Simpson") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot

######LCBD Lake q0######


#standard “Smithson & Verkuilen” (2006) transformation so that we can use zero inflated beta with zeros and 1s
n <- nrow(Nit.div.L)
Nit.div.L$lake.LCBD.q0 <- (Nit.div.L$lake.LCBD.q0 * (n - 1) + 0.5) / n
Nit.div.L$lake.LCBD.q1 <- (Nit.div.L$lake.LCBD.q1 * (n - 1) + 0.5) / n
Nit.div.L$lake.LCBD.q2 <- (Nit.div.L$lake.LCBD.q2 * (n - 1) + 0.5) / n


##building model
time.model <- brm(formula = lake.LCBD.q0 ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.L,
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

saveRDS(time.model, file = "time_model_TN_lake.LCBD.q0.L")
saveRDS(pred.a, file = "pred.A_timemodel_TN_lake.LCBD.q0.L")
saveRDS(pred.b, file = "pred.B_timemodel_TN_lake.LCBD.q0.L")

time.model <- readRDS(file = "time_model_TN_lake.LCBD.q0.L")
pred.a <- readRDS(file = "pred.A_timemodel_TN_lake.LCBD.q0.L")
pred.b <- readRDS(file = "pred.B_timemodel_TN_lake.LCBD.q0.L")


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

pred.a$x <- pred.a$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitrogen_lake.LCBD.q0.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = lake.LCBD.q0, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "LCBD (q = 0)", 
                       title = "Time Model of Total Nitrogen & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot

######LCBD Lake q1######

##building model
time.model <- brm(formula = lake.LCBD.q1 ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.L,
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

pred.a <- ggpredict(time.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("TN_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_TN_lake.LCBD.q1.L")
saveRDS(pred.a, file = "pred.A_timemodel_TN_lake.LCBD.q1.L")
saveRDS(pred.b, file = "pred.B_timemodel_TN_lake.LCBD.q1.L")

time.model <- readRDS(file = "time_model_TN_lake.LCBD.q1.L")
pred.a <- readRDS(file = "pred.A_timemodel_TN_lake.LCBD.q1.L")
pred.b <- readRDS(file = "pred.B_timemodel_TN_lake.LCBD.q1.L")


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

pred.a$x <- pred.a$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitrogen_lake.LCBD.q1.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = lake.LCBD.q1, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "LCBD (q = 1)", 
                       title = "Time Model of Total Nitrogen & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot
dev.off()
ggpred_plot

######LCBD Lake q2######

##building model
time.model <- brm(formula = lake.LCBD.q2 ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.L,
                  family = Beta(link = "cloglog"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .85))



#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 6)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

pred.a <- ggpredict(time.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("TN_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_TN_lake.LCBD.q2.L")
saveRDS(pred.a, file = "pred.A_timemodel_TN_lake.LCBD.q2.L")
saveRDS(pred.b, file = "pred.B_timemodel_TN_lake.LCBD.q2.L")

time.model <- readRDS(file = "time_model_TN_lake.LCBD.q2.L")
pred.a <- readRDS(file = "pred.A_timemodel_TN_lake.LCBD.q2.L")
pred.b <- readRDS(file = "pred.B_timemodel_TN_lake.LCBD.q2.L")


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

pred.a$x <- pred.a$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitrogen_lake.LCBD.q2.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = lake.LCBD.q2, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Nitrate Nitrite (mg/l)", y = "LCBD (q = 2)", 
                       title = "Time Model of Total Nitrogen & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot





#####Nitrogen SPACE MODEL#####
######Hill.Richness####
##building model
Nit.div.L$year <- as.factor(Nit.div.L$year)
summary(Nit.div.L)

space.model <- brm(formula = hill.rich ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.L,
                   family = Gamma(link = "log"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .8))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence

#model checks
pairs(space.model)
summary(space.model)
print(space.model, digits = 6)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)


pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_TN_rich.L")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_rich.L")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_rich.L")

space.model <- readRDS(file = "space_model_TN_rich.L")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_rich.L")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_rich.L")


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

pred.a$x <- pred.a$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_Rich.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.rich, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "Richness (as Hill numbers)", 
                       title = "Space Model of Total Nitrogen & Species Richness") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot


#######Hill.Shannon########
hist(Nit.div.L$hill.Shannon)
qqPlot(Nit.div.L$hill.Shannon)


##building model
space.model <- brm(formula = hill.Shannon ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.L,
                   family = Gamma(link = "log"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .85))
#model checks
pairs(space.model)
print(time.model, digits = 5)
summary(space.model)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)


pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_TN_shan.L")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_shan.L")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_shan.L")

space.model <- readRDS(file = "space_model_TN_shan.L")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_shan.L")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_shan.L")


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

pred.a$x <- pred.a$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_Shannon.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.Shannon, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "Shannon (as Hill numbers)", 
                       title = "Space Model of Total Nitrogen & Shannon") + 
                  theme_minimal()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot


#######hill.Inverse Simpson####

space.model <- brm(formula = hill.InverseSimpson ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.L,
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

pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]","year"), type = "re", ppd = FALSE)

summary(pred.a)
summary(pred.b)

saveRDS(space.model, file = "space_model_TN_simp.L")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_simp.L")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_simp.L")

space.model <- readRDS(file = "space_model_TN_simp.L")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_simp.L")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_simp.L")


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

pred.a$x <- pred.a$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_InverseSimpson.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.InverseSimpson, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x,predicted), linewidth = 1.5, linetype = "longdash") +
                  
                  labs(x = "Total Nitrogen (mg/l)", y = "Inverse Simpson (as Hill numbers)", 
                       title = "Space Model of Total Nitrogen & Inverse Simpson") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()

ggpred_plot

######LCBD Year q0####
##building model
Nit.div.L$year <- as.factor(Nit.div.L$year)
summary(Nit.div.L)

#standard “Smithson & Verkuilen” (2006) transformation so that we can use zero inflated beta with zeros and 1s
n <- nrow(Nit.div.L)
Nit.div.L$year.LCBD.q0 <- (Nit.div.L$year.LCBD.q0 * (n - 1) + 0.5) / n
Nit.div.L$year.LCBD.q1 <- (Nit.div.L$year.LCBD.q1 * (n - 1) + 0.5) / n
Nit.div.L$year.LCBD.q2 <- (Nit.div.L$year.LCBD.q2 * (n - 1) + 0.5) / n

space.model <- brm(formula = year.LCBD.q0 ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.L,
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


pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_TN_year.LCBD.q0.L")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_year.LCBD.q0.L")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_year.LCBD.q0.L")

space.model <- readRDS(file = "space_model_TN_year.LCBD.q0.L")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_year.LCBD.q0.L")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_year.LCBD.q0.L")


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

pred.a$x <- pred.a$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_year.LCBD.q0.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = year.LCBD.q0, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Nitrate Nitrite (mg/l)", y = "LCBD (q = 1)", 
                       title = "Space Model of Total Nitrogen & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot

######LCBD Year q1####
##building model
Nit.div.L$year <- as.factor(Nit.div.L$year)
summary(Nit.div.L)

space.model <- brm(formula = year.LCBD.q1 ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.L,
                   family = Beta(link = "cloglog"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .8, max_treedepth = 12))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence

#model checks
pairs(space.model)
print(space.model, digits = 7)
summary(space.model)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)


pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_TN_year.LCBD.q1.L")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_year.LCBD.q1.L")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_year.LCBD.q1.L")

space.model <- readRDS(file = "space_model_TN_year.LCBD.q1.L")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_year.LCBD.q1.L")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_year.LCBD.q1.L")


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

pred.a$x <- pred.a$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_year.LCBD.q1.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = year.LCBD.q1, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "LCBD (q = 1)", 
                       title = "Space Model of Total Nitrogen & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot


######LCBD Year q2####
##building model
Nit.div.L$year <- as.factor(Nit.div.L$year)
summary(Nit.div.L)

space.model <- brm(formula = year.LCBD.q2 ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.L,
                   family = Beta(link = "cloglog"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .8, max_treedepth = 12))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence

#model checks
pairs(space.model)
summary(space.model)
print(space.model, digits = 5)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)


pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_TN_year.LCBD.q2.L")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_year.LCBD.q2.L")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_year.LCBD.q2.L")

space.model <- readRDS(file = "space_model_TN_year.LCBD.q2.L")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_year.LCBD.q2.L")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_year.LCBD.q2.L")


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

pred.a$x <- pred.a$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.L$Total_Nitrogen) + mean(Nit.div.L$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_year.LCBD.q2.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = year.LCBD.q2, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "LCBD (q = 2)", 
                       title = "Space Model of Total Nitrogen & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot


#####Nitrogen TIME MODEL#####
#####LAKES on SID#####
######Richness#####
##TN Trial from grace's code using TN instead of nitrate-nitrite
Nit.div.S <- Temp.div.S %>%
  group_by(lake.date) %>%
  mutate(
    Total_Nitrogen = if_else(is.na(Nutrient.nitrogen_mg.l), Kjeldahl_nitrogen.as_N. + Nitrate.Nitrite_.as_N._mg.l, Nutrient.nitrogen_mg.l)
  ) %>%
  ungroup()  
summary(Nit.div.S)

#Calculating Z
Nit.div.S$TN_z <- scale(Nit.div.S$Total_Nitrogen, center = TRUE, scale = TRUE)
summary(Nit.div.S)



##Removing NAs
Nit.div.S <- Nit.div.S[complete.cases(Nit.div.S$TN_z), ]  
##There are no NAs in Rich
na_count <- sum(is.na(Nit.div.S$hill.rich))



##building model
time.model <- brm(formula = hill.rich ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.S,
                  family = Gamma(link = "log"),
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

pred.a <- ggpredict(time.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("TN_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_TN_rich.S")
saveRDS(pred.a, file = "pred.A_timemodel_TN_rich.S")
saveRDS(pred.b, file = "pred.B_timemodel_TN_rich.S")

time.model <- readRDS(file = "time_model_TN_rich.S")
pred.a <- readRDS(file = "pred.A_timemodel_TN_rich.S")
pred.b <- readRDS(file = "pred.B_timemodel_TN_rich.S")


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

pred.a$x <- pred.a$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitriogen_Rich.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.rich, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "Richness (as Hill numbers)", 
                       title = "Time Model of Total Nitrogen & Species Richness") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot




######Hill.Shannon#####
##building model
hist(Nit.div.S$hill.Shannon)


time.model <- brm(formula = hill.Shannon ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.S,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .95))
#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 5)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

pred.a <- ggpredict(time.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("TN_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_TN_shan.S")
saveRDS(pred.a, file = "pred.A_timemodel_TN_shan.S")
saveRDS(pred.b, file = "pred.B_timemodel_TN_shan.S")

time.model <- readRDS(file = "time_model_TN_shan.S")
pred.a <- readRDS(file = "pred.A_timemodel_TN_shan.S")
pred.b <- readRDS(file = "pred.B_timemodel_TN_shan.S")


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

pred.a$x <- pred.a$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitrogen_Hill.Shannon.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.Shannon, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "Shannon (as Hill numbers)", 
                       title = "Time Model of Total Nitrogen & Shannon") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot


######Hill.Inverse Simpson####
## building model
hist(Nit.div.S$hill.InverseSimpson)


time.model <- brm(formula = hill.InverseSimpson ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.S,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .97))
#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 5)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

pred.a <- ggpredict(time.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("TN_z [all]", "lakeID"), type = "re", ppd = FALSE)

summary(pred.a)
summary(pred.b)

saveRDS(time.model, file = "time_model_TN_simp.S")
saveRDS(pred.a, file = "pred.A_timemodel_TN_simp.S")
saveRDS(pred.b, file = "pred.B_timemodel_TN_simp.S")

time.model <- readRDS(file = "time_model_TN_simp.S")
pred.a <- readRDS(file = "pred.A_timemodel_TN_simp.S")
pred.b <- readRDS(file = "pred.B_timemodel_TN_simp.S")


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

pred.a$x <- pred.a$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitrogen_Hill.InverseSimpson.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.InverseSimpson, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  
                  labs(x = "Total Nitrogen (mg/l)", y = "Inverse Simpson (as Hill numbers)", 
                       title = "Time Model of Total Nitrogen & Inverse Simpson") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot

######LCBD Lake q0######

#standard “Smithson & Verkuilen” (2006) transformation so that we can use zero inflated beta with zeros and 1s
n <- nrow(Nit.div.S)
Nit.div.S$lake.LCBD.q0 <- (Nit.div.S$lake.LCBD.q0 * (n - 1) + 0.5) / n
Nit.div.S$lake.LCBD.q1 <- (Nit.div.S$lake.LCBD.q1 * (n - 1) + 0.5) / n
Nit.div.S$lake.LCBD.q2 <- (Nit.div.S$lake.LCBD.q2 * (n - 1) + 0.5) / n

##building model
time.model <- brm(formula = lake.LCBD.q0 ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.S,
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

saveRDS(time.model, file = "time_model_TN_lake.LCBD.q0.S")
saveRDS(pred.a, file = "pred.A_timemodel_TN_lake.LCBD.q0.S")
saveRDS(pred.b, file = "pred.B_timemodel_TN_lake.LCBD.q0.S")

time.model <- readRDS(file = "time_model_TN_lake.LCBD.q0.S")
pred.a <- readRDS(file = "pred.A_timemodel_TN_lake.LCBD.q0.S")
pred.b <- readRDS(file = "pred.B_timemodel_TN_lake.LCBD.q0.S")


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

pred.a$x <- pred.a$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitrogen_lake.LCBD.q0.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = lake.LCBD.q0, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "LCBD (q = 0)", 
                       title = "Time Model of Total Nitrogen & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot

######LCBD Lake q1######

##building model
time.model <- brm(formula = lake.LCBD.q1 ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.S,
                  family = Beta(link = "cloglog"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .85))



#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 6)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

pred.a <- ggpredict(time.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("TN_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_TN_lake.LCBD.q1.S")
saveRDS(pred.a, file = "pred.A_timemodel_TN_lake.LCBD.q1.S")
saveRDS(pred.b, file = "pred.B_timemodel_TN_lake.LCBD.q1.S")

time.model <- readRDS(file = "time_model_TN_lake.LCBD.q1.S")
pred.a <- readRDS(file = "pred.A_timemodel_TN_lake.LCBD.q1.S")
pred.b <- readRDS(file = "pred.B_timemodel_TN_lake.LCBD.q1.S")


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

pred.a$x <- pred.a$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitrogen_lake.LCBD.q1.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = lake.LCBD.q1, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "LCBD (q = 1)", 
                       title = "Time Model of Total Nitrogen & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot
dev.off()
ggpred_plot

######LCBD Lake q2######

##building model
time.model <- brm(formula = lake.LCBD.q2 ~ TN_z + (1 + TN_z | lakeID),
                  data = Nit.div.S,
                  family = Beta(link = "cloglog"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .85))



#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 6)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

pred.a <- ggpredict(time.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("TN_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_TN_lake.LCBD.q2.S")
saveRDS(pred.a, file = "pred.A_timemodel_TN_lake.LCBD.q2.S")
saveRDS(pred.b, file = "pred.B_timemodel_TN_lake.LCBD.q2.S")

time.model <- readRDS(file = "time_model_TN_lake.LCBD.q2.S")
pred.a <- readRDS(file = "pred.A_timemodel_TN_lake.LCBD.q2.S")
pred.b <- readRDS(file = "pred.B_timemodel_TN_lake.LCBD.q2.S")


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

pred.a$x <- pred.a$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)


png(file = "lakeID.RandomEffect_TotalNitrogen_lake.LCBD.q2.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = lake.LCBD.q2, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Nitrate Nitrite (mg/l)", y = "LCBD (q = 2)", 
                       title = "Time Model of Total Nitrogen & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot
dev.off()
ggpred_plot





#####Nitrogen SPACE MODEL#####
######Hill.Richness####
##building model
Nit.div.S$year <- as.factor(Nit.div.S$year)
summary(Nit.div.S)

space.model <- brm(formula = hill.rich ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.S,
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


pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_TN_rich.S")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_rich.S")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_rich.S")

space.model <- readRDS(file = "space_model_TN_rich.S")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_rich.S")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_rich.S")


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

pred.a$x <- pred.a$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_Rich.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.rich, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "Richness (as Hill numbers)", 
                       title = "Space Model of Total Nitrogen & Species Richness") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot


#######Hill.Shannon########
hist(Nit.div.S$hill.Shannon)
qqPlot(Nit.div.S$hill.Shannon)


##building model
space.model <- brm(formula = hill.Shannon ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.S,
                   family = Gamma(link = "log"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .8))
#model checks
pairs(space.model)
print(time.model, digits = 5)
summary(space.model)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)


pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_TN_shan.S")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_shan.S")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_shan.S")

space.model <- readRDS(file = "space_model_TN_shan.S")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_shan.S")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_shan.S")


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

pred.a$x <- pred.a$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_Shannon.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.Shannon, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "Shannon (as Hill numbers)", 
                       title = "Space Model of Total Nitrogen & Shannon") + 
                  theme_minimal()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot


#######hill.Inverse Simpson####

space.model <- brm(formula = hill.InverseSimpson ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.S,
                   family = Gamma(link = "log"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .87))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence
#model checks
pairs(space.model)
summary(space.model)

pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)


pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]","year"), type = "re", ppd = FALSE)

summary(pred.a)
summary(pred.b)

saveRDS(space.model, file = "space_model_TN_simp.S")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_simp.S")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_simp.S")

space.model <- readRDS(file = "space_model_TN_simp.S")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_simp.S")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_simp.S")


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

pred.a$x <- pred.a$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_InverseSimpson.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = hill.InverseSimpson, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x,predicted), linewidth = 1.5, linetype = "longdash") +
                  
                  labs(x = "Total Nitrogen (mg/l)", y = "Inverse Simpson (as Hill numbers)", 
                       title = "Space Model of Total Nitrogen & Inverse Simpson") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()

ggpred_plot

######LCBD Year q0####
##building model
Nit.div.S$year <- as.factor(Nit.div.S$year)
summary(Nit.div.S)

#standard “Smithson & Verkuilen” (2006) transformation so that we can use zero inflated beta with zeros and 1s
n <- nrow(Nit.div.S)
Nit.div.S$year.LCBD.q0 <- (Nit.div.S$year.LCBD.q0 * (n - 1) + 0.5) / n
Nit.div.S$year.LCBD.q1 <- (Nit.div.S$year.LCBD.q1 * (n - 1) + 0.5) / n
Nit.div.S$year.LCBD.q2 <- (Nit.div.S$year.LCBD.q2 * (n - 1) + 0.5) / n

space.model <- brm(formula = year.LCBD.q0 ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.S,
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


pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_TN_year.LCBD.q0.S")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_year.LCBD.q0.S")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_year.LCBD.q0.S")

space.model <- readRDS(file = "space_model_TN_year.LCBD.q0.S")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_year.LCBD.q0.S")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_year.LCBD.q0.S")


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

pred.a$x <- pred.a$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_year.LCBD.q0.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = year.LCBD.q0, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Nitrate Nitrite (mg/l)", y = "LCBD (q = 1)", 
                       title = "Space Model of Total Nitrogen & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot

######LCBD Year q1####
##building model
Nit.div.S$year <- as.factor(Nit.div.S$year)
summary(Nit.div.S)

space.model <- brm(formula = year.LCBD.q1 ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.S,
                   family = Beta(link = "cloglog"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .8, max_treedepth = 12))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence

#model checks
pairs(space.model)
print(space.model, digits = 7)
summary(space.model)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)


pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_TN_year.LCBD.q1.S")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_year.LCBD.q1.S")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_year.LCBD.q1.S")

space.model <- readRDS(file = "space_model_TN_year.LCBD.q1.S")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_year.LCBD.q1.S")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_year.LCBD.q1.S")


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

pred.a$x <- pred.a$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_year.LCBD.q1.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = year.LCBD.q1, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "LCBD (q = 1)", 
                       title = "Space Model of Total Nitrogen & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot


######LCBD Year q2####
##building model
Nit.div.S$year <- as.factor(Nit.div.S$year)
summary(Nit.div.S)

space.model <- brm(formula = year.LCBD.q2 ~ TN_z + (1 + TN_z | year),
                   data = Nit.div.S,
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


pred.a <- ggpredict(space.model, terms = c("TN_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("TN_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_TN_year.LCBD.q2.S")
saveRDS(pred.a, file = "pred.A_spacemodel_TN_year.LCBD.q2.S")
saveRDS(pred.b, file = "pred.B_spacemodel_TN_year.LCBD.q2.S")

space.model <- readRDS(file = "space_model_TN_year.LCBD.q2.S")
pred.a <- readRDS(file = "pred.A_spacemodel_TN_year.LCBD.q2.S")
pred.b <- readRDS(file = "pred.B_spacemodel_TN_year.LCBD.q2.S")


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

pred.a$x <- pred.a$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)
pred.b$x <- pred.b$x * sd(Nit.div.S$Total_Nitrogen) + mean(Nit.div.S$Total_Nitrogen)


png(file = "year.RandomEffect_TotalNitrogen_year.LCBD.q2.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = Nit.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Nitrogen, y = year.LCBD.q2, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Total Nitrogen (mg/l)", y = "LCBD (q = 2)", 
                       title = "Space Model of Total Nitrogen & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = "none")
ggpred_plot

dev.off()
ggpred_plot





#####Phosphorus TIME MODEL#####
#####Not SID LAKES#####
######Hill.Richness####

##Determine the total number of NAs within Water Temp Column
na_count <- sum(is.na(WQ.div.early$Total_Phosphorus))
##Removing NAs
TP.div.early <- WQ.div.early[complete.cases(WQ.div.early$Total_Phosphorus), ]  
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
TP.div.early$Phosphorus_z <- scale(TP.div.early$Total_Phosphorus, center = TRUE, scale = TRUE)

#####FILTER FOR ON or off the Southern Iowa drift plain
TP.div.L <- TP.div.early %>%
  filter(SouthernIowaDrift == "0")


TP.div.S <- TP.div.early %>%
  filter(SouthernIowaDrift == "1")



##building model
time.model <- brm(formula = hill.rich ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.L,
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

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_Phos_rich.L")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_rich.L")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_rich.L")

time.model <- readRDS(file = "time_model_Phos_rich.L")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_rich.L")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_rich.L")

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

pred.a$x <- pred.a$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)

png(file = "lakeID.RandomEffect_Total_Phosphorus_Rich.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.rich, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "Richness (as Hill numbers)", 
                       title = "Time Model of Phosphorus & Species Richness") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot
dev.off()
ggpred_plot


######Hill.Shannon####
##building model
hist(TP.div.L$hill.Shannon)


time.model <- brm(formula = hill.Shannon ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.L,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .92))
#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 5)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z[all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z[all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_Phos_shan.L")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_shan.L")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_shan.L")

time.model <- readRDS(file = "time_model_Phos_shan.L")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_shan.L")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_shan.L")

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

pred.a$x <- pred.a$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)

png(file = "lakeID.RandomEffect_Phosphorus_z_Hill.Shannon.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.Shannon, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  
                  labs(x = "Phosphorus (mg/l)", y = "Shannon (as Hill numbers)", 
                       title = "Time Model of Phosphorus & Shannon") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot
dev.off()
ggpred_plot


######Hill.Inverse Simpson#####
## building model
hist(TP.div.L$hill.InverseSimpson)


time.model <- brm(formula = hill.InverseSimpson ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.L,
                  family = Gamma(link = "log"),
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

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_Phos_simp.L")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_simp.L")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_simp.L")

time.model <- readRDS(file = "time_model_Phos_simp.L")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_simp.L")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_simp.L")

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

pred.a$x <- pred.a$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


png(file = "lakeID.RandomEffect_Phosphorus_z_Hill.InverseSimpson.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.InverseSimpson, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  
                  labs(x = "Phosphorus (mg/l)", y = "Inverse Simpson (as Hill numbers)", 
                       title = "Time Model of Phosphorus & Inverse Simpson") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot
dev.off()
ggpred_plot


######LCBD Lake q0####

n <- nrow(TP.div.L)
TP.div.L$lake.LCBD.q0 <- (TP.div.L$lake.LCBD.q0 * (n - 1) + 0.5) / n
TP.div.L$lake.LCBD.q1 <- (TP.div.L$lake.LCBD.q1 * (n - 1) + 0.5) / n
TP.div.L$lake.LCBD.q2 <- (TP.div.L$lake.LCBD.q2 * (n - 1) + 0.5) / n


##building model
time.model <- brm(formula = lake.LCBD.q0 ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.L,
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

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_Phos_lake.LCBD.q0.L")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_lake.LCBD.q0.L")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_lake.LCBD.q0.L")

time.model <- readRDS(file = "time_model_Phos_lake.LCBD.q0.L")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_lake.LCBD.q0.L")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_lake.LCBD.q0.L")

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

pred.a$x <- pred.a$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)

png(file = "lakeID.RandomEffect_Total_Phosphorus_lake.LCBD.q0.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = lake.LCBD.q0, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 0)", 
                       title = "Time Model of Phosphorus & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot
dev.off()
ggpred_plot


######LCBD Lake q1####

##building model
time.model <- brm(formula = lake.LCBD.q1 ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.L,
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

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_Phos_lake.LCBD.q1.L")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_lake.LCBD.q1.L")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_lake.LCBD.q1.L")

time.model <- readRDS(file = "time_model_Phos_lake.LCBD.q1.L")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_lake.LCBD.q1.L")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_lake.LCBD.q1.L")

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

pred.a$x <- pred.a$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)

png(file = "lakeID.RandomEffect_Total_Phosphorus_lake.LCBD.q1.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = lake.LCBD.q1, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 1)", 
                       title = "Time Model of Phosphorus & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot
dev.off()
ggpred_plot

######LCBD Lake q2####

##building model
time.model <- brm(formula = lake.LCBD.q2 ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.L,
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

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_Phos_lake.LCBD.q2.L")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_lake.LCBD.q2.L")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_lake.LCBD.q2.L")

time.model <- readRDS(file = "time_model_Phos_lake.LCBD.q2.L")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_lake.LCBD.q2.L")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_lake.LCBD.q2.L")

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

pred.a$x <- pred.a$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)

png(file = "lakeID.RandomEffect_Total_Phosphorus_lake.LCBD.q2.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = lake.LCBD.q2, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 2)", 
                       title = "Time Model of Phosphorus & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot
dev.off()
ggpred_plot


#####Phosphorus SPACE MODEL#####
######Richness######
##building model
TP.div.L$year <- as.factor(TP.div.L$year)
summary(TP.div.L)

space.model <- brm(formula = hill.rich ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.L,
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


pred.a <- ggpredict(space.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_rich.L")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_rich.L")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_rich.L")

space.model <- readRDS(file = "space_model_Phos_rich.L")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_rich.L")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_rich.L")

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

pred.a$x <- pred.a$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


png(file = "year.RandomEffect_Phosphorus_z_Rich.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.rich, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "Richness (as Hill numbers)", 
                       title = "Space Model of Phosphorus & Species Richness") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot

dev.off()
ggpred_plot



######Hill.Shannon####
hist(TP.div.L$hill.Shannon)
qqPlot(TP.div.L$hill.Shannon)


##building model
space.model <- brm(formula = hill.Shannon ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.L,
                   family = Gamma(link = "log"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .83))
#model checks
pairs(space.model)
summary(space.model)
print(space.model, digits = 4)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)


pred.a <- ggpredict(space.model, terms = c("Phosphorus_z"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_shan.L")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_shan.L")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_shan.L")

space.model <- readRDS(file = "space_model_Phos_shan.L")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_shan.L")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_shan.L")


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

pred.a$x <- pred.a$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


png(file = "year.RandomEffect_Phosphorus_Shannon.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.Shannon, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x,predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "Shannon (as Hill numbers)", 
                       title = "Space Model of Phosphorus & Shannon") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot

dev.off()
ggpred_plot


######Inverse Simpson####

space.model <- brm(formula = hill.InverseSimpson ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.L,
                   family = Gamma(link = "log"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .9))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence
#model checks
pairs(space.model)
summary(space.model)
print(space.model, digits = 4)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)

pred.a <- ggpredict(space.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_simp.L")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_simp.L")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_simp.L")

space.model <- readRDS(file = "space_model_Phos_simp.L")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_simp.L")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_simp.L")

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

pred.a$x <- pred.a$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)

png(file = "year.RandomEffect_Phosphorus_InverseSimpson.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.InverseSimpson, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),linewidth = 1) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  
                  labs(x = "Phosphorus (mg/l)", y = "hill.InverseSimpson", 
                       title = "Time Model of Phosphorus & Inverse Simpson") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot

dev.off()
ggpred_plot

######LCBD Year q0######
##building model
TP.div.L$year <- as.factor(TP.div.L$year)
summary(TP.div.L)

n <- nrow(TP.div.L)
TP.div.L$year.LCBD.q0 <- (TP.div.L$year.LCBD.q0 * (n - 1) + 0.5) / n
TP.div.L$year.LCBD.q1 <- (TP.div.L$year.LCBD.q1 * (n - 1) + 0.5) / n
TP.div.L$year.LCBD.q2 <- (TP.div.L$year.LCBD.q2 * (n - 1) + 0.5) / n

space.model <- brm(formula = year.LCBD.q0 ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.L,
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


pred.a <- ggpredict(space.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_year.LCBD.q0.L")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_year.LCBD.q0.L")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_year.LCBD.q0.L")

space.model <- readRDS(file = "space_model_Phos_year.LCBD.q0.L")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_year.LCBD.q0.L")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_year.LCBD.q0.L")

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

pred.a$x <- pred.a$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


png(file = "year.RandomEffect_Phosphorus_z_year.LCBD.q0.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = year.LCBD.q0, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 0)", 
                       title = "Space Model of Phosphorus & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot

dev.off()
ggpred_plot


######LCBD Year q1######
##building model
TP.div.L$year <- as.factor(TP.div.L$year)
summary(TP.div.L)

space.model <- brm(formula = year.LCBD.q1 ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.L,
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


pred.a <- ggpredict(space.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_year.LCBD.q1.L")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_year.LCBD.q1.L")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_year.LCBD.q1.L")

space.model <- readRDS(file = "space_model_Phos_year.LCBD.q1.L")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_year.LCBD.q1.L")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_year.LCBD.q1.L")

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

pred.a$x <- pred.a$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


png(file = "year.RandomEffect_Phosphorus_z_year.LCBD.q1.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = year.LCBD.q1, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 1)", 
                       title = "Space Model of Phosphorus & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot


######LCBD Year q2######
##building model
TP.div.L$year <- as.factor(TP.div.L$year)
summary(TP.div.L)

space.model <- brm(formula = year.LCBD.q2 ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.L,
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


pred.a <- ggpredict(space.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_year.LCBD.q2.L")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_year.LCBD.q2.L")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_year.LCBD.q2.L")

space.model <- readRDS(file = "space_model_Phos_year.LCBD.q2.L")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_year.LCBD.q2.L")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_year.LCBD.q2.L")

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

pred.a$x <- pred.a$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.L$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


png(file = "year.RandomEffect_Phosphorus_z_year.LCBD.q2.L.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.L,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = year.LCBD.q2, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 2)", 
                       title = "Space Model of Phosphorus & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot

#####Phosphorus Time Models#####
#####Lakes on SID LAKES#####

##building model
time.model <- brm(formula = hill.rich ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.S,
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

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_Phos_rich.S")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_rich.S")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_rich.S")

time.model <- readRDS(file = "time_model_Phos_rich.S")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_rich.S")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_rich.S")

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

pred.a$x <- pred.a$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)

png(file = "lakeID.RandomEffect_Total_Phosphorus_Rich.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.rich, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "Richness (as Hill numbers)", 
                       title = "Time Model of Phosphorus & Species Richness") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot
dev.off()
ggpred_plot


######Hill.Shannon####
##building model
hist(TP.div.S$hill.Shannon)


time.model <- brm(formula = hill.Shannon ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.S,
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

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z[all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z[all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_Phos_shan.S")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_shan.S")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_shan.S")

time.model <- readRDS(file = "time_model_Phos_shan.S")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_shan.S")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_shan.S")

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

pred.a$x <- pred.a$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)

png(file = "lakeID.RandomEffect_Phosphorus_z_Hill.Shannon.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.Shannon, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  
                  labs(x = "Phosphorus (mg/l)", y = "Shannon (as Hill numbers)", 
                       title = "Time Model of Phosphorus & Shannon") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot
dev.off()
ggpred_plot


######Hill.Inverse Simpson#####
## building model
hist(TP.div.S$hill.InverseSimpson)


time.model <- brm(formula = hill.InverseSimpson ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.S,
                  family = Gamma(link = "log"),
                  warmup = 1000,
                  iter = 3000,
                  chains = 4,
                  seed = 32,
                  control = list(adapt_delta = .85))
#Model Checks
pairs(time.model)
summary(time.model)
print(time.model, digits = 5)
pp_check(time.model, ndraws = 100)
diagnostic_posterior(time.model, diagnostic = "all", effects = c("fixed"))
plot(time.model)
bayes_R2(time.model)

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_Phos_simp.S")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_simp.S")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_simp.S")

time.model <- readRDS(file = "time_model_Phos_simp.S")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_simp.S")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_simp.S")

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

pred.a$x <- pred.a$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


png(file = "lakeID.RandomEffect_Phosphorus_z_Hill.InverseSimpson.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.InverseSimpson, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  
                  labs(x = "Phosphorus (mg/l)", y = "Inverse Simpson (as Hill numbers)", 
                       title = "Time Model of Phosphorus & Inverse Simpson") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot
dev.off()
ggpred_plot


######LCBD Lake q0####

n <- nrow(TP.div.S)
TP.div.S$lake.LCBD.q0 <- (TP.div.S$lake.LCBD.q0 * (n - 1) + 0.5) / n
TP.div.S$lake.LCBD.q1 <- (TP.div.S$lake.LCBD.q1 * (n - 1) + 0.5) / n
TP.div.S$lake.LCBD.q2 <- (TP.div.S$lake.LCBD.q2 * (n - 1) + 0.5) / n


##building model
time.model <- brm(formula = lake.LCBD.q0 ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.S,
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

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_Phos_lake.LCBD.q0.S")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_lake.LCBD.q0.S")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_lake.LCBD.q0.S")

time.model <- readRDS(file = "time_model_Phos_lake.LCBD.q0.S")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_lake.LCBD.q0.S")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_lake.LCBD.q0.S")

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

pred.a$x <- pred.a$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)

png(file = "lakeID.RandomEffect_Total_Phosphorus_lake.LCBD.q0.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = lake.LCBD.q0, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 0)", 
                       title = "Time Model of Phosphorus & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot
dev.off()
ggpred_plot


######LCBD Lake q1####

##building model
time.model <- brm(formula = lake.LCBD.q1 ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.S,
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

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_Phos_lake.LCBD.q1.S")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_lake.LCBD.q1.S")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_lake.LCBD.q1.S")

time.model <- readRDS(file = "time_model_Phos_lake.LCBD.q1.S")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_lake.LCBD.q1.S")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_lake.LCBD.q1.S")

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

pred.a$x <- pred.a$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)

png(file = "lakeID.RandomEffect_Total_Phosphorus_lake.LCBD.q1.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = lake.LCBD.q1, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 1)", 
                       title = "Time Model of Phosphorus & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot
dev.off()
ggpred_plot

######LCBD Lake q2####

##building model
time.model <- brm(formula = lake.LCBD.q2 ~ Phosphorus_z + (1 + Phosphorus_z | lakeID),
                  data = TP.div.S,
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

pred.a <- ggpredict(time.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(time.model, terms = c("Phosphorus_z [all]", "lakeID"), type = "re", ppd = FALSE)

saveRDS(time.model, file = "time_model_Phos_lake.LCBD.q2.S")
saveRDS(pred.a, file = "pred.A_timemodel_Phos_lake.LCBD.q2.S")
saveRDS(pred.b, file = "pred.B_timemodel_Phos_lake.LCBD.q2.S")

time.model <- readRDS(file = "time_model_Phos_lake.LCBD.q2.S")
pred.a <- readRDS(file = "pred.A_timemodel_Phos_lake.LCBD.q2.S")
pred.b <- readRDS(file = "pred.B_timemodel_Phos_lake.LCBD.q2.S")

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

pred.a$x <- pred.a$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)

png(file = "lakeID.RandomEffect_Total_Phosphorus_lake.LCBD.q2.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = lake.LCBD.q2, colour = lakeID),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 2)", 
                       title = "Time Model of Phosphorus & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot
dev.off()
ggpred_plot


#####Phosphorus SPACE MODEL#####
######Richness######
##building model
TP.div.S$year <- as.factor(TP.div.S$year)
summary(TP.div.S)

space.model <- brm(formula = hill.rich ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.S,
                   family = Gamma(link = "log"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .82))#0.8 is default, determines how long it will "bounce", increasing can aid in convergence

#model checks
pairs(space.model)
summary(space.model)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)


pred.a <- ggpredict(space.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_rich.S")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_rich.S")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_rich.S")

space.model <- readRDS(file = "space_model_Phos_rich.S")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_rich.S")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_rich.S")

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

pred.a$x <- pred.a$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


png(file = "year.RandomEffect_Phosphorus_z_Rich.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.rich, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "Richness (as Hill numbers)", 
                       title = "Space Model of Phosphorus & Species Richness") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot

dev.off()
ggpred_plot



######Hill.Shannon####
hist(TP.div.S$hill.Shannon)
qqPlot(TP.div.S$hill.Shannon)


##building model
space.model <- brm(formula = hill.Shannon ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.S,
                   family = Gamma(link = "log"),
                   warmup = 1000,
                   iter = 3000,
                   chains = 4,
                   seed = 32,
                   control = list(adapt_delta = .8))
#model checks
pairs(space.model)
summary(space.model)
print(space.model, digits = 4)
pp_check(space.model, ndraws = 100)
plot(space.model)
diagnostic_posterior(space.model, diagnostic = "all", effects = c("fixed"))
bayes_R2(space.model)


pred.a <- ggpredict(space.model, terms = c("Phosphorus_z"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_shan.S")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_shan.S")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_shan.S")

space.model <- readRDS(file = "space_model_Phos_shan.S")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_shan.S")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_shan.S")


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

pred.a$x <- pred.a$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


png(file = "year.RandomEffect_Phosphorus_Shannon.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.Shannon, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x,predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "Shannon (as Hill numbers)", 
                       title = "Space Model of Phosphorus & Shannon") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot

dev.off()
ggpred_plot


######Inverse Simpson####

space.model <- brm(formula = hill.InverseSimpson ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.S,
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

pred.a <- ggpredict(space.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z [all]","year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_simp.S")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_simp.S")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_simp.S")

space.model <- readRDS(file = "space_model_Phos_simp.S")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_simp.S")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_simp.S")

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

pred.a$x <- pred.a$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)

png(file = "year.RandomEffect_Phosphorus_InverseSimpson.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = hill.InverseSimpson, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group),linewidth = 1) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  
                  labs(x = "Phosphorus (mg/l)", y = "hill.InverseSimpson", 
                       title = "Time Model of Phosphorus & Inverse Simpson") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot

dev.off()
ggpred_plot

######LCBD Year q0######
##building model
TP.div.S$year <- as.factor(TP.div.S$year)
summary(TP.div.S)

n <- nrow(TP.div.S)
TP.div.S$year.LCBD.q0 <- (TP.div.S$year.LCBD.q0 * (n - 1) + 0.5) / n
TP.div.S$year.LCBD.q1 <- (TP.div.S$year.LCBD.q1 * (n - 1) + 0.5) / n
TP.div.S$year.LCBD.q2 <- (TP.div.S$year.LCBD.q2 * (n - 1) + 0.5) / n


space.model <- brm(formula = year.LCBD.q0 ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.S,
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


pred.a <- ggpredict(space.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_year.LCBD.q0.S")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_year.LCBD.q0.S")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_year.LCBD.q0.S")

space.model <- readRDS(file = "space_model_Phos_year.LCBD.q0.S")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_year.LCBD.q0.S")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_year.LCBD.q0.S")

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

pred.a$x <- pred.a$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


png(file = "year.RandomEffect_Phosphorus_z_year.LCBD.q0.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = year.LCBD.q0, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 0)", 
                       title = "Space Model of Phosphorus & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot

dev.off()
ggpred_plot


######LCBD Year q1######
##building model
TP.div.S$year <- as.factor(TP.div.S$year)
summary(TP.div.S)

space.model <- brm(formula = year.LCBD.q1 ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.S,
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


pred.a <- ggpredict(space.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_year.LCBD.q1.S")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_year.LCBD.q1.S")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_year.LCBD.q1.S")

space.model <- readRDS(file = "space_model_Phos_year.LCBD.q1.S")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_year.LCBD.q1.S")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_year.LCBD.q1.S")

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

pred.a$x <- pred.a$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


png(file = "year.RandomEffect_Phosphorus_z_year.LCBD.q1.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = year.LCBD.q1, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 1)", 
                       title = "Space Model of Phosphorus & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot


######LCBD Year q2######
##building model
TP.div.S$year <- as.factor(TP.div.S$year)
summary(TP.div.S)

space.model <- brm(formula = year.LCBD.q2 ~ Phosphorus_z + (1 + Phosphorus_z | year),
                   data = TP.div.S,
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


pred.a <- ggpredict(space.model, terms = c("Phosphorus_z [all]"), type = "re", ppd = FALSE)
pred.b <- ggpredict(space.model, terms = c("Phosphorus_z [all]", "year"), type = "re", ppd = FALSE)

saveRDS(space.model, file = "space_model_Phos_year.LCBD.q2.S")
saveRDS(pred.a, file = "pred.A_spacemodel_Phos_year.LCBD.q2.S")
saveRDS(pred.b, file = "pred.B_spacemodel_Phos_year.LCBD.q2.S")

space.model <- readRDS(file = "space_model_Phos_year.LCBD.q2.S")
pred.a <- readRDS(file = "pred.A_spacemodel_Phos_year.LCBD.q2.S")
pred.b <- readRDS(file = "pred.B_spacemodel_Phos_year.LCBD.q2.S")

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

pred.a$x <- pred.a$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)
pred.b$x <- pred.b$x * sd(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.) + mean(TP.div.S$Total_Phorphorus_Mixed_Forms_mg.l.as.P.)


png(file = "year.RandomEffect_Phosphorus_z_year.LCBD.q2.S.png")
ggpred_plot <- (ggplot() + 
                  geom_ribbon(data = pred.a, 
                              aes(ymin = conf.low, ymax = conf.high, x = x), 
                              alpha = 0.075) +  # error band
                  geom_point(data = TP.div.S,                      # adding the raw data (scaled values)
                             aes(x = Total_Phorphorus_Mixed_Forms_mg.l.as.P., y = year.LCBD.q2, colour = year),
                             alpha = 0.1) +
                  geom_line(data = pred.b, aes(x, predicted, colour = group)) +      
                  geom_line(data = pred.a, aes(x, predicted), linewidth = 1.5, linetype = "longdash") +
                  labs(x = "Phosphorus (mg/l)", y = "LCBD (q = 2)", 
                       title = "Space Model of Phosphorus & Local Contribution of Beta Diversity") + 
                  theme_bw()) +
  guides(color = FALSE)
ggpred_plot


#####Comparing models####
######temperature Rich####
time.model <- readRDS(file = "time_model_temp_rich.rds")
space.model <- readRDS(file = "space_model_temp_rich.rds")

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



#########temperature Rich Random############
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

# Extract the random slopes for 'WaterTemp_z' by lakeID
random_slopes <- ranef_values$lakeID[, ,"WaterTemp_z" ]

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


png(file = "BOX.time.RandomEffect_WaterTemp_Rich.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Temp-Richness Time Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot


#####Just looking at the slope direction (positive/negative)
random_slopes_direction <- data.frame(
  Estimate = random_slopes[, "Estimate"],
  CI_Lower = random_slopes[, "Q2.5"],
  CI_Upper = random_slopes[, "Q97.5"],
  Classification = ifelse(random_slopes[, "Estimate"] > 0, "Positive",
                          ifelse(random_slopes[, "Estimate"] < 0, "Negative",
                                 "Zero Effect")))
# Define all possible categories
categories <- c("Negative", "Positive")

# Create the table, ensuring all categories are included
randefftable <- (table(factor(random_slopes_direction$Classification, levels = categories)))

# Print the result
print(randefftable)

randefftable_df <- as.data.frame(randefftable)
colnames(randefftable_df) <- c("Classification", "Count")

png(file = "BOX.direction.time.RandomEffect_WaterTemp_Rich.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Temp - Alpha Diversity (q = 0) Time Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "gray", "Negative" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot



#####Space Model
# Extract random effects
ranef_values <- ranef(space.model)  # Returns a list

# Check structure
str(ranef_values)

# Extract the random slopes for 'WaterTemp_z' by year
random_slopes <- ranef_values$year[, ,"WaterTemp_z"]

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


png(file = "BOX.space.RandomEffect_WaterTemp_Rich.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Temp-Richness Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot

#####Just looking at the slope direction (positive/negative)
random_slopes_direction <- data.frame(
  Estimate = random_slopes[, "Estimate"],
  CI_Lower = random_slopes[, "Q2.5"],
  CI_Upper = random_slopes[, "Q97.5"],
  Classification = ifelse(random_slopes[, "Estimate"] > 0, "Positive",
                          ifelse(random_slopes[, "Estimate"] < 0, "Negative",
                                 "Zero Effect")))
# Define all possible categories
categories <- c("Negative", "Positive")

# Create the table, ensuring all categories are included
randefftable <- (table(factor(random_slopes_direction$Classification, levels = categories)))

# Print the result
print(randefftable)

randefftable_df <- as.data.frame(randefftable)
colnames(randefftable_df) <- c("Classification", "Count")

png(file = "BOX.direction.space.RandomEffect_WaterTemp_Rich.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Temp - Alpha Diversity (q = 0) Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "gray", "Negative" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot



######temperature Shannon####
time.model <- readRDS(file = "time_model_temp_shan.rds")
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


#########temperature Shan Random############
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

# Extract the random slopes for 'WaterTemp_z' by lakeID
random_slopes <- ranef_values$lakeID[, ,"WaterTemp_z" ]

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


png(file = "BOX.time.RandomEffect_WaterTemp_Shan.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Temp-Shannon Time Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot


#####Just looking at the slope direction (positive/negative)
random_slopes_direction <- data.frame(
  Estimate = random_slopes[, "Estimate"],
  CI_Lower = random_slopes[, "Q2.5"],
  CI_Upper = random_slopes[, "Q97.5"],
  Classification = ifelse(random_slopes[, "Estimate"] > 0, "Positive",
                          ifelse(random_slopes[, "Estimate"] < 0, "Negative",
                                 "Zero Effect")))
# Define all possible categories
categories <- c("Negative", "Positive")

# Create the table, ensuring all categories are included
randefftable <- (table(factor(random_slopes_direction$Classification, levels = categories)))

# Print the result
print(randefftable)

randefftable_df <- as.data.frame(randefftable)
colnames(randefftable_df) <- c("Classification", "Count")

png(file = "BOX.direction.time.RandomEffect_WaterTemp_Shan.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Temp - Alpha Diversity (q = 1) Time Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "gray", "Negative" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot


###Space Model
# Extract random effects
ranef_values <- ranef(space.model)  # Returns a list

# Check structure
str(ranef_values)

# Extract the random slopes for 'WaterTemp_z' by year
random_slopes <- ranef_values$year[, ,"WaterTemp_z"]

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


png(file = "BOX.space.RandomEffect_WaterTemp_Shan.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Temp-Shannon Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot


#####Just looking at the slope direction (positive/negative)
random_slopes_direction <- data.frame(
  Estimate = random_slopes[, "Estimate"],
  CI_Lower = random_slopes[, "Q2.5"],
  CI_Upper = random_slopes[, "Q97.5"],
  Classification = ifelse(random_slopes[, "Estimate"] > 0, "Positive",
                          ifelse(random_slopes[, "Estimate"] < 0, "Negative",
                                 "Zero Effect")))
# Define all possible categories
categories <- c("Negative", "Positive")

# Create the table, ensuring all categories are included
randefftable <- (table(factor(random_slopes_direction$Classification, levels = categories)))

# Print the result
print(randefftable)

randefftable_df <- as.data.frame(randefftable)
colnames(randefftable_df) <- c("Classification", "Count")

png(file = "BOX.direction.space.RandomEffect_WaterTemp_Shan.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Temp - Alpha Diversity (q = 1) Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "gray", "Negative" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot

######temperature Inverse Simpson####
time.model <- readRDS(file = "time_model_temp_simp.rds")
space.model <- readRDS(file = "space_model_temp_simp.rds")

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

#########temperature Inverse Simpson Random############
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

# Extract the random slopes for 'WaterTemp_z' by lakeID
random_slopes <- ranef_values$lakeID[, ,"WaterTemp_z" ]

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


png(file = "BOX.time.RandomEffect_WaterTemp_Simp.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Temp-Inverse Simpson Time Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot

#####Just looking at the slope direction (positive/negative)
random_slopes_direction <- data.frame(
  Estimate = random_slopes[, "Estimate"],
  CI_Lower = random_slopes[, "Q2.5"],
  CI_Upper = random_slopes[, "Q97.5"],
  Classification = ifelse(random_slopes[, "Estimate"] > 0, "Positive",
                          ifelse(random_slopes[, "Estimate"] < 0, "Negative",
                                 "Zero Effect")))
# Define all possible categories
categories <- c("Negative", "Positive")

# Create the table, ensuring all categories are included
randefftable <- (table(factor(random_slopes_direction$Classification, levels = categories)))

# Print the result
print(randefftable)

randefftable_df <- as.data.frame(randefftable)
colnames(randefftable_df) <- c("Classification", "Count")

png(file = "BOX.direction.time.RandomEffect_WaterTemp_Simp.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Temp - Alpha Diversity (q = 2) Time Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "gray", "Negative" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot






###Space Model
# Extract random effects
ranef_values <- ranef(space.model)  # Returns a list

# Check structure
str(ranef_values)

# Extract the random slopes for 'WaterTemp_z' by year
random_slopes <- ranef_values$year[, ,"WaterTemp_z"]

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


png(file = "BOX.space.RandomEffect_WaterTemp_Simp.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Temp-Inverse Simpson Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot


#####Just looking at the slope direction (positive/negative)
random_slopes_direction <- data.frame(
  Estimate = random_slopes[, "Estimate"],
  CI_Lower = random_slopes[, "Q2.5"],
  CI_Upper = random_slopes[, "Q97.5"],
  Classification = ifelse(random_slopes[, "Estimate"] > 0, "Positive",
                          ifelse(random_slopes[, "Estimate"] < 0, "Negative",
                                 "Zero Effect")))
# Define all possible categories
categories <- c("Negative", "Positive")

# Create the table, ensuring all categories are included
randefftable <- (table(factor(random_slopes_direction$Classification, levels = categories)))

# Print the result
print(randefftable)

randefftable_df <- as.data.frame(randefftable)
colnames(randefftable_df) <- c("Classification", "Count")

png(file = "BOX.direction.space.RandomEffect_WaterTemp_Simp.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Temp - Alpha Diversity (q = 2) Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "gray", "Negative" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot

######temperature LCBD q = 0####
time.model <- readRDS(file = "time_model_temp_LCBD.L.q0.rds")
space.model <- readRDS(file = "space_model_temp_year.LCBD.q0.rds")

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

#########temperature LCBD q = 0 Random############
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

# Extract the random slopes for 'WaterTemp_z' by lakeID
random_slopes <- ranef_values$lakeID[, ,"WaterTemp_z" ]

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


png(file = "BOX.time.RandomEffect_WaterTemp_LCBD.q0.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Temp-LCBD (q = 0) Time Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot


#####Just looking at the slope direction (positive/negative)
random_slopes_direction <- data.frame(
  Estimate = random_slopes[, "Estimate"],
  CI_Lower = random_slopes[, "Q2.5"],
  CI_Upper = random_slopes[, "Q97.5"],
  Classification = ifelse(random_slopes[, "Estimate"] > 0, "Positive",
                          ifelse(random_slopes[, "Estimate"] < 0, "Negative",
                                 "Zero Effect")))
# Define all possible categories
categories <- c("Negative", "Positive")

# Create the table, ensuring all categories are included
randefftable <- (table(factor(random_slopes_direction$Classification, levels = categories)))

# Print the result
print(randefftable)

randefftable_df <- as.data.frame(randefftable)
colnames(randefftable_df) <- c("Classification", "Count")

png(file = "BOX.direction.time.RandomEffect_WaterTemp_LCBD.q0.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Temp - LCBD Diversity (q = 0) Time Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "gray", "Negative" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot

###Space Model
# Extract random effects
ranef_values <- ranef(space.model)  # Returns a list

# Check structure
str(ranef_values)

# Extract the random slopes for 'WaterTemp_z' by year
random_slopes <- ranef_values$year[, ,"WaterTemp_z"]

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


png(file = "BOX.space.RandomEffect_WaterTemp_LCBD.q0.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Temp-LCBD (q = 0) Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot

#####Just looking at the slope direction (positive/negative)
random_slopes_direction <- data.frame(
  Estimate = random_slopes[, "Estimate"],
  CI_Lower = random_slopes[, "Q2.5"],
  CI_Upper = random_slopes[, "Q97.5"],
  Classification = ifelse(random_slopes[, "Estimate"] > 0, "Positive",
                          ifelse(random_slopes[, "Estimate"] < 0, "Negative",
                                 "Zero Effect")))
# Define all possible categories
categories <- c("Negative", "Positive")

# Create the table, ensuring all categories are included
randefftable <- (table(factor(random_slopes_direction$Classification, levels = categories)))

# Print the result
print(randefftable)

randefftable_df <- as.data.frame(randefftable)
colnames(randefftable_df) <- c("Classification", "Count")

png(file = "BOX.direction.space.RandomEffect_WaterTemp_LCBD.q0.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Temp - LCBD (q = 0) Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "gray", "Negative" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot


######temperature LCBD q = 1####
time.model <- readRDS(file = "time_model_temp_LCBD.L.q1.rds")
space.model <- readRDS(file = "space_model_temp_year.LCBD.q1.rds")

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

#########temperature LCBD q = 1 Random############
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

# Extract the random slopes for 'WaterTemp_z' by lakeID
random_slopes <- ranef_values$lakeID[, ,"WaterTemp_z" ]

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


png(file = "BOX.time.RandomEffect_WaterTemp_LCBD.q1.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Temp-LCBD (q = 1) Time Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot


#####Just looking at the slope direction (positive/negative)
random_slopes_direction <- data.frame(
  Estimate = random_slopes[, "Estimate"],
  CI_Lower = random_slopes[, "Q2.5"],
  CI_Upper = random_slopes[, "Q97.5"],
  Classification = ifelse(random_slopes[, "Estimate"] > 0, "Positive",
                          ifelse(random_slopes[, "Estimate"] < 0, "Negative",
                                 "Zero Effect")))
# Define all possible categories
categories <- c("Negative", "Positive")

# Create the table, ensuring all categories are included
randefftable <- (table(factor(random_slopes_direction$Classification, levels = categories)))

# Print the result
print(randefftable)

randefftable_df <- as.data.frame(randefftable)
colnames(randefftable_df) <- c("Classification", "Count")

png(file = "BOX.direction.time.RandomEffect_WaterTemp_LCBD.q1.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Temp - LCBD Diversity (q = 1) Time Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "gray", "Negative" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot


###Space Model
# Extract random effects
ranef_values <- ranef(space.model)  # Returns a list

# Check structure
str(ranef_values)

# Extract the random slopes for 'WaterTemp_z' by year
random_slopes <- ranef_values$year[, ,"WaterTemp_z"]

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


png(file = "BOX.space.RandomEffect_WaterTemp_LCBD.q1.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Temp-LCBD (q = 1) Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot


#####Just looking at the slope direction (positive/negative)
random_slopes_direction <- data.frame(
  Estimate = random_slopes[, "Estimate"],
  CI_Lower = random_slopes[, "Q2.5"],
  CI_Upper = random_slopes[, "Q97.5"],
  Classification = ifelse(random_slopes[, "Estimate"] > 0, "Positive",
                          ifelse(random_slopes[, "Estimate"] < 0, "Negative",
                                 "Zero Effect")))
# Define all possible categories
categories <- c("Negative", "Positive")

# Create the table, ensuring all categories are included
randefftable <- (table(factor(random_slopes_direction$Classification, levels = categories)))

# Print the result
print(randefftable)

randefftable_df <- as.data.frame(randefftable)
colnames(randefftable_df) <- c("Classification", "Count")

png(file = "BOX.direction.space.RandomEffect_WaterTemp_LCBD.q1.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Temp - LCBD Diversity (q = 1) Spcae Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Positive" = "gray", "Negative" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot

######temperature LCBD q = 2####
time.model <- readRDS(file = "time_model_temp_LCBD.L.q2.rds")
space.model <- readRDS(file = "space_model_temp_year.LCBD.q2.rds")

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

#########temperature LCBD q = 2 Random############
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

# Extract the random slopes for 'WaterTemp_z' by lakeID
random_slopes <- ranef_values$lakeID[, ,"WaterTemp_z" ]

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


png(file = "BOX.time.RandomEffect_WaterTemp_LCBD.q2.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Temp-LCBD (q = 2) Time Model") + 
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

# Extract the random slopes for 'WaterTemp_z' by year
random_slopes <- ranef_values$year[, ,"WaterTemp_z"]

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


png(file = "BOX.space.RandomEffect_WaterTemp_LCBD.q2.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Temp-LCBD (q = 2) Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot




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

#########Total Nitrogen Rich Random############
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

# Extract the random slopes for 'TN_z' by lakeID
random_slopes <- ranef_values$lakeID[, ,"TN_z" ]

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


png(file = "BOX.time.RandomEffect_TN_Rich.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Nitrogen-Richness Time Model") + 
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

# Extract the random slopes for 'WaterTemp_z' by year
random_slopes <- ranef_values$year[, ,"TN_z"]

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


png(file = "BOX.space.RandomEffect_Total Nitrogen_Rich.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Nitrogen-Richness Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot



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

#########Total Nitrogen Shannon Random############
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

# Extract the random slopes for 'NitrateNitrite_z' by lakeID
random_slopes <- ranef_values$lakeID[, ,"TN_z" ]

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


png(file = "BOX.time.RandomEffect_TotalNitrogen_Shan.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Nitrogen-Shannon Time Model") + 
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

# Extract the random slopes for 'NitrateNitrite_z' by year
random_slopes <- ranef_values$year[, ,"TN_z"]

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


png(file = "BOX.space.RandomEffect_TotalNitrogen_Shan.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Nitrogen-Shannon Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot


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

#########Total Nitrogen Inverse Simpson Random############
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

# Extract the random slopes for 'NitrateNitrite_z' by lakeID
random_slopes <- ranef_values$lakeID[, ,"TN_z" ]

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


png(file = "BOX.time.RandomEffect_TotalNitrogen_Simp.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Nitrogen-Inverse Simpson Time Model") + 
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

# Extract the random slopes for 'NitrateNitrite_z' by year
random_slopes <- ranef_values$year[, ,"TN_z"]

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


png(file = "BOX.space.RandomEffect_TotalNitrogen_Simp.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Nitrogen-Inverse Simpson Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot



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

#########TotalNitrogen LCBD q = 0 Random############
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

# Extract the random slopes for 'NitrateNitrite_z' by lakeID
random_slopes <- ranef_values$lakeID[, ,"TN_z" ]

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


png(file = "BOX.time.RandomEffect_TotalNitrogen_LCBD.q0.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Nitrogen-LCBD (q = 0) Time Model") + 
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

# Extract the random slopes for 'NitrateNitrite_z' by year
random_slopes <- ranef_values$year[, ,"TN_z"]

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


png(file = "BOX.space.RandomEffect_Total Nitrogen_LCBD.q0.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Nitrogen-LCBD (q = 0) Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot


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

#########Total Nitrogen LCBD q = 1 Random############
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

# Extract the random slopes for 'NitrateNitrite_z' by lakeID
random_slopes <- ranef_values$lakeID[, ,"NitrateNitrite_z" ]

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


png(file = "BOX.time.RandomEffect_NitrateNitrite_LCBD.q1.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Nitrogen-LCBD (q = 1) Time Model") + 
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

# Extract the random slopes for 'NitrateNitrite_z' by year
random_slopes <- ranef_values$year[, ,"NitrateNitrite_z"]

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


png(file = "BOX.space.RandomEffect_NitrateNitrite_LCBD.q1.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Nitrogen-LCBD (q = 1) Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot


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

#########Total nitrogen LCBD q = 2 Random############
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

# Extract the random slopes for 'NitrateNitrite_z' by lakeID
random_slopes <- ranef_values$lakeID[, ,"NitrateNitrite_z" ]

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


png(file = "BOX.time.RandomEffect_NitrateNitrite_LCBD.q2.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Nitrogen-LCBD (q = 2) Time Model") + 
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

# Extract the random slopes for 'NitrateNitrite_z' by year
random_slopes <- ranef_values$year[, ,"NitrateNitrite_z"]

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


png(file = "BOX.space.RandomEffect_NitrateNitrite_LCBD.q2.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Nitrogen-LCBD (q = 2) Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot




######Total Phosphorus Rich####
time.model <- readRDS(file = "time_model_Phos_rich.rds")
space.model <- readRDS(file = "space_model_Phos_rich.rds")

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

#########Total Phosphorus Rich Random############
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


png(file = "BOX.time.RandomEffect_Phosphorus_Rich.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Phosphorus-Richness Time Model") + 
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


png(file = "BOX.space.RandomEffect_Phosphorus_Rich.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Phosphorus-Richness Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot


######Total Phosphorus Shannon####
time.model <- readRDS(file = "time_model_Phos_shan.rds")
space.model <- readRDS(file = "space_model_Phos_shan.rds")

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

#########Total Phosphorus Shan Random############
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


png(file = "BOX.time.RandomEffect_Phosphorus_Shan.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Phosphorus-Shannon Time Model") + 
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


png(file = "BOX.space.RandomEffect_Phosphorus_Shan.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Phosphorus-Shannon Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot


######Total Phosphorus Inverse Simpson####
time.model <- readRDS(file = "time_model_Phos_simp.rds")
space.model <- readRDS(file = "space_model_Phos_simp.rds")

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

#########Total Phosphorus Inverse Simpson Random############
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


png(file = "BOX.time.RandomEffect_Phosphorus_Simp.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Phosphorus-Inverse Simpson Time Model") + 
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


png(file = "BOX.space.RandomEffect_Phosphorus_Simp.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Phosphorus-Inverse Simpson Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot


######Total Phosphorus LCBD q = 0 ####
time.model <- readRDS(file = "time_model_Phos_lake.LCBD.q0.rds")
space.model <- readRDS(file = "space_model_Phos_year.LCBD.q0.rds")

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

#########Total Phosphorus LCBD q = 0 Random############
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


png(file = "BOX.time.RandomEffect_Phosphorus_LCBD.q0.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Phosphorus-LCBD (q = 0) Time Model") + 
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


png(file = "BOX.space.RandomEffect_Phosphorus_LCBD.q0.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Phosphorus-LCBD (q = 0) Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot


######Total Phosphorus LCBD q = 1 ####
time.model <- readRDS(file = "time_model_Phos_lake.LCBD.q1.rds")
space.model <- readRDS(file = "space_model_Phos_year.LCBD.q1.rds")

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

#########Total Phosphorus LCBD q = 1 Random############
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


png(file = "BOX.time.RandomEffect_Phosphorus_LCBD.q1.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of Lakes", 
       title = "Direction of Random Effects of Phosphorus-LCBD (q = 1) Time Model") + 
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


png(file = "BOX.space.RandomEffect_Phosphorus_LCBD.q1.png")

ggpred_plot <- ggplot(randefftable_df, aes(x = Classification, y = Count, fill = Classification)) + 
  geom_col() +  # Use geom_col() for pre-counted data
  labs(x = "Random Effect Classification", 
       y = "Number of years", 
       title = "Direction of Random Effects of Phosphorus-LCBD (q = 1) Space Model") + 
  theme_bw() +
  scale_fill_manual(values = c("Both Positive" = "gray", "Both Negative" = "gray", "Contains Zero" = "gray")) +
  theme(legend.position = "none")  # Hide legend since colors are self-explanatory

ggpred_plot

dev.off()
ggpred_plot



######Total Phosphorus LCBD q = 2 ####
time.model <- readRDS(file = "time_model_Phos_lake.LCBD.q2.rds")
space.model <- readRDS(file = "space_model_Phos_year.LCBD.q2.rds")

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

TempPlot.A.q0 <- TempPlot.A.q0 + coord_cartesian(xlim = shared_xlim) + theme_no_x +
  theme(
    panel.border = element_rect(color = "darkviolet", fill = NA, linewidth = 2))
TempPlot.A.q1 <- TempPlot.A.q1 + coord_cartesian(xlim = shared_xlim) + theme_no_x
TempPlot.A.q2 <- TempPlot.A.q2 + coord_cartesian(xlim = shared_xlim) + theme_no_x
TempPlot.B.q0 <- TempPlot.B.q0 + coord_cartesian(xlim = shared_xlim) + theme_no_x+
  theme(
    panel.border = element_rect(color = "darkviolet", fill = NA, linewidth = 2))
TempPlot.B.q1 <- TempPlot.B.q1 + coord_cartesian(xlim = shared_xlim) + theme_no_x+
  theme(
    panel.border = element_rect(color = "darkviolet", fill = NA, linewidth = 2))
TempPlot.B.q2 <- TempPlot.B.q2 + coord_cartesian(xlim = shared_xlim) + theme_with_x+
  theme(
    panel.border = element_rect(color = "darkviolet", fill = NA, linewidth = 2))
TempPlot.A.q0.S <- TempPlot.A.q0.S + coord_cartesian(xlim = shared_xlim) + theme_no_x_or_y +
  theme(
    panel.border = element_rect(color = "darkviolet", fill = NA, linewidth = 2))
TempPlot.A.q1.S <- TempPlot.A.q1.S + coord_cartesian(xlim = shared_xlim) + theme_no_x_or_y+
  theme(
    panel.border = element_rect(color = "darkviolet", fill = NA, linewidth = 2))
TempPlot.A.q2.S <- TempPlot.A.q2.S + coord_cartesian(xlim = shared_xlim) + theme_no_x_or_y+
  theme(
    panel.border = element_rect(color = "darkviolet", fill = NA, linewidth = 2))
TempPlot.B.q0.S <- TempPlot.B.q0.S + coord_cartesian(xlim = shared_xlim) + theme_no_x_or_y
TempPlot.B.q1.S <- TempPlot.B.q1.S + coord_cartesian(xlim = shared_xlim) + theme_no_x_or_y+
  theme(
    panel.border = element_rect(color = "darkviolet", fill = NA, linewidth = 2))
TempPlot.B.q2.S <- TempPlot.B.q2.S + coord_cartesian(xlim = shared_xlim) + theme_with_x_but_no_y


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

