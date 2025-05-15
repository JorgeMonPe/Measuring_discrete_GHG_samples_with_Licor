#Library----
library(tidyverse)
library(errors)
source("https://raw.github.com/JorgeMonPe/Functions_headspace_calculation/main/Functions_gas_concentration.R")

#One extra function----
  #Select 3 injections with the lower CV for each volume-dilution treatment----
  #This step can be skip
  ##Function to find the best combination for a group----
  select_lowest_cv <- function(measurements) {
    # Generate all combinations of 3 measurements
    combos <- combn(measurements, 3, simplify = FALSE)
    
    # Calculate the CV for each combination
    cvs <- sapply(combos, function(x) sd(x) / mean(x) * 100)
    
    # Select the combination with the lowest CV
    best_combo <- combos[[which.min(cvs)]]
    return(best_combo)
  }

# ---- Directory ----
#Root
folder_root <- "/home/jorge/Documentos/Postdoctoral/Onedrive_UB/UB/NaturBPond/Hidrogeologos/Muestras_Licor" # You have to make sure this is pointing to the write folder on your local machine

# Read all files with concentration calculated in ppm-----
ppmfiles<- list.files(path = paste0(folder_root, "/Results_ppm"), pattern = "^.*ppm_samples_")
ppmdata <- data.frame()
for(i in ppmfiles){
  data <- read_csv(paste0(folder_root, "/Results_ppm/", i))
  data <- data %>% select(sample, peak_id, matches("^.*ppm")) %>% pivot_longer(-c(sample, peak_id), names_to = "Species", values_to =  "ppm")
  ppmdata <- bind_rows(ppmdata, data)
}

#Auxiliary matrix----
ppmdata <- ppmdata %>%
  separate_wider_delim(sample, delim = "-", names = c("Sample", "Replicate"))
#Add sample type (headspace or air)
ppmdata <- ppmdata %>% mutate(Sample_type = case_when(Replicate == "A" ~ "Air_reference",
                                                      TRUE ~ "Headspace"))
ppmdata <- ppmdata %>% mutate(Vol_H2O = 30, Vol_air = 30, T_Celsius = 25, Patm_eq = 1)

#Set the concentration of the air used for the headspace equilibration----
  ##Samples with air used in the headspaces----
  #Do you have samples which are the air reference used in the headspace equilibration?
  #If you took sample of air used in the headspace, you must indicate this in the column `Sample_type` of the auxiliary data frame
  #Air reference - Concentration in the air use for the headspace----
  #Create a df with just the selected injections 
  injections_selected <- ppmdata %>%
    filter(Sample_type == "Air_reference") %>% 
    group_by(Sample, Species) %>%
    reframe(
      ppm = select_lowest_cv(ppm)
    ) %>% 
    mutate(Selected = "Yes") %>%  #This is to track the best combination after join
    ungroup()
  
  #Mean with sd?
  #AirRef_ppmv <- injections_selected %>% group_by(Species) %>% summarise(Air_reference_ppm = set_errors(mean(ppm, na.rm = T), sd(ppm, na.rm = T)))

  #One atm reference per sample
  AirRef_ppmv <- injections_selected %>% group_by(Sample, Species) %>% summarise(Air_reference_ppm = mean(ppm, na.rm = T))

  #Median?
  #AirRef_ppmv <- injections_selected %>% group_by(Species) %>% summarise(Air_reference_ppm = median(ppm))
  ppmdata <- ppmdata %>% left_join(AirRef_ppmv)

##Fix values extracted from other source (air pure bottle, NOAA value references, etc.)
  # ppmdata <- ppmdata %>% mutate(Air_reference_ppm = case_when(Species == "CO2_ppm" ~ 420,
  #                                                             Species == "CH4_ppm" ~ 2.6,
  #                                                             Species == "N2O_ppm" ~ 0.35))
#Concentration in uM in the water----
data <- ppmdata %>% filter(Sample_type == "Headspace") %>% mutate(Conc_water_uM = case_when(Species == "CO2_ppm" ~ nGHG_water_uM("CO2", ppm, Vol_H2O, Vol_air, T_Celsius, Patm_eq, R = 0.08206, GHG_atm_ppmv = Air_reference_ppm),
                                                                                            Species == "CH4_ppm" ~ nGHG_water_uM("CH4", ppm, Vol_H2O, Vol_air, T_Celsius, Patm_eq, R = 0.08206, GHG_atm_ppmv = Air_reference_ppm),
                                                                                            Species == "N2O_ppm" ~ nGHG_water_uM("N2O", ppm, Vol_H2O, Vol_air, T_Celsius, Patm_eq, R = 0.08206, GHG_atm_ppmv = Air_reference_ppm)))
#Partial pressure in uatm----
#First we calculate solubility constant for each gas
data <- data %>% mutate(Kh = case_when(Species == "CO2_ppm" ~ Kh_CO2(T_Celsius),
                                       Species == "CH4_ppm" ~ Kh_CH4(T_Celsius),
                                       Species == "N2O_ppm" ~ Kh_N2O(T_Celsius)))

#Esto lo hago según Ray et al., 2023.
data <- data %>% mutate(PartialPressure_uatm = Conc_water_uM/Kh)

#Expresarlo en mg/L
data <- data %>% mutate(Conc_water_mgL = case_when(Species == "CO2_ppm" ~ Conc_water_uM*16/1000,
                                       Species == "CH4_ppm" ~ Conc_water_uM*44/1000,
                                       Species == "N2O_ppm" ~ Conc_water_uM*44/1000)) 
#Export results----
#Folder for results
folder_results<- paste0(folder_root,"/Results_headspaces")#One csv per dayofinjections will be created (auto-name from rawfile name), with individual peak parameters (label, peak_id, peaksum, peakmax, unixtime_ofmax, raw_peaksum, dayofanalysis, SNR)
if (!dir.exists(folder_results)) {
  # If it doesn't exist, create the folder
  dir.create(folder_results)
}

write_csv(data, paste0(folder_results, "/water_concentration_samples_headspaces.csv"))
