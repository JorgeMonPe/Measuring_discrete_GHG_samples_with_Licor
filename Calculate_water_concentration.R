library(tidyverse)
library(errors)
source("https://raw.github.com/JorgeMonPe/Functions_headspace_calculation/main/Functions_gas_concentration.R")

ppmfiles<- list.files(path = "Results_ppm", pattern = "^.*ppm_samples_")
ppmdata <- data.frame()
for(i in ppmfiles){
  data <- read_csv(paste0("Results_ppm/", i))
  data <- data %>% select(sample, peak_id, matches("^.*ppm")) %>% pivot_longer(-c(sample, peak_id), names_to = "Species", values_to =  "ppm")
  ppmdata <- bind_rows(ppmdata, data)
}

##Select 3 samples with the lower CV for each volume-dilution treatment
# Function to find the best combination for a group
select_lowest_cv <- function(measurements) {
  # Generate all combinations of 3 measurements
  combos <- combn(measurements, 3, simplify = FALSE)
  
  # Calculate the CV for each combination
  cvs <- sapply(combos, function(x) sd(x) / mean(x) * 100)
  
  # Select the combination with the lowest CV
  best_combo <- combos[[which.min(cvs)]]
  return(best_combo)
}
#Create a df with the combination selected
result <- ppmdata %>%
  group_by(sample, Species) %>%
  reframe(
    ppm = select_lowest_cv(ppm)
  ) %>% 
  mutate(Selected = "Yes") #This is to track the best combination after join
#Join with data
ppmdata <- ppmdata %>% left_join(result)
ppmdata <- ppmdata %>% 
  filter(Selected == "Yes")
#Split sample name into Pond and Exetainer ID
ppmdata <- ppmdata %>% 
  filter(!str_starts(sample, "6ppm")) %>%
  separate_wider_delim(sample, delim = "-", names = c("Pond", "Exetainer_ID")) %>% 
  mutate(Exetainer_ID = as.numeric(Exetainer_ID))
#Import field data
#Wate temp
Twater <- read_csv("/home/jorge/Documentos/Postdoctoral/Onedrive_UB/UB/NaturBPond/GHG/Pond_element_flux/Raw_data_to_final_matrix/December/Hobos/temp_headspace_incubations.csv")
Twater <- Twater %>%
  mutate(Datetime = strptime(Datetime, format = "%d/%m/%Y %H:%M:%S", tz = "UTC")) %>% 
  mutate(Datetime = as.POSIXct(Datetime)) %>% 
  select(Datetime, Temp) %>%
  mutate(Datetime = update(Datetime, year = year(Datetime) + 2000))

#Incubation info
Bot <- read_csv("/home/jorge/Documentos/Postdoctoral/Onedrive_UB/UB/NaturBPond/GHG/Pond_element_flux/Raw_data_to_final_matrix/December/Estadillos_campo/AUX_complete_Botanic.csv")
Bot <- Bot %>% mutate(Pond = "Botanic")
Ene <- read_csv("/home/jorge/Documentos/Postdoctoral/Onedrive_UB/UB/NaturBPond/GHG/Pond_element_flux/Raw_data_to_final_matrix/December/Estadillos_campo/AUX_complete_Enees.csv")
Ene <- Ene %>% mutate(Pond = "Enees")
Bot <- Bot %>% mutate(Final_CH4 = as.numeric(Final_CH4),
                      Final_CO2 = as.numeric(Final_CO2),
                      Initial_N2O = as.numeric(Initial_N2O))
Field_data <- Ene %>% bind_rows(Bot)

#Exetainer fieldsheet
Exetainer <- read_csv("/home/jorge/Documentos/Postdoctoral/Onedrive_UB/UB/NaturBPond/GHG/Pond_element_flux/Raw_data_to_final_matrix/December/Estadillos_campo/Template_exetainers.csv")

#Join with exetainer fieldsheet to get the Plot information
#Join with field data to get the time of incubaction of each exetainer
data <- Exetainer %>% left_join(ppmdata)
data <- data %>% select(-Strata) %>% left_join(Field_data)

#Now, from the hobo, I extract the temperature incubation for each exetainer
# Function to find the first datetime after the target datetime
find_first_after <- function(x, hobo_data){
  #x <- field %>% filter(Plot_ID == 1)
  Start_time <- x$start.time
  result <- hobo_data %>%
    filter(Datetime > Start_time) %>% #The name of the column with POSIXct must be Datetime
    slice_min(Datetime, with_ties = FALSE)
}
Temp_sel <- Field_data %>%
  group_by(Pond, Plot_ID) %>%
  group_modify(~find_first_after(.x, hobo_data = Twater)) %>% 
  ungroup()
data <- Temp_sel %>% left_join(data)
#Add air and water volumes
data <- data %>% mutate(Vol_H2O = 30, Vol_air = 30)

data <- data %>% select(Plot_ID, Strata, Exetainer_ID, Type, Pond, peak_id, Species, ppm, Vol_H2O, Vol_air, Temp)

data <- data %>% pivot_wider(names_from = "Species", values_from = "ppm")

#Concentration in the air use for the headspace
# CH4atm_ppmv <-data %>% filter(Type == "A") %>% select(CH4_ppm) %>% pull() %>% {set_errors(mean(., na.rm = T), sd(., na.rm = T))}
# CO2atm_ppmv <-data %>% filter(Type == "A") %>% select(CO2_ppm) %>% pull() %>% {set_errors(mean(., na.rm = T), sd(., na.rm = T))}
# N2Oatm_ppmv <-data %>% filter(Type == "A") %>% select(N2O_ppm) %>% pull() %>% {set_errors(mean(., na.rm = T), sd(., na.rm = T))}

Atm_ppmv <-data %>% group_by(Pond) %>% filter(Type == "A") %>% summarise(CH4_aire = set_errors(mean(CH4_ppm, na.rm = T), sd(CH4_ppm, na.rm = T)),
                                                                         CO2_aire = set_errors(mean(CO2_ppm, na.rm = T), sd(CO2_ppm, na.rm = T)),
                                                                         N2O_aire = set_errors(mean(N2O_ppm, na.rm = T), sd(N2O_ppm, na.rm = T)))
data <- data %>% left_join(Atm_ppmv)

#Caculate concentration uM in the water
data <- data %>% filter(Type != "A") %>% mutate(CH4_water_uM = nGHG_water_uM("CH4", CH4_ppm, Vol_H2O, Vol_air, T_Celsius = Temp, Patm_eq = 1, R = 0.08206, GHG_atm_ppmv = CH4_aire))

data <- data %>% filter(Type != "A") %>% mutate(N2O_water_uM = nGHG_water_uM("N2O", N2O_ppm, Vol_H2O, Vol_air, T_Celsius = Temp, Patm_eq = 1, R = 0.08206, GHG_atm_ppmv = N2O_aire))

data <- data %>% filter(Type != "A") %>% mutate(CO2_water_uM = nGHG_water_uM("CO2", CO2_ppm, Vol_H2O, Vol_air, T_Celsius = Temp, Patm_eq = 1, R = 0.08206, GHG_atm_ppmv = CO2_aire))

#Hallar la concentracion en ppmv, creo que solo hay que dividir la concentracion en uM entre la kh para cara especie

data <- data %>% mutate(Kh_CH4 = Kh_CH4(Temp), Kh_CO2 = Kh_CO2(Temp), Kh_N2O = Kh_N2O(Temp))

#Esto lo hago según Ray et al., 2023.
data <- data %>% mutate(CH4_water_uatm = CH4_water_uM/Kh_CH4, CO2_water_uatm = CO2_water_uM/Kh_CO2, N2O_water_uatm = N2O_water_uM/Kh_N2O)

#Expresarlo en mg/L
data <- data %>% mutate(CH4_water_mgL = CH4_water_uM*16/1000, CO2_water_mgL = CO2_water_uM*44/1000, N2O_water_mgL = N2O_water_uM*44/1000)

#Expresarlo en mg de C /L
data <- data %>% mutate(CH4_water_mgCL = CH4_water_uM*12/1000, CO2_water_mgCL = CO2_water_uM*12/1000, N2O_water_mgCL = N2O_water_uM*12/1000)

#Los valores negativos entiendo que no son válidos y pongo NA
data %>% filter(CO2_water_uatm < 0)
data %>% filter(CH4_water_uatm < 0)
data %>% filter(N2O_water_uatm < 0)

#uM
data_export_uM <- data %>% select(Pond, Plot_ID, Type, Strata, CH4_water_uM, CO2_water_uM, N2O_water_uM, CO2_aire, CH4_aire, N2O_aire)
write_csv(data_export_uM, "Results_water/GHGs_concentration_uM.csv")
saveRDS(data_export_uM, "Results_water/GHGs_concentration_uM.RDS")
#uatm
data_export_uatm <- data %>% select(Pond, Plot_ID, Type, Strata, CH4_water_uatm, CO2_water_uatm, N2O_water_uatm, CO2_aire, CH4_aire, N2O_aire)
write_csv(data_export_uatm, "Results_water/GHGs_concentration_uatm.csv")
saveRDS(data_export_uatm, "Results_water/GHGs_concentration_uatm.RDS")
#mgL
data_export_mgL <- data %>% select(Pond, Plot_ID, Type, Strata, CH4_water_mgL, CO2_water_mgL, N2O_water_mgL, CO2_aire, CH4_aire, N2O_aire)
write_csv(data_export_mgL, "Results_water/GHGs_concentration_mgL.csv")
saveRDS(data_export_mgL, "Results_water/GHGs_concentration_mgL.RDS")

#mgCL
data_export_mgCL <- data %>% select(Pond, Plot_ID, Type, Strata, CH4_water_mgCL, CO2_water_mgCL, N2O_water_mgCL, CO2_aire, CH4_aire, N2O_aire)
write_csv(data_export_mgCL, "GHGs_concentration_mgCL.csv")
saveRDS(data_export_mgCL, "GHGs_concentration_mgCL.RDS")

