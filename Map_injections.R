#Map injections rawdata

# ---
# This script has been modified from https://github.com/MCabreraBrufau/Licor_N2O_scripts to identify and integrate peak not only for N2O but also for CO2 and CH4
# ---

#Description: This script creates a map_injection csv for every rawdata file. These are stored in the Map_injections folder and will be manually edited to correct label text, adapt time of labels or add labels that were not written in the data at the time of collection. 

# ---- Directories ----

#Root
#Usually you will be working on your working directory
# folder_root <- dirname(rstudioapi::getSourceEditorContext()$path)
#But you can set the folder in other path
folder_root <- "/home/jorge/Documentos/Postdoctoral/Onedrive_UB/UB/NaturBPond/Hidrogeologos/Muestras_Licor" # You have to make sure this is pointing to the write folder on your local machine

#Data folders, 
folder_raw <- paste0(folder_root,"/Rawdata") #contains unedited files downloaded from licor

#First we check if the folder exist and, if not, create one
folder_mapinjections<- paste0(folder_root,"/Map_injections") #Where csv files with the injection details will be saved. 
if (!dir.exists(folder_mapinjections)) {
  # If it doesn't exist, create the folder
  dir.create(folder_mapinjections)
}

# ---- packages & functions ----
library(tidyverse)

#Import functions of repo 
repo_root <- dirname(dirname(rstudioapi::getSourceEditorContext()$path))
repo_root <- getwd()
files.sources = list.files(path = paste0(repo_root,"/functions"), full.names = T)
for (f in files.sources){source(f)}


#Check rawdata files and map_injection files

#raw_files are named as: paste("TG20-01377-", year-month-day, "-T",hourstartexport, ".data")
#maps_done will be named as raw_map_injection_rawfile.csv, the rawfile part will not end in .data

#List maps that are already created in folder_mapinjections
maps_done<- list.files(path=folder_mapinjections, pattern = "^raw_.*_map_injection")

#List raw files (for Li-7820 and Li-7810) present in folder_raw
raw_files<- list.files(path = folder_raw, pattern =  "^TG")
#Aquí debería revisar todos los archivos licor (TG10 y TG20) e identificar que gas es cada uno


#Get raw files without corresponding map injection
raw_files_withoutmap<- raw_files[!gsub(".data", "",raw_files)%in%gsub(".csv", "", gsub(pattern = "raw_.*_map_injection_","",maps_done))]


#Collect Tstart Tend and labels for all unique remarks of every raw_file_withoutmap
#Save these details in csv files named "raw_map_injection_"[rawfilename without".data"].csv 
for (i in raw_files_withoutmap){
  #Here we check if the raw file is for CO2 and CH4 or for N2O
  if (grepl("TG10", i)) {
    a<- read_Licor_TG10(paste0(folder_raw,"/",i))
    gas <- "CO2_and_CH4"
  }
  if (grepl("TG20", i)) {
    a<- read_Licor_TG20(paste0(folder_raw,"/",i))
    gas <- "N2O"
  }
  a <- a %>% 
    group_by(label) %>% 
    summarise(date=first(date),Tstart=first(UTCtime), Tend =last(UTCtime)) %>% 
    arrange(Tstart) %>% 
    mutate(rawfile=i) %>% 
    select(date, Tstart, Tend, label, rawfile) %>% 
    mutate(Tstart_correct=NA,	Tend_correct=NA,	label_correct=NA)#Add empty columns to manually correct the data
  
  #Estaría bien añadir una marca al nombre del archivo para saber si es de N2O o de CO2 aunque en realidad ya estaría con el nombre del rawfile (TG10 o TG20)
  write.csv(a,file = paste0(folder_mapinjections,"/raw_", gas, "_map_injection_", gsub(".data","",i), ".csv"),row.names = F)
}

rm(a,i, raw_files,raw_files_withoutmap,maps_done)
