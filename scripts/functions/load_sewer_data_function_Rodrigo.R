#  Load data  #

## Goal 
  ### Order metadata 
  ### Make new dataframes
  ### Load old dataframes for fast analysis

## Output 
 ### Dataset for analysis
 ### Dataset for control 
 ### T/F for generaltion of new



## Define function ad varibles for making new

make_sewer_data <- function(
    save_unclassified = F,
    new_long_df = T,
    new_core_groups = T,
    new_control_df = T,
    make_usearch_plot = F, 
    WD = "C:/Users/HD95LP/OneDrive - Aalborg Universitet/visual_studio_code/2023_sewer_microbial_communities/",
    DataPath = paste0("C:/Users/HD95LP/OneDrive - Aalborg Universitet/visual_studio_code/2023_sewer_microbial_communities/", "data/"),
    OutputPath = paste0("C:/Users/HD95LP/OneDrive - Aalborg Universitet/visual_studio_code/2023_sewer_microbial_communities/", "output/"),
    sintax_file = "ASVs_0.6_cut_off_MIDAS_5.2.R1.sintax"
){
  

source(paste0(WD, "../2023_RETHINK_sample_collection/scripts/functions/change_tax_names.R"))   # change to midas 5.1 names (2023-11-03)
  
  
# Metadata parameters
print("Import of rain and reordering rain and other envionmental parameters")  
  
# Summarise hourly rain to only n hour before sampling
n_hours_before = "1899-12-31 01:00:00"
print(paste0("For the daily rain, the rain is summrised for the ", hours(n_hours_before), 
             " hours before sampling. Meaning if hours before is 1 and the sampling was at 11.30, the rain is summerised for the hours 10 and 11"))
# Define rain events
rain_event_mm = 5
print(paste0("For calculation of *days_since_rain_event* and *time_before_rain_event*, a rain event is defined as a day with >= ", 
                          rain_event_mm, "mm rain during a day"))

# Earth temperature in 1 metr depth (https://www.bolius.dk/se-jordtemperaturen-i-danmark-maaned-for-maaned-35105)
et <- tibble(Month = seq_len(12), 
             Earth_temp_1m = c(3.6, 3.5, 2.9, 4.5, 8.8, 11.3, 15.1, 15.6, 14.7, 12.0, 8.9, 5.9))


# 
source(paste0(WD, "scripts/new_temp_by_date_function.R"))
new_temp <- new_temp_by_date()


# DMI rain 
DMI_rain_aalborg <- list.files(path = paste0(DataPath, "/DMI/Aalborg/"), 
                               pattern = ".csv", full.names = T
) %>% as.data.frame() %>% filter(str_detect(., "aalborg")) %>% unlist()
DMI_rain_hjorring <- list.files(path = paste0(DataPath, "/DMI/Hjorring/"), 
                                pattern = ".csv", full.names = T
)
DMI_rain_aalborg_daily <- list.files(path = paste0(DataPath, "/DMI/Aalborg_daily/"), 
                                     pattern = ".csv", full.names = T
)

DMI_rain_df <- vroom::vroom(DMI_rain_aalborg, delim = ";", col_names = F, skip = 1, col_types = "cc") %>% 
  rename(DateTime = X1, Rainfall = X2) %>%  
  mutate(Rainfall = as.numeric(sub(",", ".", Rainfall, fixed = TRUE)), 
         SampleDate = as.Date(DateTime, format ="%Y-%m-%d"),
         Municipalty = "Aalborg"
  ) %>% 
  full_join(., 
            vroom::vroom(DMI_rain_hjorring, delim = ";", col_names = F, skip = 1, col_types = "cc") %>%
              rename(DateTime = X1, Rainfall = X2) %>%  
              mutate(Rainfall = as.numeric(sub(",", ".", Rainfall, fixed = TRUE)), 
                     SampleDate = as.Date(DateTime, format ="%Y-%m-%d"),
                     Municipalty = "Hjorring"), by = join_by(DateTime, Rainfall, SampleDate, Municipalty)) %>%
  group_by(Municipalty) %>% 
  arrange(Municipalty, DateTime) %>% 
  mutate(index = 1:n(), 
         event = if_else(Rainfall >= rain_event_mm, true = 1, false = 0),
         tmpG = cumsum(c(FALSE, as.logical(diff(event)))), 
         tmp_a = c(0, diff(index)) * !event,
         tmp_b = c(diff(index), 0) * !event) %>%
  ungroup() %>% 
  group_by(Municipalty, tmpG) %>%
  mutate(tae = cumsum(tmp_a),
         tbe = rev(cumsum(rev(tmp_b)))
  ) %>%
  ungroup() %>%
  select(-c(tmp_a, tmp_b, tmpG, DateTime)) %>% 
  rename("days_since_rain_event" = "tae", 
         "time_before_rain_event" = "tbe") %>% 
  select(-index, -event) #

daily_rain_df <- vroom::vroom(DMI_rain_aalborg_daily, delim = ";", col_names = F, skip = 1, col_types = "cc") %>% 
  rename(DateTime = X1, Rainfall = X2) %>%  
  mutate(Rainfall_hourly = as.numeric(sub(",", ".", Rainfall, fixed = TRUE)), 
         DateTime = as.POSIXct(DateTime),
         SampleDate = format(DateTime, "%Y-%m-%d"),
         SampleDate = as.Date(DateTime, format ="%Y-%m-%d"),
         Time_rain = format(DateTime, "%H:%M:%S"),
         Municipalty = "Aalborg") %>% 
  select(-DateTime, -Rainfall)


metadata_l <- readxl::read_xlsx(paste0(DataPath, "GiuliaMeta.xlsx")) %>% 
  select(-MiSeqID, -ProjectName, -Investigator, -Supervisor, -LibraryID, -LibraryConc, -LibraryMethod, 
         -ExtractionMethod, -Primer, -Barcode, -TubeName, -SubmittedBy, -Comments, -Season) 

rain_fall_hourly <- metadata_l %>% 
  filter(!SampleContent %in% c("AS", "IWW")) %>% 
  filter(!SampleSite %in% c("Skyttegarden", "Vennebjerg", "SkagenLandevej", "AAU")) %>% 
  mutate(Municipalty = "Aalborg", 
         #SamplingTime = format(SamplingTime, "%H:%M:%S")
  ) %>% 
  distinct(SampleDate, SampleSite, SamplingTime, Municipalty) %>% 
  full_join(daily_rain_df, join_by(SampleDate, Municipalty), relationship = "many-to-many") %>% 
  mutate(Time_rain = as.POSIXct(paste("1899-12-31", Time_rain), format = "%Y-%m-%d %H:%M:%S")) %>% 
  filter(
    hours(Time_rain) >= hours(SamplingTime) - hours(n_hours_before) & 
      hours(Time_rain) < hours(SamplingTime)) %>% 
  group_by(SampleSite, SampleDate, SamplingTime, Municipalty) %>% 
  summarise(Rainfall_hourly_sum = sum(Rainfall_hourly), .groups = "drop")  

### Astronomical season funtion
getSeason <- function(input.date){
  numeric.date <- 100*month(input.date)+day(input.date)
  cuts <- base::cut(numeric.date, breaks = c(0,320,0621,0923,1222,1231))
  levels(cuts) <- c("Winter","Spring","Summer","Autumn","Winter")
  return(cuts)
}


print("Make metadata")
### Make metadata

metadata <- metadata_l %>% 
  filter(Old != "old" | is.na(Old)) %>%    # Remove samples from old homogenisation 
  mutate(
    Season = getSeason(SampleDate),
    Municipalty = case_when(SampleSite %in% 
                              c("Skyttegarden", "Vennebjerg", "SkagenLandevej", "Hjørring") ~
                              "Hjorring", .default = "Aalborg"), 
    Month = month(SampleDate), 
    Date_month = format(SampleDate, "%d-%m"), 
    #SamplingTime = format(SamplingTime, "%H:%M:%S"), #times(SamplingTime),
    SampleNote = case_when(SampleSite == "Ferslev" & SampleContent == "biofilm" ~ "suspended sand in the SWW", 
                           .default = PrimarySettler),
    PrimarySettler = ifelse(PrimarySettler == "first_flush", "after", PrimarySettler),
    SampleSite = str_replace(SampleSite, "AAW", "AalborgWest"),
    SampleContent = if_else(SampleName == "PCRpos", "control_pos", SampleContent),
    SampleContent = if_else(SampleName == "PCRneg", "control_neg", SampleContent),
    SampleContent = if_else(SampleName == "ExtractionControl", "control_neg", SampleContent),
    SampleDate_numeric = as.numeric(SampleDate), 
    SampleEnvironment = if_else(str_detect(SampleContent, "control"), "control", "mistake"), 
    SampleEnvironment = if_else(str_detect(SampleName, "biofilm") | 
                                  str_detect(SampleName, "sediment") | 
                                  str_detect(SampleName, "bio") |
                                  str_detect(SampleName, "cocio") | 
                                  str_detect(SampleName, "SWW_RG")
                                , "sewer_wet_solids", SampleEnvironment), 
    SampleName = case_when(str_detect(SampleName, "cocio") ~ "biofilm_RG",
                           str_detect(SampleName, "SWW_RG") ~ "biofilm_RG", 
                           str_detect(SampleName, "bio") & 
                             (str_detect(SampleSite, "Sejlflod") | 
                                str_detect(SampleSite, "Doctorvej")) ~ "biofilm_end_of_pressure",    
                           str_detect(SampleName, "bio") & str_detect(SampleSite, "Vennebjerg") ~ "biofilm_RG",
                           str_detect(SampleName, "bio") & str_detect(SampleSite, "Skyttegarden") ~ "biofilm_RG",
                           str_detect(SampleName, "bio") & str_detect(SampleNote, "pump") ~ "biofilm_pre_RG",
                           .default = SampleName),
    SampleEnvironment = case_when(str_detect(SampleName, "biofilm_dry") ~ "biofilm_dry",
                                  str_detect(SampleName, "AS") ~ "AS", 
                                  SampleName =="SWW" | str_detect(SampleName, "IWW") ~ "wastewater",
                                  .default = SampleEnvironment), 
    SampleContent = if_else(str_detect(SampleName, "_RG"), "biofilm_RG", SampleContent),
    SampleContent2 = case_when(
      SampleSite %in% c("Doctorvej", "Sejlflod") & str_detect(SampleContent, "biofilm") ~ "biofilm_end_pressure",
      .default = SampleContent),
    Configuration = case_when(
      SampleSite %in% c("Skyttegarden") ~ "wide_pressure_pipes", 
      SampleSite %in% c("Vennebjerg", "SkagenLandevej") ~ "small_pressure_pipes", 
      SampleSite %in% c("AAU") ~ "control", 
      SampleSite %in% c("AalborgØst", "AalborgWest", "Hjørring", "AAW") ~ "WWTP", 
      SampleSite %in% c("Doctorvej", "Sejlflod") ~ "end_of_pressure",
      .default = "gravity_pipes"),
    SewerType = case_when(SampleSite %in% c("Frejlev", "Tvaesgade", "Visse", "Roden") ~ "combined", 
                          SampleSite %in% c("Doctorvej","Stadionvej","Sejlflod") ~ "seperated",
                          SampleSite %in% c("Skolesti") ~ "combined", 
                          SampleSite %in% c("Skyttegarden") ~ "mixed", 
                          SampleSite %in% c("Vennebjerg", "SkagenLandevej") ~ "seperated", 
                          SampleSite %in% c("AalborgWest", "AalborgØst", "Hjørring") ~ "WWTP"
    ),
    SampleSite = ifelse(Configuration == "WWTP" & SampleContent == "IWW", 
                        paste0(SampleSite, "_", PrimarySettler), SampleSite),
    SampleName = ifelse(Configuration == "WWTP"& SampleContent == "IWW", 
                        paste0(SampleName, "_", PrimarySettler), SampleName),
    SampleName_Season = paste0(SampleName, "_", Season), 
    SampleSite_Season = paste0(SampleSite, "_", Season),
    SampleContent_SampleDate = paste0(SampleContent, "_", as.character(SampleDate)),
    SampleContent_SampleSite_SampleDate = paste0(SampleContent, "_", SampleSite,"_", as.character(SampleDate)),
    SampleContent_SampleSite = paste0(SampleContent, "_", SampleSite),
    SampleContent2_SampleSite = paste0(SampleContent2, "_", SampleSite),
    SampleContent2_Season = paste0(SampleContent2, "_", Season),
    Catchment_WWTP = case_when(SampleSite %in% c("Roden","Doctorvej","Stadionvej",
                                                 "Sejlflod", "Tvaesgade", "Visse", "Skolesti") ~ "AalborgEast",
                               SampleSite %in% c("Frejlev") ~ "AalborgWest"),
    SampleName = factor(SampleName, levels = c("AS", "IWW_before", "IWW_after", "SWW", 
                                               "biofilm_wet", "biofilm_end_of_pressure",
                                               "biofilm_RG", "biofilm_pre_RG", 
                                               "sediment", "biofilm_dry" ,"PCRpos", "PCRneg", "ExtractionControl")),
    SampleContent = factor(SampleContent, 
                           levels = c("AS", "IWW", "sewage", "biofilm", "biofilm_RG","sediment","biod","control_pos", "control_neg"), 
                           labels = c("AS", "IWW", "SWW", "biofilm", "biofilm_RG","sediment","biofilm_dry","control_pos", "control_neg")), 
    Season = factor(Season, levels = c("Winter", "Spring", "Summer", "Autumn")), 
    Configuration = factor(Configuration, 
                           levels = c("WWTP", "gravity_pipes", "end_of_pressure", 
                                      "wide_pressure_pipes", "small_pressure_pipes","control")),
    System = ifelse(SampleContent2 %in% c("AS", "IWW"), "WWTP", "Sewer system") 
    ) %>%
  left_join(., DMI_rain_df, 
            by = join_by(Municipalty, SampleDate), relationship = "many-to-many") %>% 
  left_join(et, by = c("Month")) %>% 
  filter(SampleSite != "Ferslev") %>%     ### Remove Ferslev 
  filter(SampleContent != "biofilm_dry") %>%        ### Remove biofilm dry
  filter(SampleID != "MQ220419-206") %>%    ### Remove biofilm dry from Roden which were wrongly labled
  filter(!SampleID %in% c("MQ220419-208", "MQ220419-210")) %>% # Remove biofilms taken other places than the pipe
  filter(!SampleID %in% c("MQ230313-117", "MQ230313-124")) %>% # Remove outliers    
  #filter(SampleContent != "biofilm_RG") %>%                    # Remove rensegris samples
  filter(SampleSite != "Hjørring") %>%                    # Remove Hjørring AS samples
  left_join(rain_fall_hourly, by = join_by(SampleSite, SampleDate, SamplingTime, Municipalty)) %>% 
  left_join(new_temp, by = c("Date_month"))
  
 
# Merge dublicates 
dublicates <- metadata %>% 
  #filter(Old == "rep") %>%  ### tabbed out to merged all replicates
  group_by(SampleSite, SampleContent, SampleDate, SampleNote) %>% 
  mutate(SampleID2 = paste0(unique(str_sub(SampleID, 1, 9)), paste0(str_sub(SampleID, 10, 12), collapse = ";"))) %>% 
  relocate(SampleID2, .after = "SampleID")

metadata <- metadata %>% 
  left_join(., dublicates, by = join_by(SampleID, SampleName, SampleContent, SampleContent2, SampleSite, SampleDate, PrimarySettler, SamplingTime, Old, Season, Municipalty, Month,
                                        SampleNote, SampleDate_numeric, SampleEnvironment, Configuration, System,SewerType, SampleName_Season, SampleSite_Season, SampleContent_SampleDate,
                                        SampleContent_SampleSite_SampleDate, SampleContent_SampleSite, SampleContent2_Season, SampleContent2_SampleSite, Catchment_WWTP, Rainfall, days_since_rain_event, time_before_rain_event, Earth_temp_1m,
                                        Rainfall_hourly_sum, temp, Date_month, day_number, cold_warm)) %>% 
  mutate(SampleID2 = if_else(is.na(SampleID2), SampleID, SampleID2))

rm(DMI_rain_aalborg, DMI_rain_hjorring, et, 
   DMI_rain_aalborg_daily, n_hours_before, rain_event_mm, 
   dublicates, DMI_rain_df, getSeason)


print("Load ampvis object and merge replicate samples") 
data <- amp_load(otutable = paste0(DataPath, "ASVtable_notax.tsv"),
                 taxonomy = paste0(DataPath, sintax_file),
                 metadata = metadata)


data <- amp_merge_replicates(data, merge_var = "SampleID2", round = "up")
data$tax <- data$tax %>% change_tax_names()

#


##################################################################################
##      REMOVE  CONTROLS                                                        ##
##################################################################################

print("Remove controls and make a seperate dataframe including the control samples")
# Make control dataframe
data_control <- data
amp_diversity_control_df <- 
  data_control %>% amp_alpha_diversity()

## REMOVE  CONTROLS                                                        ##

data <- data %>% 
  amp_subset_samples(str_detect(SampleContent, "control_",negate = T), 
                     minreads = 5000) 

amp_diversity_df <-
  data %>% amp_alpha_diversity()

##################################################################################

# amp_diversity_df %>%
#   group_by(System) %>%
#   summarise(median(Reads),
#             max(Reads),
#             quantile(Reads, 0.75))
#   
# print(paste0("AS and and IWW is rarefied to the median read counts of the sewer samples: median = 11214 reads (biofilm (e.p.))"))
# data$metadata %>% distinct(System)
# 
# data_sewer <- data %>% 
#   amp_subset_samples(System == "Sewer system")
# 
# data_wwtp <- data %>% 
#   amp_subset_samples(System == "WWTP", rarefy = 14728, removeAbsentOTUs = F)
# 
# 
# data <- amp_merge_ampvis2(data_sewer, data_wwtp, by_refseq = F)


#################################################################################

print("Make long DF")

if (new_long_df) {
  # Create new tidy_long object
  print("Make new long DF")
  source(paste0({{WD}}, "../2022_initial_microbial_analysis/general_scripts/tidy_long_function.R"))
  data_long <-
    tidy_midas_global(data, SampleID_column_name = "SampleID2")
  #saveRDS(data_long, file=paste0(OutputPath, "R_environment/", "20230523_data_long.Rda"))
  #rm(tidy_midas_global)
  
  tot_count <- data_long %>%  
    mutate(samples = 
             map(.x = samples, ~
                   distinct(., tot_count))
    ) %>% 
    unnest(samples) %>% select(SampleID, tot_count)
  
  data_long <- data_long %>% 
    mutate(samples = 
             map(.x = samples, ~
                   select(., -tot_count))
    ) %>% 
    left_join(., tot_count, by = join_by(SampleID)) %>% relocate(tot_count, .after = "SampleDate")
  
} else {
  print("Import old long DF")
  # Load old tidylong
  data_long <- readRDS(file=paste0(OutputPath, "R_environment/", "20230523_data_long.Rda"))
}


# Save unclassified ASV
if (save_unclassified) {
  print("Saving a new unclassified species txt file (for usearch)")
  unclassified_specieds <- data_long %>% 
    sample_n(1) %>% 
    unnest(samples) %>% 
    distinct(Species, OTU) %>% 
    filter(str_detect(Species, "unclassified")) %>% 
    select(OTU) %>% 
    unlist()
  write(unclassified_specieds, file=paste0(OutputPath, "tables/", "20230623_unclassified_species.txt"))
}

## GUT
  ### First I import the Species names of gut bacteria
  ### Limit of >=0.6 used for finding the gut species
  source(paste0({{WD}}, "../2023_sewer_microbial_communities/scripts/functions/species_level_classification_gut_microbes.R"))
  gut_species <- species_level_classification_gut_microbes(gut_database = "16S_no_requirement_MIDAS_5.2", ampvis_data = data)
  print(paste0(unlist(length((gut_species$Species))), " species are identified as gut species using a 16S database (latest update 2021)"))

## Classification of unclassified ASVs as gut species
  # If species within a query (ASV) are both gut and not gut, also the number of hits the different species are observed included 
  # A limit of >= 0.6 of the total hits
  source(file = paste0({{WD}}, "../2023_sewer_microbial_communities/scripts/functions/usearch_results_of_unclassified_gut_species.R"))
  gut_OTU <- including_unclassified_ASV_in_gut(gut_species_df = gut_species, 
                                               data_long_format = data_long)
##

## Growth 
 ### Raw growth groups from Miriam 
  source(file = paste0({{WD}}, "scripts/functions/growth_group_and_guild_classification_new_version_20240215.R"))
  growth_groups_list <- growth_groups_function_miriams_new_groups(
    amp_data = data)
  growth_df_raw <- growth_groups_list[[1]]
  growth_not_in_data_df <- growth_groups_list[[2]]
  #growth_all <- growth_groups_list[[3]]
  
 ### Using useach to assign more growth groups based on species classification
  source(paste0(WD, "scripts/functions/usearch_results_of_unclassified_new_groups_Miriam.R"))
  usearch_script_list <- including_unclassified_ASV_in_growth(
    blast_file = "data/usearch_results/20240215_usearch_ASVs_unclassified_0.6_MIDAS_5.2",
    WD = "C:/Users/HD95LP/OneDrive - Aalborg Universitet/visual_studio_code/2023_sewer_microbial_communities/",
    ASVs_and_species_not_in_data = growth_not_in_data_df,
    ASVs_in_data = growth_df_raw,
    make_plot = make_usearch_plot, 
    amp_data = data, 
    data_long_format = data_long
    )
  usearch_effect_plot <- usearch_script_list[[1]]
  usearch_result <- usearch_script_list[[2]]
  growth_df_usearch <-  usearch_script_list[[3]]
##


## Metabolic 
source(file = paste0(WD, "../2022_initial_microbial_analysis/general_scripts/import_MiDAS_metabolic_metadata.R"))
MIDAS_metabolic_df <- import_MiDAS_metabolic_metadata(file_name = "20240618_MiDAS_Metadata")
print(paste0("Midas metadata file: ", "20240618_MiDAS_Metadata"))



## DATAFRAMES TO LOAD 

writeLines("\n
1 -> metadata: Metadata with controls (dublicates are merged)
2 -> data: Ampvis object witout controls
3 -> data_long: Data in long and nested format without control
4 -> growth_df_usearch: DF with all OTU and growth group 
5 -> gut_OTU: DF with gut OTUs and the corresponding species (including unclassified gut species from usearch) 
6 -> empty ##sewer_core_species: DF with core species
7 -> MIDAS_metabolic_df: Metabolic classifications from MIDAS
8 -> amp_diversity_df: Alpha diversity
9 -> amp_diversity_control_df: Alpha diversity with controls
10-> data_control: Ampvis with control
11-> usearch_result: See all ASV that was classified with usearch at the species level
12-> usearch_effect_plot: Plot to show the effect of usearch 
13-> growth_df_raw: DF with all OTUs and growth group before during the usearch 
14 <- empty ##adonis    
15 <- empty ##sewer_core_genus:sewer core genera
           ")

empty <- "empty"

list(
  metadata, # Metadata with controls (dublicates are merged)
  data, # ampvis no control
  data_long, # data long format without control
  growth_df_usearch, # DF with all OTU and growth group 
  gut_OTU, # DF with gut species 
  empty, # DF with core taxa (loaded from object)
  MIDAS_metabolic_df, # Metabolic classifications from MIDAS
  amp_diversity_df, # Alpha diversity 
  amp_diversity_control_df, # Alpha diversity with controls
  data_control,  # Ampvis with control
  usearch_result, # Usearch results
  usearch_effect_plot, # Usearch effect plot
  growth_df_raw # DF with all OTUs and growth group before during the usearch  
    # Adonis function
  # sewer core genera
  ) 


}

