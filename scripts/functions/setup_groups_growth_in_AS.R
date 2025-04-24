growth_groups_function <- function(
    amp_data = data, 
    WD = "C:/Users/HD95LP/OneDrive - Aalborg Universitet/visual_studio_code/2023_sewer_microbial_communities/"){

  writeLines(
  "I merge on species level, meaning if a there is an ASV in my data not present in the growth characterization study, 
  but that ASV belongs to a species that were investigated in the growth characterization study,
  this ASV is assigned to the same groups as the species.
  Species not investigated is assigned to *unknown*")
  
growth_raw <- vroom::vroom(
  paste0(WD, "data/bacterial_groups/growth_in_AS/NEW_listgrowthgroups_MiDAS48.csv"), delim = ",", show_col_types = FALSE) %>% 
  mutate(growth_fate = if_else(is.na(growth_fate), "too_low_abundance", growth_fate)) %>% 
  filter(str_detect(Taxa, "newASV", negate = T)) %>% 
  change_tax_names()  %>%  # Use the change tax names function to update to MIDAS 5.3 
  mutate(Taxa = if_else(is.na(Species), OTU, Species))    # New column filling in ASV name for unclassified ASV


## Growth by species classification
growth_by_species <- growth_raw %>% 
  filter(!is.na(Species)) %>% 
  group_by(Species, growth_fate) %>% 
  summarise(n = n(), .groups = "drop") %>% 
  group_by(Species) %>% 
  reframe(
    growth_fate = paste0(unique(growth_fate), collapse = ";"),
    same_group = str_detect(growth_fate, ";", negate = T),
    n_list = paste0(n, collapse = ";")) %>%
  mutate(
    growth_fate1 = growth_fate,
    growth_fate = if_else(!same_group & str_detect(growth_fate, "too_low_abundance;"), 
                              str_remove(growth_fate, "too_low_abundance;"),
                              growth_fate), 
         growth_fate = if_else(!same_group & str_detect(growth_fate, ";too_low_abundance"), 
                               str_remove(growth_fate, ";too_low_abundance"),
                               growth_fate), 
         growth_fate = if_else(!same_group & str_detect(growth_fate, ";"), 
                               "variable",
                               growth_fate)
         ) %>% 
  rename( "growth_AS" = "growth_fate")


## Check if ASV within which species are inconistent (only 26!!)
  #see <- growth_by_species %>% filter(same_group == F)

## Growth by ASV classification (only for unclassified species)
growth_by_ASV <- growth_raw %>% filter(str_detect(Taxa, "ASV")) %>% 
  rename("growth_AS" = "growth_fate") %>% 
  select(OTU, growth_AS)


# Importing the relevant dataframe (ampvis object)
tidy_all_OTU <- amp_data$tax %>% as_tibble() %>% select(Species, OTU) %>% distinct(.keep_all = T)

species <- tidy_all_OTU %>% filter(str_detect(Species, "s__"))
ASVs <- tidy_all_OTU %>% filter(str_detect(Species, "s__", negate = T)) %>% select(-Species)

species_w_growth <- full_join(species %>% mutate(in_data = T), growth_by_species, by = join_by(Species)) %>% 
  mutate(growth_AS = if_else(is.na(growth_AS), "not_investigated", growth_AS),
         assigned = "by_species") %>%
  select(-same_group, -n_list)

ASVs_w_growth <- full_join(ASVs %>% mutate(in_data = T), growth_by_ASV, by = join_by(OTU)) %>% 
  mutate(growth_AS = if_else(is.na(growth_AS), "not_investigated", growth_AS), 
         assigned = "by_ASV") 

all_ASVs_growth <- full_join(species_w_growth, ASVs_w_growth, join_by(OTU, in_data, growth_AS, assigned))
ASVs_and_species_not_in_data <- all_ASVs_growth %>% filter(is.na(in_data)) %>% 
  mutate(Taxa = case_when(str_detect(Species, "s__") ~ "Species", .default = "OTU")) %>% select(-in_data) %>% 
  mutate(growth_AS = if_else(growth_AS == "too_low_abundance" |growth_AS == "not_investigated", 
                             "unknown", growth_AS))
ASVs_in_data <- all_ASVs_growth %>% filter(in_data) %>% 
  select(OTU, growth_AS, assigned) %>% 
  mutate(growth_AS = if_else(growth_AS == "too_low_abundance" |growth_AS == "not_investigated", 
                             "unknown", growth_AS))


writeLines("\n
1 dataframe: Growth groups for ASVs in the data\n
2 dataframe: Growth groups for ASVs and/or speices NOT in the data\n
3 dataframe: Growth groups for all ASVs and/or speices")
list(ASVs_in_data, ASVs_and_species_not_in_data, all_ASVs_growth)

}


