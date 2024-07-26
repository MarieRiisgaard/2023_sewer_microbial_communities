############################################################
### Species level classification of gut bacteria          ##
############################################################


# Step 1: By matching the gut databases to the Midas database 

species_level_classification_gut_microbes <- 
  function(ampvis_data, 
           WD = "C:/Users/HD95LP/OneDrive - Aalborg Universitet/visual_studio_code/2023_sewer_microbial_communities/"
           ){
    


# To classify gut bacteria, we used the 16S downloaded from these 2 studies: 
  ## Ritari, J., Salojärvi, J., Lahti, L. et al. 
  ## Improved taxonomic assignment of human intestinal 16S rRNA sequences by a dedicated reference database. 
  ## BMC Genomics 16, 1056 (2015). https://doi.org/10.1186/s12864-015-2265-y
  # AND 
  ## Kim, C.Y., Lee, M., Yang, S. et al. 
  ## Human reference gut microbiome catalog including newly assembled genomes from under-represented Asian metagenomes. 
  ## Genome Med 13, 134 (2021). https://doi.org/10.1186/s13073-021-00950-7

# Details: 

    # HITdb - Human Intestinal 16S rRNA gene reference taxonomy (2015)
      ## https://github.com/openresearchlabs/HITdb/tree/master
      ##  2473 unique prokaryotic species-like group (sequences was clustered at 97% identity OTU)

    # Human Reference Gut Microbiome (HRGM) is an updated human gut microbiome catalog by including newly assembled 29,082 genomes on 845 fecal samples collected from three under-represented Asian countries—Korea, India and Japan
      ## https://www.mbiomenet.org/HRGM/
      ## I could not find the data online anymore - Previous download used

# The database was downloaded (the 16S sequences) 
    ## Uploaded in the datafolder
    ## The Midas database + sequences for the unclassified ASV were mapped to both database using usearch
    ## See script 'reclassify_with_updated_database.sh' that resulted in the two files imported here
      ##  gut_matches_HRMG.txt
      ##  gut_matches_HitDB.txt
    ## Filtering 
      ## Identity limit of > 98.7
      ## 50 species hits were allowed. Only species constituted > 60% of these hit were included as gut species (similar to sintax 0.6)  

## Importing Midas tax 
gut_midas <- read.csv(paste0({{WD}}, "data/bacterial_groups/20240215_16S_database_HRGM_HitDB_MiDAS_5/gut_matches_HRMG.txt"), 
                      sep = "\t", header = F)
names(gut_midas) <- c("query", "DB", "percent_identity", "aligment_lenght", 
                      "mismatches", "gaps", "query_start_position", "query_end_position", 
                      "DB_start_position", "DDB_end_position", "e_value", "bit_score")
gut_midas <- gut_midas %>% separate(DB, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), ",") 
gut_midas_HRGM <- gut_midas %>% separate(query, c("v1", "v2", "v3", "V4", "V5", "V6", "V7", "V8", "V9"), "\\|") %>% 
  filter(percent_identity > 98.7) %>% 
  separate(Domain, c("OTU", "Domian"), ";tax=") %>% 
  mutate(
    Species = str_remove(Species, ";"), 
    Species = str_replace(Species, ":", "__")
  ) %>% 
  group_by(V5, Genus, Species) %>% 
  group_by(V5, Genus, Species) %>% 
  summarise(n_obs_per_species = n(), .groups = "drop") %>% 
  group_by(V5) %>% 
  reframe(n_species = n_distinct(Species), 
          sum_n_obs_per_species = sum(n_obs_per_species),
          fraction = n_obs_per_species/sum_n_obs_per_species,
          remove = if_else(fraction < 0.6, T, F),
          Species, n_obs_per_species) %>% 
  distinct(Species, remove)


gut_midas <- read.csv(paste0({{WD}}, "data/bacterial_groups/20240215_16S_database_HRGM_HitDB_MiDAS_5/gut_matches_HitDB.txt"), 
                      sep = "\t", header = F)
#names(gut_midas_HitDB) <- c(start_position", "DDB_end_position", "e_value", "bit_score")
gut_midas <- gut_midas %>% separate(V2, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), ",") 
gut_midas_HitDB <- gut_midas %>% #separate(query, c("v1", "v2", "v3", "V4", "V5", "V6", "V7", "V8", "V9"), "\\|") %>% 
  filter(V3 > 98.7) %>% 
  separate(Domain, c("OTU", "Domian"), ";tax=") %>% 
  mutate(
    Species = str_remove(Species, ";"), 
    Species = str_replace(Species, ":", "__")
    #OTU = gsub("\\..*","", OTU), 
    #OTU = str_remove(OTU, "FL")
  ) %>% 
  group_by(V1, Genus, Species) %>% 
  summarise(n_obs_per_species = n(), .groups = "drop") %>% 
  group_by(V1) %>% 
  reframe(n_species = n_distinct(Species), 
          sum_n_obs_per_species = sum(n_obs_per_species),
          fraction = n_obs_per_species/sum_n_obs_per_species,
          remove = if_else(fraction < 0.6, T, F),
          Species, n_obs_per_species) %>% 
  distinct(Species, remove)

# anti_join(gut_midas_HitDB %>% filter(remove), gut_midas_HRGM %>% filter(remove))
# anti_join(gut_midas_HRGM %>% filter(remove), gut_midas_HitDB %>% filter(remove))
# anti_join(gut_midas_HitDB %>% filter(!remove), gut_midas_HRGM %>% filter(!remove))
# anti_join(gut_midas_HRGM %>% filter(!remove), gut_midas_HitDB %>% filter(!remove))
# anti_join(gut_midas_HRGM, gut_midas_HitDB)

gut_tax <- rbind(gut_midas_HitDB, gut_midas_HitDB) %>% filter(!remove) %>%
  distinct(Species)

  }





# Step 2: By checking the unclassified ASV that are due to multiple species matches 

# Define function to read BLAST output file
read_blast_output <- function(file_path) {
  # Read the file
  blast_data <- read_delim(file_path, delim = "\t", col_names = FALSE, show_col_types = FALSE)
  
  # Rename columns
  colnames(blast_data) <- c("query", "subject", "percent_identity", "alignment_length",
                            "mismatches", "gap_opens", "query_start", "query_end",
                            "subject_start", "subject_end", "e_value", "bit_score")
  
  # Return the parsed data
  return(blast_data)
}


  
including_unclassified_ASV_in_gut <- function(
    blast_file = "data/usearch_results/20240215_usearch_ASVs_unclassified_0.6_MIDAS_5.2",
    WD = "C:/Users/HD95LP/OneDrive - Aalborg Universitet/visual_studio_code/2023_sewer_microbial_communities/",
    gut_species_df = gut_species,
    data_long_format = data_long){
  
  
# Approach 
  ## When we go down to species level sequences are often unclassified due to multiple species matches (>98.7% identity)
    ## - thus cannot disqinguish between 2 or more species
  ## Here I find all the ASVs where that is the case while ASV with no MiDAS database match are removed.
  ## Useach command: 
    ### usearch -usearch_global $query -db $MIDAS_DB -id 0.987 -blast6out 
    ###     $results -strand plus -query_cov 1.0 -maxaccepts 50
  ## I check whether the species matches are gut species.
  ## In case that > 60% of the assigned species are gut species I classify the ASV as gut 
  ## E.g. ASV10133 that maches to 3 disitnct species
    # ASV10133 s__midas_s_50099            12 98.8-100           0-3          0-0         TRUE  
    # ASV10133 s__midas_s_50946            10 98.8-100           0-3          0-0         TRUE  
    # ASV10133 s__midas_s_68495             2 99.2-99.2          2-2          0-0         FALSE 
  ## More than 60% of the n_obs are to a gut species (TRUE) --> assigned as gut.
  
  
  
  # Read the BLAST output file
  blast_results <- read_blast_output(paste0(WD, blast_file)) 
  length(unique(blast_results$query))  

  results <- blast_results %>% 
    mutate(subject = str_replace_all(subject,  ":", "__")) %>% 
    separate(subject, c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species"), ",") %>% 
    separate(Domain, c("OTU", "Domian"), ";tax=") %>% 
    mutate(
      Species = str_remove(Species, ";")) %>% 
    group_by(query, Species) %>% 
    summarise(n_obs = n(), 
              percent_identity_r = paste0(min(percent_identity),"-",max(percent_identity)),
              mismatches_r = paste0(min(mismatches),"-",max(mismatches)),
              gap_opens_r = paste0(min(gap_opens),"-",max(gap_opens)),
              .groups = "drop"
    ) 
  
  
  ## Check which of the assigned species that are also gut species 
  is_gut_df <- results %>% 
    left_join(gut_species_df %>% mutate(is_gut = T), by = join_by(Species)) %>% 
    mutate(is_gut = if_else(is.na(is_gut), F, is_gut)) 
  
  is_gut_df_ASVs <- is_gut_df %>% 
    group_by(query) %>% 
    reframe(
      Species,
      sum_observed_species = sum(n_obs), 
      n_obs, is_gut) %>% 
    ungroup %>% 
    group_by(query, is_gut) %>% 
    reframe(
      n_obs_by_gut = sum(n_obs), Species,
      frac = n_obs_by_gut/sum_observed_species) %>% 
    mutate(keep = ifelse(frac >= 0.6, T, F)) %>% 
    filter(keep & is_gut)  
 
  # Merged the gut species from the classified and unclassiifed ASV and extract both ASV and Species names for all 
  gut_OTUs <- data_long_format %>% 
    sample_n(1) %>% unnest(samples) %>% 
    filter(Species %in% unlist(gut_species_df$Species) | 
             OTU %in% unlist(is_gut_df_ASVs$query)) %>% 
    select(OTU, Species)
  

}















