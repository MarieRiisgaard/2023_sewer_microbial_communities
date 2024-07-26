#############################################
## Convert ampvis object to nested format  ##
#############################################


# This script inputs a ampvis object and outputs a tibble where the OTU table is nested for each sample
# Makes it easier to work in a the tidyverse long format 

ampvis_to_long_data <- function(data, SampleID_column_name = "SampleID"){
  
  data$tax <- data$tax %>% 
    mutate(Phylum = if_else(condition = str_detect(Phylum, 'p__'), true = Phylum, false = paste0("unclassified_", Kingdom))) %>% 
    mutate(Class = if_else(condition = str_detect(Class, 'c__'), true = Class, false = paste0("unclassified_", Phylum))) %>% 
    mutate(Order = if_else(condition = str_detect(Order, 'o__'), true = Order, false = paste0("unclassified_", Class))) %>% 
    mutate(Family = if_else(condition = str_detect(Family, 'f__'), true = Family, false = paste0("unclassified_", Order))) %>% 
    mutate(Genus = if_else(condition = str_detect(Genus, 'g__'), true = Genus, false = paste0("unclassified_", Family))) %>% 
    mutate(Species = if_else(condition = str_detect(Species, 's__'), true = Species, false = paste0("unclassified_", Genus))) %>% 
    mutate(Class = str_replace(string = Class, pattern = ".*(_[a-z]__.*)", replacement = "unclassified\\1")) %>% 
    mutate(Order = str_replace(string = Order, pattern = ".*(_[a-z]__.*)", replacement = "unclassified\\1")) %>% 
    mutate(Family = str_replace(string = Family, pattern = ".*(_[a-z]__.*)", replacement = "unclassified\\1")) %>% 
    mutate(Genus= str_replace(string = Genus, pattern = ".*(_[a-z]__.*)", replacement = "unclassified\\1")) %>% 
    mutate(Species = str_replace(string = Species, pattern = ".*(_[a-z]__.*)", replacement = "unclassified\\1")) %>% 
    mutate(Phylum = if_else(str_detect(Phylum, "unclassified"), paste0(Phylum, "_", OTU), Phylum),
           Class = if_else(str_detect(Class, "unclassified"), paste0(Class, "_", OTU), Class),
           Order = if_else(str_detect(Order, "unclassified"), paste0(Order, "_", OTU), Order),
           Family = if_else(str_detect(Family, "unclassified"), paste0(Family, "_", OTU), Family),
           Genus = if_else(str_detect(Genus, "unclassified"), paste0(Genus, "_", OTU), Genus),
           Species = if_else(str_detect(Species, "unclassified"), paste0(Species, "_", OTU), Species)
    )
  
  
  tidy_test <- data %>% 
    amp_export_long() %>%
    as_tibble() %>% 
    rename("SampleID" = SampleID_column_name) %>% 
    mutate(SampleID_n = SampleID) %>% 
    #filter(count != 0) %>%   #  <--- filtering step to make datahandling faster (8/8-22) 
    nest(samples = c(SampleID_n, Kingdom, Phylum, Class, Order, Family, Genus, Species, OTU, count
                     )) %>%
    mutate(samples = map(.x = samples, 
                         ~mutate(.x, 
                                 tot_count = sum(count)))) %>% 
    mutate(samples = map(.x = samples, 
                         ~mutate(.x, 
                                 rel_abun_ASV = count/tot_count*100))) %>% 
    mutate(samples = map(.x = samples, ~group_by(.x, Species) %>% mutate(rel_abun_species = sum(rel_abun_ASV)) %>% ungroup())) %>% 
    mutate(samples = map(.x = samples, ~group_by(.x, Genus) %>% mutate(rel_abun_genus = sum(rel_abun_ASV))%>% ungroup())) %>%
    ungroup()
    }




