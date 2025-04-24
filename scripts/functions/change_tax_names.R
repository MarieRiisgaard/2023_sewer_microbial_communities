#####################
# Change tax names  # 
#####################

## Change tax names to correspond to MiDAS 5, version 3  
## For details see MiDAS 5 change logs. 

####### USAGE #######

##  AMB_OBJECT$tax  <- AMB_OBJECT$tax %>% change_tax_names()

## Input must include columns: Family, Genus, Species

#####################

change_tax_names <- function(
  tax_data
  ){

#### Update tax names

  #  MiDAS v 5.1 
  tax_data <- 
    tax_data %>%  
  mutate(
    #g__midas_g_65
    Species = ifelse(Genus == "g__midas_g_65" & Species == "s__midas_s_177", "s__midas_s_177", Species),
    Genus = ifelse(Genus == "g__midas_g_65" & Species == "s__midas_s_177", "g__Ca_Vicinibacter", Genus),
    Species = ifelse(Genus == "g__midas_g_65" & Species == "s__midas_s_65", "s__Ca_Vicinibacter_haderslevense", Species),
    Genus = ifelse(Genus == "g__midas_g_65" & Species == "s__Ca_Vicinibacter_haderslevense", "g__Ca_Vicinibacter", Genus),
    #g__midas_g_17
    Genus = ifelse(Genus == "g__midas_g_17", "g__Ca_Opimibacter", Genus),
    #OLB
    Species = ifelse(Genus == "g__OLB8" & Species == "s__midas_s_29", "s__Ca_Brachybacter_algidus", Species),
    Genus = ifelse(Genus == "g__OLB8" & Species == "s__Ca_Brachybacter_algidus", "g__Ca_Brachybacter", Genus),
    Species = ifelse(Genus == "g__OLB8" & Species == "s__midas_s_3279", "s__Ca_Parvibacillus_calidus", Species),
    Genus = ifelse(Genus == "g__OLB8" & Species == "s__Ca_Parvibacillus_calidus", "g__Ca_Parvibacillus", Genus),
    # Accumulibactor
    Species = ifelse(Genus == "g__Ca_Accumulibacter" & Species == "s__midas_s_168", "s__Ca_Proximibacter_danicus", Species),
    Genus = ifelse(Genus == "g__Ca_Accumulibacter" & Species == "s__Ca_Proximibacter_danicus", "g__Ca_Proximibacter", Genus),
    Species = ifelse(Genus == "g__Ca_Accumulibacter" & Species == "s__midas_s_315", "s__Ca_Propionivibrio_dominans", Species),
    Genus = ifelse(Genus == "g__Ca_Accumulibacter" & Species == "s__Ca_Propionivibrio_dominans", "g__Ca_Propionivibrio", Genus),
    # Tetrespeara
    Species = ifelse(Genus == "g__Tetrasphaera" & Species == "s__midas_s_5", "s__midas_s_5", Species),
    Genus = ifelse(Genus == "g__Tetrasphaera" & Species == "s__midas_s_5", "g__Ca_Phosphoribacter", Genus),
    Species = ifelse(Genus == "g__Tetrasphaera" & Species == "s__midas_s_45", "s__Ca_Lutibacillus_vidarii", Species),
    Genus = ifelse(Genus == "g__Tetrasphaera" & Species == "s__Ca_Lutibacillus_vidarii", "g__Ca_Lutibacillus", Genus)
  )

  # Additional changes midas 5.2
  tax_data <- 
    tax_data %>%  
    mutate(
      Genus = case_when(
        Genus == "g__Ca_Vicinibacter" & Species == "s__Ca_Vicinibacter_haderslevense" ~ "g__Defluviibacterium",
        Genus == "g__Ca_Lutibacillus" & Species == "s__Ca_Lutibacillus_vidarii" ~ "g__Ca_Lutibacillus",
        Genus == "g__Ca_Propionivibrio" & Species == "s__Ca_Propionivibrio_dominans" ~ "g__Propionivibrio",
        TRUE ~ Genus
        ),
      Species = case_when(
        Species == "s__Ca_Vicinibacter_haderslevense" ~ "s__Defluviibacterium_haderslevense",
        Species == "s__Ca_Lutibacillus_vidarii" ~ "s__Ca_Lutibacillus_vidarii",
        Species == "s__Ca_Propionivibrio_dominans" ~ "s__Propionivibrio_dominans",
        TRUE ~ Species)
      )
  
  tax_data <- 
    tax_data %>%  
    mutate(
      Genus = case_when(
        Genus == "g__Dechloromonas" ~ "g__Azonexus",
        TRUE ~ Genus  # For other cases, keep the original value
      ),
      Species = case_when(
        Species == "s__Ca_Dechloromonas_phosphoritropha" ~ "s__Azonexus_phosphoritropha",
        Species == "s__Ca_Dechloromonas_phosphorivorans" ~ "s__Azonexus_phosphorivorans",
        Species == "s__Dechloromonas_hortensis" ~ "s__Azonexus_hortensis",
        Species == "s__Dechloromonas_agitata" ~ "s__Azonexus_agitata",
        TRUE ~ Species  # For other cases, keep the original value
      )
    )

tax_data <- 
    tax_data %>%  
    mutate(
      Family = case_when(
        Family == "f__midas_f_119" ~ "f__Ca_Epilineaceae",
        Family == "f__midas_f_72" ~ "f__Ca_Brachytrichaceae",
        Family == "f__A4b" ~ "f__Flexifilaceae",
        TRUE ~ Family  # For other cases, keep the original value
      ),
      Genus = case_when(
        Genus == "g__midas_g_119" ~ "g__Ca_Epilinea",
        Genus == "g__midas_g_1676" ~ "g__Ca_Avedoeria",
        Genus == "g__midas_g_72" ~ "g__Ca_Brachythrix",
        Genus == "g__UTCFX1" & Species == "s__midas_s_9708" ~ "g__Ca_Defluviilinea",
        Genus == "g__UTCFX1" & Species == "s__midas_s_12690" ~ "g__Ca_Defluviilinea",
        Genus == "g__midas_g_9708" ~ "g__Ca_Defluviilinea",
        Genus == "g__midas_g_2111" ~ "g__Ca_Hadersleviella",
        Genus == "g__midas_g_1951" ~ "g__Ca_Trichofilum",
        Genus == "g__midas_g_461" ~ "g__Ca_Leptovillus",
        Genus == "g__OLB13" & Species == "s__midas_s_19645" ~ "g__Ca_Flexicrinis",
        Genus == "g__OLB13" & Species == "s__midas_s_1091" ~ "g__Ca_Flexicrinis",
        Genus == "g__OLB13" & Species == "s__midas_s_81436" ~ "g__Ca_Flexicrinis",
        Genus == "g__OLB13" & Species == "s__midas_s_21047" ~ "g__Ca_Flexicrinis",
        Genus == "g__OLB13" & Species == "s__midas_s_17901" ~ "g__Ca_Flexicrinis",
        Genus == "g__midas_g_4871" ~ "g__Ca_Flexifilum",
        Genus == "g__midas_g_5525" ~ "g__Ca_Flexifilum",
        Genus == "g__midas_g_2265" ~ "g__Ca_Fredericiella",
        Genus == "g__midas_g_105" ~ "g__Ca_Caldilinea",
        Genus == "g__midas_g_2775" ~ "g__Ca_Ribeiella",
        Genus == "g__midas_g_731" ~ "g__Ca_Amarobacter",
        Genus == "g__midas_g_1412" ~ "g__Ca_Amarobacillus",
        Genus == "g__midas_g_391" ~ "g__Ca_Amarofilum",
        Genus == "g__midas_g_550" ~ "g__Ca_Pachofilum",
        Genus == "g__midas_g_9648" ~ "g__Ca_Tricholinea",
        Genus == "g__midas_g_169" ~ "g__Ca_Defluviifilum",
        TRUE ~ Genus
      ),
      Species = case_when(
        Species == "s__midas_s_119" ~ "s__Ca_Epilinea_brevis",
        Species == "s__midas_s_57933" ~ "s__Ca_Epilinea_brevis",
        Species == "s__midas_s_5085" ~ "s__Ca_Epilinea_brevis",
        Species == "s__midas_s_1676" ~ "s__Ca_Avedoeria_danica",
        Species == "s__midas_s_72" ~ "s__Ca_Brachythrix_odensensis",
        Species == "s__midas_s_9708" ~ "s__Ca_Defluviilinea_gracilis",
        Species == "s__midas_s_12690" ~ "s__Ca_Defluviilinea_gracilis",
        Species == "s__midas_s_16519" ~ "s__Ca_Defluviilinea_proxima",
        Species == "s__midas_s_6664" ~ "s__Ca_Villigracilis_vicinus",
        Species == "s__midas_s_15163" ~ "s__Ca_Villigracilis_adiacens",
        Species == "s__midas_s_9223" ~ "s__Ca_Villigracilis_propinquus",
        Species == "s__midas_s_2642" ~ "s__Ca_Hadersleviella_danica",
        Species == "s__midas_s_1951" ~ "s__Ca_Trichofilum_aggregatum",
        Species == "s__midas_s_176" ~ "s__Ca_Promineofilum_glycogenico",
        Species == "s__midas_s_40638" ~ "s__Ca_Promineofilum_glycogenico",
        Species == "s__midas_s_441" ~ "s__Ca_Promineofilum_glycogenico",
        Species == "s__midas_s_19645" ~ "s__Ca_Flexicrinis_affinis", 
        Species == "s__midas_s_1091" ~ "s__Ca_Flexicrinis_proximus", 
        Species == "s__midas_s_4871" ~ "s__Ca_Flexifilum_breve",
        Species == "s__midas_s_23442" ~ "s__Ca_Flexifilum_breve",
        Species == "s__midas_s_5525" ~ "s__Ca_Flexifilum_affine",
        Species == "s__midas_s_1" ~ "s__Ca_Amarolinea_dominans",
        Species == "s__midas_s_4097" ~ "s__Ca_Fredericiella_danica",
        Species == "s__midas_s_35902" ~ "s__Ca_Fredericiella_danica",
        Species == "s__midas_s_105" ~ "s__Ca_Caldilinea_saccharophila",
        Species == "s__midas_s_2775" ~ "s__Ca_Ribeiella_danica",
        Species == "s__midas_s_2308" ~ "s__Ca_Kouleothrix_ribensis",
        Species == "s__midas_s_849" ~ "s__Ca_Amarobacter_glycogenicus",
        Species == "s__midas_s_44698" ~ "s__Ca_Amarobacter_glycogenicus",
        Species == "s__midas_s_1880" ~ "s__Ca_Amarobacter_glycogenicus",
        Species == "s__midas_s_74685" ~ "s__Ca_Amarobacter_glycogenicus",
        Species == "s__midas_s_30179" ~ "s__Ca_Amarobacter_glycogenicus",
        Species == "s__midas_s_2745" ~ "s__Ca_Amarobacter_glycogenicus",
        Species == "s__midas_s_8871" ~ "s__Ca_Amarobacillus_elongatus",
        # From Jette: Arachnia is not a genus name in Midas (to be implemented later) 
        #Species == "s__Tessaracoccus_flavescens" ~ "s__Arachnia_flavescens",
        #Species == "s__Tessaracoccus_defluvii" ~ "s__Arachnia_defluvii",
        #Species == "s__Tessaracoccus_lapidicaptus" ~ "s__Arachnia_lapidicaptus",
        #Species == "s__Tessaracoccus_bendigoensis" ~ "s__Arachnia_bendigoensis",
        TRUE ~ Species
      )
    )
  
}
