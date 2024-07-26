import_MiDAS_metabolic_metadata <- function(
    path,
    file_name){

#Load data into ampvis format
dataPath <- path 

transpose_df <- function(df) {
  t_df <- data.table::transpose(df)
  colnames(t_df) <- rownames(df)
  rownames(t_df) <- colnames(df)
  t_df <- t_df %>%
    tibble::rownames_to_column(.data = .) %>%
    tibble::as_tibble(.)
  return(t_df)
}


#file_name = "20240618_MiDAS_Metadata"

meta_midas_raw <-  vroom(paste0(dataPath, {{file_name}}, ".csv"), show_col_types = FALSE, delim = ";"
                          ) 

meta_midas <- meta_midas_raw %>% select(-"Alternative names", -mesophilic, -thermophilic) %>% 
  rename("Genus" = 1) %>% 
  select(-Kingdom, -Phylum, -Class, -Order, -Family, - influent, -sludge) %>% 
  mutate(Genus = str_replace_all(Genus, " ", "_"), 
         Genus = str_c( "g__", Genus)) %>% 
  #column_to_rownames("Genus") %>% #na.omit() %>% 
  pivot_longer(cols = -Genus, names_to = "metabolism", values_to = "value") %>%
  mutate(value = if_else(is.na(value), "na", value)) %>% 
  pivot_wider(names_from = metabolism, values_from = value)

sulfate_red <- meta_midas %>% select(Genus, `Sulfate reduction:Other`, `Sulfate reduction:In situ`) %>% 
  mutate("Sulfate red." = if_else(`Sulfate reduction:Other` %in% c("neg") | `Sulfate reduction:In situ` %in% c("neg"), "negative", "unknown"),
         "Sulfate red." = if_else(`Sulfate reduction:Other` %in% c("var") | `Sulfate reduction:In situ` %in% c("var"), "varible", `Sulfate red.`),
         "Sulfate red." = if_else(`Sulfate reduction:Other` %in% c("pos") | `Sulfate reduction:In situ` %in% c("pos"), "positive", `Sulfate red.`),
  ) %>% select(Genus, "Sulfate red.")

fermentation <- meta_midas %>% select(Genus, `Fermentation:In situ`, `Fermentation:Other`) %>% 
  mutate("Fermentation" = if_else(`Fermentation:In situ` %in% c("neg") | `Fermentation:Other` %in% c("neg"), "negative", "Unknown"),
         "Fermentation" = if_else(`Fermentation:In situ` %in% c("var") | `Fermentation:Other` %in% c("var"), "varible", Fermentation),
         "Fermentation" = if_else(`Fermentation:In situ` %in% c("pos") | `Fermentation:Other` %in% c("pos"), "positive", Fermentation),
         ) %>% 
  select(Genus, "Fermentation")

Aero_hetero <- meta_midas %>% select(Genus, `Aerobic heterotroph:In situ`, `Aerobic heterotroph:Other`) %>% 
  mutate("Aerobic heterotroph" = if_else(`Aerobic heterotroph:In situ` %in% c("neg") | `Aerobic heterotroph:Other` %in% c("neg"), "negative", "unknown"),
         "Aerobic heterotroph" = if_else(`Aerobic heterotroph:In situ` %in% c("var") | `Aerobic heterotroph:Other` %in% c("var"), "varible", `Aerobic heterotroph`),
         "Aerobic heterotroph" = if_else(`Aerobic heterotroph:In situ` %in% c("pos") | `Aerobic heterotroph:Other` %in% c("pos"), "positive", `Aerobic heterotroph`),
  ) %>% select(Genus, "Aerobic heterotroph")


chemo_auto <- meta_midas %>% select(Genus, `Chemoautotroph/mixotroph:In situ`, `Chemoautotroph/mixotroph:Other`) %>% 
  mutate("Chemoauto-/mixotroph" = if_else(`Chemoautotroph/mixotroph:In situ` %in% c("neg") | `Chemoautotroph/mixotroph:Other` %in% c("neg"), "negative", "unknown"),
         "Chemoauto-/mixotroph" = if_else(`Chemoautotroph/mixotroph:In situ` %in% c("var") | `Chemoautotroph/mixotroph:Other` %in% c("var"), "varible", `Chemoauto-/mixotroph`),
         "Chemoauto-/mixotroph" = if_else(`Chemoautotroph/mixotroph:In situ` %in% c("pos") | `Chemoautotroph/mixotroph:Other` %in% c("pos"), "positive", `Chemoauto-/mixotroph`)) %>% 
  select(Genus, "Chemoauto-/mixotroph")


## Changed 

pao <- meta_midas %>% select(Genus, `PAO:In situ`, `PAO:Other`) %>% 
  mutate("PAO" = if_else(`PAO:In situ` %in% c("neg") | `PAO:Other` %in% c("neg"), "negative", "unknown"),
         "PAO" = if_else(`PAO:In situ` %in% c("var") | `PAO:Other` %in% c("var"), "varible", `PAO`),
         "PAO" = if_else(`PAO:In situ` %in% c("pos") | `PAO:Other` %in% c("pos"), "positive", `PAO`)) %>% 
  select(Genus, "PAO")

gao <- meta_midas %>% select(Genus, `GAO:In situ`, `GAO:Other`) %>% 
  mutate("GAO" = if_else(`GAO:In situ` %in% c("neg") | `GAO:Other` %in% c("neg"), "negative", "unknown"),
         "GAO" = if_else(`GAO:In situ` %in% c("var") | `GAO:Other` %in% c("var"), "varible", `GAO`),
         "GAO" = if_else(`GAO:In situ` %in% c("pos") | `GAO:Other` %in% c("pos"), "positive", `GAO`)) %>% 
  select(Genus, "GAO")

Nitrite <- meta_midas %>% select(Genus, `Nitrite reduction:In situ`, `Nitrite reduction:Other`) %>% 
  mutate("Nitrite" = if_else(`Nitrite reduction:In situ` %in% c("neg") | `Nitrite reduction:Other` %in% c("neg"), "negative", "unknown"),
         "Nitrite" = if_else(`Nitrite reduction:In situ` %in% c("var") | `Nitrite reduction:Other` %in% c("var"), "varible", `Nitrite`),
         "Nitrite" = if_else(`Nitrite reduction:In situ` %in% c("pos") | `Nitrite reduction:Other` %in% c("pos"), "positive", `Nitrite`)) %>% 
  select(Genus, Nitrite) %>% 
  rename("Nitrite_reduction" = "Nitrite")

Filamentous <- meta_midas %>% select(Genus, `Filamentous:In situ`, `Filamentous:Other`) %>% 
  mutate("Filamentous" = if_else(`Filamentous:In situ` %in% c("neg") | `Filamentous:Other` %in% c("neg"), "negative", "unknown"),
         "Filamentous" = if_else(`Filamentous:In situ` %in% c("var") | `Filamentous:Other` %in% c("var"), "varible", `Filamentous`),
         "Filamentous" = if_else(`Filamentous:In situ` %in% c("pos") | `Filamentous:Other` %in% c("pos"), "positive", `Filamentous`)) %>% 
  select(Genus, Filamentous) 

AOB <-meta_midas %>% select(Genus, `AOB:In situ`, `AOB:Other`) %>% 
  mutate("AOB" = if_else(`AOB:In situ` %in% c("neg") | `AOB:Other` %in% c("neg"), "negative", "unknown"),
         "AOB" = if_else(`AOB:In situ` %in% c("var") | `AOB:Other` %in% c("var"), "varible", `AOB`),
         "AOB" = if_else(`AOB:In situ` %in% c("pos") | `AOB:Other` %in% c("pos"), "positive", `AOB`)) %>% 
  select(Genus, "AOB")

NOB <-meta_midas %>% select(Genus, `NOB:In situ`, `NOB:Other`) %>% 
  mutate("NOB" = if_else(`NOB:In situ` %in% c("neg") | `NOB:Other` %in% c("neg"), "negative", "unknown"),
         "NOB" = if_else(`NOB:In situ` %in% c("var") | `NOB:Other` %in% c("var"), "varible", `NOB`),
         "NOB" = if_else(`NOB:In situ` %in% c("pos") | `NOB:Other` %in% c("pos"), "positive", `NOB`)) %>% 
  select(Genus, "NOB")

# anti_join(NOB, AOB)
# see <- inner_join(NOB, AOB) %>% mutate(E = if_else(NOB == AOB, T, F)) %>% arrange()

SOB_vec <- c("g__Thiobacillus", "g__Thiothrix", "g__Sulfolobus", "g__Beggiatoa", "g__Thermothrix")

SOB <- meta_midas %>% 
  select(Genus) %>% 
  mutate(SOB = 
           if_else(Genus %in% SOB_vec, "positive", "unknown")) %>% 
  arrange(SOB)

metabolic_by_genus <- 
  full_join(sulfate_red, fermentation, by = c("Genus")) %>% 
  full_join(., Aero_hetero, by = c("Genus")) %>% 
  full_join(., chemo_auto, by = c("Genus")) %>% 
  full_join(., NOB, by = c("Genus")) %>% 
  full_join(., AOB, by = c("Genus")) %>% 
  full_join(., Filamentous, by = c("Genus")) %>% 
  full_join(., Nitrite, by = c("Genus")) %>% 
  full_join(., gao, by = c("Genus")) %>% 
  full_join(., pao, by = c("Genus")) %>% 
  full_join(., SOB, by = c("Genus")) %>% 
  pivot_longer(-Genus, names_to = "meta") %>% 
  mutate(
    meta = if_else(is.na(meta), "unknown", meta),
    value = if_else(is.na(value), "unknown", value))


metabolic_by_genus

}
