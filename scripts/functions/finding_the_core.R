#####################################
##   Make sewer core groups
#####################################


sewer_core <- function(
    long_df = data_long,
    tax_level = Species,
    tax_rel_abun = rel_abun_species,
    abun_limit = 0.01
  ){
  core_with_species_names <- 
    long_df %>% 
    #sample_n(10) %>% 
    group_by(SampleContent2) %>% 
    mutate(n_samples = n()) %>% 
    ungroup() %>% 
    mutate(samples = 
             map(.x = samples, ~
                   ungroup(.) %>% 
                   distinct(across(c({{tax_level}}, {{tax_rel_abun}}))) %>% 
                   filter({{tax_rel_abun}} >= abun_limit) 
             )) %>% 
    unnest(samples) %>%
    distinct() %>% 
    group_by(across(c({{tax_level}}, SampleContent2))) %>%
    reframe(n_obs_tax = n(),
              is_core = ifelse(n_obs_tax >= 0.8*n_samples, T, F), 
            core_type = case_when(n_obs_tax >= 0.8*n_samples ~ "strict core",
                                  n_obs_tax >= 0.5*n_samples ~ "general core",
                                  n_obs_tax >= 0.2*n_samples & 0.2*n_samples > 1 ~ "loose core", 
                                  n_obs_tax < 0.2*n_samples ~ "detected", 
                                  .default = "not detected"),
            n_samples)  %>% 
    distinct()
    
  # 80% (strict core), 50% (general core) and 20% (loose core)
  
  summarised_core <- 
    core_with_species_names %>% 
    group_by(SampleContent2, is_core, n_samples) %>% 
    summarise(n_core = n(), .groups = "drop")
  
  
  list(core_with_species_names, summarised_core)
}

upset_of_core <- function(
    core_taxa_df = sewer_core_species_names,
    tax_level_name = "Species",
    tax_level = Species
    ){
  
  matrix_w_tax <- 
    core_taxa_df %>% #sample_n(40) %>% 
    select(-n_samples, -n_obs_tax) %>% 
    filter(!(str_detect({{tax_level}}, "unclass") & tax_level_name == "Genus")) %>% 
    # group_by(across(c(group, Species))) %>% 
    # summarise(
    #   #group_name = paste0(unique(SampleContent), collapse = "_"), 
    #   mean = mean(rel_abun_species), .groups = "drop") %>% 
    #        #group_name = as.factor(group_name)
    #        ) %>% 
    #select(-tax_level) %>% 
    #pivot_wider(names_from = tax_level, values_from = tax_name)
    complete({{tax_level}}, SampleContent2, fill = list(is_core = F)) %>% 
  # %>% 
  #   mutate(present = if_else(!is.na({{tax_level}}), T, F)) %>%  
  #   select(#-{{tax_level}},
  #          -is_core, -n_obs_tax, -n_samples) %>% 
   pivot_wider(names_from = SampleContent2, values_from = is_core)
  # # 
  
  matrix_w_tax

  matrix_w_tax_filterd <- matrix_w_tax[rowSums(matrix_w_tax[,-1]) >= 1, ]
  matrix_wo_tax <-  matrix_w_tax %>% select(-{{tax_level}})
  matrix_wo_tax <- matrix_wo_tax[rowSums(matrix_wo_tax) >= 1, ]
  names <- colnames(matrix_wo_tax)
  n_groups = length(names)
  
  # Make uset plot
  plot <- ComplexUpset::upset(matrix_wo_tax, names, min_size = 5,
                              base_annotations=list('Intersection size'= 
                                                      ComplexUpset::intersection_size(
                                counts=T, text = list(size = 4, color = "black"),
                                text_colors = c("black", "black", "black")
                                # ,
                                # mapping=aes(fill=core_type)
                              )),
                              width_ratio=0.3, height_ratio = 0.3, # 0.5
                              sort_sets=F, guides='over',
                              themes=ComplexUpset::upset_default_themes(text=
                                                                          element_text(size = 12, color = "black"), 
                                                                        title = element_text(color = "black"), 
                                                                        axis.text = element_text(color = "black")))
  
  shared <- matrix_w_tax_filterd %>%
    mutate(sum_logical =  rowSums(across(where(is.logical)))) %>% 
    mutate(across(where(is.logical), ~ifelse(. == TRUE, sum_logical, .))) %>% 
    pivot_longer(cols = c(2:4), names_to = "group_name", values_to = "shared_with") %>% 
    select({{tax_level}}, group_name, shared_with) %>% 
    mutate(core_group = 
             case_when(shared_with == 3 ~ "all types", 
                       shared_with == 2 ~ "two types", 
                       shared_with == 1 ~ paste0("unique ", str_sub(group_name, start = 0, end = 8)), 
                       shared_with == 0 ~ "not found")
    )
  
  
  
  list(plot, shared)
  
}











