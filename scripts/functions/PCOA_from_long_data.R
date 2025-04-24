##############

# My own PCA 

# 
# 
# ####### Species ################
# 
#   df_tax_level <- data_long_merged %>% 
#   mutate(samples = 
#            map(.x = samples, ~
#                  ungroup(.) %>% 
#                  distinct(Species, rel_abun_species) %>% 
#                  filter(!str_detect(Species, "unclassified")) %>% 
#                  mutate(over_limit = if_else(rel_abun_species >= 0.01, T, F))
#                  )) %>% 
#   unnest(samples) %>% 
#   rename(
#     "tax_level" = "Species",
#     "tax_abun" = "rel_abun_species"
#   )
# 
# 
# ####### GENUS ################
# 
# ## Genus level: 
# df_tax_level <- data_long_merged %>% 
#   mutate(samples = 
#            map(.x = samples, ~
#                  ungroup(.) %>% 
#                  distinct(Genus, rel_abun_genus) %>% 
#                  filter(!str_detect(Genus, "unclassified")) %>% 
#                  mutate(over_limit = if_else(rel_abun_genus >= 0.01, T, F))
#            )) %>% 
#   unnest(samples) %>% 
#   rename(
#     "tax_level" = "Genus",
#     "tax_abun" = "rel_abun_genus"
#   )
# 
# 
# #######
# 
# ## Filter low abundant Species
# tax_over_limit <- 
#   df_tax_level %>% 
#   group_by(tax_level) %>% 
#   reframe(sum = sum(over_limit)) %>% 
#   filter(sum > 0) %>% 
#   distinct(tax_level) %>% unlist()
# 
# df_tax_level_filtered <- 
#   df_tax_level %>% 
#   filter(tax_level %in% tax_over_limit)
# 
# 
# ####### GENUS END ################
# 
# 
# ## Calculate distance matrix 
# abun_matrix <- df_tax_level_filtered %>% 
#   select(SampleID, tax_level, tax_abun) %>% 
#   pivot_wider(names_from = tax_level, values_from = tax_abun, values_fill = 0) %>% 
#   column_to_rownames("SampleID") %>% 
#   as.matrix()
# 
# bray_curtis_dist <- vegan::vegdist((abun_matrix), method = "bray")
# 
# 
# ## Making PCoA plot 
#  ### Form Code club
#  ### https://github.com/riffomonas/distances/blob/f1d31bf29443ee4986e16ea045e35946a3e06bf5/code/pcoa.R
# pcoa <- cmdscale(bray_curtis_dist, eig=TRUE, add=TRUE)
# 
# positions <- pcoa$points
# colnames(positions) <- c("pcoa1", "pcoa2")
# 
# percent_explained <- 100 * pcoa$eig / sum(pcoa$eig)
# 
# pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)
# 
# labels <- c(glue::glue("PCo Axis 1 ({pretty_pe[1]}%)"),
#             glue::glue("PCo Axis 2 ({pretty_pe[2]}%)"))
# 
# positions %>%
#   as_tibble(rownames = "SampleID") %>%
#   left_join(df_tax_level %>% distinct(SampleID, country, WWTP_ID, ProjectName)) %>%
#   mutate(xx = paste0(country, ProjectName)) %>% 
#   filter(country != "United States" & country != "USA") %>% 
#   ggplot(aes(x=pcoa1, y=pcoa2, color = country)) +
#   geom_point() +
#   geom_text(aes(label = ProjectName)) + 
#   labs(x=labels[1], y=labels[2])
# 
# 
# df_tax_level %>% distinct(SampleID, country, WWTP_ID, ProjectName) %>% 
#   filter(str_detect(country, "Argen")) %>% 
#   print(n=100)
# 
# positions %>%
#   as_tibble(rownames = "SampleID") %>%
#   left_join(df_tax_level %>% distinct(SampleID, country, WWTP_ID, ProjectName)) %>%
#   filter(WWTP_ID %in% WWTP_ID_common | ProjectName == "RETHINK" | country == "Germany") %>% 
#   mutate(xx = paste0(country, ProjectName)) %>% 
#   filter(country != "United States" & country != "USA") %>% 
#   ggplot(aes(x=pcoa1, y=pcoa2, color = country)) +
#   geom_point() +
#   geom_text(aes(label = ProjectName)) + 
#   labs(x=labels[1], y=labels[2])
# 

# 
# tibble(pe = cumsum(percent_explained),
#        axis = 1:length(percent_explained)) %>%
#   ggplot(aes(x=axis, y=pe)) +
#   geom_line() +
#   coord_cartesian(xlim = c(1, 10), ylim=c(0, 50)) +
#   scale_x_continuous(breaks=1:10)

###


####################################

# PCoA from long dataframe

PCOA_from_long_df <- function(
  data = data_long,
  tax_level = "Species", 
  abun_limit = 0.1, 
  color_by = "country", 
  extra_column = "ProjectName",
  type = "PCOA"
  ){
  
  
  
  if (tax_level == "Species") {
    df_tax_level <- {{data}} %>% 
      mutate(samples = 
               map(.x = samples, ~
                     ungroup(.) %>% 
                     distinct(Species, rel_abun_species) %>% 
                     filter(!str_detect(Species, "unclassified")) %>% 
                     mutate(over_limit = if_else(rel_abun_species >= abun_limit, T, F))
               )) %>% 
      unnest(samples) %>% 
      rename(
        "tax_level" = "Species",
        "tax_abun" = "rel_abun_species"
      )
  } else if (tax_level == "Genus") {
    df_tax_level <- {{data}} %>% 
      mutate(samples = 
               map(.x = samples, ~
                     ungroup(.) %>% 
                     distinct(Genus, rel_abun_genus) %>% 
                     filter(!str_detect(Genus, "unclassified")) %>% 
                     mutate(over_limit = if_else(rel_abun_genus >= abun_limit, T, F))
               )) %>% 
      unnest(samples) %>% 
      rename(
        "tax_level" = "Genus",
        "tax_abun" = "rel_abun_genus"
      )
  } else if (tax_level == "Family") {

    df_tax_level <- {{data}} %>% 
      ungroup() %>% 
      mutate(samples = map(.x = samples, ~group_by(.x, Family) %>% mutate(rel_abun_family = sum(rel_abun_ASV)) %>% ungroup())) %>% 
      mutate(samples = 
               map(.x = samples, ~
                     ungroup(.) %>% 
                     distinct(Family, rel_abun_family) %>% 
                     filter(!str_detect(Family, "unclassified")) %>% 
                     mutate(over_limit = if_else(rel_abun_family >= abun_limit, T, F))
               )) %>% 
      unnest(samples) %>% 
      rename(
        "tax_level" = "Family",
        "tax_abun" = "rel_abun_family"
      )
  } else if (tax_level == "Order") {
    
    df_tax_level <- {{data}} %>% 
      ungroup() %>% 
      mutate(samples = map(.x = samples, ~group_by(.x, Order) %>% mutate(rel_abun_order = sum(rel_abun_ASV)) %>% ungroup())) %>% 
      mutate(samples = 
               map(.x = samples, ~
                     ungroup(.) %>% 
                     distinct(Order, rel_abun_order) %>% 
                     filter(!str_detect(Order, "unclassified")) %>% 
                     mutate(over_limit = if_else(rel_abun_order >= abun_limit, T, F))
               )) %>% 
      unnest(samples) %>% 
      rename(
        "tax_level" = "Order",
        "tax_abun" = "rel_abun_order"
      )
  } else if (tax_level == "Class") {
    
    df_tax_level <- {{data}} %>% 
      ungroup() %>% 
      mutate(samples = map(.x = samples, ~group_by(.x, Class) %>% mutate(rel_abun_class = sum(rel_abun_ASV)) %>% ungroup())) %>% 
      mutate(samples = 
               map(.x = samples, ~
                     ungroup(.) %>% 
                     distinct(Class, rel_abun_class) %>% 
                     filter(!str_detect(Class, "unclassified")) %>% 
                     mutate(over_limit = if_else(rel_abun_class >= abun_limit, T, F))
               )) %>% 
      unnest(samples) %>% 
      rename(
        "tax_level" = "Class",
        "tax_abun" = "rel_abun_class"
      )
  } else if (tax_level == "Phylum") {
    
    df_tax_level <- {{data}} %>% 
      ungroup() %>% 
      mutate(samples = map(.x = samples, ~group_by(.x, Phylum) %>% mutate(rel_abun_phylum = sum(rel_abun_ASV)) %>% ungroup())) %>% 
      mutate(samples = 
               map(.x = samples, ~
                     ungroup(.) %>% 
                     distinct(Phylum, rel_abun_phylum) %>% 
                     filter(!str_detect(Phylum, "unclassified")) %>% 
                     mutate(over_limit = if_else(rel_abun_phylum >= abun_limit, T, F))
               )) %>% 
      unnest(samples) %>% 
      rename(
        "tax_level" = "Phylum",
        "tax_abun" = "rel_abun_phylum"
      )
  } else if (tax_level == "OTU") {
      df_tax_level <- {{data}} %>% 
        ungroup() %>% 
        mutate(samples = 
                 map(.x = samples, ~
                       ungroup(.) %>% 
                       distinct(OTU, rel_abun_ASV) %>% 
                       mutate(over_limit = if_else(rel_abun_ASV >= abun_limit, T, F))
                 )) %>% 
        unnest(samples) %>% 
        rename(
          "tax_level" = "OTU",
          "tax_abun" = "rel_abun_ASV"
        )
  } else {
    print("Check spelling of tax_level ")
  }
  
  
## Filter low abundant Species
tax_over_limit <- 
  df_tax_level %>% 
  group_by(tax_level) %>% 
  reframe(sum = sum(over_limit)) %>% 
  filter(sum > 0) %>% 
  distinct(tax_level) %>% unlist()

df_tax_level_filtered <- 
  df_tax_level %>% 
  filter(tax_level %in% tax_over_limit) %>% 
  rename("color" = {{color_by}}, 
         "extra_column_name" = {{extra_column}})
  

## Calculate distance matrix 
abun_matrix <- df_tax_level_filtered %>%
  distinct(SampleID, tax_level, tax_abun) %>% 
  pivot_wider(names_from = tax_level, values_from = tax_abun, values_fill = 0) %>% 
  column_to_rownames("SampleID") %>% 
  as.matrix()

bray_curtis_dist <- vegan::vegdist((abun_matrix), method = "bray")



## Making PCoA plot 
### Form Code club
### https://github.com/riffomonas/distances/blob/f1d31bf29443ee4986e16ea045e35946a3e06bf5/code/pcoa.R
pcoa <- cmdscale(bray_curtis_dist, eig=TRUE, add=TRUE)

positions <- pcoa$points #%>% as_tibble()

colnames(positions) <- c("pcoa1", "pcoa2")

percent_explained <- 100 * pcoa$eig / sum(pcoa$eig)

pretty_pe <- format(round(percent_explained, digits =1), nsmall=1, trim=TRUE)

labels <- c(glue::glue("PCo1 ({pretty_pe[1]}%)"),
            glue::glue("PCo2 ({pretty_pe[2]}%)"))

plot <- positions %>%
  as_tibble(rownames = "SampleID") %>% 
  distinct() %>% 
  left_join(df_tax_level_filtered %>% distinct(SampleID, color, extra_column_name) 
            )%>%
  ggplot(aes(x=pcoa1, y=pcoa2, color = color)) +
  geom_point() +
  labs(x=labels[1], y=labels[2])

plot



if(type == "RDA"){
  #type = "PCOA"
  
  abun.helliger.transformed <- decostand(abun_matrix, "hellinger")
  
  spe.rda <- rda(abun.helliger.transformed ~ ., 
                 data = {{data}} %>% 
                   column_to_rownames("SampleID") %>% 
                   rename("color" = {{color_by}}) %>% 
                   select(color#, time_before_rain_event#, SampleContent2
                          ))
  
  # summary(spe.rda)
  # plot(spe.rda)
  # RsquareAdj(spe.rda) 
  # anova.cca(spe.rda, permutations = 999)
  #anova.cca(spe.rda, by = "axis")  # Test which axis are significant
  #anova.cca(spe.rda, by = "terms")  # Test which terms are significant
  
}

if(type == "db-RDA"){
  #type = "PCOA"
  
  #abun.helliger.transformed <- decostand(abun_matrix, "hellinger")
  
  db.rda <- dbrda(formula = abun_matrix ~ ., 
        data = {{data}} %>% 
          column_to_rownames("SampleID") %>% 
          rename("color" = {{color_by}}) %>% 
          select(color#, time_before_rain_event#, SampleContent2
          ), 
        distance = "bray")
  
  # spe.rda <- rda(abun.helliger.transformed ~ ., 
  #                data = {{data}} %>% 
  #                  column_to_rownames("SampleID") %>% 
  #                  rename("color" = {{color_by}}) %>% 
  #                  select(color#, time_before_rain_event#, SampleContent2
  #                  ))
  
  # summary(spe.rda)
  # plot(spe.rda)
  # RsquareAdj(spe.rda) 
  # anova.cca(spe.rda, permutations = 999)
  #anova.cca(spe.rda, by = "axis")  # Test which axis are significant
  #anova.cca(spe.rda, by = "terms")  # Test which terms are significant
  
}

if(type == "db-RDA_and_season"){
  #type = "PCOA"
  
  #abun.helliger.transformed <- decostand(abun_matrix, "hellinger")
  
  db.rda <- dbrda(formula = abun_matrix ~ ., 
                  data = {{data}} %>% 
                    column_to_rownames("SampleID") %>% 
                    rename("color" = {{color_by}}) %>% 
                    select(color, Season#, time_before_rain_event#, SampleContent2
                    ), 
                  distance = "bray")
  
 
}





if(type == "RDA") {
  spe.rda
} else if(type == "db-RDA"| type == "db-RDA_and_season") {
  db.rda
} else {
  plot
}


}

# PCOA_from_long_df(
#   data = data_long_merged_all %>% filter(country != "Hong Kong"),
#     tax_level = "Species",
#     abun_limit = 0.1,
#     color_by = "country",
#     extra_column = "ProjectName")

  
