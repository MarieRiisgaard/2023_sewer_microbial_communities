

#### Make the growing/dying segment plots


# library(ampvis2)
# library(rmarkdown)
# library(patchwork)
# library(tidyverse)
# library(readxl)
# library(scales)
# library(reshape2)


growing_dying_bar_plot <- 
  function(data = data_long, 
         include = c("Biofilm (g.)", "Biofilm (e.p.)", "Sediment", "SWW","IWW", "AS"), 
         #gut_bacteria_as_group = T, 
         #gut_df = gut_OTU, 
         abun_limit = 0.6,
         make_label = T
){
  
  
  n = length({{include}})
  
 
    growth_levels = c("unknown",'disappearing','surviving/variable','growing') 
    growth_labels = c("Unknown",'Disappearing', 'Surviving/variable','Growing')
    
    col <- c("grey", "#2E7052", "#6D9E77", "#A7C7A3")
    
    # define colors
    #col = c("grey",'#ea9999','#ffe599',"lightblue",'#b6d7a8')
    #col = c("#BADDCF","#BADDCF",'#BADDCF','#BADDCF',"#BADDCF",'#BADDCF')
    
  growth_df_usearch <- growth_df  
  
df_groups <- {{data}} %>% 
  filter(SampleContent2 %in% {{include}}) %>% 
  droplevels() %>% 
  mutate(samples = 
           map(.x = samples, ~
                 ungroup(.) %>% 
                 left_join(growth_df_usearch, by = c("OTU")) %>% 
                 mutate(growth_AS = if_else(growth_AS %in% c("surviving", "variable"), "surviving/variable", growth_AS)) %>% 
                 #mutate(growth_AS = if_else(OTU %in% unlist(gut_OTU$OTU) & gut_bacteria_as_group, "Gut bacteria", growth_AS)) %>% 
                 #mutate(growth_AS = as.character(growth_AS), growth_AS = "Inconclusive") %>% # Changed
                 group_by(growth_AS) %>% 
                 filter(rel_abun_ASV != 0) %>% 
                 summarise(sum_growth = sum(rel_abun_ASV), 
                           n_species = n_distinct(Species),
                           .groups = "drop")
           )) %>% 
  unnest(samples) %>% 
  group_by(SampleContent2, growth_AS#, Configuration
           ) %>% 
  summarise(mean_abun = mean(sum_growth), 
            mean_n_species = mean(n_species), 
            .groups = "drop") %>% 
  mutate(growth_AS = factor(growth_AS,
                            levels = growth_levels, ## changed
                            labels = growth_labels  ## Changed
                            ),
         #SampleContent2 = factor(SampleContent2, c("biofilm","biofilm_end_pressure", "sediment","sewage","IWW", "AS"), 
         #                       labels = c("Biofilm (g)", "Biofilm (p)", "Sediment", "SWW","IWW", "AS")),
  )

# Desription of the dataframe: df_groups; 
# It is in a long format suitable for ggplot. 
# The relative abundance are summarised based on growth in AS (with mean for samples across sampletypes)  


# I make a dublicate dataframe
# IMPORTANT: I make the SampleContent numeric!  
make_dub <- df_groups %>% 
  group_by(SampleContent2#, Configuration
           ) %>% 
  arrange(SampleContent2,desc(growth_AS)) %>% 
  mutate(mean_abun_sum = cumsum(mean_abun), 
         SampleContent_num = as.numeric(SampleContent2), 
         x = "x")

make_dub1 <- mutate(make_dub, SampleContent_num = SampleContent_num - 0.7)

df_for_plot <- 
  rbind(make_dub, make_dub1)


#indcludes <- c("biofilm", "biofilm_end_pressure", "sediment","sewage","IWW", "AS")
#labels <-  c("Biofilm (g.)", "Biofilm (e.p.)", "Sediment", "SWW","IWW", "AS")

xintercept_vec <- c(1, 1.3, 2, 2.3, 3, 3.3, 4, 4.3, 5, 5.3,6, 6.3,7, 7.3,8, 8.3,9, 9.3,
                    10, 10.3,11, 11.3,12, 12.3,13, 13.3,14, 14.3,15, 15.3)[1:(n*2-2)]
breaks_vec = c(0.7, 1.7, 2.7, 3.7, 4.7, 5.7, 6.7, 7.7, 8.7, 9.7, 10.7, 11.7, 12.7, 13.7, 14.7, 15.7, 16.7)[1:(n)]

plot <- df_for_plot %>% 
  ggplot(aes(x=SampleContent_num, y=mean_abun, fill=growth_AS)) + 
  geom_area(linewidth = 0.5, color = "grey30", alpha = 0.8) + 
  geom_vline(xintercept = xintercept_vec, linetype = "dashed", color = "grey30", linewidth = 0.3) +
  geom_label(
    data = make_dub %>% 
      filter(mean_abun > abun_limit | growth_AS == "Growing") %>% 
      filter(SampleContent2 %in% {{include}}),
    aes(
      y = mean_abun_sum - mean_abun*0.6, 
      # y = case_when(mean_abun > 2 | SampleContent2 == "AS" ~ mean_abun_sum - mean_abun*0.6, 
      #               
      #                   mean_abun <= 3 & growth_AS == "Growing" #& SampleContent2 != "AS"
      #               ~ 0,
      #                   mean_abun <= 3 & growth_AS == "Inconclusive" #& SampleContent2 != "AS" 
      #               ~ mean_abun_sum - 2,
      #                   mean_abun <= 2 & growth_AS == "Surviving" #& SampleContent2 != "AS"
      #               ~ mean_abun_sum + 2,
      # #                  mean_abun <= 2 & growth_AS == "Gut bacteria" & SampleContent2 != "AS"~ mean_abun_sum + 2,
      # #                  mean_abun <= 2 & growth_AS == "Disappearing" & SampleContent2 != "AS"~ mean_abun_sum + 2,
      #                .default = mean_abun_sum - mean_abun*0.8
    #),
    x = SampleContent_num - 0.3,
    label = ifelse(((mean_abun > abun_limit & make_label) | (growth_AS == "Growing" & mean_abun >= 0.1)), 
                   paste0(round(mean_abun, 1), "%(", round(mean_n_species, 0), ")"), NA),
    color = case_when(str_detect(growth_AS, "Gut")~"white", 
                      str_detect(growth_AS, "isappearing")~"white", 
                      str_detect(growth_AS, "urviving")~"white", 
                      .default = "black"),
    ), 
    show.legend = F, vjust = 0, size = 4, hjust = 0.6,
    label.padding = unit(0.2, "mm")) +
  scale_fill_manual(values = col) +
  scale_color_manual(values = c("black","grey90")) + 
  scale_y_continuous(expand = expansion(mult = c(0, .03)), 
                                      ) +
  ylab("Cummulative relative abundance [%]") +
  scale_x_continuous(breaks = breaks_vec,labels = {{include}}, expand = c(0,0)) + 
  theme(
    axis.title.x = element_blank(), 
    plot.background = element_rect(fill = "white"), 
    panel.background = element_rect(fill = "NA"), 
    axis.text.x = element_text(size = 15, color = "black"), 
    legend.text = element_text(size = 15, color = "black"),
    axis.text.y = element_text(size = 15, color = "black"), 
    axis.title.y = element_text(size = 15, color = "black"),
    axis.ticks.x = element_blank(), #legend.title = element_blank()
  )

plot

}




