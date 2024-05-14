# 08 - Maps and Project Description Table
# Katy Dynarski

# 0 - Load data ----
project <- read.csv(here("data_processed", "05_project_data.csv"))

# 1 - Map of all DSP4SH projects ----
# Download map of USA from NaturalEarth as basemap
usa_ne <- ne_states(country="united states of america")

# Make dataframe for adding project data annotation to map
# What I want - just one point for each project label combination, text label describing project
# Many projects have more than one type of treatment within a label - look at the data to be able to write a description that generally covers things: 
project_distinct <- project %>%
  distinct(project, label, lu, till, trt)
flextable(project_distinct)

# Write tibble of annotation for each project/label based on project_distinct
annotation <- tibble(project=c(rep("Illinois", 3), 
                               rep("KansasState", 3),
                               rep("NCState", 3),
                               rep("OregonState", 3),
                               rep("TexasA&MPt-1", 3),
                               rep("TexasA&MPt-2", 3),
                               rep("UConn", 3),
                               rep("UTRGV", 2),
                               rep("UnivOfMinnesota", 3),
                               rep("WashingtonState", 3)),
                     label=c(rep(c("BAU", "Ref", "SHM"), 7), "BAU", "Ref", rep(c("BAU", "Ref", "SHM"), 2)),
                     annotation=c("BAU: Conventional and organic corn-soybean rotation",
                                  "Ref: Forest restored in 1990",
                                  "SHM: No-till corn-soybean rotation",
                                  "BAU: Diverse crops with conventional tillage",
                                  "Ref: Native rangeland",
                                  "SHM: No-till diverse crops and cover cropping",
                                  "BAU: Wheat and corn with conventional tillage",
                                  "Ref: Forest",
                                  "SHM: Hayed perennial grass",
                                  "BAU: Ryegrass and Christmas tree with tillage",
                                  "Ref: Timber forest and hazelnut orchard",
                                  "SHM: No-till grass and vineyard",
                                  "BAU: Cropping with conventional tillage",
                                  "Ref: Native rangeland",
                                  "SHM: No-till rye and mixed crops",
                                  "BAU: Wheat and sorghum with conventional tillage",
                                  "Ref: Native perennial forage",
                                  "SHM: No-till wheat and sorghum",
                                  "BAU: Corn with conventional tillage",
                                  "Ref: Forest",
                                  "SHM: No-till corn and hay",
                                  "BAU: Diverse crops with conventional tillage",
                                  "Ref: Native forest",
                                  "BAU: Soybean with conventional tillage",
                                  "Ref: Native rangeland",
                                  "SHM: No-till soybean",
                                  "BAU: Cropping with conventional tillage",
                                  "Ref: Perennial grassland",
                                  "SHM: No-till cropping"))

# Get unique soil series
project_soil <- project %>%
  distinct(project, soil) %>%
  group_by(project) %>%
  mutate(count = paste("soil",seq(n()), sep="_")) %>%
  pivot_wider(names_from=count, values_from=soil) %>%
  unite("soils", soil_1:soil_3, sep=", ", na.rm=TRUE)
flextable(project_soil)

# Make annotation dataframe that includes project, label, soil, and xy coords
project_annotate <- project %>% 
  group_by(project, label) %>%
  mutate(avg_lat = mean(pedon_y, na.rm=TRUE),
         avg_long = mean(pedon_x, na.rm=TRUE)) %>%
  distinct(project, label, avg_lat, avg_long) %>%
  left_join(annotation, by=c("project", "label")) %>%
  left_join(project_soil, by="project")

# Make a dataframe with only one point per project to make the map labels
project_labs <- project_annotate %>%
  group_by(project) %>%
  filter(label=="BAU") %>%
  mutate(name_long = ifelse(project=="Illinois", "University of Illinois",
                            ifelse(project=="OregonState", "Oregon State University",
                                   ifelse(project=="UConn", "University of Connecticut",
                                          ifelse(project=="UTRGV", "University of Texas - Rio Grande Valley",
                                                 ifelse(project=="KansasState", "Kansas State University",
                                                        ifelse(project=="NCState", "North Carolina State University",
                                                               ifelse(project=="TexasA&MPt-1", "Texas A&M - 1",
                                                                      ifelse(project=="TexasA&MPt-2", "Texas A&M - 2",
                                                                             ifelse(project=="UnivOfMinnesota", "University of Minnesota", "Washington State University"))))))))))

# Make map
ggplot(data=usa_ne) +
  geom_sf(fill=NA) +
  coord_sf(xlim=c(-125.0, -66.93457), ylim=c(23.5, 49.384358)) + # set bounding box around CONUS
  annotation_north_arrow(location="bl", which_north="true", height=unit(.25, "in"), width=unit(.25, "in"),
                         pad_x = unit(0.4, "in"), pad_y = unit(0.25, "in"), style=north_arrow_fancy_orienteering) + # add north arrow
  annotation_scale(location = "bl") +
  geom_point(data=project_annotate, aes(x=avg_long, y=avg_lat, color=project)) +
  geom_label_repel(data=project_labs, aes(x=avg_long, y=avg_lat, label=name_long),
                   min.segment.length = 0, seed = 42, box.padding = 0.5) +
  scale_color_viridis(discrete=TRUE) +
  theme_classic() +
  easy_remove_axes() +
  theme(legend.position="none")
ggsave(here("figs", "project_map.png"), height=5, width=8, units="in", dpi=400)

# make version of the map with no annotation - so I can manually annotate with project information in Powerpoint
ggplot(data=usa_ne) +
  geom_sf(fill=NA) +
  coord_sf(xlim=c(-125.0, -66.93457), ylim=c(23.5, 49.384358)) + # set bounding box around CONUS
  annotation_north_arrow(location="bl", which_north="true", height=unit(.25, "in"), width=unit(.25, "in"),
                         pad_x = unit(0.4, "in"), pad_y = unit(0.25, "in"), style=north_arrow_fancy_orienteering) + # add north arrow
  annotation_scale(location = "bl") +
  geom_point(data=project_labs, aes(x=avg_long, y=avg_lat, color=project)) +
  scale_color_viridis(discrete=TRUE) +
  theme_classic() +
  easy_remove_axes() +
  theme(legend.position="none")
ggsave(here("figs", "project_map_no_labs.png"), height=5, width=8, units="in", dpi=600)

# 2 - Table of project information ----
# Average climate data
site_clim_sum <- project %>%
  group_by(project) %>%
  summarize(across(mat:map, ~ mean(.x, na.rm = TRUE))) %>%
  mutate(mat = round(mat, 1),
         map = round(map, 0))

project_table <- project_annotate %>%
  group_by(project) %>%
  select(project, soils, annotation) %>%
  separate_wider_delim(annotation, delim=": ", names=c("label", "description")) %>%
  ungroup() %>%
  left_join(site_clim_sum, by="project") %>%
  arrange(project, label)
flextable(project_table)
write_csv(project_table, here("figs", "project_descriptions.csv"))

# alternative version of project information that has details on each treatment- maybe for supplement?
treatment_table <- meta_df %>%
  select(project, label, lu, till, trt, explanation) %>%
  distinct() %>%
  arrange(project)
write_csv(treatment_table, here("figs", "project_treatment_table.csv"))
