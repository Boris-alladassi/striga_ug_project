setwd("C:/Users/aboris/Box/Postdoc_UIUC/Truc_perso/Striga_analysis")
rm(list = ls())

library(tidyverse)
library(data.table)
library(ggpubfigs)

boris_theme_rect <- theme_bw() + 
  theme(axis.text = element_text(color = "black", size = 12),
        axis.title = element_text(colour = "black", face = "bold", size = 13),
        #axis.line = element_line(size = 1, colour = "black"),
        panel.grid = element_blank(),
        panel.border = element_rect(linewidth = 1.2, color = "black"),
        axis.ticks = element_line(linewidth = 1, colour = "black"),
        legend.text = element_text(colour = "black", size = 12),
        legend.position = "bottom", legend.direction = "horizontal",
        legend.title = element_blank(),
        text = element_text(colour = "black", size = 12))

transpoz <- function(data, new_col1_name = NULL){
  require(tidyverse)
  name <- names(data)[1]
  data <- column_to_rownames(data, var = name[1])
  data <- as.data.frame(t(as.matrix(data)))
  if(!(is.null(new_col1_name))){name <- new_col1_name}
  data <- rownames_to_column(data, var = name[1])
  return(data)
}

difference <- function(a,b){
  res <- c(setdiff(a,b), setdiff(b,a))
  return(res)
}

### A function to convert the marker data from 0,1,2 to ==> -1,0,1 format.
snp_dt_012_to_101 <- function(x){
  y = rep(0,length(x))
  y[which(x == 0)] <- -1
  y[which(x == 1)] <- 0
  y[which(x == 2)] <- 1
  y = matrix(y, nrow = nrow(x), ncol = ncol(x), byrow = F)
  rownames(y) <- rownames(x); colnames(y) <- colnames(x)
  return(y)
}

############################################################################
##### Phenotypic data analysis ============
pheno.dt <- readxl::read_xlsx("./data/Phenotypic_data.xlsx")
dim(pheno.dt)
head(pheno.dt)
colnames(pheno.dt) <- c("Location", "IID", "Flower_col", "Leaf_col", "Leaf_pigm", "Leaf_texture", "Leaf_size")
pheno.dt <- pheno.dt %>%
  mutate(Location = recode(Location, KUM = "Kumi", NGO = "Ngora", BUK = "Bukedea", PAL = "Pallisa",
                           TOR = "Tororo", BUS = "Busia", NAM = "Namutumba", IGA = "Iganga",
                           KAL = "Kalaki", KAB = "Kaberamaido", SER = "Serere"))


write.csv(pheno.dt, "./data/pheno_dt_final.csv", row.names = F)

# "#648FFF" "#785EF0" "#DC267F" "#FE6100" "#FFB000"
# "#332288" "#117733" "#CC6677" "#88CCEE" "#999933" "#882255" "#44AA99" "#DDCC77" "#AA4499"
loc_colors <- c(ggpubfigs::friendly_pal("muted_nine"), "#FE6100", "grey50")
# loc_order <- pheno.dt %>% count(Location) %>% arrange(desc(n)) %>% pull(Location)
loc_order <- c("Iganga", "Namutumba", "Ngora", "Serere", "Kaberamaido", "Pallisa",
                "Bukedea", "Kumi", "Tororo", "Busia", "Kalaki")

locs_summary <- pheno.dt %>% count(Location) %>% 
  mutate(Location = factor(Location, levels = loc_order))

p_loc <- ggplot(locs_summary, aes(x = reorder(Location, -n), y = n, fill = Location)) + 
  geom_bar(stat = "identity", color = "black") + 
  geom_text(aes(label = n), vjust = -0.2) + 
  labs(x = "Location", y = "Number of samples") + boris_theme_rect +
  scale_fill_manual(values = loc_colors)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")

ggsave("./output/Location_summary_stats.jpg", width = 4, 
       height = 5, dpi = 600, units = "in", plot = p_loc)

#############################################################
#### Creating a map with the locations ============
library(sf)
# Reading Uganda district shapefile
uga <- st_read("./data/gadm41_UGA_shp/gadm41_UGA_4.shp")

locations <- c("Kumi","Ngora","Bukedea","Pallisa","Tororo","Busia",
               "Kalika", "Namutumba","Iganga","Kaberamaido","Serere")

uga_sel_1 <- uga[uga$NAME_1 %in% loc_order, ]
sel1 <- unique(uga$NAME_1[which(uga$NAME_1 %in% loc_order)])
uga_sel_2 <- uga[uga$NAME_2 %in% loc_order, ]
uga_sel_3 <- uga[uga$NAME_3 %in% loc_order, ]
uga_sel_4 <- uga[uga$NAME_4 %in% loc_order, ]
# uga_sel_2 <- uga[uga$NAME_4 %in% loc_order[!(locations %in% unique(uga_sel_4$NAME_4))], ]
# uga_sel <- rbind(uga[uga$NAME_2 %in% locations, ], uga[uga$NAME_4 %in% locations, ])
col1 <- loc_colors[which(loc_order %in% unique(uga_sel_1$NAME_1))]
col2 <- loc_colors[which(loc_order %in% unique(uga_sel_2$NAME_2))]
col3 <- loc_colors[which(loc_order %in% unique(uga_sel_3$NAME_3))]
col4 <- loc_colors[which(loc_order %in% unique(uga_sel_4$NAME_4))]
new_col <- loc_colors[c(which(loc_order %in% col1),
                        which(loc_order %in% col2),
                        which(loc_order %in% col3),
                        which(loc_order %in% col4))]
# map <- ggplot(uga) +
#   geom_sf(color = "#00000000") +
#   geom_sf(data = uga_sel_1, aes(fill = factor(NAME_1, levels = loc_order)), color="#00000000") +
#   geom_sf(data = uga_sel_2, aes(fill = factor(NAME_2, levels = loc_order)), color="#00000000") +
#   geom_sf(data = uga_sel_3, aes(fill = factor(NAME_3, levels = loc_order)), color="#00000000") +
#   geom_sf(data = uga_sel_4, aes(fill = factor(NAME_4, levels = loc_order)), color="#00000000") +
#   theme_bw() + scale_fill_manual(values = new_col)+
#   # scale_color_manual(values = rep("#00000000", length(unique(uga$NAME_1)))) +
#   theme(axis.text = element_text(color = 'black', size = 12),
#         axis.title = element_text(color = 'black', size = 12, face = "bold"),
#         legend.title = element_blank())+ labs(x = "Longitude", y = "Latitude")
library(ggnewscale)
map <- ggplot(uga) +
  geom_sf(color = "#00000000") +
  geom_sf(data = uga_sel_1, aes(fill = factor(NAME_1, levels = loc_order)), color="#00000000") +
  scale_fill_manual(values = col1) + new_scale_fill() +
  geom_sf(data = uga_sel_2, aes(fill = factor(NAME_2, levels = loc_order)), color="#00000000") +
  scale_fill_manual(values = col2) + new_scale_fill() +
  geom_sf(data = uga_sel_3, aes(fill = factor(NAME_3, levels = loc_order)), color="#00000000") +
  scale_fill_manual(values = col3) + new_scale_fill() +
  geom_sf(data = uga_sel_4, aes(fill = factor(NAME_4, levels = loc_order)), color="#00000000") +
  theme_bw() + 
  scale_fill_manual(values = col4) +
  theme(axis.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 12, face = "bold"),
        legend.position = "none") + 
  labs(x = "Longitude", y = "Latitude")

ggsave("./output/Locations_map.jpg", width = 4, 
       height = 4, dpi = 600, units = "in", plot = map)

## Map2 below helps to confirm location color on map.
map2 <- ggplot(uga) +
  geom_sf(color = "#00000000") +
  geom_sf(data = uga_sel_1, aes(fill = factor(NAME_1, levels = loc_order)), color="#00000000") +
  scale_fill_manual(values = col1) + new_scale_fill() +
  geom_sf(data = uga_sel_2, aes(fill = factor(NAME_2, levels = loc_order)), color="#00000000") +
  scale_fill_manual(values = col2) + new_scale_fill() +
  geom_sf(data = uga_sel_3, aes(fill = factor(NAME_3, levels = loc_order)), color="#00000000") +
  scale_fill_manual(values = col3) + new_scale_fill() +
  geom_sf(data = uga_sel_4, aes(fill = factor(NAME_4, levels = loc_order)), color="#00000000") +
  theme_bw() + 
  scale_fill_manual(values = col4) +
  theme(axis.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 12, face = "bold")) + 
  labs(x = "Longitude", y = "Latitude")

ggsave("./output/Locations_map2.jpg", width = 8, 
       height = 8, dpi = 600, units = "in", plot = map2)

# cropping out some areas to Zoom in on the map.
xlow = 32.9 ; xhigh = 34.5
ylow = 0.2 ; yhigh = 2.1
uga_filtered <- uga %>%
  mutate(centroid = st_centroid(geometry)) %>%
  filter(st_coordinates(centroid)[,1] > xlow & 
           st_coordinates(centroid)[,1] < xhigh) %>% 
  filter(st_coordinates(centroid)[,2] > ylow & 
           st_coordinates(centroid)[,2] < yhigh)

saveRDS(uga_filtered, "./data/UG_map_cropped_data.RDS")
map2 <- ggplot(uga_filtered) +
  geom_sf(color = "#00000000") +
  geom_sf(data = uga_sel_1, aes(fill = factor(NAME_1, levels = loc_order)), color="#00000000") +
  scale_fill_manual(values = col1) + new_scale_fill() +
  geom_sf(data = uga_sel_2, aes(fill = factor(NAME_2, levels = loc_order)), color="#00000000") +
  scale_fill_manual(values = col2) + new_scale_fill() +
  geom_sf(data = uga_sel_3, aes(fill = factor(NAME_3, levels = loc_order)), color="#00000000") +
  scale_fill_manual(values = col3) + new_scale_fill() +
  geom_sf(data = uga_sel_4, aes(fill = factor(NAME_4, levels = loc_order)), color="#00000000") +
  scale_fill_manual(values = col4) + new_scale_fill() +
  theme_bw() + scale_fill_manual(values = loc_colors)+
  theme(axis.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 12, face = "bold"),
        legend.position = "none") + 
  scale_x_continuous(name = "Longitude", limits = c(xlow, xhigh),
                     breaks = seq(round(xlow,0), xhigh, 0.5),
                     labels = function(x) paste0(x, "°E"), expand = c(0, 0)) +
  scale_y_continuous(name = "Latitude", limits = c(ylow, yhigh),
                     labels = function(y) paste0(y, "°N"), expand = c(0, 0))

ggsave("./output/Locations_map_cropped2.jpg", width = 5, 
       height = 5, dpi = 600, units = "in", plot = map2)

####### Plotting the samples as dots on the map
uga_filtered <- readRDS("./data/UG_map_cropped_data.RDS")
pheno.dt_sub <- read.csv("./data/pheno_dt_final.csv")[,1:2]
coordinates <- read.csv("./data/Sample_coordinates.csv")
admix_k3 <- read.csv("./data/Admixture_results_k3_final.csv")
loc_coord <- dplyr::left_join(pheno.dt_sub, coordinates, by = c("IID"="Sample.name")) |> 
  dplyr::rename(Longitude = longitude, Latitude = Latitude)

locs_and_struct <- left_join(admix_k3, pheno.dt_sub, by = "IID") %>%
  na.omit() %>%
  group_by(Location, Cluster) %>%
  summarize(count = n(), .groups = "drop") |> 
  tidyr::pivot_wider(names_from = "Cluster", values_from = "count") |> 
  dplyr::mutate(total = rowSums(across(starts_with("Subpopulation")), na.rm = T),
                   `Subpopulation 1` = `Subpopulation 1`/total,
                   `Subpopulation 2` = `Subpopulation 2`/total,
                   `Subpopulation 3` = `Subpopulation 3`/total) |> 
  dplyr::mutate(across(starts_with("Subpopulation"), ~ replace_na(., 0))) |> 
  as.data.frame()

loc_coord_summarized <- loc_coord |> 
  filter(IID %in% admix_k3$IID) |> 
  group_by(Location) |> 
  summarise(n_samples = n(), Longitude = mean(Longitude),
            Latitude = mean(Latitude), .groups = "drop")
loc_admix_coord <- locs_and_struct |> 
  left_join(loc_coord_summarized) 
  # st_as_sf(coords = c("Longitude", "Latitude"), 
  #          crs = st_crs(4326), remove = F)
write.csv(loc_admix_coord, "./data/final_data_scatterpie.csv", row.names = F)
  
sf_coord <- st_as_sf(loc_coord[,c(2,4,5)], 
                     coords = c("Longitude", "Latitude"), crs = 4326)


map3 <- ggplot(uga_filtered) +
  geom_sf(color = "#00000000") +
  geom_sf(data = uga_sel_1, aes(fill = factor(NAME_1, levels = loc_order)), color="#00000000") +
  scale_fill_manual(values = col1) + new_scale_fill() +
  geom_sf(data = uga_sel_2, aes(fill = factor(NAME_2, levels = loc_order)), color="#00000000") +
  scale_fill_manual(values = col2) + new_scale_fill() +
  geom_sf(data = uga_sel_3, aes(fill = factor(NAME_3, levels = loc_order)), color="#00000000") +
  scale_fill_manual(values = col3) + new_scale_fill() +
  geom_sf(data = uga_sel_4, aes(fill = factor(NAME_4, levels = loc_order)), color="#00000000") +
  scale_fill_manual(values = col4) + new_scale_fill() +
  theme_bw() + scale_fill_manual(values = loc_colors)+
  # geom_point(data = loc_coord_summarized, aes(x = Longitude, y = Latitude, size = n_samples),
  #            color = 'black', alpha = 0.5) +
  geom_point(data = loc_coord[,4:5], aes(x = Longitude, y = Latitude),
             color = 'black', size = 1.5, alpha = 0.5) +
  theme(axis.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 12, face = "bold"),
        legend.position = "none")+  
  scale_x_continuous(name = "Longitude", limits = c(xlow, xhigh),
                     breaks = seq(round(xlow,0), xhigh, 0.5),
                     labels = function(x) paste0(x, "°E"), expand = c(0, 0)) +
  scale_y_continuous(name = "Latitude", limits = c(ylow, yhigh),
                     labels = function(y) paste0(y, "°N"), expand = c(0, 0))

ggsave("./output/Locations_map_cropped_samples.jpg", width = 4.5, 
       height = 4.5, dpi = 600, units = "in", plot = map3)

library(scatterpie) ### Needed to plot pies
loc_admix_coord <- read.csv("./data/final_data_scatterpie.csv")
map4 <- ggplot() +
  geom_sf(data = uga_filtered, color = "#00000000", fill = "grey90") +
  geom_sf(data = uga_sel_1, aes(fill = factor(NAME_1, levels = loc_order)), color="#00000000", fill = "grey70") +
  geom_sf(data = uga_sel_2, aes(fill = factor(NAME_2, levels = loc_order)), color="#00000000", fill = "grey70") +
  geom_sf(data = uga_sel_3, aes(fill = factor(NAME_3, levels = loc_order)), color="#00000000", fill = "grey70") +
  geom_sf(data = uga_sel_4, aes(fill = factor(NAME_4, levels = loc_order)), color="#00000000", fill = "grey70") +
  theme_bw() + 
  ## Adding the pies
  geom_scatterpie(data = loc_admix_coord,
          aes(x = Longitude, y = Latitude, r = sqrt(total)/50, group = Location),
          cols = c("Subpopulation.1", "Subpopulation.2", "Subpopulation.3")) +
  scale_fill_manual(values = friendly_pal("bright_seven")[c(2,3,1)]) +
  theme_bw() +
  theme(axis.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 12, face = "bold"),
        legend.position = "none") +  
  scale_x_continuous(name = "Longitude", limits = c(xlow, xhigh),
                     breaks = seq(round(xlow,0), xhigh, 0.5),
                     labels = function(x) paste0(x, "°E"), expand = c(0, 0)) +
  scale_y_continuous(name = "Latitude", limits = c(ylow, yhigh),
                     labels = function(y) paste0(y, "°N"), expand = c(0, 0))

ggsave("./output/Locations_map_cropped_piechart.jpg", width = 4, 
       height = 4, dpi = 600, units = "in", plot = map4)

##################################################################
## Import and process SNP dosage data set ===============
dosage_dt <- read.table("./data/dosage_output.raw", header = TRUE)
colnames(dosage_dt) <- c(colnames(dosage_dt)[1:6], paste0("SNP", 1:(ncol(dosage_dt) - 6)))

fwrite(dosage_dt, "./data/dosage_new.csv", row.names = F)

dosage_dt <- fread("./data/dosage_new.csv")
indiv_missing <- apply(dosage_dt, 1, function(x){100*length(which(is.na(x)))/84814})
index <- which(indiv_missing > 25) ## Individuals with >= 40% missing data
length(index) ### Index is empty
dosage_dt2 <- dosage_dt#[-index,]

write.table(dosage_dt2[,-c(1,3:6)], "./data/dosage_dt_final.txt", 
            row.names = F, col.names = T, quote = F, sep = "\t")

### Create a matrix for analysis
geno_mat <- as.matrix(dosage_dt2[, -c(1:6)])
rownames(geno_mat) <- dosage_dt2$IID
100*length(which(is.na(geno_mat)))/length(which(!is.na(geno_mat))) # Percent of missing values
geno_mat[which(is.na(geno_mat))] <- 1

write.table(geno_mat, "./data/geno_mat_final.txt", 
            row.names = T, col.names = T, quote = F, sep = "\t")


###############################################################
## Principal coordinate analysis ============
pheno.dt <- read.csv("./data/pheno_dt_final.csv")
geno_mat <- as.matrix(read.table("./data/geno_mat_final.txt"))

## Converting character variables to factors
str(pheno.dt)
pheno.dt[,-c(1:2)] <- lapply(pheno.dt[,-c(1:2)], function(x) if (is.character(x)) factor(x) else x)

## Retain only the 170 accessions in snp data set
pheno.dt <- filter(pheno.dt, IID %in% rownames(geno_mat))

library(cluster)
d <- daisy(pheno.dt[,-c(1:2)], metric = "gower")

# Run PCoA
library(ape)
pcoa_res <- pcoa(as.matrix(d))

# k3_results <- read.csv("./data/Admixture_results_k2_final.csv")
k3_results <- read.csv("./data/Admixture_results_k3_final.csv")

pcoa_dt <- as.data.frame(pcoa_res$vectors) %>% 
  mutate(., IID = pheno.dt$IID) %>% 
  full_join(., k3_results) %>% na.omit() %>% 
  mutate(., Cluster = as.character(Cluster))

pcoa_plt <- ggplot(pcoa_dt, aes(x = Axis.1, y = Axis.2, color = Cluster)) + 
  geom_point(size = 2) +
  scale_color_manual(values = friendly_pal("bright_seven")[c(2,3,1)])+
  labs(x = "PCoA1", y = "PCoA2") + boris_theme_rect +
  geom_jitter(width = .05, height = .05) +
  theme(legend.position = "right", legend.margin = margin(0,0,0,0),
        legend.direction = "vertical", legend.text = element_text(size = 8))

ggsave("./output/Phenotype_PCoA_biplot.jpg", width = 5, height = 4, 
       dpi = 600, units = "in", plot = pcoa_plt)

############################################################
## Admixture results processing =============
final_IDs <- fread("./data/dosage_dt_final.txt")[["IID"]]
samples <- read.table("./data/ld_pruned/data_pruned_final_numeric.fam")[["V2"]]
results_k2 <- read.table("./data/admixture_results3/data_pruned_final_numeric.2.Q") %>% 
  mutate(., IID = samples) %>% filter(., `IID` %in% final_IDs) %>% 
  rename(., Subpopulation1 = V1, Subpopulation2 = V2) %>% 
  arrange(., desc(Subpopulation1))

results_k3_raw <- read.table("./data/admixture_results3/data_pruned_final_numeric.3.Q") %>% 
  mutate(., IID = samples) %>% filter(., `IID` %in% final_IDs) %>% 
  rename(., Subpopulation1 = V1, Subpopulation2 = V2, Subpopulation3 = V3) %>% 
  arrange(., desc(Subpopulation1), desc(Subpopulation2))

### Creating the barplot of the structure analysis
geno_lvl <- results_k2$IID
results_k2 <- results_k2 %>% 
  mutate(IID = factor(IID, levels = geno_lvl)) %>% 
  pivot_longer(., 1:2, names_to = "Cluster", values_to = "Clust_prob")

(vcp = ggplot(data = results_k2, aes(x = IID, y = Clust_prob, fill = Cluster)) + 
    geom_col(position = "fill") + 
    scale_fill_manual(values = friendly_pal("bright_seven")[c(2,3,1)])+
    labs(y = "Probability", x = "") + theme_bw() +
    scale_y_continuous(limits = c(0.00,1.00), breaks = seq(0.00, 1.00, 0.20)) +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank()))
ggsave("./output/admixture_K2.jpg", width = 6, height = 2, dpi = 600, units = "in", plot = vcp)

## Barplot for K=3
geno_lvl <- results_k3_raw$IID
results_k3 <- results_k3_raw %>%
  mutate(IID = factor(IID, levels = geno_lvl)) %>%
  pivot_longer(., 1:3, names_to = "Cluster", values_to = "Clust_prob")

#"#4477AA" "#228833" "#AA3377" "#BBBBBB" "#66CCEE" "#CCBB44" "#EE6677"

(vcp2 = ggplot(data = results_k3, aes(x = IID, y = Clust_prob, fill = Cluster)) +
    geom_col(position = "fill") +
    scale_fill_manual(values = friendly_pal("bright_seven")[c(2:3,1)])+
    labs(y = "Probability", x = "") + theme_bw() +
    scale_y_continuous(limits = c(0.00,1.00), breaks = seq(0.00, 1.00, 0.20)) +
    theme(axis.text.x = element_blank(), axis.ticks = element_blank()))
ggsave("./output/admixture_K3.jpg", width = 6, height = 2, dpi = 600, units = "in", plot = vcp2)

### Geographical Location barplot with admixture results, K=3
clust1 <- filter(results_k3_raw, Subpopulation1 > .5)
clust2 <- filter(results_k3_raw, Subpopulation2 > .5)
clust3 <- filter(results_k3_raw, Subpopulation3 > .5)
clust_mixed <- filter(results_k3_raw, !(IID %in% c(clust1$IID, clust2$IID, clust3$IID)))
struct_results_k3 <- rbind(clust1, clust2, clust3)

admix_k3 <- rbind(mutate(clust1, Cluster = "Subpopulation 1"),
                   mutate(clust2, Cluster = "Subpopulation 2"),
                   mutate(clust3, Cluster = "Subpopulation 3"),
                  mutate(clust_mixed, Cluster = c("Subpopulation 2", "Subpopulation 3", "Subpopulation 2"))) |> 
  select(IID, Cluster)
write.csv(admix_k3, "./data/Admixture_results_k3_final.csv", row.names = F)
admix_k3 <- read.csv("./data/Admixture_results_k3_final.csv")
locs_and_struct <- left_join(admix_k3, pheno.dt[,1:2], by = "IID") %>%
  na.omit() %>%
  group_by(Location, Cluster) %>%
  summarize(count = n(), .groups = "drop")

ploc_admix <- ggplot(locs_and_struct, aes(x = factor(Location, levels = loc_order), y = count, fill = Cluster)) + 
  geom_bar(stat = "identity", position = "stack", color = "black") + 
  scale_fill_manual(values = friendly_pal("bright_seven")[c(2:3,1)])+
  labs(y = "Number of samples", x = "Location") + boris_theme_rect +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "right", legend.margin = margin(0,0,0,0),
        legend.direction = "vertical")
ggsave("./output/Location_admixture_K3.jpg", width = 6, 
       height = 4, dpi = 600, units = "in", plot = ploc_admix)

xtable <- locs_and_struct %>% 
  pivot_wider(names_from = Cluster, values_from = count, values_fill = 0) %>%
  mutate(Total = `Subpopulation 1` + `Subpopulation 2` + `Subpopulation 3`) %>%
  arrange(desc(Total))
write.csv(xtable, "./output/Location_admixture_K3_table.csv", row.names = F)

### ADMIXTURE results, optimum K value
cv_err <- read.table("./data/admixture_results3/cv_error.txt")
cv_error_dt <- data.frame(K_value = c(1,10, 2:9), CV_error = cv_err$V4) |> 
  arrange(K_value) |> as.data.frame()
cv_plot <- ggplot(cv_error_dt, aes(x = K_value, y = CV_error)) + geom_point(size = 2) + 
  geom_line(linewidth = 1) + labs(x = "K", y = "Cross-validation error") + boris_theme_rect +
  scale_x_continuous(limits = c(1,10), breaks = seq(1,10,1))

ggsave("./output/admixture_cv_error_plot.jpg",width = 3, height = 3, 
       dpi = 600, units = "in", plot = cv_plot)

# ## Investigating the proportion of admixed indiv for k =2
# struc_clust_admixed <- results_k2 %>% 
#   mutate(., Cluster = ifelse(Clust_prob > 0.5 & Clust_prob <= 0.52,"Admixed", Cluster))
# table(struc_clust_admixed$Cluster)
# write.csv(struc_clust_admixed, "./output/structure_cluster_and_admixed.csv", row.names = F)
# 
# (vcp = ggplot(data = results_k2, aes(x = IID, y = Clust_prob, fill = Cluster)) + 
#     geom_col(position = "fill") + 
#     scale_fill_manual(values = friendly_pal("bright_seven")[c(2,3,1)])+
#     labs(y = "Probability", x = "") + theme_bw() +
#     scale_y_continuous(limits = c(0.00,1.00), breaks = seq(0.00, 1.00, 0.20)) +
#     theme(axis.text.x = element_blank(), axis.ticks = element_blank()))
# ggsave("./output/structure_K2.jpg", width = 6, height = 2, dpi = 600, units = "in", plot = vcp)

#################################################
## Principal component analysis ==========
pc <- prcomp(geno_mat)

plot((pc$sdev^2)[1:50], type = "l", xlab = "Number of axes", ylab = "Component variance")
plot((100*(pc$sdev^2)/sum(pc$sdev^2))[1:50], type = "l", 
     xlab = "Number of axes", ylab = "PC variance contribution (%)")
summary(pc)
## Saving the scores of first 3 PCs
scor <- as.data.frame(pc$x[,(1:3)])
scor$Label <- rownames(geno_mat)
write.csv(scor, "./data/PCA_3_scores.csv", row.names = F)

### Plotting pca scores
library(tidyverse)
scor <- read.csv("./data/PCA_3_scores.csv")

## PCA Colored by pop structure results ========
scor <- read.csv("./data/PCA_3_scores.csv")
pc1_lab <- "PC1 (6%)"; pc2_lab = "PC2 (2%)"; pc3_lab = "PC3 (1%)"

k3_results <- read.csv("./data/Admixture_results_k3_final.csv")
scor_k3 <-left_join(k3_results, scor, by = c("IID"= "Label")) %>% 
  mutate(., Cluster = as.character(Cluster))

pc12 = ggplot(scor_k3, aes(x = PC1, y = PC2, color = Cluster)) + geom_point(size =3) + 
  labs(x = pc1_lab, y = pc2_lab) + boris_theme_rect +
  scale_color_manual(values = friendly_pal("bright_seven")[c(2,3,1)]) +
  theme(legend.position = "right", legend.margin = margin(0,0,0,0),
        legend.direction = "vertical", legend.text = element_text(size = 8))
pc13 = ggplot(scor_k3, aes(x = PC1, y = PC3, color = Cluster)) + geom_point(size =3) + 
  labs(x = pc1_lab, y = pc3_lab) + boris_theme_rect +
  scale_color_manual(values = friendly_pal("bright_seven")[c(2,3,1)]) +
  theme(legend.position = "right", legend.margin = margin(0,0,0,0),
        legend.direction = "vertical", legend.text = element_text(size = 8))
pc23 = ggplot(scor_k3, aes(x = PC2, y = PC3, color = Cluster)) + geom_point(size =3) + 
  labs(x = pc2_lab, y = pc3_lab) + boris_theme_rect +
  scale_color_manual(values = friendly_pal("bright_seven")[c(2,3,1)]) +
  theme(legend.position = "right", legend.margin = margin(0,0,0,0),
        legend.direction = "vertical", legend.text = element_text(size = 8))

ggsave("./output/Structure_pc12_biplot.jpg", width = 5, height = 4, 
       dpi = 600, units = "in", plot = pc12)
ggsave("./output/Structure_pc13_biplot.jpg", width = 5, height = 4, 
       dpi = 600, units = "in", plot = pc13)
ggsave("./output/Structure_pc23_biplot.jpg", width = 5, height = 4, 
       dpi = 600, units = "in", plot = pc23)

## Color by location
pheno.dt <- read.csv("./data/pheno_dt_final.csv")
locs <- pheno.dt[,1:2] %>% as.data.frame()
pc1_lab <- "PC1 (6%)"; pc2_lab = "PC2 (2%)"; pc3_lab = "PC3 (1%)"
scor <- read.csv("./data/PCA_3_scores.csv")
scor <-left_join(locs, scor, by = c("IID" ="Label")) %>% 
  mutate(Location = factor(Location, levels = loc_order))

pc12 = ggplot(scor, aes(x = PC1, y = PC2, color = Location)) + geom_point(size =2) + 
  labs(x = pc1_lab, y = pc2_lab) + boris_theme_rect +
  scale_color_manual(values = loc_colors) +
  theme(legend.position = "right", legend.margin = margin(0,0,0,0),
        legend.direction = "vertical", legend.text = element_text(size = 8))
pc13 = ggplot(scor, aes(x = PC1, y = PC3, color = Location)) + geom_point(size =2) + 
  labs(x = pc1_lab, y = pc3_lab) + boris_theme_rect +
  scale_color_manual(values = loc_colors) +
  theme(legend.position = "right", legend.margin = margin(0,0,0,0),
        legend.direction = "vertical", legend.text = element_text(size = 8))
pc23 = ggplot(scor, aes(x = PC2, y = PC3, color = Location)) + geom_point(size =2) + 
  labs(x = pc2_lab, y = pc3_lab) + boris_theme_rect +
  scale_color_manual(values = loc_colors) +
  theme(legend.position = "right", legend.margin = margin(0,0,0,0),
        legend.direction = "vertical", legend.text = element_text(size = 8))

ggsave("./output/Location_pc12_biplot.jpg", width = 5, height = 4, 
       dpi = 600, units = "in", plot = pc12)
ggsave("./output/Location_pc13_biplot.jpg", width = 5, height = 4, 
       dpi = 600, units = "in", plot = pc13)
ggsave("./output/Location_pc23_biplot.jpg", width = 5, height = 4, 
       dpi = 600, units = "in", plot = pc23)

### 3D plotting of PCA results
library(scatterplot3d)
locations <- unique(locs$Location)
col_map <- setNames(loc_colors, loc_order)

jpeg("./output/3D_PCA_location_plot.jpg", width = 10, height = 6, res = 600, units = "in")
scatterplot3d(x = scor$PC1, y = scor$PC3, z = scor$PC2,
              color = col_map[scor$Location], pch = 16, cex.symbols = 1,
              xlab = pc1_lab, ylab = pc3_lab, zlab = pc2_lab, angle = 55)

legend("topright", inset = c(-0.01, 0), legend = as.character(loc_order), 
       col = loc_colors, pch = 19, cex = 0.5)
dev.off()



#################################################
## Cluster analysis using neighbor-joining tree ===========
library(adegenet)
library(data.table)

clus_mat = geno_mat
gen_ind_object <- df2genind(X = clus_mat, ploidy = 2, ncode = 1)
gen_pop_obj <- as.genpop(gen_ind_object$tab)
#cschord_dist <- dist.genpop(gen_ind_object , method = "CSChord")
cschord_dist <- dist.genpop(gen_pop_obj , method = 2)

k3_results <- read.csv("./data/Admixture_results_k3_final.csv")

library(ape)
nj_clust <- nj(cschord_dist)
rm(list= list(gen_ind_object, gen_pop_obj, dat, clus_mat))

labels_tree2 <- nj_clust[["tip.label"]]

jpeg("./output/cluster_dendrogram.jpg", 
     width = 8, height = 8, res = 600, units = "in")
plot(nj_clust, cex = 0.6)
dev.off()

jpeg("./output/cluster_dendrogram_fan_shape.jpg", 
     width = 8, height = 8, res = 600, units = "in")
plot(nj_clust, cex = 0.6, type = "fan")
dev.off()

## Attribute colors to subpopulation
subpop1 <- k3_results$IID[which(k3_results$Cluster == "Subpopulation 1")]
subpop2 <- k3_results$IID[which(k3_results$Cluster == "Subpopulation 2")]
subpop3 <- k3_results$IID[which(k3_results$Cluster == "Subpopulation 3")]

labels.color <- rep("x",length(labels_tree2))

labels.color[which(labels_tree2 %in% subpop1)] <- friendly_pal("bright_seven")[2]
labels.color[which(labels_tree2 %in% subpop2)] <- friendly_pal("bright_seven")[3]
labels.color[which(labels_tree2 %in% subpop3)] <- friendly_pal("bright_seven")[1]
labels.color[which(labels.color == "x")] <- "white"

edge_colors <- rep("x", nrow(nj_clust$edge))
edges <- which(nj_clust$edge[,2] <= length(nj_clust$tip.label))
edge_colors[edges] <- labels.color[nj_clust$edge[edges,2]]

jpeg("./output/cluster_fan_Subpopulation_color_code.jpg", width = 8, height = 8, res = 600, units = "in")
plot(nj_clust, type = "fan", cex = .7, tip.color = labels.color)
dev.off()

jpeg("./output/cluster_unrooted_Subpopulation_color_code.jpg", width = 8, height = 8, res = 600, units = "in")
plot(nj_clust, type = "unrooted", cex = .5, tip.color = labels.color)
dev.off()

##### DAPC analysis #######
## Reference tutorial for adegenet DAPC functions
#chrome-extension://efaidnbmnnnibpcajpcglclefindmkaj/https://adegenet.r-forge.r-project.org/files/tutorial-dapc.pdf
geno_mat <- as.matrix(read.table("./data/geno_mat_final.txt"))

library(adegenet)
class(geno_mat)
clusters <- find.clusters(geno_mat, max.n.clust = 40, n.pca = 50)
dapc_K3 = dapc(x = geno_mat, grp = clusters$grp, n.pca = 50)
dapc_loadings <- dapc_K3$ind.coord |> as.data.frame() |>  
  rownames_to_column(var = "IID")

dapc_orig <- adegenet::scatter.dapc(dapc_K3, scree.da = TRUE, posi.da = "bottomright",
                                    col = c("darkred", "darkblue", "darkgreen", "purple"))
jpeg("./output/DAPC_plot_default.jpg", width = 4, height = 4, res = 300, units = "in")
adegenet::scatter.dapc(dapc_K3, scree.da = TRUE, posi.da = "bottomright",
                                    col = c("darkred", "darkblue", "darkgreen", "purple"))
dev.off()

k3_results <- read.csv("./data/Admixture_results_k3_final.csv")
admix_labels <- k3_results$Cluster
names(admix_labels) <- k3_results$IID
admix_labels <- factor(admix_labels[match(names(dapc_K3$grp), names(admix_labels))])
cols <- friendly_pal("bright_seven")[c(2,3,1)]

jpeg("./output/DAPC_plot_admixture_res.jpg", width = 4.5, height = 4.5, res = 300, units = "in")
adegenet::scatter.dapc(dapc_K3, grp = admix_labels,
                       scree.da = TRUE, posi.da = "bottomright",
                       col = cols)
dev.off()

############################################################
##### FST computation ====================================
library(adegenet)
library(hierfstat)
library(tidyverse)

k3_results <- read.csv("./data/Admixture_results_k3_final.csv")
geno_mat <- read.table("./data/geno_mat_final.txt")

ind <- row.names(geno_mat)
pop <- left_join(data.frame(IID = ind), k3_results) %>% pull(Cluster)
ssdata <- data.frame(population = pop, taxa = ind)
gd_dt <- df2genind(X = geno_mat, ind.names = ind, pop = pop, 
                   ploidy = 2, type = "codom", ncode = 1)
saveRDS(gd_dt, "./data/genind_object_final.rds")

fst <- wc(gd_dt)
saveRDS(fst, "./data/wc_fst_results.rds")
global_fst <- fst$FST ## Fst = 0.03
loc_fst <- fst$per.loc
write.csv(loc_fst, "./data/per_locus_fst_values.csv", row.names = F)

fst_plt <- ggplot(loc_fst, aes(x = FST)) + 
  geom_histogram(fill = "grey40", color = "white") +
  geom_histogram(data = filter(loc_fst, FST >= 0.2), 
                 aes(x = FST), fill = "darkred", color = "white") +
  labs(y = "Frequency", x = "Fst value per SNP") + boris_theme_rect +
  scale_x_continuous(limits = c(0,0.5), breaks = seq(0,0.5,0.1)) +
  scale_y_continuous(limits = c(0,4200), breaks = seq(0,4000,2000))
ggsave("./output/Fst_value_per_locus.jpg", width = 4, height = 4, 
       dpi = 600, units = "in", plot = fst_plt)

##### Compute the pairwise FSt values of the population
gd_dt <- readRDS("./data/genind_object_final.rds")
paire_fst <- pairwise.WCfst(gd_dt)

### FST Manhathan plot
fst <- readRDS("./data/wc_fst_results.rds")
fst_dt <- fst$per.loc |> as.data.frame() |> 
  rownames_to_column(var = "snp") %>%
  mutate(., snp_index = as.numeric(gsub("SNP", "", snp)))


fst_manhattan <- ggplot(data = fst_dt, aes(x= snp_index, y = FST)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "Fst value") + boris_theme_rect +
  scale_y_continuous(limits = c(0,0.5), breaks = seq(0,0.5,0.1)) +
  scale_x_continuous(limits = c(0,23000), breaks = seq(0,22000,5000)) +
  geom_hline(yintercept = 0.2, color = "darkred", linetype = "dashed", linewidth = 1)

ggsave("./output/Fst_mahattan_plot.jpg", width = 8, height = 4, 
       dpi = 600, units = "in", plot = fst_manhattan)

####################################################################
### Computing expected He
He_function <- function(mkr_dt, ID_col_name,
                        homo_code1 = "0", homo_code2 = "2", hetero_code = "1"){
  # name_ind <- which(ID_col_name %in% colnames(mkr_dt))
  mkr_names <- setdiff(colnames(mkr_dt), ID_col_name)
  #Obtain genotypic counts
  snp_dt <- mkr_dt |>
    tidyr::pivot_longer(cols = all_of(mkr_names),
                        names_to = "snp", values_to = "dosage") |>
    dplyr::group_by(snp, dosage) |>
    dplyr::summarise(counts = dplyr::n(),
                     .groups = "drop") |>
    tidyr::pivot_wider(names_from = dosage, values_from = counts)
  ##Compute allele frequencies and He
  if(is.null(hetero_code)){
    freq_he <- snp_dt |>
      dplyr::mutate(across(where(is.numeric), ~tidyr::replace_na(., 0))) |>
      dplyr::mutate(total = rowSums(dplyr::across(where(is.numeric)), na.rm = T),
                    p = .data[[homo_code1]]/total,
                    q = 1-p,
                    He = 1 - (p^2 + q^2))
  }else{
    freq_he <- snp_dt |>
      dplyr::mutate(across(where(is.numeric), ~tidyr::replace_na(., 0))) |>
      dplyr::mutate(total = rowSums(dplyr::across(where(is.numeric)), na.rm = T),
                    p = (2*(.data[[homo_code1]]) + .data[[hetero_code]])/(2*total),
                    q = 1-p,
                    He = 1 - (p^2 + q^2))
  }
  return(freq_he)
}

He_values <- He_function(mkr_dt = geno_dt[,-22449], ID_col_name = "IID") |>
  dplyr::mutate(Cluster = "Whole population")
He_pop1 <- He_function(mkr_dt = filter(geno_dt, Cluster == "Subpopulation 1")[,-22449],
                       ID_col_name = "IID") |> dplyr::mutate(Cluster = "Subpopulation 1")
He_pop2 <- He_function(mkr_dt = filter(geno_dt, Cluster == "Subpopulation 2")[,-22449],
                       ID_col_name = "IID") |> dplyr::mutate(Cluster = "Subpopulation 1")
He_pop3 <- He_function(mkr_dt = filter(geno_dt, Cluster == "Subpopulation 3")[,-22449],
                       ID_col_name = "IID") |> dplyr::mutate(Cluster = "Subpopulation 1")
He_all <- rbind(He_values, He_pop1, He_pop2, He_pop3)

gd_dt <- readRDS("./data/genind_object_final.rds")
stats <- hierfstat::basic.stats(gd_dt, diploid = TRUE)
saveRDS(stats, "./data/Diversity_basic_stats_final.rds")
diversity <- stats$perloc |> rownames_to_column(var = "snp")
summary(diversity)
### Hs is not right for He
# pop_He <- stats$Hs |> as.data.frame() |> rownames_to_column(var = "snp")
pop_Ho <- stats$Ho |> as.data.frame() |> rownames_to_column(var = "snp")

#### For some reason, some code went missing
### Saving their outputs for backup
He_list <- list(He_all = He_all,
                He_Ho = He_Ho,
                He_Fst = He_Fst)
saveRDS(He_list, "./data/He_list.RDS")
### Heterozygosity ---------------------
# He <- full_join(pop_He, select(diversity,snp,Hs)) |> 
#   pivot_longer(cols = 2:5, names_to = "Population", values_to = "He") |> 
#   mutate(Population = recode(Population, "Hs" = "Whole population"))

Ho <- full_join(pop_Ho, select(diversity,snp,Ho)) |> 
  pivot_longer(cols = 2:5, names_to = "Population", values_to = "Ho") |> 
  mutate(Population = recode(Population, "Ho" = "Whole population"))
He_Ho <- full_join(He, Ho, by = c("snp", "Population"))

p_He <- ggplot(He_all, aes(x = He, fill = Cluster)) +
  geom_histogram(position = "dodge", bins = 10, color = "white") +
  labs(x = "Expected heterozygosity", y = "Frequency") +
  scale_fill_manual(values = friendly_pal("bright_seven")[c(2,3,1,4)]) +
  scale_x_continuous(breaks = seq(0,0.5,0.1)) +
  boris_theme_rect + theme(legend.position = "none",
                           legend.direction = "vertical")
ggsave("./output/Expect_hetero_distribution.jpg", width = 4, height = 3, 
       dpi = 600, units = "in", plot = p_He)

pb_He <- ggplot(He_all, aes(y = He,  x = Cluster, fill = Cluster)) +
  geom_boxplot()+
  labs(y = "Expected heterozygosity", x = "Population") +
  scale_fill_manual(values = friendly_pal("bright_seven")[c(2,3,1,4)]) + 
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  boris_theme_rect + 
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("./output/Expect_hetero_distribution_boxplot.jpg", width = 4, height = 4, 
       dpi = 600, units = "in", plot = pb_He)
p_Ho <- ggplot(He_Ho, aes(x = Ho, fill = Population)) + 
  geom_histogram(position = "dodge", bins = 10, color = "white") +
  labs(x = "Observed heterozygosity", y = "Frequency") + 
  scale_fill_manual(values = friendly_pal("bright_seven")[c(2,3,1,4)]) + 
  # scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  boris_theme_rect + theme(legend.position = "none")
ggsave("./output/Observed_hetero_distribution.jpg", width = 4, height = 3, 
       dpi = 600, units = "in", plot = p_Ho)

pb_Ho <- ggplot(He_Ho, aes(y = Ho, x = Cluster, fill = Cluster)) + 
  geom_boxplot()+
  labs(y = "Observed heterozygosity", x = "Population") +
  scale_fill_manual(values = friendly_pal("bright_seven")[c(2,3,1,4)]) + 
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  boris_theme_rect + 
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("./output/Observed_hetero_distribution_boxplot.jpg", width = 4, height = 4, 
       dpi = 600, units = "in", plot = pb_Ho)

p_Ho_He <- ggplot(He_Ho, aes(x = Ho, y = He, color = Population)) + 
  geom_point(size = 1) + 
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "black") +
  labs(x = "Observed heterozygosity", y = "Expected heterozygosity") + 
  scale_color_manual(values = friendly_pal("bright_seven")[c(2,3,1,4)]) + 
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.25))+
  boris_theme_rect + theme(legend.position = "none")+
  facet_wrap(~Population, nrow =2, ncol =2)

ggsave("./output/Observed_Expect_hetero_scatterplot.jpg", width = 5, height = 5, 
       dpi = 600, units = "in", plot = p_Ho_He)

### MAF --------------------------------
maf_function <- function(vec){
  tb <- table(vec)
  nallele <- length(tb)
  
  if(nallele > 2){
    af <- (2*tb['0']+ tb['1'])/(2*sum(tb))
  }else if(nallele == 2 & ('1' %in% names(tb))){
    af <- (2*tb[setdiff(names(tb), '1')] + tb['1'])/(2*sum(tb))
    
  }else if(nallele == 2 & !('1' %in% names(tb))){
    af <- (tb['0'])/(sum(tb))
  }else if(nallele == 1 & (names(tb) == '1')){
    af <- 0.5
  }else if(nallele == 1 & (names(tb) != '1')){
    af <- 0
  }
  ## Make sure MAF is <= 0.5
  if(af > 0.5){
    af <- (1 - af)
  }
  return(af)
}

maf_gmat <- apply(geno_mat, 2, maf_function)

geno_pop_dt <- data.frame(pop = pop, geno_mat)
pop1_geno_mat <- filter(geno_pop_dt, pop == "Subpopulation 1") %>% select(-pop) %>% as.matrix()
pop2_geno_mat <- filter(geno_pop_dt, pop == "Subpopulation 2") %>% select(-pop) %>% as.matrix()
pop3_geno_mat <- filter(geno_pop_dt, pop == "Subpopulation 3") %>% select(-pop) %>% as.matrix()
maf_pop1 <- apply(pop1_geno_mat, 2, maf_function)
maf_pop2 <- apply(pop2_geno_mat, 2, maf_function)
maf_pop3 <- apply(pop3_geno_mat, 2, maf_function)

maf_summary <- data.frame(snp = names(maf_gmat),
                     maf_whole = maf_gmat,
                     maf_subpop1 = maf_pop1,
                     maf_subpop2 = maf_pop2,
                     maf_subpop3 = maf_pop3)
write.csv(maf_summary, "./data/MAF_per_locus_and_population.csv", row.names = F)
maf_summary <- read.csv("./data/MAF_per_locus_and_population.csv")
maf_summary2 <- maf_summary |> 
  pivot_longer(cols = starts_with("maf"), 
               names_to = "Population", values_to = "MAF") |> 
  as.data.frame() |> 
  dplyr::mutate(Population = recode(Population, 
         "maf_whole" = "Whole population",
         "maf_subpop1" = "Subpopulation 1",
         "maf_subpop2" = "Subpopulation 2",
         "maf_subpop3" = "Subpopulation 3"))

p_maf <- ggplot(maf_summary2, aes(x = MAF, fill = Population)) + 
  geom_histogram(position = "dodge", bins = 10, color = "white") +
  labs(x = "Minor allele frequency (MAF)", y = "Frequency") + 
  scale_fill_manual(values = friendly_pal("bright_seven")[c(2,3,1,4)]) + 
  scale_x_continuous(limits = c(0,0.5), breaks = seq(0,0.5,0.1)) +
  boris_theme_rect + theme(legend.position = "right", 
                           legend.direction = "vertical")
ggsave("./output/MAF_distribution_histo.jpg", width = 6, height = 3, 
       dpi = 600, units = "in", plot = p_maf)

pb_maf <- ggplot(maf_summary2, aes(y = MAF, x = Population, fill = Population)) + 
  geom_boxplot() +
  labs(y = "Minor allele frequency (MAF)", y = "Population") + 
  scale_fill_manual(values = friendly_pal("bright_seven")[c(2,3,1,4)]) + 
  scale_y_continuous(limits = c(0,0.5), breaks = seq(0,0.5,0.1)) +
  boris_theme_rect + 
  theme(legend.position = "none",axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("./output/MAF_distribution_boxplot.jpg", width = 4, height = 4, 
       dpi = 600, units = "in", plot = pb_maf)


### FST & He -------------------------------
loc_fst <- read.csv("./data/per_locus_fst_values.csv")
He_Fst <- full_join(He_values, loc_fst)

loc_fst2 <- fst$per.loc |> as.data.frame() |> 
  rownames_to_column(var = "snp") |> 
  full_join(select(diversity, snp, Hs,Ho, Fst, Fstp), by = "snp")
cor(loc_fst2$Fstp, loc_fst2$FST)

p_He_Fst <- ggplot(He_Fst, aes(x = He, y = FST)) +
  geom_point(color = "grey40") +
  geom_point(data = filter(He_Fst, FST > 0.2), aes(x = He, y = FST), color = "darkred") +
  labs(x = "Expected heterozygosity", y = "Fst value") +
  scale_x_continuous(breaks = seq(0,0.5,0.1)) +
  boris_theme_rect

ggsave("./output/Expect_hetero_vs_Fst_scatterplot.jpg", width = 4, height = 3, 
       dpi = 600, units = "in", plot = p_He_Fst)

p_Ho_Fst <- ggplot(diversity, aes(x = Ho, y = Fstp)) + 
  geom_point(color = "grey40") +
  geom_point(data = filter(diversity, Fstp > 0.2), aes(x = Ho, y = Fstp), color = "darkred") +
  labs(x = "Observed heterozygosity", y="Fst value") +
  scale_x_continuous(limits = c(0,1), breaks = seq(0,1,0.25)) +
  boris_theme_rect

ggsave("./output/Observed_hetero_vs_Fst_scatterplot.jpg", width = 4, height = 3, 
       dpi = 600, units = "in", plot = p_Ho_Fst)


### Perform Amova
library(pegas)
dt_loci <- pegas::genind2loci(gd_dt)
strata_df <- data.frame(Pop = pop(dt_loci))
geno_mat <- as.matrix(geno_mat)
dist_mat <- dist(dt_loci)
amova_res <- pegas::amova(geno_mat ~ population, data = ssdata, nperm = 999)
amova_res2 <- pegas::amova(dist_mat ~ population, data = ssdata, nperm = 999)

saveRDS(amova_res2, "./data/Amova_results_final.rds")
pairewise_fst <- pairwise.WCfst(gd_dt)

### .. ################
bim <- read.table("./data/ld_pruned/data_pruned_final_numeric.bim")
tr_dosage <- transpoz(dosage_dt, new_col1_name = "snp")
#Investigate #snps per scaffold/contig
length(unique(bim$V1)) ## 4093 unique contigs
suma <- table(bim$V1) |> as.data.frame()
suma_500 <- suma |> filter(Freq >= 500) ## 44 contigs, 21 contigs have >= 500, 1000 snps, respectively.

##### Morphological data analysis
traits <- c(rep("Flower color", 3), rep("Leaf color", 9),
            rep("Leaf pigmentation", 2), rep("Leaf size", 2))
values <- c("Dark purple","Purple","Pale purple",
            "Dark green", "Green",  "Pale green", 
            "Yellowish green", "Green purple", "Yellowish green purple", 
            "Purple", "Dark purple", "Pink","Pigmented", "Non-pigmented","Narrow", "Broad")
values_f <- c("Pink", "Dark purple","Purple","Pale purple", 
              "Green purple", "Yellowish green purple","Dark green", 
              "Green",  "Pale green", "Yellowish green",
              "Pigmented", "Non-pigmented", "Narrow", "Broad")
proportion <- c(84.5, 0.5, 15.0,2.1, 28.0, 5.2, 18.7,
                26.9, 5.2, 11.9, 0.5, 1.6,59.6, 40.4, 47.2, 52.8)
summary_dt <- data.frame(Trait = traits, Value = values, Proportion = proportion) |>
  mutate(Value = factor(Value, levels = values_f))
View(summary_dt)

trait_colors <- c(
  "Dark purple"              = "#4B0082",
  "Purple"                   = "#9933FF",
  "Pale purple"              = "#D8BFD8",
  "Dark green"               = "#005C00",
  "Green"                    = "#228B22",
  "Pale green"               = "#98FB98",
  "Yellowish green"          = "#9ACD32",
  "Green purple"             = "#6A5ACD",
  "Yellowish green purple"   = "#C3B091",
  "Purple"                   = "#9933FF",
  "Dark purple"              = "#4B0082",
  "Pink"                     = "#FF69B4",
  "Pigmented"                = "#005C00",
  "Non-pigmented"            = "white",
  "Narrow"                   = "grey60",
  "Broad"                    = "grey20"
)

plt_Fcl <- ggplot(data = filter(summary_dt, Trait == "Flower color"),
                  aes(x = Value, y = Proportion, fill = Value)) +
geom_bar(stat = "identity", color = "black") +
  labs(x = "Flower color", y = "Proportion of samples (%)") + 
  boris_theme_rect + scale_fill_manual(values = trait_colors) +
  geom_text(aes(label = paste0( Value, "\n", Proportion, "%")), 
            position = position_stack(vjust = 0.6), size = 4) +
  theme(legend.position = "none", legend.margin = margin(0,0,0,0),
        legend.direction = "vertical", legend.text = element_text(size = 12))

export::graph2ppt(plt_Fcl, "./output/Leaf_color_proportions.pptx", width = 4, height = 4, append = T)

plt_Lpm <- ggplot(data = filter(summary_dt, Trait == "Leaf pigmentation"),
                  aes(x = Value, y = Proportion, fill = Value)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "Leaf pigmentation", y = "Proportion of samples (%)") + 
  boris_theme_rect + scale_fill_manual(values = trait_colors) +
  geom_text(aes(label = paste0( Value, "\n", Proportion, "%")), 
            position = position_stack(vjust = 0.6), size = 4) +
  theme(legend.position = "none", legend.margin = margin(0,0,0,0),
        legend.direction = "vertical", legend.text = element_text(size = 12))

export::graph2ppt(plt_Lpm, "./output/Leaf_color_proportions.pptx", width = 6, height = 5, append = T)
plt_Lsz <- ggplot(data = filter(summary_dt, Trait == "Leaf size"),
                  aes(x = Value, y = Proportion, fill = Value)) +
  geom_bar(stat = "identity", color = "black") +
  labs(x = "Leaf size", y = "Proportion of samples (%)") + 
  boris_theme_rect + scale_fill_manual(values = trait_colors) +
  geom_text(aes(label = paste0( Value, "\n", Proportion, "%")), 
            position = position_stack(vjust = 0.6), size = 4) +
  theme(legend.position = "none", legend.margin = margin(0,0,0,0),
        legend.direction = "vertical", legend.text = element_text(size = 12))

export::graph2ppt(plt_Lsz, "./output/Leaf_color_proportions.pptx", width = 6, height = 5, append = T)


######## GWAS analysis ============== ####################

bim <- read.table("./data/ld_pruned/data_pruned_final_numeric.bim")
colnames(bim) <- c("snp", "chr", "void", "pos", "allele1", "allele2")
head(bim)
bim2 <- bim |> 
  dplyr::mutate(chr = sub("\\..*", "", chr), snp = paste0("SNP",1:nrow(bim))) |> 
  dplyr::select(-void)
head(bim2)
bim_summary <- bim2$chr |> table() |> as.data.frame() |>
  dplyr::arrange(desc(Freq)) |> dplyr::slice_head(n = 40) #### Check if not missing those with high FST values
saveRDS(bim_summary, "./data/bim_summary.RDS")
bim_final <- dplyr::filter(bim2, chr %in% bim_summary$Var1) |> 
  as.data.frame()
saveRDS(bim_final, "./data/genetic_map_final.RDS")

dosage_dt <- read.table("./data/dosage_output.raw", header = TRUE)
snps <- sub("\\..*", "", colnames(dosage_dt[,-c(1:6)]))

gwas_snps <- which(snps %in% bim_summary$Var1)
gwas_GD <- geno_mat[,gwas_snps] |> as.data.frame() |> 
  tibble::rownames_to_column(var = "ID")
saveRDS(gwas_GD, "./data/gwas_genomic_data.RDS")

### GWAS analysis with milorgwas
library(milorGWAS)
help(package = "milorGWAS")
gwas_GD <- readRDS("./data/gwas_genomic_data.RDS")
t_gwa_gd <- transpoz(gwas_GD, new_col1_name = "snp")
gm <- readRDS("./data/genetic_map_final.RDS")
colnames(gm) <- c("id", "chr", "pos", "A1", "A2")
gwas_dosage <- dplyr::left_join(gm, t_gwa_gd, by = c("id" = "snp"))
write.table(gwas_dosage, "./data/gwas_dosage_dt.txt", row.names = F)

gwas_bim <- data.frame(gwas_dosage[,c(2,1)], dist = rep(0, nrow(gwas_dosage)), gwas_dosage[,c(3:5)])
gwas_fam <- read.table("./data/ld_pruned/data_pruned_final_numeric.fam")
colnames(gwas_fam) <- c("famid", "id", "father", "mother", "sex", "pheno")


pheno.dt <- read.csv("./data/pheno_dt_final.csv")
pheno.dt2 <- dplyr::mutate(pheno.dt, )
pheno.dt2 <- filter(pheno.dt, IID %in% gwas_fam$id) |> 
  mutate(Flower_col2 = recode(Flower_col, "pale purple" = 0, "purple" = 1, "dark purple" = 1),
         Leaf_pigm2 = recode(Leaf_pigm,"Non-Pigmented" = 0,"Pigmented" = 1),
         Leaf_size2 = recode(Leaf_size, "narrow" = 0, "broad" = 1),
         Leaf_col2 = ifelse(Leaf_col %in% c("pink", "dark purple", "green purple", "purple green", 
                                            "purple", "redish purple"), 0, 1))
saveRDS(pheno.dt2, "./data/gwas_pheno_dt.RDS")
pheno.dt2 <- readRDS("./data/gwas_pheno_dt.RDS")

gwas_gen <- transpoz(gwas_dosage[,-c(2:5)], new_col1_name = "id") |> 
  tibble::column_to_rownames(var = "id") |> as.matrix()

which(!(rownames(gwas_gen) %in% pheno.dt2$IID)) #### 9 which is Results is "A5"
gwas_gen <-  gwas_gen[-9,]
gwas_fam <- filter(gwas_fam, id != "A5")

#### Create a bed matrix, a required format
gwas_bed <- as.bed.matrix(gwas_gen, gwas_fam, gwas_bim)

#### Create Genomic relationship matrix
Amat_format <- snp_dt_012_to_101(gwas_gen)
grm = rrBLUP::A.mat(Amat_format)

#### Import the Q matrix form ADMIXTURE
pop_str <- read.table("./data/admixture_results3/data_pruned_final_numeric.3.Q")
pop_str <- as.matrix(pop_str[-9,])

#### Save all the necessary objects and reload them when needed.
gwas_list <- list(gwas_gen = gwas_gen, gwas_fam = gwas_fam, gwas_bim = gwas_bim,
                  grm = grm, pop_str = pop_str, pheno.dt2 = pheno.dt2)
saveRDS(gwas_list, "./data/gwas_list_40_contigs.RDS")

gwas_imputs <- readRDS("./data/gwas_list_40_contigs.RDS")
list2env(gwas_imputs, envir = .GlobalEnv)
gwas_bed <- as.bed.matrix(gwas_gen, gwas_fam, gwas_bim)

#### Computation of simple M significanct threshold
cor_r <- cor(gwas_gen)
eigen_values <- eigen(cor_r)
saveRDS(eigen_values, "./data/Eigen_values_genomic_cor_matrix.RDS")
eigen_values_dt <- data.frame(eigen_values = 100*sort(eigen_values$values, decreasing = T)/sum(eigen_values$values))
eigen_values_dt$cumulative <- cumsum(eigen_values_dt$eigen_values)
nrow(eigen_values_dt) - length(eigen_values_dt$cumulative[eigen_values_dt$cumulative > 99.990]) ## Results 167
### Meff = 167 => threshold = -log10(0.05/167)
threshold = -log10(0.05/167) # = 3.52 so let us use 4

#### Finally GWAS analysis 40 contigs !!!!!!!!!!!!!!!! #####
gwas_imputs <- readRDS("./data/gwas_list_RDS")
list2env(gwas_imputs, envir = .GlobalEnv)
gwas_bed <- as.bed.matrix(gwas_gen, gwas_fam, gwas_bim)

nmarkers <- nrow(gwas_bim) # =11383
bonferroni = -log10(0.05/nmarkers)
bonferroni = -log10(0.0001)

# The Error: external pointer is not valid is due to the gwas_bed object

### Leaf pigmentation
pigm_lmm <- association.test(x = gwas_bed, Y = pheno.dt2$Leaf_pigm2, 
                             K = grm, X = pop_str, method = "lmm", response = "bin") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

pigm_off <- association.test.logistic(x = gwas_bed, Y = pheno.dt2$Leaf_pigm2, 
                                      K = grm, X = pop_str, algorithm = "offset") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

pigm_amle <- association.test.logistic(x = gwas_bed, Y = pheno.dt2$Leaf_pigm2, 
                                       K = grm, X = pop_str,  algorithm = "amle") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

pig_lmm_plt <- ggplot(data = pigm_lmm, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Pigment LMM") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)
pig_off_plt <- ggplot(data = pigm_off, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Pigment OFF") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)

pig_amle_plt <- ggplot(data = pigm_amle, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Pigment AMLE") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)

pig_comb <- ggpubr::ggarrange(pig_lmm_plt, pig_off_plt, pig_amle_plt, nrow = 3)
ggsave(plot = pig_comb, filename = "./output/pigment_gwas_com.jpg",
       width = 6, height = 6.2, dpi = 300)

### Leaf size
size_lmm <- association.test(x = gwas_bed, Y = pheno.dt2$Leaf_size2, 
                             K = grm, X = pop_str, method = "lmm", response = "bin") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

size_off <- association.test.logistic(x = gwas_bed, Y = pheno.dt2$Leaf_size2, 
                                      K = grm, X = pop_str, algorithm = "offset") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

size_amle <- association.test.logistic(x = gwas_bed, Y = pheno.dt2$Leaf_size2, 
                                       K = grm, X = pop_str,  algorithm = "amle") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

siz_lmm_plt <- ggplot(data = size_lmm, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Leaf size LMM") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)
siz_off_plt <- ggplot(data = size_off, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Leaf size OFF") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)

siz_amle_plt <- ggplot(data = size_amle, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Leaf size AMLE") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)
siz_comb <- ggpubr::ggarrange(siz_lmm_plt, siz_off_plt, siz_amle_plt, nrow = 3)
ggsave(plot = siz_comb, filename = "./output/leaf_size_gwas_comb.jpg",
       width = 6, height = 6.2, dpi = 300)

### Leaf color
lcol_lmm <- association.test(x = gwas_bed, Y = pheno.dt2$Leaf_col2, 
                             K = grm, X = pop_str, method = "lmm", response = "bin") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

lcol_off <- association.test.logistic(x = gwas_bed, Y = pheno.dt2$Leaf_col2, 
                                      K = grm, X = pop_str, algorithm = "offset") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

lcol_amle <- association.test.logistic(x = gwas_bed, Y = pheno.dt2$Leaf_col2, 
                                       K = grm, X = pop_str,  algorithm = "amle") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

lcol_lmm_plt <- ggplot(data = lcol_lmm, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Leaf color LMM") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)
lcol_off_plt <- ggplot(data = lcol_off, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Leaf color OFF") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)

lcol_amle_plt <- ggplot(data = lcol_amle, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Leaf color AMLE") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)
lcol_comb <- ggpubr::ggarrange(lcol_lmm_plt, lcol_off_plt, lcol_amle_plt, nrow = 3)
ggsave(plot = lcol_comb, filename = "./output/leaf_color_gwas_comb.jpg",
       width = 6, height = 6.2, dpi = 300)

## Flower color
fcol_lmm <- association.test(x = gwas_bed, Y = pheno.dt2$Flower_col2, 
                             K = grm, X = pop_str, method = "lmm", response = "bin") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

fcol_off <- association.test.logistic(x = gwas_bed, Y = pheno.dt2$Flower_col2, 
                                      K = grm, X = pop_str, algorithm = "offset") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

fcol_amle <- association.test.logistic(x = gwas_bed, Y = pheno.dt2$Flower_col2, 
                                       K = grm, X = pop_str,  algorithm = "amle") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

fcol_lmm_plt <- ggplot(data = fcol_lmm, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Flower color LMM") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)
fcol_off_plt <- ggplot(data = fcol_off, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Flower color OFF") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)

fcol_amle_plt <- ggplot(data = fcol_amle, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Flower color AMLE") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)
fcol_comb <- ggpubr::ggarrange(fcol_lmm_plt, fcol_off_plt, fcol_amle_plt, nrow = 3)
ggsave(plot = fcol_comb, filename = "./output/flower_color_gwas_comb.jpg",
       width = 6, height = 6.2, dpi = 300)


gwas_res_list <- list(fcol_amle = fcol_amle, fcol_lmm = fcol_lmm, fcol_off = fcol_off,
                  pigm_amle = pigm_amle, pigm_lmm = pigm_lmm, pigm_off = pigm_off,
                  lcol_amle = lcol_amle, lcol_lmm = lcol_lmm, lcol_off = lcol_off,
                  size_amle = size_amle, size_lmm = size_lmm, size_off = size_off)
saveRDS(gwas_res_list, "./data/GWAS_results_40_contigs.RDS")
openxlsx::write.xlsx(gwas_list, "./data/GWAS_results_40contigs.xlsx")


#### Saving only the signficant snps and add trait and model in sep columns
library(openxlsx)
gwas_result_list <- lapply(getSheetNames("./data/GWAS_results_40contigs.xlsx"), function(s) {
  read.xlsx("./data/GWAS_results_40contigs.xlsx", sheet = s)
})
names(gwas_result_list) <- getSheetNames("./data/GWAS_results_40contigs.xlsx")

gwas_significant <- list()
for(i in seq_along(names(gwas_result_list))){
  obj <- names(gwas_result_list)[i]
  ftmp <- gwas_result_list[[obj]] |> 
    dplyr::filter(logp >=4) |> 
    dplyr::mutate(Trait = strsplit(obj, split = "_")[[1]][1],
                  Model = strsplit(obj, split = "_")[[1]][2])
  gwas_significant[[i]] <- ftmp
}

gwas_significant_dt <- data.table::rbindlist(gwas_significant[c(1,10)]) |> 
  as.data.frame()

write.csv(gwas_significant_dt, "./data/significant_snps_trait_gwas_40contigs.csv", row.names = F)


################# ----------- Creating a customized Manhattan plots --------- ###########

#### Preprocessing the data
library(openxlsx)
gwas_result_list_40 <- lapply(getSheetNames("./data/GWAS_results_40contigs.xlsx"), function(s) {
  read.xlsx("./data/GWAS_results_40contigs.xlsx", sheet = s)
})
data_40 <- gwas_result_list_40[[1]][c("chr", "pos", "id", "logp")]
dt <- data_40 |>
  dplyr::mutate(chr = factor(chr, levels = bim_summary$Var1),
                chr_num = factor(as.numeric(chr))) |> 
  dplyr::arrange(chr_num) |> 
  dplyr::mutate(pos_cum = 1:length(chr))

centers <- dt |> dplyr::group_by(chr_num) |>
  summarise_at("pos_cum", mean)

gwas_significant_dt_40 <- read.csv("./data/significant_snps_trait_gwas_40contigs.csv")
dt_sign_40 <- gwas_significant_dt_40[c("chr", "pos", "id", "logp", "Trait", "Model")]
dt_sign_40 <- dplyr::left_join(dt_sign_40, dt[,c(3,5,6)]) |> 
  dplyr::mutate(Trait = recode(Trait, fcol = "Flower color", size = "Leaf size"))

dt_40 <- dplyr::filter(dt, logp <4) ### First remove significant snps

manhplot_40 <- ggplot(dt_40, aes(x = pos_cum, y = logp)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,ymin = threshold, 
                max = Inf), fill = "#F0DAFF", alpha = 0.2) +
  geom_hline(yintercept = threshold, color = "darkred", linetype = "dashed") +
  geom_point(aes(color = factor(chr_num)), shape = 19) +
  scale_color_manual(values = rep(c("#274AB3", "#7CA4E1"), 
                                  length(unique(dt$chr_num))), guide = "none") +
  
  ggnewscale::new_scale_color() +
  
  geom_point(data = dt_sign_40, aes(x = pos_cum, y = logp, color = Trait), size = 3) + 
  # scale_color_brewer(palette = "Set1", name = "Trait") +
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2")[4:3], name = "Trait") +
  scale_x_continuous(labels = centers$chr_num, breaks = centers$pos_cum) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
  labs(x = "Contig", y = expression(bold(-log[10](p)))) +
  boris_theme_rect +
  theme(legend.position = c(0.1, 0.85), legend.direction = "vertical",
        axis.text.x = element_text(angle = 90, vjust = 0.2, size = 10))
manhplot_40
ggsave("./output/Manhattan_plot_multivar_40.jpg", width = 10, height = 3.5, 
       plot = manhplot_40, dpi = 400, units = "in")

######### Creating Q-Q plot
library(openxlsx)
gwas_result_list_40 <- lapply(getSheetNames("./data/GWAS_results_40contigs.xlsx"), function(s) {
  read.xlsx("./data/GWAS_results_40contigs.xlsx", sheet = s)
})
names(gwas_result_list_40) <- getSheetNames("./data/GWAS_results_40contigs.xlsx")
gwas_list2 <- list()
for(i in seq_along(names(gwas_result_list_40))){
  obj <- names(gwas_result_list_40)[i]
  ftmp <- gwas_result_list_40[[obj]] |> 
    dplyr::arrange(desc(logp)) |> 
    dplyr::mutate(Trait = strsplit(obj, split = "_")[[1]][1],
                  Model = strsplit(obj, split = "_")[[1]][2])
  gwas_list2[[i]] <- ftmp
}
names(gwas_list2) <- names(gwas_result_list_40)

gwas_amle <- data.table::rbindlist(gwas_list2[grep("amle", names(gwas_list2))]) |> 
  as.data.frame() |> 
  dplyr::mutate(exptp = rep(-log10(1:nrow(gwas_list2[[1]])/nrow(gwas_list2[[1]])), 4),
                Trait = recode(Trait, fcol = "Flower color", size = "Leaf size",
                                       lcol = "Leaf color", pigm = "Leaf pigmentation"))

jpeg("./output/QQplot_GWAS_AMLE.jpg", width = 4, height = 4, res = 400, units = "in")
ggplot(gwas_amle, aes(x = exptp, y = logp, color = Trait)) +
  geom_point(alpha = 0.5, shape = 19, size = 2.5) + boris_theme_rect +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = expression(bold(Expected -log[10](p))), 
       y = expression(bold(Observed -log[10](p)))) +
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2")[c(4,1,2,3)], name = "Trait") +
  theme(legend.position = c(0.3, 0.8), 
        legend.direction = "vertical")
dev.off()

gwas_off <- data.table::rbindlist(gwas_list2[grep("off", names(gwas_list2))]) |> 
  as.data.frame() |> 
  dplyr::mutate(exptp = rep(-log10(1:nrow(gwas_list2[[1]])/nrow(gwas_list2[[1]])), 4),
                Trait = recode(Trait, fcol = "Flower color", size = "Leaf size",
                               lcol = "Leaf color", pigm = "Leaf pigmentation"))

jpeg("./output/QQplot_GWAS_OFF.jpg", width = 4, height = 4, res = 400, units = "in")
ggplot(gwas_off, aes(x = exptp, y = logp, color = Trait)) +
  geom_point(alpha = 0.5, shape = 19, size = 2.5) + boris_theme_rect +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = expression(bold(Expected -log[10](p))), 
       y = expression(bold(Observed -log[10](p)))) +
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2")[c(4,1,2,3)], name = "Trait") +
  theme(legend.position = c(0.25, 0.8), 
        legend.direction = "vertical")
dev.off()

gwas_lmm <- data.table::rbindlist(gwas_list2[grep("lmm", names(gwas_list2))]) |> 
  as.data.frame() |> 
  dplyr::mutate(exptp = rep(-log10(1:nrow(gwas_list2[[1]])/nrow(gwas_list2[[1]])), 4),
                Trait = recode(Trait, fcol = "Flower color", size = "Leaf size",
                               lcol = "Leaf color", pigm = "Leaf pigmentation"))

jpeg("./output/QQplot_GWAS_LMM.jpg", width = 4, height = 4, res = 400, units = "in")
ggplot(gwas_lmm, aes(x = exptp, y = logp, color = Trait)) +
  geom_point(alpha = 0.5, shape = 19, size = 2.5) + boris_theme_rect +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = expression(bold(Expected -log[10](p))), 
       y = expression(bold(Observed -log[10](p)))) +
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2")[c(4,1,2,3)], name = "Trait") +
  theme(legend.position = c(0.3, 0.8), 
        legend.direction = "vertical")
dev.off()

######## GWAS analysis with whole data set ============== ####################

bim <- read.table("./data/ld_pruned/data_pruned_final_numeric.bim")
colnames(bim) <- c("snp", "chr", "void", "pos", "allele1", "allele2")
head(bim)
bim2 <- bim |> 
  dplyr::mutate(chr = sub("\\..*", "", chr), snp = paste0("SNP",1:nrow(bim))) |> 
  dplyr::select(-void)
head(bim2)
bim_summary <- bim2$chr |> table() |> as.data.frame() |>
  dplyr::arrange(desc(Freq)) |> 
  dplyr::mutate(Category = ifelse(Freq >= 24, "Selected", "Grouped"))#### Check if not missing those with high FST values
head(bim_summary)

bim_final <- bim2
saveRDS(bim_final, "./data/genetic_map_final_whole_set.RDS")

##Preparing GWAS dosage data
dosage_dt <- read.table("./data/dosage_output.raw", header = TRUE)
geno_mat <- as.matrix(read.table("./data/geno_mat_final.txt"))
snps <- sub("\\..*", "", colnames(dosage_dt[,-c(1:6)]))
gwas_snps <- which(snps %in% bim_summary$Var1)
gwas_GD <- geno_mat[,gwas_snps] |> as.data.frame() |> 
  tibble::rownames_to_column(var = "ID")
saveRDS(gwas_GD, "./data/gwas_genomic_data_whole_set.RDS")
rm(dosage_dt)

### GWAS analysis with milorgwas
library(milorGWAS)
help(package = "milorGWAS")
gwas_GD <- readRDS("./data/gwas_genomic_data_whole_set.RDS")
t_gwa_gd <- transpoz(gwas_GD, new_col1_name = "snp")
gm <- readRDS("./data/genetic_map_final_whole_set.RDS")
colnames(gm) <- c("id", "chr", "pos", "A1", "A2")
gwas_dosage <- dplyr::left_join(gm, t_gwa_gd, by = c("id" = "snp"))
write.table(gwas_dosage, "./data/gwas_dosage_dt_whole_set.txt", row.names = F)

gwas_bim <- data.frame(gwas_dosage[,c(2,1)], 
                       dist = rep(0, nrow(gwas_dosage)), 
                       gwas_dosage[,c(3:5)])
saveRDS(gwas_bim, "./data/gwas_bim_whole_set.RDS")
gwas_fam <- read.table("./data/ld_pruned/data_pruned_final_numeric.fam")
colnames(gwas_fam) <- c("famid", "id", "father", "mother", "sex", "pheno")

##Import alread formatted phenotypic data
pheno.dt2 <- readRDS("./data/gwas_pheno_dt.RDS")

gwas_gen <- transpoz(gwas_dosage[,-c(2:5)], new_col1_name = "id") |> 
  tibble::column_to_rownames(var = "id") |> as.matrix()

which(!(rownames(gwas_gen) %in% pheno.dt2$IID)) #### 9 which is Results is "A5"
gwas_gen <-  gwas_gen[-9,]
gwas_fam <- filter(gwas_fam, id != "A5")
saveRDS(gwas_gen, "./data/gwas_gen_whole_set.RDS")
saveRDS(gwas_fam, "./data/gwas_fam_final_whole_set.RDS")

#### Create a bed matrix, a required format for milorGWAS
gwas_bim <- readRDS("./data/gwas_bim_whole_set.RDS")
gwas_gen <- readRDS("./data/gwas_gen_whole_set.RDS")
gwas_fam <- readRDS("./data/gwas_fam_final_whole_set.RDS")
gwas_bed <- gaston::as.bed.matrix(gwas_gen, gwas_fam, gwas_bim)

#### Create Genomic relationship matrix
Amat_format <- snp_dt_012_to_101(gwas_gen)
grm = rrBLUP::A.mat(Amat_format)

#### Import the Q matrix form ADMIXTURE
pop_str <- read.table("./data/admixture_results3/data_pruned_final_numeric.3.Q")
pop_str <- as.matrix(pop_str[-9,])

gwas_list <- list(gwas_bed = gwas_bed, grm = grm, Qmat = pop_str)
saveRDS(gwas_list, "./data/gwas_list_whole_set.RDS")

gwas_imputs <- readRDS("./data/gwas_list_whole_set.RDS")
list2env(gwas_imputs, envir = .GlobalEnv)

#### Computation of simple M significanct threshold
cor_r <- cor(gwas_gen)
eigen_values <- eigen(cor_r)
saveRDS(eigen_values, "./data/Eigen_values_genomic_cor_matrix.RDS")
eigen_values_dt <- data.frame(eigen_values = 100*sort(eigen_values$values, decreasing = T)/sum(eigen_values$values))
eigen_values_dt$cumulative <- cumsum(eigen_values_dt$eigen_values)
nrow(eigen_values_dt) - length(eigen_values_dt$cumulative[eigen_values_dt$cumulative > 99.990]) ## Results 167
### Meff = 167 => threshold = -log10(0.05/167)
threshold = -log10(0.05/167) # = 3.52 so let us use 4

#### Finally GWAS analysis whole set!!!!!!!!!!!!!!!! #####
library(milorGWAS)
nmarkers <- nrow(gwas_bim) # =22447
bonferroni = -log10(0.05/nmarkers) #5.65
bonferroni = -log10(0.0001) #4.00, SimpleM

### Leaf pigmentation
pigm_lmm <- association.test(x = gwas_bed, Y = pheno.dt2$Leaf_pigm2, 
                             K = grm, X = pop_str, method = "lmm", response = "bin") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

pigm_off <- association.test.logistic(x = gwas_bed, Y = pheno.dt2$Leaf_pigm2, 
                                      K = grm, X = pop_str, algorithm = "offset") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

pigm_amle <- association.test.logistic(x = gwas_bed, Y = pheno.dt2$Leaf_pigm2, 
                                       K = grm, X = pop_str,  algorithm = "amle") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

pig_lmm_plt <- ggplot(data = pigm_lmm, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Pigment LMM") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)
pig_off_plt <- ggplot(data = pigm_off, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Pigment OFF") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)

pig_amle_plt <- ggplot(data = pigm_amle, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Pigment AMLE") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)

pig_comb <- ggpubr::ggarrange(pig_lmm_plt, pig_off_plt, pig_amle_plt, nrow = 3)
ggsave(plot = pig_comb, filename = "./output/pigment_gwas_comb_whole_set.jpg",
       width = 6, height = 6.2, dpi = 300)

### Leaf size
size_lmm <- association.test(x = gwas_bed, Y = pheno.dt2$Leaf_size2, 
                             K = grm, X = pop_str, method = "lmm", response = "bin") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

size_off <- association.test.logistic(x = gwas_bed, Y = pheno.dt2$Leaf_size2, 
                                      K = grm, X = pop_str, algorithm = "offset") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

size_amle <- association.test.logistic(x = gwas_bed, Y = pheno.dt2$Leaf_size2, 
                                       K = grm, X = pop_str,  algorithm = "amle") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

siz_lmm_plt <- ggplot(data = size_lmm, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Leaf size LMM") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)
siz_off_plt <- ggplot(data = size_off, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Leaf size OFF") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)

siz_amle_plt <- ggplot(data = size_amle, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Leaf size AMLE") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)

siz_comb <- ggpubr::ggarrange(siz_lmm_plt, siz_off_plt, siz_amle_plt, nrow = 3)
ggsave(plot = siz_comb, filename = "./output/leaf_size_gwas_comb_whole_set.jpg",
       width = 6, height = 6.2, dpi = 300)

### Leaf color
lcol_lmm <- association.test(x = gwas_bed, Y = pheno.dt2$Leaf_col2, 
                             K = grm, X = pop_str, method = "lmm", response = "bin") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

lcol_off <- association.test.logistic(x = gwas_bed, Y = pheno.dt2$Leaf_col2, 
                                      K = grm, X = pop_str, algorithm = "offset") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

lcol_amle <- association.test.logistic(x = gwas_bed, Y = pheno.dt2$Leaf_col2, 
                                       K = grm, X = pop_str,  algorithm = "amle") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

lcol_lmm_plt <- ggplot(data = lcol_lmm, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Leaf color LMM") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)
lcol_off_plt <- ggplot(data = lcol_off, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Leaf color OFF") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)

lcol_amle_plt <- ggplot(data = lcol_amle, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Leaf color AMLE") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)
lcol_comb <- ggpubr::ggarrange(lcol_lmm_plt, lcol_off_plt, lcol_amle_plt, nrow = 3)
ggsave(plot = lcol_comb, filename = "./output/leaf_color_gwas_comb_whole_set.jpg",
       width = 6, height = 6.2, dpi = 300)

## Flower color
fcol_lmm <- association.test(x = gwas_bed, Y = pheno.dt2$Flower_col2, 
                             K = grm, X = pop_str, method = "lmm", response = "bin") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

fcol_off <- association.test.logistic(x = gwas_bed, Y = pheno.dt2$Flower_col2, 
                                      K = grm, X = pop_str, algorithm = "offset") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

fcol_amle <- association.test.logistic(x = gwas_bed, Y = pheno.dt2$Flower_col2, 
                                       K = grm, X = pop_str,  algorithm = "amle") |> 
  as.data.frame() |>  dplyr::mutate(logp = -log10(p))

fcol_lmm_plt <- ggplot(data = fcol_lmm, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Flower color LMM") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)
fcol_off_plt <- ggplot(data = fcol_off, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Flower color OFF") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)

fcol_amle_plt <- ggplot(data = fcol_amle, aes(x= pos, y = logp)) + 
  geom_point(color = "grey40", size = 1.5) +
  labs(x = "SNP marker", y = "-log10(p-values)", title = "Flower color AMLE") + boris_theme_rect +
  scale_y_continuous(limits = c(0,8), breaks = seq(0,8, 2)) +
  scale_x_continuous(limits = c(0,34e6), breaks = seq(0,30e6,5e6)) +
  geom_hline(yintercept = bonferroni, color = "darkred", linetype = "dashed", linewidth = 1)

fcol_comb <- ggpubr::ggarrange(fcol_lmm_plt, fcol_off_plt, fcol_amle_plt, nrow = 3)
ggsave(plot = fcol_comb, filename = "./output/flower_color_gwas_comb_whole_set.jpg",
       width = 6, height = 6.2, dpi = 300)


gwas_result_list <- list(fcol_amle = fcol_amle, fcol_lmm = fcol_lmm, fcol_off = fcol_off,
                  pigm_amle = pigm_amle, pigm_lmm = pigm_lmm, pigm_off = pigm_off,
                  lcol_amle = lcol_amle, lcol_lmm = lcol_lmm, lcol_off = lcol_off,
                  size_amle = size_amle, size_lmm = size_lmm, size_off = size_off)
openxlsx::write.xlsx(gwas_list, "./data/GWAS_results_whole_data_set.xlsx")

#### Saving only the signficant snps and add trait and model in sep columns
#### Saving only the signficant snps and add trait and model in sep columns
library(openxlsx)
gwas_result_list <- lapply(getSheetNames("./data/GWAS_results_whole_data_set.xlsx"), function(s) {
  read.xlsx("./data/GWAS_results_whole_data_set.xlsx", sheet = s)
})
names(gwas_result_list) <- getSheetNames("./data/GWAS_results_whole_data_set.xlsx")

gwas_significant <- list()
for(i in seq_along(names(gwas_result_list))){
  obj <- names(gwas_result_list)[i]
  ftmp <- gwas_result_list[[obj]] |> 
    dplyr::filter(logp >=4) |> 
    dplyr::mutate(Trait = strsplit(obj, split = "_")[[1]][1],
                  Model = strsplit(obj, split = "_")[[1]][2])
  gwas_significant[[i]] <- ftmp
}

gwas_significant_dt <- data.table::rbindlist(gwas_significant[c(1,4,7,10)]) |> 
  as.data.frame()

write.csv(gwas_significant_dt, "./data/significant_snps_trait_gwas_whole_data.csv", row.names = F)

## Plotting across_traits
p_amle_comb <- ggpubr::ggarrange(pig_amle_plt, fcol_amle_plt, lcol_amle_plt, 
                                 siz_amle_plt, nrow = 4)
ggsave(plot = p_amle_comb, filename = "./output/amle_gwas_comb_whole_set.jpg",
       width = 6, height = 8.2, dpi = 300)

p_lmm_comb <- ggpubr::ggarrange(pig_lmm_plt, fcol_lmm_plt, lcol_lmm_plt,
                                siz_lmm_plt, nrow = 4)
ggsave(plot = p_lmm_comb, filename = "./output/LMM_gwas_comb_whole_set.jpg",
       width = 6, height = 8.2, dpi = 300)

p_off_comb <- ggpubr::ggarrange(pig_off_plt, fcol_off_plt, lcol_off_plt,
                                siz_off_plt, nrow = 4)
ggsave(plot = p_off_comb, filename = "./output/OFF_gwas_comb_whole_set.jpg",
       width = 6, height = 8.2, dpi = 300)

###### Creating Manhattan plot for whole dataset ------------
library(openxlsx)

##First loading stuff of 40 contigs to keep same numbering of contigs
gwas_result_list_40 <- lapply(getSheetNames("./data/GWAS_results_40contigs.xlsx"), function(s) {
  read.xlsx("./data/GWAS_results_40contigs.xlsx", sheet = s)
})
data_40 <- gwas_result_list_40[[1]][c("chr", "pos", "id", "logp")]

dt <- dplyr::arrange(data_40, chr) |>
  dplyr::mutate(chr_num = factor(as.numeric(factor(chr))),
                pos_cum = 1:length(chr))

##Now loading stuff of whole data set
gwas_result_list_whole <- lapply(getSheetNames("./data/GWAS_results_whole_data_set.xlsx"), function(s) {
  read.xlsx("./data/GWAS_results_whole_data_set.xlsx", sheet = s)
})

data_whole <- gwas_result_list_whole[[1]][c("chr", "pos", "id", "logp")]
dt_whole <-  left_join(data_whole, dt[,c(3,5,6)]) |> 
  arrange(pos_cum) |> 
  dplyr::mutate(chr_num = ifelse(is.na(chr_num), "Other", chr_num),
                pos_cum = 1:length(chr))

centers <- dt_whole |> dplyr::group_by(chr_num) |>
  summarise_at("pos_cum", mean)

gwas_significant_dt <- read.csv("./data/significant_snps_trait_gwas_whole_data.csv")
dt_sign <- gwas_significant_dt[c("chr", "pos", "id", "logp", "Trait", "Model")]
dt_sign <- dplyr::left_join(dt_sign, dt_whole[,c(3,5,6)])|> 
  dplyr::mutate(Trait = recode(Trait, fcol = "Flower color", size = "Leaf size",
                               lcol = "Leaf color", pigm = "Leaf pigmentation"))

threshold = 4
ymax = ceiling(max(dt_sign$logp, na.rm = T))

library(ggnewscale)

jpeg("./output/Manhattan_plot_GWAS_whole_set.jpg", 
     width = 12, height = 3.5, res = 400, units = "in")

ggplot(dt_whole, aes(x = pos_cum, y = logp)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,ymin = threshold, 
                max = Inf), fill = "#F0DAFF", alpha = 0.2) +
  geom_hline(yintercept = threshold, color = "darkred", linetype = "dashed") +
  geom_point(aes(color = factor(chr_num)), shape = 19) +
  scale_color_manual(values = rep(c("#274AB3", "#7CA4E1"), 
                                  length(unique(dt$chr_num))), guide = "none") +
  
  ggnewscale::new_scale_color() +
  
  geom_point(data = dt_sign, aes(x = pos_cum, y = logp, color = Trait), size = 3) + 
  # scale_color_brewer(palette = "Set1", name = "Trait") +
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2")[c(4,1,2,3)], name = "Trait") +
  # scale_color_manual(values = c("purple", "darkred")) +
  scale_x_continuous(labels = centers$chr_num, breaks = centers$pos_cum) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
  labs(x = "Contig", y = expression(bold(-log[10](p)))) +
  boris_theme_rect +
  theme(legend.position = c(0.07, 0.45), legend.direction = "vertical",
        axis.text.x = element_text(angle = 90, vjust = 0.5),
        legend.margin = margin(0,0,0,0), 
        legend.text = element_text(size =10))
dev.off()

######### Creating Q-Q plot
library(openxlsx)
gwas_result_list <- lapply(getSheetNames("./data/GWAS_results_whole_data_set.xlsx"), function(s) {
  read.xlsx("./data/GWAS_results_whole_data_set.xlsx", sheet = s)
})
names(gwas_result_list) <- getSheetNames("./data/GWAS_results_whole_data_set.xlsx")

gwas_list2 <- list()
for(i in seq_along(names(gwas_result_list))){
  obj <- names(gwas_result_list)[i]
  ftmp <- gwas_result_list[[obj]] |>
    dplyr::arrange(desc(logp)) |>
    dplyr::mutate(Trait = strsplit(obj, split = "_")[[1]][1],
                  Model = strsplit(obj, split = "_")[[1]][2])
  gwas_list2[[i]] <- ftmp
}
names(gwas_list2) <- names(gwas_result_list)

gwas_amle <- data.table::rbindlist(gwas_list2[grep("amle", names(gwas_list2))]) |> 
  as.data.frame() |> 
  dplyr::mutate(exptp = rep(-log10(1:nrow(gwas_list2[[1]])/nrow(gwas_list2[[1]])), 4),
                Trait = recode(Trait, fcol = "Flower color", size = "Leaf size",
                               lcol = "Leaf color", pigm = "Leaf pigmentation"))

jpeg("./output/QQplot_GWAS_AMLE_whole_data.jpg", width = 4, height = 4, res = 400, units = "in")
ggplot(gwas_amle, aes(x = exptp, y = logp, color = Trait)) +
  geom_point(alpha = 0.5, shape = 19, size = 2.5) + boris_theme_rect +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = expression(bold(Expected -log[10](p))), 
       y = expression(bold(Observed -log[10](p)))) +
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2")[c(4,1,2,3)], name = "Trait") +
  theme(legend.position = c(0.3, 0.8), 
        legend.direction = "vertical")
dev.off()

gwas_off <- data.table::rbindlist(gwas_list2[grep("off", names(gwas_list2))]) |> 
  as.data.frame() |> 
  dplyr::mutate(exptp = rep(-log10(1:nrow(gwas_list2[[1]])/nrow(gwas_list2[[1]])), 4),
                Trait = recode(Trait, fcol = "Flower color", size = "Leaf size",
                               lcol = "Leaf color", pigm = "Leaf pigmentation"))

jpeg("./output/QQplot_GWAS_OFF_whole_data.jpg", width = 4, height = 4, res = 400, units = "in")
ggplot(gwas_off, aes(x = exptp, y = logp, color = Trait)) +
  geom_point(alpha = 0.5, shape = 19, size = 2.5) + boris_theme_rect +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = expression(bold(Expected -log[10](p))), 
       y = expression(bold(Observed -log[10](p)))) +
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2")[c(4,1,2,3)], name = "Trait") +
  theme(legend.position = c(0.25, 0.8), 
        legend.direction = "vertical")
dev.off()

gwas_lmm <- data.table::rbindlist(gwas_list2[grep("lmm", names(gwas_list2))]) |> 
  as.data.frame() |> 
  dplyr::mutate(exptp = rep(-log10(1:nrow(gwas_list2[[1]])/nrow(gwas_list2[[1]])), 4),
                Trait = recode(Trait, fcol = "Flower color", size = "Leaf size",
                               lcol = "Leaf color", pigm = "Leaf pigmentation"))

jpeg("./output/QQplot_GWAS_LMM_whole_data.jpg", width = 4, height = 4, res = 400, units = "in")
ggplot(gwas_lmm, aes(x = exptp, y = logp, color = Trait)) +
  geom_point(alpha = 0.5, shape = 19, size = 2.5) + boris_theme_rect +
  geom_abline(slope = 1, intercept = 0) +
  labs(x = expression(bold(Expected -log[10](p))), 
       y = expression(bold(Observed -log[10](p)))) +
  scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2")[c(4,1,2,3)], name = "Trait") +
  theme(legend.position = c(0.3, 0.8), 
        legend.direction = "vertical")
dev.off()


####### Genome-Environment association with climatic data ###############

####1. Download weather data
### The samples were collected between Nov 15-20 2023
coordinates <- read.csv("./data/Sample_coordinates.csv")
unique_coordinates <- coordinates[,-c(1:2)] |> 
  dplyr::distinct() |> na.omit()
hist(unique_coordinates$longitude)
hist(unique_coordinates$Latitude)

#### Using the nasapower in apsimx package
install.packages("apsimx")
install.packages("nasapower")

library(apsimx)
library(nasapower)

data_list <- vector("list", nrow(unique_coordinates))

for (i in seq_len(nrow(unique_coordinates))) {
  print(paste0("Starting to download", i, "th data"))
  lonlat <- c(unique_coordinates$longitude[i],
              unique_coordinates$Latitude[i])
  
  data_list[[i]] <- get_power_apsim_met(lonlat = lonlat,
                                        dates = c("1993-01-01", "2023-12-31")
                                        )
  print(paste0("Finished downloading ", i, "th data"))
}

names(data_list) <- paste0("Coord_", unique_coordinates$longitude,
                           "_", unique_coordinates$Latitude)

saveRDS(data_list, "./data/Weather_data_30_years.RDS")
openxlsx::write.xlsx(data_list, "./data/Weather_data_30_years.xlsx")

### Process and summarize weather data --------------------
comb_coord <- paste0("Coord_", unique_coordinates$longitude,
                    "_", unique_coordinates$Latitude)
data_list <- readRDS( "./data/Weather_data_30_years.RDS")
for (j in seq_along(comb_coord))  {
  data_list[[j]] <- dplyr::mutate(as.data.frame(data_list[[j]]), coordinates = comb_coord[j])
}

pheno.dt_sub <- read.csv("./data/pheno_dt_final.csv")[,1:2]
loc_coord <- dplyr::left_join(pheno.dt_sub, coordinates,
                       by = c("IID"="Sample.name")) |> 
  dplyr::mutate(comb_coord = paste0("Coord_", longitude,"_", Latitude))

data_list_df <- data.table::rbindlist(data_list)
#Season 1: days 74 - 196 (Mach 15th - July 15th)
#season 2: days 227 349 (Aug 15th - Dec 15th)
#Here, we are computing the annual averages using 2 seasons data

data_list_df <- data_list_df |> 
  dplyr::mutate(Season = ifelse(day > 73 & day < 197, "season 1",
                    ifelse(day > 226 & day < 350, "season 2", "remove"))) |> 
  dplyr::filter(Season != "remove") |> 
  dplyr::group_by(coordinates, year) |> 
  dplyr::summarise(dplyr::across(radn:windspeed, ~ mean(.x, na.rm = TRUE)), .groups = "drop")

## Add location names to it. There was a duplicate between rows 19 & 20, so filtered out 19
uniq_loc_coord <- dplyr::distinct(dplyr::select(loc_coord, Location, comb_coord))
data_list_df_plt <- dplyr::left_join(data_list_df, 
                                 uniq_loc_coord[-19,],
                                 by = c("coordinates" = "comb_coord")) |> 
  dplyr::filter(!is.na(Location))

## Let us explore the data
loc_colors <- c(ggpubfigs::friendly_pal("muted_nine"), "#FE6100", "grey50")
# loc_order <- pheno.dt %>% count(Location) %>% arrange(desc(n)) %>% pull(Location)
loc_order <- c("Iganga", "Namutumba", "Ngora", "Serere", "Kaberamaido", "Pallisa",
               "Bukedea", "Kumi", "Tororo", "Busia", "Kalaki")
data_list_df_plt <- dplyr::mutate(data_list_df_plt, 
                              Location = factor(Location, levels = loc_order))
prh <- ggplot(data_list_df_plt, aes(x = year, y = rh, color = Location)) +
  geom_line(aes(group = coordinates), linewidth = 0.8, alpha = 0.5) +
  scale_color_manual(values = loc_colors) +
  labs(x = "Year", y = "Relative humidity (%)") +
  scale_x_continuous(limits = c(1993, 2023), breaks = seq(1993, 2023, 10))+
  scale_y_continuous(limits = c(0,83), breaks = seq(0,80,20)) +
  theme_bw() + boris_theme_rect +
  theme(legend.direction = "vertical", legend.position = "right")

pmaxt <- ggplot(data_list_df_plt, aes(x = year, y = maxt, color = Location)) +
  geom_line(aes(group = coordinates), linewidth = 0.8, alpha = 0.5) +
  scale_color_manual(values = loc_colors) +
  labs(x = "Year", y = expression(bold("Maximum temperature (" * degree * "C)"))) +
  scale_x_continuous(limits = c(1993, 2023), breaks = seq(1993, 2023, 10))+
  scale_y_continuous(limits = c(0,32), breaks = seq(0,30,5)) +
  theme_bw() + boris_theme_rect +
  theme(legend.direction = "vertical", legend.position = "right")

pmint <- ggplot(data_list_df_plt, aes(x = year, y = mint, color = Location)) +
  geom_line(aes(group = coordinates), linewidth = 0.8, alpha = 0.5) +
  scale_color_manual(values = loc_colors) +
  labs(x = "Year", y = expression(bold("Minimum temperature (" * degree * "C)"))) +
  scale_x_continuous(limits = c(1993, 2023), breaks = seq(1993, 2023, 10))+
  scale_y_continuous(limits = c(0,32), breaks = seq(0,30,5)) +
  theme_bw() + boris_theme_rect +
  theme(legend.direction = "vertical", legend.position = "right")

prain <- ggplot(data_list_df_plt, aes(x = year, y = 25.4*rain, color = Location)) +
  geom_line(aes(group = coordinates), linewidth = 0.8, alpha = 0.5) +
  scale_color_manual(values = loc_colors) +
  labs(x = "Year", y = "Precipitation (mm)") +
  scale_x_continuous(limits = c(1993, 2023), breaks = seq(1993, 2023, 10))+
  scale_y_continuous(limits = c(0,261), breaks = seq(0,250,50)) +
  theme_bw() + boris_theme_rect +
  theme(legend.direction = "vertical", legend.position = "right")

pradn <- ggplot(data_list_df_plt, aes(x = year, y = radn, color = Location)) +
  geom_line(aes(group = coordinates), linewidth = 0.8, alpha = 0.5) +
  scale_color_manual(values = loc_colors)+
  labs(x = "Year", y = expression(bold("Radiation (" *W/m^2* ")"))) +
  scale_x_continuous(limits = c(1993, 2023), breaks = seq(1993, 2023, 10))+
  scale_y_continuous(limits = c(0,23), breaks = seq(0,20,5)) +
  theme_bw() + boris_theme_rect +
  theme(legend.direction = "vertical", legend.position = "right")

pwsp <- ggplot(data_list_df_plt, aes(x = year, y = windspeed , color = Location)) +
  geom_line(aes(group = coordinates), linewidth = 0.8, alpha = 0.5) +
  scale_color_manual(values = loc_colors)+
  labs(x = "Year", y = "Windspeed (m/s)") +
  scale_x_continuous(limits = c(1993, 2023), breaks = seq(1993, 2023, 10))+
  scale_y_continuous(limits = c(0,2.5), breaks = seq(0,2.5,0.5)) +
  theme_bw() + boris_theme_rect + 
  theme(legend.direction = "vertical", legend.position = "right")


pcomb <- ggpubr::ggarrange(pmaxt,pmint, prain,pradn, prh,pwsp,
                           nrow = 3, ncol = 2, common.legend = T,
                           legend = "right")

ggsave("./output/Weather_data_profile.jpg", plot = pcomb, 
       width = 7, height = 8, dpi = 300, units = "in")

### Pivoting wider to have #years times columns per variables
data_list_df2 <- data_list_df |> 
  tidyr::pivot_wider(names_from = c(year), 
                     values_from = c(radn:windspeed),
                     names_glue = "{.value}_{year}")
weather_mat <- data_list_df2 |> 
  tibble::column_to_rownames(var = "coordinates") |> 
  as.matrix()

### Save weather_mat in case
saveRDS(weather_mat, "./data/Matrix_of_weather_data.RDS")

#### create separate matrices for the variables
weather_mat <- readRDS("./data/Matrix_of_weather_data.RDS")
radn_mat <- weather_mat[,grep("radn", colnames(weather_mat))]
rain_mat <- weather_mat[,grep("rain", colnames(weather_mat))]
mint_mat <- weather_mat[,grep("mint", colnames(weather_mat))]
maxt_mat <- weather_mat[,grep("maxt", colnames(weather_mat))]
rh_mat <- weather_mat[,grep("rh", colnames(weather_mat))]
windspeed_mat <- weather_mat[,grep("windspeed", colnames(weather_mat))]


### run PCA of full data set
env_pca <- prcomp(weather_mat, center = T, scale. = T)
screeplot(env_pca)
biplot(env_pca)
env_scores <- env_pca$x[,1:2] |> as.data.frame() |> 
  dplyr::rename_all(~paste0("env", .)) |> 
  tibble::rownames_to_column(var = "comb_coord")

## Plotting pca scores for full data set
env_pca <- prcomp(weather_mat, center = T, scale. = T)
screeplot(env_pca)
biplot(env_pca)

plot((env_pca$sdev^2)[1:50], type = "l", xlab = "Number of axes", ylab = "Component variance")
plot((100*(env_pca$sdev^2)/sum(env_pca$sdev^2))[1:50], type = "l", 
     xlab = "Number of axes", ylab = "env_pca variance contribution (%)")
summary(env_pca)
loadings <- abs(env_pca$rotation[,1:3]) |> as.data.frame() |> 
  dplyr::arrange(desc(PC1))
View(loadings)
cor(weather_mat)
## Saving the scores of first 3 env_pcas
scor <- as.data.frame(env_pca$x[,(1:3)])
scor$Label <- rownames(weather_mat)
write.csv(scor, "./data/env_pca_3_scores.csv", row.names = F)

### Plotting env_pca scores full data
library(tidyverse)

env_pca1_lab <- "env_pca1 (72.28%)"; env_pca2_lab = "env_pca2 (12.83%)"; env_pca3_lab = "env_pca3 (6.82%)"

ggplot(scor, aes(x = PC1, y = PC2)) + geom_point(size =3) + 
  labs(x = env_pca1_lab, y = env_pca2_lab) + boris_theme_rect

### Running PCA on individual variables
#### Radiation 
radn_pca <- prcomp(radn_mat, center = T, scale. = T)
radn_scores <- radn_pca$x[,1:2] |> as.data.frame() |> 
  dplyr::rename_all(~paste0("radn_", .)) |> 
  tibble::rownames_to_column(var = "comb_coord")

#### Precipitation
rain_pca <- prcomp(rain_mat, center = T, scale. = T)
rain_scores <- rain_pca$x[,1:2] |> as.data.frame() |> 
  dplyr::rename_all(~paste0("rain_", .)) |> 
  tibble::rownames_to_column(var = "comb_coord")

### Minimum temperature
mint_pca <- prcomp(mint_mat, center = T, scale. = T)
mint_scores <- mint_pca$x[,1:2] |> as.data.frame() |> 
  dplyr::rename_all(~paste0("mint_", .)) |> 
  tibble::rownames_to_column(var = "comb_coord")

##Maximum temperature
maxt_pca <- prcomp(maxt_mat, center = T, scale. = T)
maxt_scores <- maxt_pca$x[,1:2] |> as.data.frame() |> 
  dplyr::rename_all(~paste0("maxt_", .)) |> 
  tibble::rownames_to_column(var = "comb_coord")

## Relative humidity
rh_pca <- prcomp(rh_mat, center = T, scale. = T)
rh_scores <- rh_pca$x[,1:2] |> as.data.frame() |> 
  dplyr::rename_all(~paste0("rh_", .)) |> 
  tibble::rownames_to_column(var = "comb_coord")

##Wind speed
windspeed_pca <- prcomp(windspeed_mat, center = T, scale. = T)
wind_scores <- windspeed_pca$x[,1:2] |> as.data.frame() |>
  dplyr::rename_all(~ paste0("wind_", .)) |> 
  tibble::rownames_to_column(var = "comb_coord")

hist(radn_pca$rotation[,1])
hist(rain_pca$rotation[,1])
hist(maxt_pca$rotation[,1])
hist(mint_pca$rotation[,1])
hist(rh_pca$rotation[,1])
hist(windspeed_pca$rotation[,1])
biplot(rain_pca)
for (ind in 1:6) {
  
}
### Creating simple combined screeplot
pct_screeplot <- function(pc_obj, num_pcs = 2, color = "grey60", main_label = NULL){
  stopifnot("You must select at least two PC axes"= num_pcs >= 2)
  dt <- data.frame("Variance" = (pc_obj$sdev^2)[1:num_pcs],
                   "Percent_variance" = (100*(pc_obj$sdev^2)/sum(pc_obj$sdev^2))[1:num_pcs],
                   "Axes" = paste0("PC", 1:num_pcs))
  plt <- barplot(Percent_variance~Axes, data = dt, col = color,
                 xlab = "Axes", ylab = "Variance explained (%)")
  usr <- par("usr")
  text(x= 0.8*usr[2], y = 0.9*usr[4], main_label)
}

jpeg("./output/PCA_screeplots.jpg", width = 7.5, height = 9, res = 300, units = "in")
par(mfrow = c(3,2), mar = c(4,4,1,2))
pct_screeplot(pc_obj = maxt_pca, num_pcs = 5, main_label = "Maximum temperature")
pct_screeplot(pc_obj = mint_pca, num_pcs = 5, main_label = "Minimum temperature")
pct_screeplot(pc_obj = rain_pca, num_pcs = 5, main_label = "Precipitation")
pct_screeplot(pc_obj = radn_pca, num_pcs = 5, main_label = "Radiation")
pct_screeplot(pc_obj = rh_pca, num_pcs = 5, main_label = "Relative humidity")
pct_screeplot(pc_obj = windspeed_pca, num_pcs = 5, main_label = "Wind speed")
dev.off()

### Create a color coded line screeplot
compute_pct_var <- function(pc_obj, num_pcs = 2, color = "grey60", main_label = NULL){
  stopifnot("You must select at least two PC axes"= num_pcs >= 2)
  dt <- data.frame("Variance" = (pc_obj$sdev^2)[1:num_pcs],
                   "Percent_variance" = (100*(pc_obj$sdev^2)/sum(pc_obj$sdev^2))[1:num_pcs],
                   "Axes" = paste0("PC", 1:num_pcs),
                   "Variable" = rep(main_label, num_pcs))
  return(dt)
}
pct_var <- data.table::rbindlist(list(
  Maxt = compute_pct_var(pc_obj = maxt_pca, num_pcs = 5, main_label = "Maximum temperature"),
  Mint = compute_pct_var(pc_obj = mint_pca, num_pcs = 5, main_label = "Minimum temperature"),
  Rain = compute_pct_var(pc_obj = rain_pca, num_pcs = 5, main_label = "Precipitation"),
  Rad = compute_pct_var(pc_obj = radn_pca, num_pcs = 5, main_label = "Radiation"),
  Rh = compute_pct_var(pc_obj = rh_pca, num_pcs = 5, main_label = "Relative humidity"),
  WSp = compute_pct_var(pc_obj = windspeed_pca, num_pcs = 5, main_label = "Wind speed")))

jpeg("./output/Screeplot_combined_bar.jpg", width = 4, height = 4, res = 400, units = "in")
ggplot(pct_var, aes(x = Axes, y = Percent_variance, fill = Variable)) +
  geom_col(position = "dodge", alpha=0.8) +
  scale_fill_brewer(palette = "Set1", name = "Variable") +
  theme(legend.position = c(0.3, 0.75), legend.direction = "vertical") +
  labs(x = "Axes", y = "Variance explained (%)") + boris_theme_rect +
  theme(legend.position = c(0.7, 0.75), legend.direction = "vertical",
        legend.text = element_text(size = 10))
dev.off()

jpeg("./output/Screeplot_combined_line.jpg", width = 4, height = 4, res = 400, units = "in")
ggplot(pct_var, aes(x = Axes, y = Percent_variance, fill = Variable)) +
  geom_line(aes(group = Variable, color = Variable), 
            position = position_dodge(width = 0.5), linewidth = 1.4) +
  scale_color_brewer(palette = "Set1", name = "Variable") +
  theme(legend.position = c(0.3, 0.75), legend.direction = "vertical") +
  labs(x = "Axes", y = "Variance explained (%)") + boris_theme_rect +
  theme(legend.position = c(0.7, 0.75), legend.direction = "vertical",
        legend.text = element_text(size = 10))
dev.off()

### Combining all scores
comb_scores <- list(radn_scores, rain_scores,
                    maxt_scores, mint_scores, rh_scores, 
                    wind_scores) |> purrr::reduce(full_join)

### Plotting pca scores combined as faceted plot--------------
uniq_comb_scores <- left_join(uniq_loc_coord, comb_scores)
uniq_comb_scores_long <- pivot_longer(uniq_comb_scores,
                                      cols = radn_PC1:wind_PC2,
                                      names_to = c("name", ".value"),
                                      names_sep = "_") |> 
  dplyr::mutate(variable = recode(name,
                                  "radn" = "Radiation",
                                  "rain" = "Precipitation",
                                  "maxt" ="Max temprature",
                                  "mint" = "Min temperature",
                                  "rh" = "Relative humidity",
                                  "wind" = "Wind speed"),
                Location = factor(Location, levels = loc_order))
scores_plt <- ggplot(uniq_comb_scores_long, aes(x = PC1, y = PC2, color = Location)) + 
  geom_point(size =3, alpha = 0.5) + scale_color_manual(values = loc_colors) +
  facet_wrap(~variable, ncol = 2) + 
  boris_theme_rect + geom_jitter(width = 3, height = 2, alpha = 0.5, size = 3) + 
  theme(legend.direction = "vertical", legend.position = "right")

ggsave("./output/pca_scores_environment.jpg", width = 7.5, height = 9, 
       dpi = 300, plot = scores_plt)

### Prepare genomic data for GEA analysis -------------
gwas_fam <- readRDS("./data/gwas_fam_final_whole_set.RDS")
genos <- unique(gwas_fam$id)
geno_scores_expanded <- inner_join(filter(loc_coord, IID %in% genos), comb_scores)
saveRDS(geno_scores_expanded, "./data/samples_coord_pca.RDS")
geno_scores <- geno_scores_expanded[,-c(1,3:7, seq(8,20,2))] |> 
  dplyr::arrange(IID) |> 
  tibble::column_to_rownames(var = "IID") |> 
  as.matrix() 

geno_matrix <- geno_mat[-which(rownames(geno_mat) == setdiff(rownames(geno_mat), geno_scores$IID)),] ## Result is row 9
geno_matrix <- geno_matrix[order(rownames(geno_matrix)),]
pheno_matrix <- geno

### Plotting PC1 scores on the map as dots on the map ---------------
library(sf)
uga_filtered <- readRDS("./data/UG_map_cropped_data.RDS")
pca_coord <- readRDS("./data/samples_coord_pca.RDS")[,c(1:5, seq(7,17,2))] |> 
  dplyr::rename(Longitude = longitude, Precipitation = rain_PC1, `Max temperature` = maxt_PC1, Radiation = radn_PC1,
                `Min temperature` = mint_PC1, `Relative humidity` = rh_PC1, `Wind speed` = wind_PC1)
##Check the correlation among PCs
pairs(pca_coord[,6:11])
pca_coord_long <- pca_coord |> pivot_longer(Radiation:`Wind speed`, names_to = "Variable", values_to = "PC1_scores")

pca_coord_long_scale <- pca_coord_long |> group_by(Variable) |> 
  dplyr::mutate(PC1_scaled = scales::rescale(PC1_scores))

### Radiation on its own
pca_map_rain <- ggplot() +
  geom_sf(data = uga_filtered, color = "#00000000", fill = "grey90") +
  geom_sf(data = uga_sel_1, aes(fill = factor(NAME_1, levels = loc_order)), color="#00000000", fill = "grey70") +
  geom_sf(data = uga_sel_2, aes(fill = factor(NAME_2, levels = loc_order)), color="#00000000", fill = "grey70") +
  # scale_fill_manual(values = col2) + new_scale_fill() +
  geom_sf(data = uga_sel_3, aes(fill = factor(NAME_3, levels = loc_order)), color="#00000000", fill = "grey70") +
  # scale_fill_manual(values = col3) + new_scale_fill() +
  geom_sf(data = uga_sel_4, aes(fill = factor(NAME_4, levels = loc_order)), color="#00000000", fill = "grey70") +
  scale_fill_manual(values = col4) + new_scale_fill() +
  theme_bw() + 
  scale_fill_manual(values = loc_colors)+ new_scale_fill() +
  geom_point(data = pca_coord, aes(x = Longitude, y = Latitude, color = Precipitation), size = 2, alpha = 0.5) +
  scale_color_gradient2(low = "#EBBE4D", high = "#0044CD", mid = "lightblue", midpoint = 0)+
  theme(axis.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 12, face = "bold"),
        legend.position = "right", legend.title = element_blank())+  
  scale_x_continuous(name = "Longitude", limits = c(xlow, xhigh),
                     breaks = seq(round(xlow,0), xhigh, 0.5),
                     labels = function(x) paste0(x, "°E"), expand = c(0, 0)) +
  scale_y_continuous(name = "Latitude", limits = c(ylow, yhigh),
                     labels = function(y) paste0(y, "°N"), expand = c(0, 0))

ggsave("./output/pca_scores_map_rain.jpg", width = 4.2, 
       height = 4, dpi = 600, units = "in", plot = pca_map_rain)

pca_map_maxt <- ggplot() +
  geom_sf(data = uga_filtered, color = "#00000000", fill = "grey90") +
  geom_sf(data = uga_sel_1, aes(fill = factor(NAME_1, levels = loc_order)), color="#00000000", fill = "grey70") +
  geom_sf(data = uga_sel_2, aes(fill = factor(NAME_2, levels = loc_order)), color="#00000000", fill = "grey70") +
  # scale_fill_manual(values = col2) + new_scale_fill() +
  geom_sf(data = uga_sel_3, aes(fill = factor(NAME_3, levels = loc_order)), color="#00000000", fill = "grey70") +
  # scale_fill_manual(values = col3) + new_scale_fill() +
  geom_sf(data = uga_sel_4, aes(fill = factor(NAME_4, levels = loc_order)), color="#00000000", fill = "grey70") +
  scale_fill_manual(values = col4) + new_scale_fill() +
  theme_bw() + 
  scale_fill_manual(values = loc_colors)+ new_scale_fill() +
  geom_point(data = pca_coord, aes(x = Longitude, y = Latitude, color = `Max temperature`), 
             size = 2, alpha = 0.5) +
  scale_color_gradient2(low = "green3", high = "red3", mid = "yellow2", midpoint = 0)+
  theme(axis.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 12, face = "bold"),
        legend.position = "right", legend.title = element_blank())+  
  scale_x_continuous(name = "Longitude", limits = c(xlow, xhigh),
                     breaks = seq(round(xlow,0), xhigh, 0.5),
                     labels = function(x) paste0(x, "°E"), expand = c(0, 0)) +
  scale_y_continuous(name = "Latitude", limits = c(ylow, yhigh),
                     labels = function(y) paste0(y, "°N"), expand = c(0, 0))

ggsave("./output/pca_scores_map_maxt.jpg", width = 4.5, 
       height = 4, dpi = 600, units = "in", plot = pca_map_maxt)

#### Radiation on its own
pca_map_radn <- ggplot() +
  geom_sf(data = uga_filtered, color = "#00000000", fill = "grey90") +
  geom_sf(data = uga_sel_1, aes(fill = factor(NAME_1, levels = loc_order)), color="#00000000", fill = "grey70") +
  geom_sf(data = uga_sel_2, aes(fill = factor(NAME_2, levels = loc_order)), color="#00000000", fill = "grey70") +
  # scale_fill_manual(values = col2) + new_scale_fill() +
  geom_sf(data = uga_sel_3, aes(fill = factor(NAME_3, levels = loc_order)), color="#00000000", fill = "grey70") +
  # scale_fill_manual(values = col3) + new_scale_fill() +
  geom_sf(data = uga_sel_4, aes(fill = factor(NAME_4, levels = loc_order)), color="#00000000", fill = "grey70") +
  scale_fill_manual(values = col4) + new_scale_fill() +
  theme_bw() + 
  scale_fill_manual(values = loc_colors)+ new_scale_fill() +
  geom_point(data = pca_coord, aes(x = Longitude, y = Latitude, color = Radiation), 
             size = 2, alpha = 0.5) +
  scale_color_gradient2(low = "green2", high = "red3", mid = "yellow2", midpoint = 0)+
  theme(axis.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 12, face = "bold"),
        legend.position = "right", legend.title = element_blank())+  
  scale_x_continuous(name = "Longitude", limits = c(xlow, xhigh),
                     breaks = seq(round(xlow,0), xhigh, 0.5),
                     labels = function(x) paste0(x, "°E"), expand = c(0, 0)) +
  scale_y_continuous(name = "Latitude", limits = c(ylow, yhigh),
                     labels = function(y) paste0(y, "°N"), expand = c(0, 0))

ggsave("./output/pca_scores_map_radn.jpg", width = 4.5, 
       height = 4, dpi = 600, units = "in", plot = pca_map_radn)

### Creating a faceted plot
pca_map_fct <- ggplot() +
  geom_sf(data = uga_filtered, color = "#00000000", fill = "grey90") +
  geom_sf(data = uga_sel_1, aes(fill = factor(NAME_1, levels = loc_order)), color="#00000000", fill = "grey70") +
  geom_sf(data = uga_sel_2, aes(fill = factor(NAME_2, levels = loc_order)), color="#00000000", fill = "grey70") +
  # scale_fill_manual(values = col2) + new_scale_fill() +
  geom_sf(data = uga_sel_3, aes(fill = factor(NAME_3, levels = loc_order)), color="#00000000", fill = "grey70") +
  # scale_fill_manual(values = col3) + new_scale_fill() +
  geom_sf(data = uga_sel_4, aes(fill = factor(NAME_4, levels = loc_order)), color="#00000000", fill = "grey70") +
  scale_fill_manual(values = col4) + new_scale_fill() +
  theme_bw() + scale_fill_manual(values = loc_colors)+ 
  new_scale_color() +
  geom_point(data = filter(pca_coord_long, Variable == "Radiation"), 
             aes(x = Longitude, y = Latitude, color = PC1_scores), size = 2, alpha = 0.5) +
  scale_color_gradient2(low = "green4", high = "red4", mid = "yellow3", midpoint = 0, name = "Radiation")+
  new_scale_color() +
  geom_point(data = filter(pca_coord_long, Variable == "Precipitation"), 
             aes(x = Longitude, y = Latitude, color = PC1_scores), size = 2, alpha = 0.5) +
  scale_color_gradient2(low = "blue4", high = "red4", mid = "white", midpoint = 0, name = "Precipitation") +
  new_scale_color() +
  geom_point(data = filter(pca_coord_long, Variable == "Min temperature"), 
             aes(x = Longitude, y = Latitude, color = PC1_scores), size = 2, alpha = 0.5) +
  scale_color_gradient2(low = "blue4", high = "red4", mid = "white", midpoint = 0, name = "Min temperature") +
  new_scale_color() +
  geom_point(data = filter(pca_coord_long, Variable == "Max temperature"), 
             aes(x = Longitude, y = Latitude, color = PC1_scores), size = 2, alpha = 0.5) +
  scale_color_gradient2(low = "blue4", high = "red4", mid = "white", midpoint = 0, name = "Max temperature") +
  new_scale_color() +
  geom_point(data = filter(pca_coord_long, Variable == "Relative humidity"), 
             aes(x = Longitude, y = Latitude, color = PC1_scores), size = 2, alpha = 0.5) +
  scale_color_gradient2(low = "blue4", high = "red4", mid = "white", midpoint = 0, name = "Relative humidity") +
  new_scale_color() +
  geom_point(data = filter(pca_coord_long, Variable == "Wind speed"), 
             aes(x = Longitude, y = Latitude, color = PC1_scores), size = 2, alpha = 0.5) +
  scale_color_gradient2(low = "blue4", high = "red4", mid = "white", midpoint = 0, name = "Wind speed") +
  facet_wrap(~Variable, ncol = 2, nrow = 3)+
  theme(axis.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 12, face = "bold"),
        legend.position = "right", panel.spacing.x = unit(1, "cm"))+  
  scale_x_continuous(name = "Longitude", limits = c(xlow, xhigh),
                     breaks = seq(round(xlow,0), xhigh, 0.5),
                     labels = function(x) paste0(x, "°E"), expand = c(0, 0)) +
  scale_y_continuous(name = "Latitude", limits = c(ylow, yhigh),
                     labels = function(y) paste0(y, "°N"), expand = c(0, 0))

ggsave("./output/pca_scores_map_facet.jpg", width = 8, 
       height = 10.5, dpi = 400, units = "in", plot = pca_map_fct)
### Creating a faceted plot
pca_map_fct_scaled <- ggplot() +
  geom_sf(data = uga_filtered, color = "#00000000", fill = "grey90") +
  geom_sf(data = uga_sel_1, aes(fill = factor(NAME_1, levels = loc_order)), color="#00000000", fill = "grey70") +
  geom_sf(data = uga_sel_2, aes(fill = factor(NAME_2, levels = loc_order)), color="#00000000", fill = "grey70") +
  # scale_fill_manual(values = col2) + new_scale_fill() +
  geom_sf(data = uga_sel_3, aes(fill = factor(NAME_3, levels = loc_order)), color="#00000000", fill = "grey70") +
  # scale_fill_manual(values = col3) + new_scale_fill() +
  geom_sf(data = uga_sel_4, aes(fill = factor(NAME_4, levels = loc_order)), color="#00000000", fill = "grey70") +
  scale_fill_manual(values = col4) + new_scale_fill() +
  theme_bw() + scale_fill_manual(values = loc_colors)+ 
  new_scale_color() +
  geom_point(data = pca_coord_long_scale, 
             aes(x = Longitude, y = Latitude, color = PC1_scaled), size = 2, alpha = 0.5) +
  scale_color_gradient2(low = "blue4", high = "red4", mid = "white", midpoint = 0.5, name = "Low vs. high") +
  facet_wrap(~Variable, ncol = 2, nrow = 3)+
  theme(axis.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 12, face = "bold"),
        legend.position = "right", panel.spacing.x = unit(1, "cm"))+  
  scale_x_continuous(name = "Longitude", limits = c(xlow, xhigh),
                     breaks = seq(round(xlow,0), xhigh, 0.5),
                     labels = function(x) paste0(x, "°E"), expand = c(0, 0)) +
  scale_y_continuous(name = "Latitude", limits = c(ylow, yhigh),
                     labels = function(y) paste0(y, "°N"), expand = c(0, 0))

ggsave("./output/pca_scores_map_facet_scaled.jpg", width = 8, 
       height = 9, dpi = 400, units = "in", plot = pca_map_fct_scaled)

## Finally Install the R package for Genome-environment association analyses ------------------
# install.packages("BiocManager")
# BiocManager::install("LEA")
##Reference for package LEA: https://besjournals.onlinelibrary.wiley.com/doi/epdf/10.1111/2041-210X.12382

library(LEA)
##Function lfmm (Bayesian) vs. lfmm2(Frequentist) approaches. The latter faster.
#input a genotypic matrix. The genotypic matrix must be in the lfmm{lfmm_format} format without missing values (9 or NA)
#env is a matrix of environmental covariates
#K integer = number of latent factors. could be estimated from the elbow point in the PCA screeplot for the genotype matrix.
modl <- LEA::lfmm2(input = geno_matrix, env = geno_scores, K=3, effect.sizes = T)

pv <- lfmm2.test(object = modl, input = geno_matrix, env = geno_scores, linear = TRUE)

cor(t(pv$pvalues))
jpeg("./output/genome_environment_association.jpg", width = 9, height = 12, res = 300, units = "in")
par(mfrow = c(6,1), mar = c(1,4,1,2))
plot(-log10(pv$pvalues[1,]), col = "grey50", cex = 1, pch = 19, ylim = c(0,15), ylab = "-log10(p.value)")
text(500, 14, expression(bold("Radiation")))
abline(h = 6, lty = "dashed", col = "darkred")
plot(-log10(pv$pvalues[2,]), col = "grey50", cex = 1, pch = 19, ylim = c(0,15), ylab = "-log10(p.value)")
text(500, 14, expression(bold("Precipitation")))
abline(h = 6, lty = "dashed", col = "darkred")
plot(-log10(pv$pvalues[3,]), col = "grey50", cex = 1, pch = 19, ylim = c(0,15), ylab = "-log10(p.value)")
text(1000, 14, expression(bold("Maximum temperature")))
abline(h = 6, lty = "dashed", col = "darkred")
plot(-log10(pv$pvalues[4,]), col = "grey50", cex = 1, pch = 19, ylim = c(0,15), ylab = "-log10(p.value)")
text(1000, 14, expression(bold("Minimum temperature")))
abline(h = 6, lty = "dashed", col = "darkred")
plot(-log10(pv$pvalues[5,]), col = "grey50", cex = 1, pch = 19, ylim = c(0,15), ylab = "-log10(p.value)")
text(500, 14, expression(bold("Relative humidity")))
abline(h = 6, lty = "dashed", col = "darkred")
plot(-log10(pv$pvalues[6,]), col = "grey50", cex = 1, pch = 19, ylim = c(0,15), ylab = "-log10(p.value)")
text(500, 14, expression(bold("Wind speed")))
abline(h = 6, lty = "dashed", col = "darkred")
dev.off()
# plot(-log10(pv$pvalues[2,]), col = "grey", cex = .8, pch = 19)

#### Format results' pvalues in a df and save significant snps
pv_dt <- t(pv$pvalues) |> as.data.frame() |> 
  tibble::rownames_to_column(var = "snp") |> 
  rename(Radiation = radnPC1, Precipitation = rainPC1,
          Max_temprature = maxtPC1, Min_temperature = mintPC1,
          Relative_humidity = rhPC1, Wind_speed = wind_PC1) |> 
  pivot_longer(cols = Radiation:Wind_speed, names_to = "Variable", values_to = "log10pvalue") |> 
  dplyr::mutate(log10pvalue = -log10(log10pvalue), snp = gsub("Response ", "", snp))
  
pv_dt_significant <- pv_dt |> 
  dplyr::filter(log10pvalue >= 6) |> 
  left_join(bim_final)

write.csv(pv_dt_significant, "./data/significant_snps_environ_gwas.csv", row.names = F)


### Save inputs and outputs objects a RDS file
genome_env_assoc <- list(geno = geno_matrix, pheno = geno_scores,
                         model = modl, p_values = pv)
saveRDS(genome_env_assoc, "./data/genome_environment_association.RDS")

################# ----------- Creating a customized Manhattan plots --------- ###########
genome_env_assoc <- readRDS("./data/genome_environment_association.RDS")
list2env(genome_env_assoc, envir = .GlobalEnv)
bim_final <- readRDS("./data/genetic_map_final_whole_set.RDS")

pv_dt <- t(p_values$pvalues) |> as.data.frame() |> 
  tibble::rownames_to_column(var = "snp") |> 
  rename(Radiation = radnPC1, Precipitation = rainPC1,
         Max_temprature = maxtPC1, Min_temperature = mintPC1,
         Relative_humidity = rhPC1, Wind_speed = wind_PC1) |> 
  pivot_longer(cols = Radiation:Wind_speed, names_to = "Variable", values_to = "log10pvalue") |> 
  dplyr::mutate(log10pvalue = -log10(log10pvalue), snp = gsub("Response ", "", snp)) |> 
  left_join(bim_final)

pv_dt2 <- left_join(pv_dt, dt[,c(3,5,6)], by = c("snp" = "id")) |> 
  tidyr::pivot_wider(names_from = Variable, values_from = log10pvalue) |> 
  arrange(pos_cum) |> 
  dplyr::mutate(chr_num = ifelse(is.na(chr_num), "Other", chr_num), pos_cum = 1:length(chr)) |> 
  tidyr::pivot_longer(cols = Radiation:Wind_speed, names_to = "Variable", values_to = "log10pvalue")

pv_dt_sign <- filter(pv_dt2, log10pvalue >= 6)
pv_dt_backgrd <- filter(pv_dt2, log10pvalue < 6 & Variable == "Radiation")
centers <- pv_dt2 |> dplyr::group_by(chr_num) |>
  summarise_at("pos_cum", mean)
ymax = 18; threshold = 6

manhplot_env <- ggplot(pv_dt_backgrd, aes(x = pos_cum, y = log10pvalue)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,ymin = threshold, ymax = Inf), 
            fill = "#F0DAFF", alpha = 0.2) +
  geom_hline(yintercept = threshold, color = "darkred", linetype = "dashed") +
  geom_point(aes(color = factor(chr_num)), shape = 19) +
  scale_color_manual(values = rep(c("#274AB3", "#7CA4E1"), 
                                  length(unique(dt$chr_num))), guide = "none") +
  
  ggnewscale::new_scale_color() +
  
  geom_point(data = pv_dt_sign, aes(x = pos_cum, y = log10pvalue, color = Variable), size = 3) + 
  scale_color_brewer(palette = "Set1", name = "Variable") +
  # scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2"), name = "Variable") +
  scale_x_continuous(labels = centers$chr_num, breaks = centers$pos_cum) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
  labs(x = "Contig", y = expression(bold(-log[10](p)))) +
  boris_theme_rect +
  theme(legend.position = c(0.1, 0.8), legend.direction = "vertical",
        axis.text.x = element_text(angle = 90, vjust = 0.2, size = 10))
manhplot_env
ggsave("./output/Manhattan_plot_Env_GWAS_multivar.jpg", width = 12, height = 3.5, 
       plot = manhplot_env, dpi = 400, units = "in")

#### OR run below if the above is too slow
jpeg("./output/Manhattan_plot_Env_GWAS_multivar.jpg", width = 12, height = 3.5, 
     res = 400, units = "in")
ggplot(pv_dt_backgrd, aes(x = pos_cum, y = log10pvalue)) +
  geom_rect(aes(xmin = -Inf, xmax = Inf,ymin = threshold, ymax = Inf), 
            fill = "#F0DAFF", alpha = 0.2) +
  geom_hline(yintercept = threshold, color = "darkred", linetype = "dashed") +
  geom_point(aes(color = factor(chr_num)), shape = 19) +
  scale_color_manual(values = rep(c("#274AB3", "#7CA4E1"), 
                                  length(unique(dt$chr_num))), guide = "none") +
  
  ggnewscale::new_scale_color() +
  
  geom_point(data = pv_dt_sign, aes(x = pos_cum, y = log10pvalue, color = Variable), size = 3) + 
  scale_color_brewer(palette = "Set1", name = "Variable") +
  # scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2"), name = "Variable") +
  scale_x_continuous(labels = centers$chr_num, breaks = centers$pos_cum) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, ymax)) +
  labs(x = "Contig", y = expression(bold(-log[10](p)))) +
  boris_theme_rect +
  theme(legend.position = c(0.1, 0.7), legend.direction = "vertical",
        axis.text.x = element_text(angle = 90, vjust = 0.2, size = 8))
dev.off()

######### Creating Q-Q plot
genome_env_assoc <- readRDS("./data/genome_environment_association.RDS")
list2env(genome_env_assoc, envir = .GlobalEnv)
bim_final <- readRDS("./data/genetic_map_final_whole_set.RDS")

pv_dt <- t(p_values$pvalues) |> as.data.frame() |> 
  tibble::rownames_to_column(var = "snp") |> 
  rename(Radiation = radnPC1, Precipitation = rainPC1,
         Max_temprature = maxtPC1, Min_temperature = mintPC1,
         Relative_humidity = rhPC1, Wind_speed = wind_PC1) |> 
  pivot_longer(cols = Radiation:Wind_speed, names_to = "Variable", values_to = "log10pvalue") |> 
  dplyr::mutate(log10pvalue = -log10(log10pvalue), snp = gsub("Response ", "", snp)) |> 
  dplyr::arrange(Variable, log10pvalue) |> 
  dplyr::mutate(exptp = rep(sort(-log10(seq(6/nrow(pv_dt), 1, 6/nrow(pv_dt)))), 6),
                Variable = recode(Variable, "Max_temprature" = "Max temprature", "Min_temprature" = "Min temprature",
                                  "Precipitation" = "Precipitation", "Radiation" = "Radiation", 
                                  "Relative_humidity" = "Relative humidity", "Wind_speed" = "Wind speed"))


jpeg("./output/QQplot_GEA_whole_data.jpg", width = 4, height = 4, res = 400, units = "in")
ggplot(pv_dt, aes(x = exptp, y = log10pvalue, color = Variable)) +
  geom_point(alpha = 0.5, shape = 19, size = 2.5) + boris_theme_rect +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = expression(bold(Expected -log[10](p))), 
       y = expression(bold(Observed -log[10](p)))) +
  scale_color_brewer(palette = "Set1", name = "Variable") +
  # scale_color_manual(values = RColorBrewer::brewer.pal(8, "Dark2")[c(4,1,2,3)], name = "Trait") +
  theme(legend.position = c(0.3, 0.75), legend.direction = "vertical")
dev.off()

### Plotting snps distribution on the map as dots on the map ---------------
library(sf)
uga_filtered <- readRDS("./data/UG_map_cropped_data.RDS")
snp_signif_gea <- read.csv("./data/significant_snps_environ_gwas.csv") |> 
  dplyr::arrange(desc(log10pvalue))
geno_scores <- geno |> as.data.frame() |> 
  tibble::rownames_to_column(var = "IID") |> 
  dplyr::select_at(c("IID", snp_signif_gea$snp))
geno_loc <- readRDS("./data/samples_coord_pca.RDS")[,1:5] |> 
  left_join(geno_scores)

geno_loc_long <- geno_loc |> pivot_longer(SNP10073:SNP13901, names_to = "snp", values_to = "genotype") |> 
  group_by(Location, snp, genotype) |> summarise(count = n())

library(scatterpie) ### Needed to plot pies
loc_admix_coord <- read.csv("./data/final_data_scatterpie.csv")

loc_admix_coord_snp <- left_join(loc_admix_coord, geno_loc_long) |> 
  dplyr::mutate(genotype = factor(paste0("G", genotype))) |> 
  pivot_wider(names_from = genotype, values_from = count) |> 
  dplyr::select(-c(Subpopulation.1:Subpopulation.3)) |> 
  dplyr::mutate_all(~replace_na(., 0))

# pca_coord_long <- pca_coord |> 
#   pivot_longer(Radiation:`Wind speed`, names_to = "Variable", values_to = "PC1_scores")
# 
# pca_coord_long_scale <- pca_coord_long |> group_by(Variable) |> 
#   dplyr::mutate(PC1_scaled = scales::rescale(PC1_scores))

cur_snp <- "SNP13901"
snp_map <- ggplot() +
  geom_sf(data = uga_filtered, color = "#00000000", fill = "grey90") +
  geom_sf(data = uga_sel_1, aes(fill = factor(NAME_1, levels = loc_order)), color="#00000000", fill = "grey70") +
  geom_sf(data = uga_sel_2, aes(fill = factor(NAME_2, levels = loc_order)), color="#00000000", fill = "grey70") +
  geom_sf(data = uga_sel_3, aes(fill = factor(NAME_3, levels = loc_order)), color="#00000000", fill = "grey70") +
  geom_sf(data = uga_sel_4, aes(fill = factor(NAME_4, levels = loc_order)), color="#00000000", fill = "grey70") +
  theme_bw() + 
  ## Adding the pies
  geom_scatterpie(data = filter(loc_admix_coord_snp, snp == cur_snp),
                  aes(x = Longitude, y = Latitude, r = sqrt(total)/50),
                  cols = c("G0", "G1", "G2")) +
  # scale_fill_manual(values = friendly_pal("bright_seven")[c(2,3,1)]) +
  annotate("text", x = xhigh-0.3, y = yhigh-0.1, label = cur_snp)+
  theme_bw() +
  theme(axis.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 12, face = "bold"),
        legend.position = "none") +  
  scale_x_continuous(name = "Longitude", limits = c(xlow, xhigh),
                     breaks = seq(round(xlow,0), xhigh, 0.5),
                     labels = function(x) paste0(x, "°E"), expand = c(0, 0)) +
  scale_y_continuous(name = "Latitude", limits = c(ylow, yhigh),
                     labels = function(y) paste0(y, "°N"), expand = c(0, 0))
ggsave(paste0("./output/pca_scores_map_radn_", cur_snp,"_.jpg"), width = 3, 
       height = 3, dpi = 600, units = "in", plot = snp_map)
### Rain on its own
pca_map_rain <- ggplot() +
  geom_sf(data = uga_filtered, color = "#00000000", fill = "grey90") +
  geom_sf(data = uga_sel_1, aes(fill = factor(NAME_1, levels = loc_order)), color="#00000000", fill = "grey70") +
  geom_sf(data = uga_sel_2, aes(fill = factor(NAME_2, levels = loc_order)), color="#00000000", fill = "grey70") +
  geom_sf(data = uga_sel_3, aes(fill = factor(NAME_3, levels = loc_order)), color="#00000000", fill = "grey70") +
  geom_sf(data = uga_sel_4, aes(fill = factor(NAME_4, levels = loc_order)), color="#00000000", fill = "grey70") +
  theme_bw() + 
  geom_point(data = pca_coord, aes(x = Longitude, y = Latitude, color = Precipitation), size = 2, alpha = 0.5) +
  scale_color_gradient2(low = "#EBBE4D", high = "#0044CD", mid = "lightblue", midpoint = 0)+
  theme(axis.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 12, face = "bold"),
        legend.position = "right", legend.title = element_blank())+  
  scale_x_continuous(name = "Longitude", limits = c(xlow, xhigh),
                     breaks = seq(round(xlow,0), xhigh, 0.5),
                     labels = function(x) paste0(x, "°E"), expand = c(0, 0)) +
  scale_y_continuous(name = "Latitude", limits = c(ylow, yhigh),
                     labels = function(y) paste0(y, "°N"), expand = c(0, 0))

ggsave("./output/pca_scores_map_rain.jpg", width = 4.2, 
       height = 4, dpi = 600, units = "in", plot = pca_map_rain)

### Creating a faceted plot
snp_map_fct <- ggplot() +
  geom_sf(data = uga_filtered, color = "#00000000", fill = "grey90") +
  geom_sf(data = uga_sel_1, aes(fill = factor(NAME_1, levels = loc_order)), color="#00000000", fill = "grey70") +
  geom_sf(data = uga_sel_2, aes(fill = factor(NAME_2, levels = loc_order)), color="#00000000", fill = "grey70") +
  # scale_fill_manual(values = col2) + new_scale_fill() +
  geom_sf(data = uga_sel_3, aes(fill = factor(NAME_3, levels = loc_order)), color="#00000000", fill = "grey70") +
  # scale_fill_manual(values = col3) + new_scale_fill() +
  geom_sf(data = uga_sel_4, aes(fill = factor(NAME_4, levels = loc_order)), color="#00000000", fill = "grey70") +
  scale_fill_manual(values = col4) + new_scale_fill() +
  theme_bw() + scale_fill_manual(values = loc_colors)+ 
  new_scale_color() +
  geom_point(data = pca_coord_long_scale, 
             aes(x = Longitude, y = Latitude, color = PC1_scaled), size = 2, alpha = 0.5) +
  scale_color_gradient2(low = "blue4", high = "red4", mid = "white", midpoint = 0.5, name = "Low vs. high") +
  facet_wrap(~Variable, ncol = 2, nrow = 3)+
  theme(axis.text = element_text(color = 'black', size = 12),
        axis.title = element_text(color = 'black', size = 12, face = "bold"),
        legend.position = "right", panel.spacing.x = unit(1, "cm"))+  
  scale_x_continuous(name = "Longitude", limits = c(xlow, xhigh),
                     breaks = seq(round(xlow,0), xhigh, 0.5),
                     labels = function(x) paste0(x, "°E"), expand = c(0, 0)) +
  scale_y_continuous(name = "Latitude", limits = c(ylow, yhigh),
                     labels = function(y) paste0(y, "°N"), expand = c(0, 0))

ggsave("./output/pca_scores_map_facet_scaled.jpg", width = 8, 
       height = 9, dpi = 400, units = "in", plot = pca_map_fct_scaled)
