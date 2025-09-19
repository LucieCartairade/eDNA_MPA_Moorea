These R codes show how to plot figures from this paper : 
"TITLE" doi : 

## Figure 1
<p align="center">
  <img src="Figures/Figure1.png" alt="Figure 1" class="center" width="50%"/>
</p>

Shapes files available on : https://www.tefenua.data.gov.pf
https://www.tefenua.data.gov.pf/datasets/ef2bdc8e55f049318a3888f8134349b0_0/explore?location=-17.535984%2C-149.840113%2C11.48
```r
library(sf)

### Reads Files
my_sf <- read_sf(paste0(Datas_path,"Moorea/WGS84_TraitDeCote.shp"), options = "ENCODING=UTF-8")
sf_Moorea <- my_sf %>% filter(ILE == "MOOREA")

sf_Moorea <- sf_Moorea %>%
  mutate(Carto = recode(Carto,
                        "Terre"   = "Land",
                        "Récif"   = "Reef",
                        "Lagon"   = "Lagoon",
                        "Tombant" = "Outer reef"
  ))

sf_Moorea$Carto <- factor(sf_Moorea$Carto, levels = c("Land", "Reef", "Lagoon", "Outer reef"))

sf_AMP <- read_sf(paste0(Datas_path,"PGEM_AMP/PGEM_Aires_marines_prot%C3%A9g%C3%A9es_(TeFenua).shp"))

st_crs(sf_AMP)
st_crs(sf_Moorea)

### Color palettes
palette <- c(
  Land = "gray80",
  Reef = "lightblue",
  Lagoon = "turquoise",
  `Outer reef` = "skyblue", 
  `Marine Protected Area` ="red"
)
palette2 <- c(
  Control = "blue", 
  MPA = "black"
)

### Coordinates
coord <- read.csv(file = "MPA_GPS-points.csv", header = T, sep = ";", dec = ",", na.strings = "NA", fileEncoding = "ISO-8859-1")
coord <- as.data.frame(coord)
coord <- coord %>% filter(Nom.de.l.Aire.Marine == "AHI" | Nom.de.l.Aire.Marine == "AFAREAITU"  | Nom.de.l.Aire.Marine == "ENTRE 2 BAIES"  | 
                            Nom.de.l.Aire.Marine == "HAAPITI" | Nom.de.l.Aire.Marine == "PIHAENA"  | Nom.de.l.Aire.Marine == "TAOTAHA" )
coord$Latitude <-paste0("-",coord$Latitude)
coord$Longitude <-paste0("-",coord$Longitude)
coord$Type.d.Aire.Marine <- recode(coord$Type.d.Aire.Marine, "AMT" = "Control", "AMP" = "MPA")
coord <- coord[!is.na(coord$Nom),]

coord$Latitude_d <- parzer::parse_lat(coord$Latitude)
coord$Longitude_d <- parzer::parse_lon(coord$Longitude)

coord <- st_as_sf(coord, coords = c("Longitude_d", "Latitude_d"),  crs = 4326)

#### Map
ggplot() +
  geom_sf(data = sf_Moorea, aes(fill = Carto), color = "black") +  
  scale_fill_manual(values = palette, name = "Zones") + 
  geom_sf(data = sf_AMP, aes(fill = "Marine Protected Area"), alpha = 0.4) +
  geom_sf(data = coord, aes(color = Type.d.Aire.Marine), size = 3, shape = 18) + 
  scale_color_manual(values = palette2, name = "Status") + 
  theme_minimal() + 
  scale_x_continuous(
    #name = "Longitude",
    labels = "149°50′W",
    breaks = -149.83333
  ) +
  scale_y_continuous(
    #name = "Latitude",
    labels = "17°32′S",
    breaks = -17.5333312
  ) + ggspatial::annotation_scale(location = 'tl')

ggsave(path = Images_path, file = "Figure1.pdf", height = 7, width = 7)  
```
## Figure 2
<p align="center">
  <img src="Figures/Figure2.png" alt="Figure 2" class="center" width="50%"/>
</p>

```r
pdf(file = paste0(Images_path, "Figure2.pdf"), width = 10, height = 5)
set.seed(10)

Tab_Euler_1 <- reshape2::acast(subset(Tax_melt_wVC, Sampling.Site != "Control" & Replica != "B" & Family != "unknown"), value.var = "relative_biomass", Taxon ~ Sample.Type, fill = 0, fun.aggregate = sum)
Tab_Euler_1 <- ifelse(Tab_Euler_1 == 0, FALSE, TRUE)
p1 <- plot(eulerr::euler(Tab_Euler_1, shape = "ellipse"), quantities = TRUE, main = "A. Species", return_grob = TRUE)

Tab_Euler_2 <- reshape2::acast(subset(Tax_melt_wVC, Sampling.Site != "Control" & Replica != "B" & Family != "unknown"), value.var = "relative_biomass", Family ~ Sample.Type, fill = 0, fun.aggregate = sum)
Tab_Euler_2 <- ifelse(Tab_Euler_2 == 0, FALSE, TRUE)
p2 <- plot(eulerr::euler(Tab_Euler_2, shape = "ellipse"), quantities = TRUE, main = "B. Family", return_grob = TRUE)

gridExtra::grid.arrange(p1, p2, ncol = 2, padding = unit(1, "cm"), top = grid::textGrob(""), bottom = grid::textGrob(""), left = grid::textGrob(""), right = grid::textGrob(""))

dev.off()
```

# Figure 3
<p align="center">
  <img src="Figures/Figure3.png" alt="Figure 3" class="center" width="50%"/>
</p>

```r
df_count <- Tax_melt_wVC %>%
  filter(!is.na(Family), !is.na(Taxon), Taxon != "unknown") %>%
  distinct(Sample.Type, Family, Taxon) %>%
  group_by(Sample.Type, Family) %>%
  summarise(count = n(), .groups = "drop")

family_order <- df_count %>%
  group_by(Family) %>%
  summarise(
    sum_count = sum(count, na.rm = TRUE),
    both_methods = as.integer(n_distinct(Sample.Type[count > 0]) == 2),
    .groups = "drop"
  ) %>%
  arrange(sum_count, desc(both_methods)) %>% 
  pull(Family)

df_count <- df_count %>%
  tidyr::complete(Sample.Type, Family, fill = list(count = 0)) %>%
  mutate(Family = factor(Family, levels = family_order),
         Sample.Type = factor(Sample.Type, levels = c("eDNA", "UVC")),
         count_plot = ifelse(Sample.Type == "eDNA", -count, count)) 

ggplot(df_count, aes(x = Family, y = count_plot, fill = Sample.Type)) +
  geom_bar(stat = "identity", width = 0.5) +
  scale_fill_manual(values = c("eDNA"= "#8f226e", "UVC" = "#f18055")) +
  scale_y_continuous(labels = abs) +
  labs(x = "Family", y = "Species count", fill = "Method") +
  coord_flip() +
  theme(
    axis.text.y = element_text(size = 10),
    axis.text.x = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "top",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.4)
  )

ggsave(path = Images_path, filename = "Figure3.pdf", width = 10, height = 7)
```
## Figure 4
```r
df_area <- Tax_melt_wVC %>%
  filter(!is.na(Family), !is.na(relative_biomass)) %>%
  group_by(Sample.Type, Family) %>%
  summarise(total = sum(relative_biomass), .groups = "drop") %>%
  group_by(Sample.Type) %>%
  mutate(prop = total / sum(total)) %>%
  ungroup()

df_complete <- df_area %>%
  tidyr::complete(Sample.Type, Family, fill = list(total = 0, prop = 0))
df_complete

family_order <- c(unlist(rev(unique(df_area[order(df_area$prop, decreasing = F), "Family"]))))
df_area$Family <- factor(df_area$Family, levels = unique(c("unknown",family_order)))

family_order <- c(unlist(rev(unique(df_complete[order(df_complete$prop, decreasing = F), "Family"]))))
df_complete$Family <- factor(df_complete$Family, levels = unique(c("unknown",family_order)))

set.seed(100)
getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Paired"))
my_Palette = getPalette(length(unique(Tax_melt_wVC$Family)))
names(my_Palette) = unique(Tax_melt_wVC$Family)
my_Palette["Hominidae"] <- "#444444"
my_Palette["unknown"] <- "#999999"

# Plot
ggplot(df_complete, aes(x = Sample.Type, y = prop, fill = Family, group = Family)) +
  geom_area(stat = "identity", position = "stack") +
  scale_x_discrete(limits = c("eDNA", "UVC")) +
  scale_fill_manual(values= my_Palette) + 
  labs(x = "Method",y = "Family relative biomass") + 
  guides(fill=guide_legend(ncol=3))
  
ggsave(path = Images_path, file = "Figure4.pdf", width = 6.5, height = 5)
```
## Figure 5 and 6
<p align="center">
  <img src="Figures/Figure5.png" alt="Figure 5" class="center" width="50%"/>
</p>

<p align="center">
  <img src="Figures/Figure6.png" alt="Figure 6" class="center" width="50%"/>
</p>

```r
my_CAP <- function(dist, n_arrows) {
  cap <- vegan::capscale(dist  ~ Sample.Type + Marine.Area + Sampling.Site + Sampling.Site:Habitat,
                         data = subset(sample, Sampling.Site != "Control" & Replica != "B"))
  cap
  summary(cap)
  anova(cap, permutations = 999, by = "terms")
  
  ef <- vegan::envfit(cap, vegan::decostand(t(Tax_table), method = "hellinger"), perm = 999)
  
  eig <- cap$CCA$eig
  eig_perc <- eig / sum(eig) * 100
  eig_perc[1:2]
  
  vec <- ef$vectors$arrows
  pvals <- ef$vectors$pvals
  
  vec_scaled <- vec * sqrt(ef$vectors$r)
  arrow_lengths <- sqrt(rowSums(vec_scaled^2))
  
  pvals <- pvals[pvals < 0.05]
  arrow_lengths <- arrow_lengths[names(arrow_lengths) %in% names(pvals)]
  top_idx <- names(sort(arrow_lengths, decreasing = TRUE))[1:15]
  
  ef_subset <- ef
  ef_subset$vectors$arrows <- ef$vectors$arrows[top_idx, , drop = FALSE]
  ef_subset$vectors$r <- ef$vectors$r[top_idx]
  ef_subset$vectors$pvals <- ef$vectors$pvals[top_idx]
  
  site_df <- as.data.frame(vegan::scores(cap, display = "sites"))
  site_df <- merge(site_df, subset(sample, Sampling.Site != "Control" & Replica != "B"), by = "row.names")
  rownames(site_df) <- site_df$Row.names; site_df$Row.names <- NULL
  
  arrows <- as.data.frame(ef$vectors$arrows)
  arrows$species <- rownames(arrows)
  arrows$length <- ef$vectors$r
  arrows$pval <- ef$vectors$pvals
  
  arrows <- subset(arrows, pval < 0.05)
  arrows <- arrows[order(arrows$length, decreasing = TRUE),]
  arrows_select <- arrows[1:n_arrows,]
  
  
  p <- ggplot(site_df, aes(x = CAP1, y = CAP2)) +
    stat_ellipse(aes(color = Sample.Type, linetype = Habitat),
                 type = "norm", level = 0.92, size = 0.75) + 
    scale_color_viridis_d(option = "magma", begin = 0.45, end = 0.75) +
    ggnewscale::new_scale_color() +
    geom_point(aes(color = Coast, shape = Marine.Area), size = 2)  + 
    scale_color_viridis_d(option = "viridis", begin = 0.4) +
    geom_segment(data = arrows_select, aes(x = 0, y = 0, xend = CAP1, yend = CAP2), 
                 arrow = arrow(length = unit(0.2, "cm")), color = "gray", linewidth = 0.5) +
    ggrepel::geom_text_repel(data = arrows_select,
                             aes(x = CAP1 * 1.1, y = CAP2 * 1.1, 
                                 label = paste0("italic('", species, "')")), 
                             color = "gray50", size = 3, 
                             max.overlaps = 1000, segment.color = NA,
                             parse = TRUE) + 
    scale_linetype_manual(values = c(3, 2, 1)) + 
    guides(shape = guide_legend(order = 1)) + 
    labs(
      x = paste0("CAP1 (", round(eig_perc[1], 1), "%)"),
      y = paste0("CAP2 (", round(eig_perc[2], 1), "%)")
    )

  return(p)
}

p <- my_CAP(dist = dist.jc.both$beta.jac, n_arrows = 15)
p
ggsave(path = paste0(Images_path, "CAP/"), file = "CAP_Jaccard_withSpecies_top15.pdf", height = 6, width = 7 )

p10 <- my_CAP(dist = dist.jc.both$beta.jac, n_arrows = 10) ; p10
ggsave(path = paste0(Images_path, "CAP/"), file = "CAP_Jaccard_withSpecies_top10.pdf", height = 6, width = 7 )

p20 <- my_CAP(dist = dist.jc.both$beta.jac, n_arrows = 20) ; p20
ggsave(path = paste0(Images_path, "CAP/"), file = "CAP_Jaccard_withSpecies_top20.pdf", height = 6, width = 7 )


p <- my_CAP(dist = dist.bc.both, n_arrows = 15)
p
ggsave(path = paste0(Images_path, "CAP/"), file = "CAP_BrayCurtis_withSpecies_top15.pdf", height = 6, width = 7 )
```
