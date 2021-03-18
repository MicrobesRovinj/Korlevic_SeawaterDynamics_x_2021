#################################################################################################################
# plot_pcoa_dbrda.R
# 
# A script to generate the PCoA and db-RDA figures.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared
#               data/raw/environmental_data.csv
#               data/raw/metadata.csv
# Produces: results/figures/pcoa_dbrda_figure.jpg
#
#################################################################################################################

#######################
# PCoA
#######################

# Loading OTU/sample data
shared <- read_tsv("data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared")

# Generating the random rarefied community data
rarefied <- shared %>%
  select(-label, -Group, -numOtus) %>%
  rrarefy(., min(rowSums(.))) %>%
  as_tibble(.name_repair="unique") %>%
  add_column("Group"=shared$Group, .before=TRUE) %>%
  select_if(list(~!is.numeric(.) || sum(.)!=0))

# Loading metadata 
metadata <- read_tsv("data/raw/metadata.csv")

# Joining metadata and OTU/sample data
rarefied_metadata <- inner_join(metadata, rarefied, by=c("ID"="Group"))

# Function for setting the number of decimal places
scaleFUN <- function(x) sprintf("%.2f", x)

# Generating a common theme for plots
theme <- theme(text=element_text(family="Times"), line=element_line(color="black"),
               panel.border=element_rect(fill=NA), panel.background=element_blank(),
               panel.grid=element_blank(), axis.line=element_blank(),
               axis.text=element_text(size=14, color="black"), axis.title=element_text(size=18, color="black"),
               plot.margin=unit(c(5.5, 16.5, 16.5, 16.5), "pt"),
               legend.text=element_text(size=18, margin=margin(r=0.2, unit="cm")), 
               legend.key.width=unit(0, "cm"), legend.key.height=unit(1.0, "cm"),
               legend.key=element_rect(fill="white"), legend.text.align=0, 
               legend.margin=margin(-5, 0, 0, 0), legend.justification=c("left", "top"),
               plot.title=element_text(size=14, hjust=0.5))

# Generating PCoA data
spe.bray <- vegdist(select(rarefied_metadata, starts_with("Otu")))
spe.b.pcoa <- cmdscale(spe.bray, k=nrow(rarefied_metadata)-1, eig=TRUE, add=FALSE)
coordinates <- spe.b.pcoa$points[,1:2] %>%
  as_tibble(.name_repair= ~c("A1", "A2")) %>%
  add_column("Group"=rarefied_metadata$label, .before=TRUE)
coordinates <- inner_join(metadata, coordinates, by=c("label"="Group"))

# PCoA plot generation
p.pcoa <- ggplot() +
  geom_point(data=coordinates, aes(x=A1, y=A2, fill=season, shape=station), size=5, stroke=0.5) +
  scale_fill_manual(name=NULL,
                    breaks=c("Spring", "Summer", "Autumn", "Winter"),
                    values=c("Spring"="#377EB8", "Summer"="#4DAF4A",
                             "Autumn"="#E41A1C", "Winter"="#984EA3")) +
  scale_shape_manual(name=NULL,
                     values=c(21, 23),
                     breaks=c("F", "S"),
                     labels=c("Funtana", "Saline")) +
  labs(x=paste0("PCoA I (",format(round(spe.b.pcoa$eig[1]/sum(spe.b.pcoa$eig)*100, digits=2), nsmall=2)," %)"), 
       y=paste0("PCoA II (",format(round(spe.b.pcoa$eig[2]/sum(spe.b.pcoa$eig)*100, digits=2), nsmall=2)," %)")) +
  scale_x_continuous(labels=scaleFUN) +
  scale_y_continuous(labels=scaleFUN) +
  theme +
  theme(legend.spacing.y=unit(0.1, "cm"),
        legend.position="", legend.justification=c("left", "bottom")) +
  guides(fill=guide_legend(override.aes=list(shape=21)))

#######################
# db-RDA
#######################

# The DNA sample originating from Funtana sampled on 24 April 2018 (ID, 23) was sequences twice.
# Sequences from both samples are summed.
shared <- shared %>%
  mutate(Group=str_replace(Group, "23_.", "23")) %>%
  group_by(Group) %>%
  summarise(across(starts_with("Otu"), list(sum), .names="{.col}"), .groups="drop")

# Generating the random rarefied community data
rarefied <- shared %>%
  select(-Group) %>%
  rrarefy(., min(rowSums(.))) %>%
  as_tibble(.name_repair="unique") %>%
  add_column("Group"=shared$Group, .before=TRUE) %>%
  select_if(list(~!is.numeric(.) || sum(.)!=0))

# Loading metadata 
metadata <- read_tsv("data/raw/metadata.csv") %>%
  mutate(ID=str_replace(ID, "23_1", "23"))

# Joining metadata and OTU/sample data
rarefied_metadata <- inner_join(metadata, rarefied, by=c("ID"="Group"))

# Loading environmental variables data
env <- read_tsv("data/raw/environmental_data.csv") %>%
  mutate(ID=as.character(ID)) %>%
  filter(!is.na(ID)) %>%
  mutate(across(!ID, list(as.numeric), .names="{.col}"))

# Joining metadata and environmental variables data
env_metadata <- inner_join(metadata, env, by=c("ID"="ID"))

# Generating db-RDA data
spe.bray <- vegdist(select(rarefied_metadata, starts_with("Otu")), method="bray")
spe.b.pcoa <- cmdscale(spe.bray, k=nrow(select(rarefied_metadata, starts_with("Otu")))-1, eig=TRUE,
                       add=TRUE)
spe.scores <- spe.b.pcoa$points
(spe.dbrda <- rda(spe.scores ~ `T` + S + PO4 + NH4 + NO2 + NO3 + SiO4 + PM + Chla + PA, env_metadata))
(R2adj <- RsquareAdj(spe.dbrda)$adj.r.squared)

# Global test of the db-RDA result
(anova.global <- anova.cca(spe.dbrda, permuatations=how(nperm=999)))

# Test of all canonical axes
(anova.axes <- anova(spe.dbrda, by="axis", permutations=how(nperm=999)))

# Variance inflation factors (VIF)
vif.cca(spe.dbrda)

# Combining and saving statistic results
statistic.dbrda <- list(R2adj=R2adj, anova.global=anova.global, anova.axes=anova.axes)
save(statistic.dbrda, file="results/numerical/statistic.dbrda.Rdata")

# Extracting data to generate the db-RDA plot
(summary.dbrda <- summary(spe.dbrda, scaling=1))
coordinates_sites <- as_tibble(summary.dbrda$sites) %>%
  select(RDA1, RDA2) %>%
  add_column("Group"=rarefied_metadata$label, .before=TRUE)
coordinates_sites <- inner_join(metadata, coordinates_sites, by=c("label"="Group"))
coordinates_variables <- as.data.frame(summary.dbrda$biplot) %>%
  select(RDA1, RDA2) %>%
  rownames_to_column("variable") %>%
  mutate(variable=case_when(variable=="PO4" ~ "bold('PO'['4']^{' 3–'})",
                            variable=="NH4" ~ "bold('NH'['4']^{' +'})",
                            variable=="NO2" ~ "bold('NO'['2']^{' –'})",
                            variable=="NO3" ~ "bold('NO'['3']^{' –'})",
                            variable=="SiO4" ~ "bold('SiO'['4']^{' 4–'})",
                            variable=="Chla" ~ "bold('Chl')~bolditalic('a')",
                            TRUE ~ paste0("bold(", variable, ")"))) %>%
  column_to_rownames("variable")
eigenvalues_constrained <- as.data.frame(summary.dbrda$concont$importance)

# db-RDA plot generation
p.dbrda <- ggplot() +
  geom_point(data=coordinates_sites, aes(x=RDA1, y=RDA2, fill=season, shape=station), size=5, stroke=0.5) +
  scale_fill_manual(name=NULL,
                    breaks=c("Spring", "Summer", "Autumn", "Winter"),
                    values=c("Spring"="#377EB8", "Summer"="#4DAF4A",
                             "Autumn"="#E41A1C", "Winter"="#984EA3")) +
  scale_shape_manual(name=NULL,
                     values=c(21, 23),
                     breaks=c("F", "S"),
                     labels=c("Funtana", "Saline")) +
  geom_segment(data=coordinates_variables, aes(x=0, y=0, xend=RDA1, yend=RDA2),
               arrow=arrow(length=unit(0.25, "cm")),  size=0.5, color="gray50") +
  geom_text(data=coordinates_variables,
            aes(x=RDA1+0.04*cos(atan(RDA2/RDA1))*sign(RDA1),
                y=RDA2+0.04*sin(atan(RDA2/RDA1))*sign(RDA1),
                label=rownames(coordinates_variables), family="Times"),
            col="gray50", size=5, parse=TRUE) +
  annotate("text", x=-0.48, y=-0.2,
           label=paste0("bolditalic('R'['a']^{bold('2')})~bold('=')~",
                 "bold('", format(round(R2adj*100, digits=2), nsmall=2), "')", "~bold('%')"),
           family="Times", parse=TRUE, col="gray50", size=9) + 
  labs(x=paste0("db-RDA I (", format(round(eigenvalues_constrained$RDA1[2]*R2adj*100, digits=2), nsmall=2), " %)"), 
       y=paste0("db-RDA II (", format(round(eigenvalues_constrained$RDA2[2]*R2adj*100, digits=2), nsmall=2), " %)")) +
  scale_x_continuous(labels=scaleFUN) +
  scale_y_continuous(labels=scaleFUN) +
  theme +
  theme(legend.spacing.y=unit(0.1, "cm"),
        legend.position="", legend.justification=c("left", "bottom")) +
  guides(fill=guide_legend(override.aes=list(shape=21)))

# Generating a plot to extract a common legend
p <- ggplot() +
  geom_point(data=coordinates_sites, aes(x=RDA1, y=RDA2, fill=season, shape=station), size=5, stroke=0.5) +
  scale_fill_manual(name=NULL,
                    breaks=c("Spring", "Summer", "Autumn", "Winter"),
                    values=c("Spring"="#377EB8", "Summer"="#4DAF4A",
                             "Autumn"="#E41A1C", "Winter"="#984EA3")) +
  scale_shape_manual(name=NULL,
                     values=c(21, 23),
                     breaks=c("F", "S"),
                     labels=c("Funtana", "Saline")) +
  geom_segment(data=coordinates_variables, aes(x=0, y=0, xend=RDA1, yend=RDA2),
               arrow=arrow(length=unit(0.25, "cm")),  size=0.5, color="gray50") +
  geom_text(data=coordinates_variables,
            aes(x=RDA1+0.04*cos(atan(RDA2/RDA1))*sign(RDA1),
                y=RDA2+0.04*sin(atan(RDA2/RDA1))*sign(RDA1),
                label=rownames(coordinates_variables), family="Times"),
            col="gray50", size=5, parse=TRUE) +
  labs(x=paste0("db-RDA I (", format(round(eigenvalues_constrained$RDA1[2]*R2adj*100, digits=2), nsmall=2), " %)"), 
       y=paste0("db-RDA II (", format(round(eigenvalues_constrained$RDA2[2]*R2adj*100, digits=2), nsmall=2), " %)")) +
  scale_x_continuous(labels=scaleFUN) +
  scale_y_continuous(labels=scaleFUN) +
  theme +
  theme(legend.spacing.y=unit(0.1, "cm"),
        legend.position=c(0, 0), legend.justification=c("left", "bottom")) +
  guides(fill=guide_legend(override.aes=list(shape=21)))
legend <- cowplot::get_legend(p)

# Combining plots and saving
p <- cowplot::plot_grid(p.pcoa, p.dbrda, labels=c("A", "B"), label_size=34, label_fontfamily="Times",
                        hjust=-0.2)
p <- cowplot::ggdraw() +
  cowplot::draw_plot(p, x=0, y=0, width=0.93, height=1) +
  cowplot::draw_plot(legend, x=0.93, y=0.1, width=0, height=1)
ggsave("results/figures/pcoa_dbrda_figure.jpg", p, width=2*297, height=210, units="mm")
