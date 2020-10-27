#################################################################################################################
# plot_pcoa.R
# 
# A script to generate the PCoA figure.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared
#               data/raw/metadata.csv
# Produces: results/figures/pcoa_figure.jpg
#
#################################################################################################################

# Loading OTU/sample data
shared <- read_tsv("data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared")

# Generation of random rarefied community data
rarefied <- shared %>%
  select(-label, -Group, -numOtus) %>%
  rrarefy(., min(rowSums(.))) %>%
  as_tibble(.name_repair="unique") %>%
  add_column("Group"=shared$Group, .before=TRUE) %>%
  select_if(list(~!is.numeric(.) || sum(.)!=0))

# Loading metadata 
metadata <- read_tsv("data/raw/metadata.csv")

# Joining metadata and point coordinates
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
         legend.margin=margin(-5, 0, 0, 0), legend.position=c(0.01, 0.99),
         legend.justification=c("left", "top"), plot.title=element_text(size=14, hjust=0.5))

# Plots generation
rarefied_metadata_select <- rarefied_metadata

# Generation of PCoA data
spe.bray <- vegdist(select(rarefied_metadata, starts_with("Otu")))
spe.b.pcoa <- cmdscale(spe.bray, k=(nrow(rarefied_metadata_select)-1), eig=TRUE)
coordinates <- spe.b.pcoa$points[,1:2] %>%
  as_tibble(.name_repair= ~c("A1", "A2")) %>%
  add_column("Group"=rarefied_metadata_select$label, .before=TRUE)
coordinates <- inner_join(metadata, coordinates, by=c("label"="Group"))

# PCoA plot generation
p <- ggplot() +
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
        legend.position=c(0.04, 0.64), legend.justification=c("left", "bottom")) +
  guides(fill=guide_legend(override.aes=list(shape=21)))

# Saving plot
ggsave("results/figures/pcoa_figure.jpg", p, width=297, height=210, units="mm")
