#################################################################################################################
# plot_calculators.R
# 
# A script to plot richness and diversity calculators.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared
#               data/raw/metadata.csv
# Produces: results/figures/calculators.jpg
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

# Copying the sample labels to the rows (input for library vegan)
row.names(rarefied) <- rarefied$Group

# Removing column containing sample labels
rarefied <- rarefied %>%
  select(-Group)

# Calculating Chao1 and ACE species estimators
estimators <- estimateR(rarefied)
estimators <- as_tibble(t(estimators)) %>%
  add_column("Group"=colnames(estimators), .before=TRUE)

# Calculating diversity indices
shannon <- diversity(rarefied, index="shannon")
shannon <- tibble(Group=names(shannon), shannon)
invsimpson <- diversity(rarefied, index="invsimpson")
invsimpson <- tibble(Group=names(invsimpson), invsimpson)

# Transforming the Shannon entropy to the effective number of species
# (http://www.loujost.com/Statistics%20and%20Physics/Diversity%20and%20Similarity/EffectiveNumberOfSpecies.htm)
shannon <- mutate(shannon, shannon=exp(shannon)) %>%
  rename(eshannon=shannon)

# Joining together estimators and indices
estimators_indices <- inner_join(estimators, shannon, by="Group") %>%
  inner_join(., invsimpson, by="Group")

# Loading metadata 
metadata <- read_tsv("data/raw/metadata.csv")

# Joining metadata with estimators and indices
Sys.setlocale(locale="en_GB.utf8")
estimators_indices_metadata <- inner_join(metadata, estimators_indices, by=c("ID"="Group")) %>%
  mutate(date=as.Date(date, "%d.%m.%Y"))

# Generating a common theme for plots
theme <- theme(text=element_text(family="Times"), line=element_line(color="black"),
               panel.border=element_rect(fill=NA), panel.background=element_blank(),
               panel.grid=element_blank(), axis.line=element_blank(),
               axis.text=element_text(size=12, color="black"), axis.text.x=element_text(angle=90, hjust=0.95, vjust=0.5),
               axis.title=element_text(size=14, color="black"),
               plot.margin=unit(c(5.5, 5.5, 5.5, 5.5), "pt"), legend.position="none",
               plot.title=element_text(size=16, hjust=0.5))
cowplot_theme <-cowplot::draw_label("Number of OTUs", x=0.07, y=0.5,
                    vjust=-0.5, angle=90, fontfamily="Times", size=14)

estimators_indices_metadata <- estimators_indices_metadata %>%
  gather("S.obs", "S.chao1", "S.ACE", "eshannon", "invsimpson", key="estimator_index", value="value")

# Defining line types, dot shapes and fill dot colors
lines_p1 <- c("S.obs"="dotted", "S.chao1"="solid", "S.ACE"="dotted")
shapes_p1 <- c("S.obs"=21, "S.chao1"=23, "S.ACE"=25)
fills_p1 <- c("S.obs"="white", "S.chao1"="black", "S.ACE"="white")

lines_p2 <- c(c("eshannon"="solid", "invsimpson"="dotted"))
shapes_p2 <- c("eshannon"=21, "invsimpson"=24)
fills_p2 <- c("eshannon"="black", "invsimpson"="white")

# Plots generation
# Seawater samples
p1 <- filter(estimators_indices_metadata, station=="F") %>%
  filter(estimator_index=="S.obs" | estimator_index=="S.chao1" | estimator_index=="S.ACE") %>%
  ggplot(aes(x=date, y=value, linetype=estimator_index, shape=estimator_index,
             fill=estimator_index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=lines_p1) +
  scale_shape_manual(values=shapes_p1) +
  scale_fill_manual(values=fills_p1) +
  scale_y_continuous(limits=c(300, 3000)) +
  labs(x="", y="") +
  ggtitle(parse(text="bold('Seawater')")) +
  theme +
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())

p2 <- filter(estimators_indices_metadata, station=="F") %>%
  filter(estimator_index=="eshannon" | estimator_index=="invsimpson") %>%
  ggplot(aes(x=date, y=value, linetype=estimator_index, shape=estimator_index,
             fill=estimator_index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=lines_p2) +
  scale_shape_manual(values=shapes_p2) +
  scale_fill_manual(values=fills_p2) +
  scale_y_continuous(limits=c(0, 150)) +
  scale_x_date(date_break ="months" , date_labels="%b %Y") +
  labs(x="", y="") +
  theme

f <- cowplot::plot_grid(p1, p2, nrow=2, ncol=1, rel_heights=c(1,0.75), align="v") +
  cowplot_theme

# Plots generation
# Cymodocea nodosa samples
p1 <- filter(estimators_indices_metadata, station=="FCyM") %>%
  filter(estimator_index=="S.obs" | estimator_index=="S.chao1" | estimator_index=="S.ACE") %>%
  ggplot(aes(x=date, y=value, linetype=estimator_index, shape=estimator_index,
             fill=estimator_index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=lines_p1) +
  scale_shape_manual(values=shapes_p1) +
  scale_fill_manual(values=fills_p1) +
  scale_y_continuous(limits=c(300, 2500)) +
  labs(x="", y="") +
  ggtitle(parse(text="bolditalic('Cymodocea nodosa')~bold('(Mixed)')")) +
  theme +
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())

p2 <- filter(estimators_indices_metadata, station=="FCyM") %>%
  filter(estimator_index=="eshannon" | estimator_index=="invsimpson") %>%
  ggplot(aes(x=date, y=value, linetype=estimator_index, shape=estimator_index,
             fill=estimator_index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=lines_p2) +
  scale_shape_manual(values=shapes_p2) +
  scale_fill_manual(values=fills_p2) +
  scale_y_continuous(limits=c(0, 450)) +
  scale_x_date(date_break="months" , date_labels="%b %Y") +
  labs(x="", y="") +
  theme

fcym <- cowplot::plot_grid(p1, p2, nrow=2, ncol=1, rel_heights=c(1,0.75), align="v")

# Plots generation
# Caulerpa cylindracea (Mixed) samples
p1 <- filter(estimators_indices_metadata, station=="FCaM") %>%
  filter(estimator_index=="S.obs" | estimator_index=="S.chao1" | estimator_index=="S.ACE") %>%
  ggplot(aes(x=date, y=value, linetype=estimator_index, shape=estimator_index,
             fill=estimator_index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=lines_p1) +
  scale_shape_manual(values=shapes_p1) +
  scale_fill_manual(values=fills_p1) +
  scale_y_continuous(limits=c(1300, 4300)) +
  labs(x="", y="") +
  ggtitle(parse(text="bolditalic('Caulerpa cylindracea')~bold('(Mixed)')")) +
  theme +
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())

p2 <- filter(estimators_indices_metadata, station=="FCaM") %>%
  filter(estimator_index=="eshannon" | estimator_index=="invsimpson") %>%
  ggplot(aes(x=date, y=value, linetype=estimator_index, shape=estimator_index,
             fill=estimator_index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=lines_p2) +
  scale_shape_manual(values=shapes_p2) +
  scale_fill_manual(values=fills_p2) +
  scale_y_continuous(limits=c(0, 700)) +
  scale_x_date(date_break ="months" , date_labels="%b %Y") +
  labs(x="Date", y="") +
  theme

fcam <- cowplot::plot_grid(p1, p2, nrow=2, ncol=1, rel_heights=c(1,0.75), align="v") +
  cowplot_theme

# Plots generation
# Caulerpa cylindracea (Monospecific) samples
p1 <- filter(estimators_indices_metadata, station=="FCa") %>%
  filter(estimator_index=="S.obs" | estimator_index=="S.chao1" | estimator_index=="S.ACE") %>%
  ggplot(aes(x=date, y=value, linetype=estimator_index, shape=estimator_index,
             fill=estimator_index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=lines_p1) +
  scale_shape_manual(values=shapes_p1) +
  scale_fill_manual(values=fills_p1) +
  scale_y_continuous(limits=c(1300, 4300)) +
  labs(x="", y="") +
  ggtitle(parse(text="bolditalic('Caulerpa cylindracea')~bold('(Monospecific)')")) +
  theme +
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())

p2 <- filter(estimators_indices_metadata, station=="FCa") %>%
  filter(estimator_index=="eshannon" | estimator_index=="invsimpson") %>%
  ggplot(aes(x=date, y=value, linetype=estimator_index, shape=estimator_index,
             fill=estimator_index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=lines_p2) +
  scale_shape_manual(values=shapes_p2) +
  scale_fill_manual(values=fills_p2) +
  scale_y_continuous(limits=c(0, 700)) +
  scale_x_date(date_break ="months" , date_labels="%b %Y") +
  labs(x="Date", y="") +
  theme

fca <- cowplot::plot_grid(p1, p2, nrow=2, ncol=1, rel_heights=c(1, 0.75), align="v")

# Generating a plot to extract a common legend
labels <- c("Observed Number of OTUs", "Chao1", "ACE", "Exponential Shannon", "Inverse Simpson")
p1 <- filter(estimators_indices_metadata, station=="F") %>%
  ggplot(aes(x=date, y=value, linetype=estimator_index, shape=estimator_index,
             fill=estimator_index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=c(lines_p1, lines_p2),
                        breaks=c(names(lines_p1), names(lines_p2)),
                        labels=labels) +
  scale_shape_manual(values=c(shapes_p1, shapes_p2),
                     breaks=c(names(lines_p1), names(lines_p2)),
                     labels=labels) +
  scale_fill_manual(values=c(fills_p1, fills_p2),
                    breaks=c(names(lines_p1), names(lines_p2)),
                    labels=labels) +
  labs(x="", y="") +
  theme +
  theme(legend.position="bottom", legend.title=element_blank(),
        legend.text=element_text(size=10, margin=margin(r=0.2, unit="cm")),
        legend.key.width=unit(1.4, "cm"), legend.key.height=unit(0.5, "cm"),
        legend.key=element_rect(fill="white"), legend.justification=c("top"),
        legend.text.align=0) +
  guides(linetype=guide_legend(ncol=2))
legend <- cowplot::get_legend(p1)

# Combining plots together and saving
plots <- cowplot::plot_grid(f, fcym, fcam, fca, ncol=2, nrow=2)
p <- cowplot::plot_grid(plots, legend, ncol=1, nrow=2, rel_heights=c(4.6, 0.4))
ggsave("results/figures/calculators.jpg", p, width=210, height=297, units="mm")
