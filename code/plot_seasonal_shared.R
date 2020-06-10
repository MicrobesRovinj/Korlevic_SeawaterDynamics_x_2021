#################################################################################################################
# code/plot_seasonal_shared.R
# 
# A script to plot the Bray-Curtis and Jaccard's Similarity Coefficients of subsequent sampling dates.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared
#               data/raw/metadata.csv
# Produces: results/figures/temporal_shared.jpg
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
  select_if(list(~ !is.numeric(.) || sum(.)!=0))

# Copying the sample labels to the rows (input for library vegan)
row.names(rarefied) <- rarefied$Group

# Removing column containing sample labels
rarefied <- rarefied %>%
  select(-Group)

# Calculating dissimilarity indices
jaccard <- vegdist(rarefied, method="jaccard", binary=T)
jaccard <- as_tibble(data.frame(t(combn(rownames(rarefied), 2)), as.numeric(jaccard)),
                     .name_repair= ~c("V1", "V2", "jaccard"))

bray <- vegdist(rarefied, method="bray", binary=F)
bray <- as_tibble(data.frame(t(combn(rownames(rarefied), 2)), as.numeric(bray)),
                  .name_repair= ~c("V1", "V2", "bray"))

distance <- inner_join(jaccard, bray, by=c("V1"="V1", "V2"="V2")) %>%
  mutate_at(c("V1", "V2"), list(~ as.character(.)))

# Loading metadata 
metadata <- read_tsv("data/raw/metadata.csv")

# Joining metadata with dissimilarity indices data
Sys.setlocale(locale="en_GB.utf8")
distance_metadata <- inner_join(distance, metadata, by=c("V1"="ID")) %>%
  inner_join(., metadata, by=c("V2"="ID")) %>%
  mutate(date.x=as.Date(date.x, "%d.%m.%Y")) %>%
  mutate(date.y=as.Date(date.y, "%d.%m.%Y")) %>%
  mutate(jaccard=(1-jaccard)*100) %>%
  mutate(bray=(1-bray)*100) %>%
  gather(key="index", value="value", jaccard, bray)

# Generating a common theme for plots
theme <- theme(text=element_text(family="Times"), line=element_line(color="black"),
               panel.border=element_rect(fill=NA), panel.background=element_blank(),
               panel.grid=element_blank(), axis.line=element_blank(),
               axis.text=element_text(size=12, color="black"), axis.text.x=element_text(angle=90, hjust=0.95, vjust=0.5),
               axis.title=element_text(size=14, color="black"),
               plot.margin=unit(c(5.5, 5.5, 5.5, 5.5), "pt"), legend.position="none",
               plot.title=element_text(size=16, hjust=0.5))

# Defining line types, dot shapes and fill dot colors
lines <- c("jaccard"="dotted", "bray"="solid")
shapes <- c("jaccard"=21, "bray"=23)
fills <- c("jaccard"="white", "bray"="black")

# Plots generation
# Seawater samples
data <- filter(distance_metadata, station.x=="F" & station.y=="F") %>%
  filter((date.x=="2017-07-13" & date.y=="2017-07-27") |
         (date.x=="2017-07-27" & date.y=="2017-08-10") |
         (date.x=="2017-08-10" & date.y=="2017-08-24") |
         (date.x=="2017-08-24" & date.y=="2017-09-19") |
         (date.x=="2017-09-19" & date.y=="2017-10-12") |
         (date.x=="2017-10-12" & date.y=="2017-11-23") |
         (date.x=="2017-11-23" & date.y=="2017-12-14") |
         (date.x=="2017-12-14" & date.y=="2018-02-13") |
         (date.x=="2018-02-13" & date.y=="2018-03-27") |
         (date.x=="2018-03-27" & date.y=="2018-04-24") |
         (date.x=="2018-04-24" & date.y=="2018-05-22") |
         (date.x=="2018-05-22" & date.y=="2018-06-19") |
         (date.x=="2018-06-19" & date.y=="2018-07-10") |
         (date.x=="2018-07-10" & date.y=="2018-08-09") |
         (date.x=="2018-08-09" & date.y=="2018-09-04") |
         (date.x=="2018-04-09" & date.y=="2018-10-05"))

p1 <- filter(data, index=="bray") %>%
  ggplot(aes(x=date.y, y=value, linetype=index, shape=index, fill=index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=lines) +
  scale_shape_manual(values=shapes) +
  scale_fill_manual(values=fills) +
  scale_y_continuous(limits=c(0, 100)) +
  labs(x="", y="%") +
  ggtitle(parse(text="bold('Seawater')")) +
  theme +
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())

p2 <- filter(data, index=="jaccard") %>%
  ggplot(aes(x=date.y, y=value, linetype=index, shape=index, fill=index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=lines) +
  scale_shape_manual(values=shapes) +
  scale_fill_manual(values=fills) +
  scale_x_date(date_break ="months" , date_labels="%b %Y") +
  scale_y_continuous(limits=c(10, 30)) +
  labs(x="", y="%") +
  theme

f <- cowplot::plot_grid(p1, p2, nrow=2, ncol=1, rel_heights=c(1, 0.75), align="v")

# Plots generation
# Cymodocea nodosa samples
data <- filter(distance_metadata, station.x=="FCyM" & station.y=="FCyM") %>%
  filter((date.x=="2017-11-23" & date.y=="2017-12-04") |
         (date.x=="2017-12-04" & date.y=="2017-12-14") |
         (date.x=="2017-12-14" & date.y=="2018-02-13") |
         (date.x=="2018-02-13" & date.y=="2018-03-27") |
         (date.x=="2018-03-27" & date.y=="2018-04-24") |
         (date.x=="2018-04-24" & date.y=="2018-05-22") |
         (date.x=="2018-05-22" & date.y=="2018-06-19") |
         (date.x=="2018-06-19" & date.y=="2018-07-10") |
         (date.x=="2018-07-10" & date.y=="2018-08-09") |
         (date.x=="2018-08-09" & date.y=="2018-09-04") |
         (date.x=="2018-09-04" & date.y=="2018-10-05"))

p1 <- filter(data, index=="bray") %>%
  ggplot(aes(x=date.y, y=value, linetype=index, shape=index, fill=index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=lines) +
  scale_shape_manual(values=shapes) +
  scale_fill_manual(values=fills) +
  scale_y_continuous(limits=c(0, 100)) +
  labs(x="", y="") +
  ggtitle(parse(text="bolditalic('Cymodocea nodosa')~bold('(Mixed)')")) +
  theme +
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())

p2 <- filter(data, index=="jaccard") %>%
  ggplot(aes(x=date.y, y=value, linetype=index, shape=index, fill=index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=lines) +
  scale_shape_manual(values=shapes) +
  scale_fill_manual(values=fills) +
  scale_x_date(date_break ="months" , date_labels="%b %Y") +
  scale_y_continuous(limits=c(20, 50)) +
  labs(x="", y="") +
  theme

fcym <- cowplot::plot_grid(p1, p2, nrow=2, ncol=1, rel_heights=c(1, 0.75), align="v")

# Plots generation
# Caulerpa cylindracea (Mixed) samples
data <- filter(distance_metadata, station.x=="FCaM" & station.y=="FCaM") %>%
  filter((date.x=="2017-11-23" & date.y=="2017-12-04") |
         (date.x=="2017-12-04" & date.y=="2017-12-14") |
         (date.x=="2017-12-14" & date.y=="2018-02-13") |
         (date.x=="2018-02-13" & date.y=="2018-03-27") |
         (date.x=="2018-03-27" & date.y=="2018-04-24") |
         (date.x=="2018-04-24" & date.y=="2018-05-22") |
         (date.x=="2018-05-22" & date.y=="2018-06-19") |
         (date.x=="2018-06-19" & date.y=="2018-07-10") |
         (date.x=="2018-07-10" & date.y=="2018-08-09") |
         (date.x=="2018-08-09" & date.y=="2018-09-04") |
         (date.x=="2018-09-04" & date.y=="2018-10-05"))

p1 <- filter(data, index=="bray") %>%
  ggplot(aes(x=date.y, y=value, linetype=index, shape=index, fill=index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=lines) +
  scale_shape_manual(values=shapes) +
  scale_fill_manual(values=fills) +
  scale_y_continuous(limits=c(0, 100)) +
  labs(x="", y="%") +
  ggtitle(parse(text="bolditalic('Caulerpa cylindracea')~bold('(Mixed)')")) +
  theme +
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())

p2 <- filter(data, index=="jaccard") %>%
  ggplot(aes(x=date.y, y=value, linetype=index, shape=index, fill=index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=lines) +
  scale_shape_manual(values=shapes) +
  scale_fill_manual(values=fills) +
  scale_x_date(date_break ="months" , date_labels="%b %Y") +
  scale_y_continuous(limits=c(20, 40)) +
  labs(x="Date", y="%") +
  theme

fcam <- cowplot::plot_grid(p1, p2, nrow=2, ncol=1, rel_heights=c(1, 0.75), align="v")

# Plots generation
# Caulerpa cylindracea (Monospecific) samples
data <- filter(distance_metadata, station.x=="FCa" & station.y=="FCa") %>%
  filter((date.x=="2017-11-23" & date.y=="2017-12-04") |
         (date.x=="2017-12-04" & date.y=="2017-12-14") |
         (date.x=="2017-12-14" & date.y=="2018-02-13") |
         (date.x=="2018-02-13" & date.y=="2018-03-27") |
         (date.x=="2018-03-27" & date.y=="2018-04-24") |
         (date.x=="2018-04-24" & date.y=="2018-05-22") |
         (date.x=="2018-05-22" & date.y=="2018-06-19") |
         (date.x=="2018-06-19" & date.y=="2018-07-10") |
         (date.x=="2018-07-10" & date.y=="2018-08-09") |
         (date.x=="2018-08-09" & date.y=="2018-09-04") |
         (date.x=="2018-09-04" & date.y=="2018-10-05"))

p1 <- filter(data, index=="bray") %>%
  ggplot(aes(x=date.y, y=value, linetype=index, shape=index, fill=index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=lines) +
  scale_shape_manual(values=shapes) +
  scale_fill_manual(values=fills) +
  scale_y_continuous(limits=c(0, 100)) +
  labs(x="", y="") +
  ggtitle(parse(text="bolditalic('Caulerpa cylindracea')~bold('(Monospecific)')")) +
  theme +
  theme(axis.ticks.x=element_blank(), axis.text.x=element_blank())

p2 <- filter(data, index=="jaccard") %>%
  ggplot(aes(x=date.y, y=value, linetype=index, shape=index, fill=index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=lines) +
  scale_shape_manual(values=shapes) +
  scale_fill_manual(values=fills) +
  scale_x_date(date_break ="months" , date_labels="%b %Y") +
  scale_y_continuous(limits=c(20, 40)) +
  labs(x="Date", y="") +
  theme

fca <- cowplot::plot_grid(p1, p2, nrow=2, ncol=1, rel_heights=c(1, 0.75), align="v")

# Generating a plot to extract a common legend
labels <- c("Bray-Curtis Similarity Coefficient", "Jaccard's Similarity Coefficient")
p1 <- distance_metadata %>%
  ggplot(aes(x=date.y, y=value, linetype=index, shape=index, fill=index)) +
  geom_line() +
  geom_point(size=3) +
  scale_linetype_manual(values=lines, labels=labels) +
  scale_shape_manual(values=shapes, labels=labels) +
  scale_fill_manual(values=fills, labels=labels) +
  labs(x="", y="") +
  theme +
  theme(legend.position="bottom", legend.title=element_blank(),
        legend.text=element_text(size=10, margin=margin(r=0.2, unit="cm")),
        legend.key.width=unit(1.4, "cm"), legend.key.height=unit(0.5, "cm"),
        legend.key=element_rect(fill="white"), legend.justification=c("top"),
        legend.text.align=0) +
  guides(linetype=guide_legend(ncol=1))
legend <- cowplot::get_legend(p1)

# Combining plots together and saving
plots <- cowplot::plot_grid(f, fcym, fcam, fca, ncol=2, nrow=2)
p <- cowplot::plot_grid(plots, legend, ncol=1, nrow=2, rel_heights=c(4.6, 0.4))
ggsave("results/figures/seasonal_shared.jpg", p, width=210, height=297, units="mm")
