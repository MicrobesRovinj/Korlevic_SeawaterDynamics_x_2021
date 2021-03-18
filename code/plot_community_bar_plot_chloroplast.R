#################################################################################################################
# plot_community_bar_plot_chloroplast.R
# 
# A script to plot chloroplast sequence abundances of each sample.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v138.wang.tax.summary
#               data/raw/metadata.csv
#               data/raw/group_colors.csv
# Produces: results/figures/chloroplast_bar_plot.jpg
#
#################################################################################################################

# Loading input data containing sequence abundances and subsequent input data customization
community <- read_tsv("data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v138.wang.tax.summary") %>%
  filter(!str_detect(taxon, "^Eukaryota")) %>%
  filter(taxon!="Root")

# Removing mitochondrial sequences and subtracting their number from higher taxonomic levels to which
# they belong
mitochondria <- filter(community, str_detect(taxon, "^Mitochondria$"))$rankID
community <- mutate_at(community, 5:ncol(community), list(~case_when(
    rankID==str_extract(mitochondria, "(\\d+\\.){4}\\d+") ~ . - .[taxon=="Mitochondria"],
    rankID==str_extract(mitochondria, "(\\d+\\.){3}\\d+") ~ . - .[taxon=="Mitochondria"],
    rankID==str_extract(mitochondria, "(\\d+\\.){2}\\d+") ~ . - .[taxon=="Mitochondria"],
    rankID==str_extract(mitochondria, "(\\d+\\.){1}\\d+") ~ . - .[taxon=="Mitochondria"],
    TRUE ~ .))) %>%
  filter(!str_detect(taxon, "^Mitochondria")) %>%
  select(-ATCC_1, -ATCC_5, -NC_1) %>%
  mutate(`23`=`23_1`+`23_2`) %>%
  select(-`23_1`, -`23_2`) %>%
  group_by(taxlevel) %>%
  mutate_at(5:ncol(.), list(~. / sum(.) * 100)) %>%
  ungroup()

# Selecting groups for plotting
plot <- community %>%
  filter(taxon=="Chloroplast") %>%
  bind_rows(summarise_all(., list(~ifelse(is.numeric(.), 100-sum(.), paste("Other")))))

# Loading colors for each group on the plot
color <- read_tsv("data/raw/group_colors.csv") %>%
  select(-Taxlevel) %>%
  deframe()

# Generating italic names for groups
names <- parse(text=case_when(plot$taxon=="Chloroplast" ~ paste0("plain('", plot$taxon,  "')"),
                              plot$taxon=="Bacteria_unclassified" ~ "italic('Bacteria')~plain('(No Relative)')",
                              plot$taxon=="Marinimicrobia_(SAR406_clade)" ~ "italic('Marinimicrobia')",
                              plot$taxon=="Other" ~ paste0("plain('", plot$taxon, "')"),
                              plot$taxon=="No Relative" ~ paste0("plain('", plot$taxon, "')"),
                              TRUE ~ paste0("italic('", plot$taxon, "')")))

# Tidying the sequence abundance data
plot <- gather(plot, key="Group", value="abundance", 6:(ncol(plot)))

# Loading metadata
metadata <- read_tsv("data/raw/metadata.csv") %>%
  filter(ID!="23_1") %>%
  mutate(ID=str_replace(ID, "23_2", "23")) %>%
  mutate(label=str_replace(label, "24/4-18 F-2", "24/4-18 F"))

# Joining sequence abundances data and metadata
Sys.setlocale(locale="en_GB.utf8")
plot <- inner_join(metadata, plot, by=c("ID"="Group")) %>%
  mutate(taxon=factor(taxon, levels=rev(unique(plot$taxon)))) %>%
  mutate(label=factor(label, levels=metadata$label)) %>%
  mutate(station=factor(station, levels=c("S", "F"))) %>%
  mutate(date=as.Date(date, "%d.%m.%Y"))

# Generating a common theme for plots
theme <- theme(text=element_text(family="Times"), line=element_line(color="black"),
               panel.grid=element_blank(), 
               axis.line.x=element_line(), axis.line.y=element_line(),
               axis.ticks.x=element_line(), axis.ticks.y.left=element_line(),
               axis.text.y=element_text(size=12, color="black"),
               axis.text.x=element_blank(), axis.title.y=element_text(size=14, color="black", vjust=-0.75),
               panel.background=element_blank(), plot.margin=unit(c(5.5, 5.5, 5.5, 5,5), "pt"), legend.position="none",
               plot.title=element_text(size=16, hjust=0.5))

# Generating plot
p <- ggplot(plot) +
  geom_bar(aes(x=date, y=abundance, fill=taxon), stat="identity", colour="black",
           size=0.3, width=8) +
  scale_fill_manual(values=color, labels=names, breaks=c("Chloroplast", "Other")) + 
  labs(x="Date", y="%", tag="%") +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0), breaks=seq(0, 100, by=10)) +
  theme(text=element_text(family="Times"),
        line=element_line(color="black"),
        panel.grid=element_blank(),
        panel.background=element_blank(),
        plot.margin=unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        plot.title=element_text(size=16, hjust=0.5),
        plot.tag.position=c(0.01, 0.312),
        plot.tag=element_text(size=14, color="black", angle=90),
        axis.line.x=element_line(),
        axis.line.y=element_line(),
        axis.ticks.x=element_line(),
        axis.ticks.y.left=element_line(),
        axis.text.x=element_text(size=12, color="black", angle=90, hjust=1,
                                 vjust=3, margin=margin(t=5.5, unit="pt")),
        axis.text.y=element_text(size=12, color="black"),
        axis.title.x=element_text(size=14, color="black"),
        axis.title.y=element_text(size=14, color="black", vjust=0, hjust=0.778),
        legend.position="right",
        legend.title=element_blank(),
        legend.text=element_text(size=12),
        legend.spacing.x=unit(0.2, "cm"),
        legend.justification=c("left", "bottom"),
        legend.box.margin=unit(c(0, 0, 0, 0), "pt"),
        legend.text.align=0,
        legend.key.size=unit(0.5, "cm")) +
  facet_rep_wrap(station ~ ., strip.position="top",
                 nrow=2, ncol=1,
                 labeller=as_labeller(c(`S`="Bay of Saline", `F`="Bay of Funtana"))) +
  theme(strip.background=element_blank(),
        strip.placement="outside",
        panel.spacing=unit(11, "pt"),
        strip.text=element_text(face="bold", size=22, hjust=0.5))

# Saving
ggsave("results/figures/chloroplast_bar_plot.jpg", p, width=297, height=210, units="mm")