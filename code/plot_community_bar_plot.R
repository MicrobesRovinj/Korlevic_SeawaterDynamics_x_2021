#################################################################################################################
# plot_community_bar_plot.R
# 
# A script to plot the community structure of each sample.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v138.wang.tax.summary
#               data/raw/metadata.csv
#               data/raw/group_colors.csv
# Produces: results/figures/community_bar_plot.jpg
#
#################################################################################################################

# Loading input data containing sequence abundances and subsequent input data customization
community <- read_tsv("data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v138.wang.tax.summary") %>%
  filter(!str_detect(taxon, "^Eukaryota")) %>%
  filter(taxon!="Root")

# Remove chloroplast and mitochondrial sequences and subtract their number from higher taxonomic levels to which
# they belong
chloroplast <- filter(community, str_detect(taxon, "^Chloroplast$"))$rankID
mitochondria <- filter(community, str_detect(taxon, "^Mitochondria$"))$rankID
community <- mutate_at(community, 5:ncol(community), list(~case_when(
    rankID==str_extract(chloroplast, "(\\d+\\.){3}\\d+") ~ . - .[taxon=="Chloroplast"],
    rankID==str_extract(chloroplast, "(\\d+\\.){2}\\d+") ~ . - .[taxon=="Chloroplast"],
    rankID==str_extract(chloroplast, "(\\d+\\.){1}\\d+") ~ . - .[taxon=="Chloroplast"],
  TRUE ~ .))) %>%
  filter(!str_detect(taxon, "^Chloroplast")) %>%
  mutate_at(5:ncol(.), list(~case_when(
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
  mutate_at(5:ncol(.), list(~ . / sum(.) * 100)) %>%
  ungroup()

# Selection of groups for plotting
plot <- filter(community, taxlevel==2 |
                 (taxlevel==4 & str_detect(taxon, "^Chloroplast$")) |
                 (taxlevel==3 & str_detect(rankID, filter(community, str_detect(taxon, "^Proteobacteria$"))$rankID))) %>%
  filter_at(6:ncol(.), any_vars(. >= 1)) %>%
  mutate_at(5:ncol(.), list(~case_when(taxon=="Proteobacteria" ~ . - sum(.[taxlevel==3 &
    str_detect(rankID, filter(community, str_detect(taxon, "^Proteobacteria$"))$rankID)]), TRUE ~ .))) %>%
  mutate(taxon=str_replace(taxon, "Proteobacteria", "Other Proteobacteria")) %>%
  mutate(taxon=str_replace_all(taxon, c("unknown_unclassified"="No_Relative", "unknown"="No_Relative"))) %>%
  filter_at(6:ncol(.), any_vars(. >= 1)) %>%
  bind_rows(summarise_all(., list(~ifelse(is.numeric(.), 100-sum(.), paste("Other"))))) %>%
  arrange(taxon %in% "No_Relative")

# Loading colors for each group on the plot
color <- read_tsv("data/raw/group_colors.csv") %>%
  select(-Taxlevel) %>%
  deframe()

# Generation of italic names for groups
names <- parse(text=case_when(plot$taxon=="Chloroplast" ~ paste0("plain('", plot$taxon,  "')"),
                              plot$taxon=="Bacteria_unclassified" ~ "italic('Bacteria')~plain('(NR)')",
                              plot$taxon=="Marinimicrobia_(SAR406_clade)" ~ "italic('Marinimicrobia')",
                              plot$taxon=="Other" ~ paste0("plain('", plot$taxon, "')"),
                              plot$taxon=="No_Relative" ~ "plain('No Relative')",
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
  mutate(taxon=factor(taxon, levels=unique(plot$taxon))) %>%
  mutate(label=factor(label, levels=metadata$label)) %>%
  mutate(date=as.Date(date, "%d.%m.%Y"))

# Generating a common theme for plots
theme <- theme(text=element_text(family="Times"), line=element_line(color="black"),
               panel.grid=element_blank(), 
               axis.line.x=element_line(), axis.line.y=element_line(),
               axis.ticks.x=element_line(), axis.ticks.y.left=element_line(),
               axis.text.y=element_text(size=12, color="black"),
               axis.text.x=element_blank(), axis.title.y=element_text(size=14, color="black", vjust=-0.75),
               panel.background=element_blank(), plot.margin=unit(c(5.5, 5.5, 5.5, -11), "pt"), legend.position="none",
               plot.title=element_text(size=16, hjust=0.13))

# Plots generation
f <- filter(plot, station=="F") %>%
  ggplot() +
  geom_bar(aes(x=date, y=abundance, fill=taxon), stat="identity", colour="black", size=0.3, width=8) +
  scale_fill_manual(values=color, labels=names) + 
  labs(x=NULL, y="%") +
  ggtitle(parse(text="bold('Seawater')")) +
  scale_x_date(date_break="months" , date_labels="%b %Y", limits=as.Date(c("2017-07-01", "2018-10-15")),
               expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0), breaks=seq(0, 100, by=10)) +
  theme +
  theme(axis.text.x=element_text(size=12, color="black", angle=90,
                                 hjust=1, vjust=0.5, margin=margin(t=5.5, unit="pt")),
        axis.title.x=element_text(size=14, color="black"), plot.margin=unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        plot.title=element_text(size=16, hjust=0.5))

fcym <- filter(plot, station=="FCyM") %>%
  ggplot() +
  geom_bar(aes(x=date, y=abundance, fill=taxon), stat="identity", colour="black", size=0.3, width=8) +
  scale_fill_manual(values=color, labels=names) + 
  labs(x=NULL, y="%") +
  ggtitle(parse(text="bolditalic('Cymodocea nodosa')~bold('(Mixed)')")) +
  scale_x_date(date_break="months" , date_labels="%b %Y", limits=as.Date(c("2017-11-01", "2018-10-15")),
               expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0), breaks=seq(0, 100, by=10)) +
  theme

fcam <- filter(plot, station=="FCaM") %>%
  ggplot() +
  geom_bar(aes(x=date, y=abundance, fill=taxon), stat="identity", colour="black", size=0.3, width=8) +
  scale_fill_manual(values=color, labels=names) + 
  labs(x=NULL, y="%") +
  ggtitle(parse(text="bolditalic('Caulerpa cylindracea')~bold('(Mixed)')")) +
  scale_x_date(date_break="months" , date_labels="%b %Y", limits=as.Date(c("2017-11-01", "2018-10-15")),
               expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0), breaks=seq(0, 100, by=10)) +
  theme

fca <- filter(plot, station=="FCa") %>%
  ggplot() +
  geom_bar(aes(x=date, y=abundance, fill=taxon), stat="identity", colour="black", size=0.3, width=8) +
  scale_fill_manual(values=color, labels=names) + 
  labs(x="Date", y="%") +
  ggtitle(parse(text="bolditalic('Caulerpa cylindracea')~bold('(Monospecific)')")) +
  scale_x_date(date_break="months" , date_labels="%b %Y", limits=as.Date(c("2017-11-01", "2018-10-15")),
               expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0), breaks=seq(0, 100, by=10)) +
  theme +
  theme(axis.text.x=element_text(size=12, color="black", angle=90, hjust=1, vjust=0.5,
                                 margin=margin(t=5.5, unit="pt")),
        axis.title.x=element_text(size=14, color="black", margin=margin(t=11, unit="pt")))

# Generating a plot to extract a common legend
p1 <- ggplot(plot) +
  geom_bar(aes(x=date, y=abundance, fill=taxon), stat="identity", colour="black", size=0.3, width=8) +
  scale_fill_manual(values=color, labels=names) + 
  labs(x=NULL, y="%") +
  scale_x_date(date_break="months" , date_labels="%b %Y", limits=as.Date(c("2017-07-01", "2018-10-15")),
               expand=c(0, 0)) +
  scale_y_continuous(expand=c(0, 0), breaks=seq(0, 100, by=10)) +
  theme +
  theme(legend.position="right", 
        legend.title=element_blank(), legend.text=element_text(size=12),
        legend.spacing.x=unit(0.2, "cm"), legend.justification=c("bottom"),
        legend.box.margin=margin(0, 0, 65, -20), legend.text.align=0,
        legend.key.size=unit(0.5, "cm"))
legend <- cowplot::get_legend(p1)

# Combining plots together and saving
p <- cowplot::ggdraw() +
  cowplot::draw_plot(f, x=0, y=0.708, width=1, height=0.288) +
  cowplot::draw_plot(fcym, x=0.267, y=0.545, width=0.733, height=0.228) +
  cowplot::draw_plot(fcam, x=0.267, y=0.317, width=0.733, height=0.228) +
  cowplot::draw_plot(fca, x=0.267, y=0, width=0.733, height=0.317) +
  cowplot::draw_plot(legend, x=0.1, y=0.28, width=0.1, height=0.2)
ggsave("results/figures/community_bar_plot.jpg", p, width=210, height=297, units="mm")
