#################################################################################################################
# plot_community_bar_plot_alphaproteobacteria.R
# 
# A script to plot sequence abundances of groups from the class Alphaproteobacteria of each sample.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v138.wang.tax.summary
#               data/raw/metadata.csv
#               data/raw/group_colors.csv
# Produces: results/figures/alphaproteobacteria_bar_plot.jpg
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
  mutate_at(5:ncol(.), list(~. / sum(.) * 100)) %>%
  ungroup()

# Selection of groups for plotting
select <- filter(community,
               taxlevel==6 &
               str_detect(rankID, filter(community, str_detect(taxon, "^Alphaproteobacteria$"))$rankID)) %>%
  filter_at(6:ncol(.), any_vars(. >= 2))

plot <- filter(community,
               taxlevel==6 &
               str_detect(rankID, filter(community, str_detect(taxon, "^Alphaproteobacteria$"))$rankID)) %>%
  mutate_at(5:ncol(.), list(~. / sum(.) * 100)) %>%
  ungroup() %>%
  filter(rankID %in% select$rankID) %>%
  bind_rows(summarise_all(., list(~ifelse(is.numeric(.), 100-sum(.), paste("Other_Alphaproteobacteria"))))) %>%
# Remove last digit from rankID so for the next step
  mutate(rankID=if_else(taxon=="uncultured", str_replace(rankID, "\\.\\d+$", ""), rankID))

# Adding data to "uncultured" taxa describing the higher taxonomic level to which they belong
uncultured <- select(filter(community, rankID %in% filter(plot, taxon=="uncultured")$rankID), rankID, taxon) %>%
  rename(rankID_uncultured=rankID, taxon_uncultured=taxon)

plot <- left_join(plot, uncultured, by=c("rankID"="rankID_uncultured")) %>%
  mutate(taxon=if_else(taxon=="uncultured", paste0(taxon, "_", taxon_uncultured), taxon)) %>%
  select(-taxon_uncultured)

# Loading colors for each group on the plot
color <- read_tsv("data/raw/group_colors.csv") %>%
  select(-Taxlevel) %>%
  deframe()

# Generation of italic names for taxonomic groups
names <- parse(text=case_when(str_detect(plot$taxon, "uncultured") ~ paste0("plain('Uncultured')~italic('", str_remove(plot$taxon, "uncultured_"), "')"),
                              str_detect(plot$taxon, "unclassified") ~ paste0("italic('", str_remove(plot$taxon, "_unclassified"), "')~plain('(NR)')"),
                              plot$taxon=="OCS116_clade_ge" ~ "plain('OCS116 Clade')",
                              plot$taxon=="Candidatus_Puniceispirillum" ~ "italic('\"Candidatus')~plain('Puniceispirillum\"')",
                              plot$taxon=="SAR116_clade_ge" ~ "plain('SAR116 Clade')",
                              plot$taxon=="Stappiaceae_ge" ~ "italic('Stappiaceae')",
                              plot$taxon=="HIMB11" ~ "plain('HIMB11')",
                              plot$taxon=="AEGEAN-169_marine_group_ge" ~ "plain('AEGEAN-169 Marine Group')",
                              plot$taxon=="Clade_Ia" ~ "plain('Clade Ia')",
                              plot$taxon=="Clade_Ib" ~ "plain('Clade Ib')",
                              plot$taxon=="Clade_II_ge" ~ "plain('Clade II')",
                              plot$taxon=="Clade_III_ge" ~ "plain('Clade III')",
                              plot$taxon=="Other_Alphaproteobacteria" ~ "plain('Other')~italic('Alphaproteobacteria')",
                              TRUE ~ paste0("italic('", plot$taxon, "')")))

# Tidying the sequence abundance data
plot <- gather(plot, key="Group", value="abundance", 6:ncol(plot))

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

# Selecting the relative abundance of the targeted group in the whole community, tidying the obtained data and
# joining with metadata
whole <- filter(community,
                taxlevel==3) %>%
  gather(key="Group", value="abundance", 6:ncol(.)) %>%
  filter(taxon=="Alphaproteobacteria")

whole <- inner_join(metadata, whole, by=c("ID"="Group")) %>%
  mutate(date=as.Date(date, "%d.%m.%Y"))

# Generating a common theme for plots
theme <- theme(text=element_text(family="Times"), line=element_line(color="black"),
               panel.grid=element_blank(), 
               axis.line.x=element_line(), axis.line.y=element_blank(),
               axis.ticks.x=element_line(), axis.ticks.y.left=element_line(),
               axis.text.y=element_text(size=12, color="black"),
               axis.text.x=element_blank(), axis.title.y=element_text(size=14, color="black", vjust=-0.75, hjust=0.365),
               panel.background=element_blank(), plot.margin=unit(c(5.5, 5.5, 5.5, -11), "pt"), legend.position="none",
               plot.title=element_text(size=16, hjust=0.13))

# Plots generation
f <- filter(plot, station=="F") %>%
  ggplot() +
  geom_bar(aes(x=date, y=abundance, fill=taxon), stat="identity", colour="black", size=0.3, width=8) +
  scale_fill_manual(values=color, labels=names) + 
  geom_text(data=filter(whole, station=="F"), aes(x=date, y=75, label=paste(format(round(abundance, 1), nsmall=1), "%")),
            family="Times", fontface="bold", angle=90, hjust=-1) +
  labs(x=NULL, y="%") +
  ggtitle(parse(text="bold('Seawater')")) +
  scale_x_date(date_break="months" , date_labels="%b %Y", limits=as.Date(c("2017-07-01", "2018-10-15")),
               expand=c(0, 0)) +
  scale_y_continuous(expand=expand_scale(mult=c(0, 0.33)), breaks=seq(0, 100, by=10)) +
  annotate("segment", x=as.Date("2017-07-01"), y=0, xend=as.Date("2017-07-01"), yend=100.5, color="black", size=0.75) +
  theme +
  theme(axis.text.x=element_text(size=12, color="black", angle=90,
                                 hjust=1, vjust=0.5, margin=margin(t=5.5, unit="pt")),
        axis.title.x=element_text(size=14, color="black"), plot.margin=unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
        plot.title=element_text(size=16, hjust=0.5))

fcym <- filter(plot, station=="FCyM") %>%
  ggplot() +
  geom_bar(aes(x=date, y=abundance, fill=taxon), stat="identity", colour="black", size=0.3, width=8) +
  scale_fill_manual(values=color, labels=names) + 
  geom_text(data=filter(whole, station=="FCyM"), aes(x=date, y=75, label=paste(format(round(abundance, 1), nsmall=1), "%")),
            family="Times", fontface="bold", angle=90, hjust=-1) +
  labs(x=NULL, y="%") +
  ggtitle(parse(text="bolditalic('Cymodocea nodosa')~bold('(Mixed)')")) +
  scale_x_date(date_break="months" , date_labels="%b %Y", limits=as.Date(c("2017-11-01", "2018-10-15")),
               expand=c(0, 0)) +
  scale_y_continuous(expand=expand_scale(mult=c(0, 0.33)), breaks=seq(0, 100, by=10)) +
  annotate("segment", x=as.Date("2017-11-01"), y=0, xend=as.Date("2017-11-01"), yend=100.5, color="black", size=0.75) +
  theme

fcam <- filter(plot, station=="FCaM") %>%
  ggplot() +
  geom_bar(aes(x=date, y=abundance, fill=taxon), stat="identity", colour="black", size=0.3, width=8) +
  scale_fill_manual(values=color, labels=names) + 
  geom_text(data=filter(whole, station=="FCaM"), aes(x=date, y=75, label=paste(format(round(abundance, 1), nsmall=1), "%")),
            family="Times", fontface="bold", angle=90, hjust=-1) +
  labs(x=NULL, y="%") +
  ggtitle(parse(text="bolditalic('Caulerpa cylindracea')~bold('(Mixed)')")) +
  scale_x_date(date_break="months" , date_labels="%b %Y", limits=as.Date(c("2017-11-01", "2018-10-15")),
               expand=c(0, 0)) +
  scale_y_continuous(expand=expand_scale(mult=c(0, 0.33)), breaks=seq(0, 100, by=10)) +
  annotate("segment", x=as.Date("2017-11-01"), y=0, xend=as.Date("2017-11-01"), yend=100.5, color="black", size=0.75) +
  theme

fca <- filter(plot, station=="FCa") %>%
  ggplot() +
  geom_bar(aes(x=date, y=abundance, fill=taxon), stat="identity", colour="black", size=0.3, width=8) +
  scale_fill_manual(values=color, labels=names) + 
  geom_text(data=filter(whole, station=="FCa"), aes(x=date, y=75, label=paste(format(round(abundance, 1), nsmall=1), "%")),
            family="Times", fontface="bold", angle=90, hjust=-1) +
  labs(x="Date", y="%") +
  ggtitle(parse(text="bolditalic('Caulerpa cylindracea')~bold('(Monospecific)')")) +
  scale_x_date(date_break="months" , date_labels="%b %Y", limits=as.Date(c("2017-11-01", "2018-10-15")),
               expand=c(0, 0)) +
  scale_y_continuous(expand=expand_scale(mult=c(0, 0.33)), breaks=seq(0, 100, by=10)) +
  annotate("segment", x=as.Date("2017-11-01"), y=0, xend=as.Date("2017-11-01"), yend=100.5, color="black", size=0.75) +
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
        legend.key.size=unit(0.5, "cm")) +
  guides(fill=guide_legend(ncol=1))
legend <- cowplot::get_legend(p1)

# Combining plots together and saving
p <- cowplot::ggdraw() +
  cowplot::draw_plot(f, x=0.068, y=0.708, width=0.932, height=0.288) +
  cowplot::draw_plot(fcym, x=0.317, y=0.545, width=0.683, height=0.228) +
  cowplot::draw_plot(fcam, x=0.317, y=0.317, width=0.683, height=0.228) +
  cowplot::draw_plot(fca, x=0.317, y=0, width=0.683, height=0.317) +
  cowplot::draw_plot(legend, x=0.120, y=0.280, width=0.100, height=0.200) +
  cowplot::draw_line(x=c(0.135, 0.135), y=c(0.285, 0.344), size=0.5) +
  cowplot::draw_label("SAR11 Clade", x=0.145, y=0.314, hjust=0,  fontfamily="Times", size=12)
ggsave("results/figures/alphaproteobacteria_bar_plot.jpg", p, width=210, height=297, units="mm")
