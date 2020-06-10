#################################################################################################################
# plot_rarefaction.R
# 
# A script to plot the rarefaction curve of each sample.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.groups.rarefaction
#               data/raw/metadata.csv
# Produces: results/figures/rarefaction.jpg
#
#################################################################################################################

# Loading input data and selection of values for plotting
rarefaction <- read_tsv(file="data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.groups.rarefaction") %>%
  select(-contains("lci-"), -contains("hci-")) %>%
  gather(-numsampled, key=sample, value=sobs) %>%
  mutate(sample=str_replace_all(sample, pattern="0.03-", replacement="")) %>%
  drop_na()

# Loading metadata
Sys.setlocale(locale="en_GB.utf8")
metadata <- read_tsv("data/raw/metadata.csv") %>%
  mutate(date=as.Date(date, "%d.%m.%Y")) %>%
  arrange(date) %>%
  mutate(date=format(date, "%d %B %Y")) %>%
  mutate(date=str_replace(date, "^0", "")) %>%
  mutate(date=if_else(label=="24/4-18 F-2", paste0(date, "a"), date))

# Joining metadata and input data
metadata_rarefaction <- inner_join(metadata, rarefaction, by=c("ID"="sample"))

# Defining line color and type for each sampling date
color_type <- tribble(~date, ~color, ~line_type,
                      "13 July 2017", "#A6CEE3", "solid",
                      "27 July 2017", "#1F78B4", "solid",
                      "10 August 2017", "#B2DF8A", "solid",
                      "24 August 2017", "#33A02C", "solid",
                      "19 September 2017", "#FB9A99", "solid",
                      "12 October 2017", "#E31A1C", "solid",
                      "23 November 2017", "#FDBF6F", "solid",
                      "4 December 2017", "#FF7F00", "solid",
                      "14 December 2017", "#CAB2D6", "solid",
                      "13 February 2018", "#6A3D9A", "solid",
                      "27 March 2018", "#A6CEE3", "dotted",
                      "24 April 2018", "#1F78B4", "dotted",
                      "24 April 2018a", "#B2DF8A", "dotted",
                      "22 May 2018", "#33A02C", "dotted",
                      "19 June 2018", "#FB9A99", "dotted",
                      "10 July 2018", "#E31A1C", "dotted",
                      "9 August 2018", "#FDBF6F", "dotted",
                      "4 September 2018", "#FF7F00", "dotted",
                      "5 October 2018", "#CAB2D6", "dotted")

# Generating a common theme for plots
theme <- theme(text=element_text(family="Times"), line=element_line(color="black"),
               panel.border=element_rect(fill=NA), panel.background=element_blank(),
               panel.grid=element_blank(), axis.line=element_blank(),
               axis.text=element_text(size=12, color="black"), axis.title=element_text(size=14, color="black"),
               plot.margin=unit(c(5.5, 16.5, 5.5, 5.5), "pt"), legend.position="none",
               plot.title=element_text(size=16, hjust=0.5))

# Plots generation
p1 <- filter(metadata_rarefaction, station=="F") %>%
  ggplot(aes(x=numsampled, y=sobs, group=ID, color=factor(date, levels=unique(date)),
             linetype=factor(date, levels=unique(date)))) +
  geom_line(size=1.5) +
  scale_colour_manual(values=set_names(color_type$color, color_type$date)) +
  scale_linetype_manual(values=set_names(color_type$line_type, color_type$date)) +
  labs(x="", y="Number of OTUs") +
  ggtitle(parse(text="bold('Seawater')")) +
  theme

p2 <- filter(metadata_rarefaction, str_detect(label, "FCyM$")) %>%
  ggplot(aes(x=numsampled, y=sobs, group=ID, color=factor(date, levels=unique(date)),
             linetype=factor(date, levels=unique(date)))) +
  geom_line(size=1.5) +
  scale_colour_manual(values=set_names(color_type$color, color_type$date)) +
  scale_linetype_manual(values=set_names(color_type$line_type, color_type$date)) +
  labs(x="", y="") +
  ggtitle(parse(text="bolditalic('Cymodocea nodosa')~bold('(Mixed)')")) +
  theme

p3 <- filter(metadata_rarefaction, str_detect(label, "FCaM$")) %>%
  ggplot(aes(x=numsampled, y=sobs, group=ID, color=factor(date, levels=unique(date)),
             linetype=factor(date, levels=unique(date)))) +
  geom_line(size=1.5) +
  scale_colour_manual(values=set_names(color_type$color, color_type$date)) +
  scale_linetype_manual(values=set_names(color_type$line_type, color_type$date)) +
  scale_x_continuous(limits=c(0, 40000)) +
  scale_y_continuous(limits=c(0, 4000)) +
  labs(x="Number of Sequences", y="Number of OTUs") +
  ggtitle(parse(text="bolditalic('Caulerpa cylindracea')~bold('(Mixed)')")) +
  theme

p4 <- filter(metadata_rarefaction, str_detect(label, "FCa$")) %>%
  ggplot(aes(x=numsampled, y=sobs, group=ID, color=factor(date, levels=unique(date)),
             linetype=factor(date, levels=unique(date)))) +
  geom_line(size=1.5) +
  scale_colour_manual(values=set_names(color_type$color, color_type$date)) +
  scale_linetype_manual(values=set_names(color_type$line_type, color_type$date)) +
  ggtitle(parse(text="bolditalic('Caulerpa cylindracea')~bold('(Monospecific)')")) +
  labs(x="Number of Sequences", y="") +
  theme

# Generating a plot with all samples to extract a common legend
p <- ggplot(metadata_rarefaction, aes(x=numsampled, y=sobs, group=ID,
            color=factor(date, levels=unique(date)), linetype=factor(date, levels=unique(date)))) +
  geom_line(size=1.5) +
  scale_colour_manual(values=set_names(color_type$color, color_type$date)) +
  scale_linetype_manual(values=set_names(color_type$line_type, color_type$date)) +
  scale_x_continuous(limits=c(0, max(metadata_rarefaction$numsampled))) +
  scale_y_continuous(limits=c(0, max(metadata_rarefaction$sobs))) +
  labs(x="Number of Sequences", y="") +
  theme +
  theme(legend.position="bottom", legend.title=element_blank(),
        legend.text=element_text(size=10, margin=margin(r=0.2, unit="cm")),
        legend.key.width=unit(1.4, "cm"), legend.key.height=unit(0.5, "cm"),
        legend.key=element_rect(fill="white"), legend.justification=c("top"),
        legend.text.align=0, legend.spacing.x=unit(0, "cm"))

# Common legend extraction
legend <- cowplot::get_legend(p)

# Combining plots together and save
plots <- cowplot::plot_grid(p1, p2, p3, p4, ncol=2, nrow=2, align="h")
p <- cowplot::plot_grid(plots, legend, ncol=1, nrow=2, rel_heights=c(4.6, 0.4))
ggsave("results/figures/rarefaction.jpg", p, width=210, height=297, units="mm")