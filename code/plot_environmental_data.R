#################################################################################################################
# plot_environmental_data.R
# 
# A script to plot environmental data.
# Dependencies: data/raw/environmental_data.csv
#               data/raw/metadata.csv
# Produces: results/figures/environmental_data.jpg
#
#################################################################################################################

# Loading input environmental data
env <- read_tsv("data/raw/environmental_data.csv") %>%
  mutate(ID=as.character(ID))

# Loading metadata
metadata <- read_tsv("data/raw/metadata.csv") %>%
  filter(ID!="23_1") %>%
  mutate(ID=str_replace(ID, "23_2", "23")) %>%
  mutate(label=str_replace(label, "24/4-18 F-2", "24/4-18 F"))

# Joining environmental data and metadata
Sys.setlocale(locale="en_GB.utf8")
env_metadata <- inner_join(metadata, env, by=c("ID"="ID")) %>%
  mutate(label=factor(label, levels=metadata$label)) %>%
  mutate(date=as.Date(date, "%d.%m.%Y"))

# Generating a common theme for plots
theme <- theme(text=element_text(family="Times"), line=element_line(color="black"),
               panel.border=element_blank(), panel.background=element_blank(),
               panel.grid=element_blank(), axis.line=element_line(color="gray60", size=0.2),
               axis.ticks=element_line(color="gray60", size=0.2),
               axis.text=element_text(size=12, color="black"), axis.text.x=element_text(angle=90, hjust=0.95, vjust=1.5),
               axis.title=element_text(size=12, color="black"),
               plot.margin=unit(c(0, 5.5, 0, 5.5), "pt"), legend.position="none",
               plot.title=element_text(size=16, hjust=0.5))

# Defining line and point size
line <- 0.3
point <- 1.2

# Generating plots
# Generating plot no. 1
# Temperature
temp_saline <- filter(env_metadata, station=="S") %>%
  ggplot(aes(x=date, y=`T`)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(8, 28), breaks=c(seq(8, 28, by=4)), expand=c(0, 0)) +
  labs(x="", y="Temperature (°C)") +
  ggtitle(parse(text="bold('Bay of Saline')")) +
  theme +
  theme(axis.text.x=element_blank(),
        plot.margin=unit(c(5.5, 5.5, 0, 5.5), "pt"))

temp_funtana <- filter(env_metadata, station=="F") %>%
  ggplot(aes(x=date, y=`T`)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(8, 28), breaks=c(seq(8, 28, by=4)), expand=c(0, 0)) +
  labs(x="", y="") +
  ggtitle(parse(text="bold('Bay of Funtana')")) +
  theme +
  theme(axis.text.x=element_blank(),
        plot.margin=unit(c(5.5, 5.5, 0, -5.5), "pt"))

# Salinity
sal_saline <- filter(env_metadata, station=="S") %>%
  ggplot(aes(x=date, y=S)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(34, 39), breaks=c(seq(34, 39, by=1)), expand=c(0, 0)) +
  labs(x="", y="Salinity") +
  ggtitle(parse(text="")) +
  theme +
  theme(axis.text.x=element_blank())

sal_funtana <- filter(env_metadata, station=="F") %>%
  ggplot(aes(x=date, y=S)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(34, 39), breaks=c(seq(34, 39, by=1)), expand=c(0, 0)) +
  labs(x="", y="") +
  ggtitle(parse(text="")) +
  theme +
  theme(axis.text.x=element_blank(),
        plot.margin=unit(c(0, 5.5, 0, -5.5), "pt"))

# Chlorophyll a
chla_saline <- filter(env_metadata, station=="S") %>%
  ggplot(aes(x=date, y=Chla)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(0.2, 0.9), breaks=c(seq(0.2, 0.9, by=0.1)), expand=c(0, 0)) +
  labs(x="", y=expression(paste("Chl ", italic("a "), "(", mu * g, " ", l^-1, ")"))) +
  ggtitle(parse(text="")) +
  theme +
  theme(axis.text.x=element_blank())

chla_funtana <- filter(env_metadata, station=="F") %>%
  ggplot(aes(x=date, y=Chla)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(0.2, 0.9), breaks=c(seq(0.2, 0.9, by=0.1)), expand=c(0, 0)) +
  labs(x="", y="") +
  ggtitle(parse(text="")) +
  theme +
  theme(axis.text.x=element_blank(),
        plot.margin=unit(c(0, 5.5, 0, -5.5), "pt"))

# Particulate matter
pm_saline <- filter(env_metadata, station=="S") %>%
  ggplot(aes(x=date, y=PM)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(2, 16), breaks=c(seq(2, 16, by=2)), expand=c(0, 0)) +
  labs(x="", y=expression(paste("Particulate matter ", "(", m * g, " ", l^-1, ")"))) +
  ggtitle(parse(text="")) +
  theme +
  theme(axis.text.x=element_blank())

pm_funtana <- filter(env_metadata, station=="F") %>%
  ggplot(aes(x=date, y=PM)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(2, 16), breaks=c(seq(2, 16, by=2)), expand=c(0, 0)) +
  labs(x="", y="") +
  ggtitle(parse(text="")) +
  theme +
  theme(axis.text.x=element_blank(),
        plot.margin=unit(c(0, 5.5, 0, -5.5), "pt"))

# Prokaryotic abundance
pa_saline <- filter(env_metadata, station=="S") %>%
  ggplot(aes(x=date, y=PA)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(2, 12), breaks=c(seq(2, 12, by=2)), expand=c(0, 0)) +
  labs(x="Date", y=expression(atop(paste("Prokaryotic abundance"), paste("(× ", 10^5, " cells ", ml^-1, ")")))) +
  ggtitle(parse(text="")) +
  theme

pa_funtana <- filter(env_metadata, station=="F") %>%
  ggplot(aes(x=date, y=PA)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(2, 12), breaks=c(seq(2, 12, by=2)), expand=c(0, 0)) +
  labs(x="Date", y="") +
  ggtitle(parse(text="")) +
  theme +
  theme(plot.margin=unit(c(0, 5.5, 0, -5.5), "pt"))

# Combining plots together and saving
saline <- cowplot::plot_grid(temp_saline,
                         sal_saline,
                         chla_saline,
                         pm_saline,
                         pa_saline,
                         ncol=1, nrow=5,
                         align="v", rel_heights=c(1.14, 1, 1, 1, 1.34))
funtana <- cowplot::plot_grid(temp_funtana,
                         sal_funtana,
                         chla_funtana,
                         pm_funtana,
                         pa_funtana,
                         ncol=1, nrow=5,
                         align="v", rel_heights=c(1.15, 1, 1, 1, 1.34))
p <- cowplot::plot_grid(saline, funtana, ncol=2, nrow=1, rel_widths=c(1.11, 1))
ggsave("results/figures/environmental_data_1.jpg", p, width=210, height=297, units="mm")

# Generating plot no. 2
# PO4
PO4_saline <- filter(env_metadata, station=="S") %>%
  ggplot(aes(x=date, y=PO4)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(0, 0.2), breaks=c(seq(0, 0.2, by=0.04)), expand=c(0, 0)) +
  labs(x="", y=expression(paste(PO[4]^{"3-"}, " (", mu, "M)"))) +
  ggtitle(parse(text="bold('Bay of Saline')")) +
  theme +
  theme(axis.text.x=element_blank(),
        plot.margin=unit(c(5.5, 5.5, 0, 5.5), "pt"))

PO4_funtana <- filter(env_metadata, station=="F") %>%
  ggplot(aes(x=date, y=PO4)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(0, 0.2), breaks=c(seq(0, 0.2, by=0.04)), expand=c(0, 0)) +
  labs(x="", y="") +
  ggtitle(parse(text="bold('Bay of Funtana')")) +
  theme +
  theme(axis.text.x=element_blank(),
        plot.margin=unit(c(5.5, 5.5, 0, -5.5), "pt"))

# NH4
NH4_saline <- filter(env_metadata, station=="S") %>%
  ggplot(aes(x=date, y=NH4)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(0, 2), breaks=c(seq(0, 2, by=0.4)), expand=c(0, 0)) +
  labs(x="", y=expression(paste(NH[4]^{"+"}, " (", mu, "M)"))) +
  ggtitle(parse(text="")) +
  theme +
  theme(axis.text.x=element_blank())

NH4_funtana <- filter(env_metadata, station=="F") %>%
  ggplot(aes(x=date, y=NH4)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(0, 2), breaks=c(seq(0, 2, by=0.4)), expand=c(0, 0)) +
  labs(x="", y="") +
  ggtitle(parse(text="")) +
  theme +
  theme(axis.text.x=element_blank(),
        plot.margin=unit(c(0, 5.5, 0, -5.5), "pt"))

# NO2
NO2_saline <- filter(env_metadata, station=="S") %>%
  ggplot(aes(x=date, y=NO2)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(0, 0.16), breaks=c(seq(0, 0.16, by=0.04)), expand=c(0, 0)) +
  labs(x="", y=expression(paste(NO[2]^textstyle("-"), " (", mu, "M)"))) +
  ggtitle(parse(text="")) +
  theme +
  theme(axis.text.x=element_blank())

NO2_funtana <- filter(env_metadata, station=="F") %>%
  ggplot(aes(x=date, y=NO2)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(0, 0.16), breaks=c(seq(0, 0.16, by=0.04)), expand=c(0, 0)) +
  labs(x="", y="") +
  ggtitle(parse(text="")) +
  theme +
  theme( axis.text.x=element_blank(),
        plot.margin=unit(c(0, 5.5, 0, -5.5), "pt"))

# NO3
NO3_saline <- filter(env_metadata, station=="S") %>%
  ggplot(aes(x=date, y=NO3)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(0, 8), breaks=c(seq(0, 8, by=2)), expand=c(0, 0)) +
  labs(x="", y=expression(paste(NO[3]^{textstyle("-")}, " (", mu, "M)"))) +
  ggtitle(parse(text="")) +
  theme +
  theme(axis.text.x=element_blank())

NO3_funtana <- filter(env_metadata, station=="F") %>%
  ggplot(aes(x=date, y=NO3)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(0, 8), breaks=c(seq(0, 8, by=2)), expand=c(0, 0)) +
  labs(x="", y="") +
  ggtitle(parse(text="")) +
  theme +
  theme(axis.text.x=element_blank(),
        plot.margin=unit(c(0, 5.5, 0, -5.5), "pt"))

# TIN
TIN_saline <- filter(env_metadata, station=="S") %>%
  ggplot(aes(x=date, y=TIN)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(0, 10), breaks=c(seq(0, 10, by=2)), expand=c(0, 0)) +
  labs(x="", y=expression(paste("TIN", " (", mu, "M)"))) +
  ggtitle(parse(text="")) +
  theme +
  theme(axis.text.x=element_blank())

TIN_funtana <- filter(env_metadata, station=="F") %>%
  ggplot(aes(x=date, y=TIN)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(0, 10), breaks=c(seq(0, 10, by=2)), expand=c(0, 0)) +
  labs(x="", y="") +
  ggtitle(parse(text="")) +
  theme +
  theme(axis.text.x=element_blank(),
        plot.margin=unit(c(0, 5.5, 0, -5.5), "pt"))

# SiO4
SiO4_saline <- filter(env_metadata, station=="S") %>%
  ggplot(aes(x=date, y=SiO4)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(0, 10), breaks=c(seq(0, 10, by=2)), expand=c(0, 0)) +
  labs(x="Date", y=expression(paste(SiO[4]^{"4-"}, " (", mu, "M)"))) +
  ggtitle(parse(text="")) +
  theme

SiO4_funtana <- filter(env_metadata, station=="F") %>%
  ggplot(aes(x=date, y=SiO4)) +
  geom_line(size=line) +
  geom_point(size=point) +
  scale_x_date(breaks=seq(as.Date("2017-07-01"), as.Date("2018-11-01"), "months"),
               labels=c("Jul 2017", "Aug 2017", "Sep 2017", "Oct 2017", "Nov 2017", "Dec 2017",
                        "Jan 2018", "Feb 2018", "Mar 2018", "Apr 2018", "May 2018", "Jun 2018",
                        "Jul 2018", "Aug 2018", "Sep 2018", "Oct 2018", ""),
               limits=as.Date(c("2017-07-01", "2018-11-01")),
               expand=c(0, 0)) +
  scale_y_continuous(limits=c(0, 10), breaks=c(seq(0, 10, by=2)), expand=c(0, 0)) +
  labs(x="Date", y="") +
  ggtitle(parse(text="")) +
  theme +
  theme(plot.margin=unit(c(0, 5.5, 0, -5.5), "pt"))

# Combining plots together and saving
saline <- cowplot::plot_grid(PO4_saline,
                         NH4_saline,
                         NO2_saline,
                         NO3_saline,
                         TIN_saline,
                         SiO4_saline,
                         ncol=1, nrow=6,
                         align="v", rel_heights=c(1.14, 1, 1, 1, 1, 1.4))
funtana <- cowplot::plot_grid(PO4_funtana,
                         NH4_funtana,
                         NO2_funtana,
                         NO3_funtana,
                         TIN_funtana,
                         SiO4_funtana,
                         ncol=1, nrow=6,
                         align="v", rel_heights=c(1.14, 1, 1, 1, 1, 1.4))
p <- cowplot::plot_grid(saline, funtana, ncol=2, nrow=1, rel_widths=c(1.06, 1))
ggsave("results/figures/environmental_data_2.jpg", p, width=210, height=297, units="mm")