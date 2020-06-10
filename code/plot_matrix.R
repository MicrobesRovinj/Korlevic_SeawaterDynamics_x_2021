#################################################################################################################
# code/plot_matrix.R
# 
# A script to plot a matrix of similarity indices.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.pick.pick.opti_mcc.shared
#               data/raw/metadata.csv
# Produces: results/figures/matrix.jpg
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

# Loading metadata 
metadata <- read_tsv("data/raw/metadata.csv")

# Joining metadata with OTU/sample data and summing sequences from each environment
rarefied_metadata <- inner_join(rarefied, metadata, by=c("Group"="ID")) %>%
  mutate(date=as.Date(date, "%d.%m.%Y")) %>%
  filter(date >= "2017-11-01") %>%
  select_if(funs(!is.numeric(.) || sum(.)!=0)) %>%
  group_by(station) %>%
  summarise_at(vars(starts_with("Otu")), sum)

# Copying the sample labels to the rows (input for library vegan)
row.names(rarefied_metadata) <- rarefied_metadata$station

# Removing column containing sample labels
rarefied_metadata <- rarefied_metadata %>%
  select(-station)

# Calculating matrices of similarity indices
jaccard <- vegdist(rarefied_metadata, method="jaccard", binary=T)
jaccard <- as_tibble(data.frame(t(combn(rownames(rarefied_metadata), 2)), as.numeric(jaccard)),
                     .name_repair= ~c("V1", "V2", "jaccard"))

bray <- vegdist(rarefied_metadata, method="bray", binary=F)
bray <- as_tibble(data.frame(t(combn(rownames(rarefied_metadata), 2)), as.numeric(bray)),
                  .name_repair= ~c("V1", "V2", "bray"))

similarity <- inner_join(jaccard, bray, by=c("V1"="V1", "V2"="V2")) %>%
  mutate(jaccard=1-jaccard) %>%
  mutate(bray=1-bray)

# Arrange factors for plotting
similarity$V1 <- factor(similarity$V1, levels=c("F", "FCyM", "FCa", "FCaM"))
similarity$V2 <- factor(similarity$V2, levels=c("FCyM", "FCaM", "FCa"))

# Function for setting the number of decimal places
scaleFUN <- function(x) sprintf("%.2f", x)

# Generate a common theme
theme <- theme(text=element_text(family="Times"), line=element_line(color="black"),
         panel.grid=element_blank(),
         axis.line.x=element_blank(), axis.line.y=element_blank(),
         axis.ticks.x=element_blank(), axis.ticks.y.left=element_blank(),
         axis.text.y=element_text(size=15, color="black", hjust=0.5, vjust=0.5, margin=margin(r=-0.2, unit="cm")),
         axis.text.x=element_text(size=15, color="black", hjust=0.5, vjust=0.5, margin=margin(t=-0.2, unit="cm")),
         panel.background=element_blank(), plot.margin=unit(c(5.5, 5.5, 5.5, 5.5), "pt"),
         legend.position="none", plot.title=element_text(face="bold", size=20))

# Plot generation
group_name <- c("F"=expression(bold("Seawater")),
                "FCyM"=expression(atop(bolditalic("Cymodocea nodosa"), bold("(Mixed)"))),
                "FCaM"=expression(atop(bolditalic("Caulerpa cylindracea"), bold("(Mixed)"))),
                "FCa"=expression(atop(bolditalic("Caulerpa cylindracea"), bold("(Monospecific)"))))

jaccard <- ggplot(similarity, aes(x=V2, y=V1, fill=jaccard)) + 
  geom_tile(colour="black", size=0.5) +
  geom_text(aes(V2, V1, label=scaleFUN(jaccard)), color="black", size=7, family="Times", fontface="bold") +
  labs(x=NULL, y=NULL) +
  ggtitle("Jaccard's Similarity Coefficient") +
  scale_fill_gradient2(low="#F7FCB9", high= "#ADDD8E", mid= "#31A354", 
                       midpoint=0.5, limit=c(0, 1)) +
  scale_x_discrete(labels=NULL) +
  scale_y_discrete(labels=group_name) +
  theme +
  theme(plot.margin=unit(c(5.5, 5.5, 33, 5.5), "pt"))

bray <- ggplot(similarity, aes(x=V2, y=V1, fill=bray)) + 
  geom_tile(colour="black", size=0.5) +
  geom_text(aes(V2, V1, label=scaleFUN(bray)), color="black", size=7, family="Times", fontface="bold") +
  labs(x=NULL, y=NULL) +
  ggtitle("Bray-Curtis Similarity Coefficient") +
  scale_fill_gradient2(low="#F7FCB9", high= "#ADDD8E", mid= "#31A354", 
                       midpoint=0.5, limit=c(0, 1)) +
  scale_x_discrete(labels=group_name) +
  scale_y_discrete(labels=group_name) +
  theme

# Combining plots together and saving
plots <- cowplot::plot_grid(jaccard, bray, ncol=1, nrow=2, rel_heights=c(1, 1))
ggsave("results/figures/matrix.jpg", plots, width=210, height=297, units="mm")