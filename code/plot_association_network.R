#################################################################################################################
# plot_association_network.R
# 
# A script to generate the association network figure.
# Dependencies: data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v138.wang.tax.summary
#               data/raw/environmental_data.csv
# Produces: results/figures/association_network.jpg
#
#################################################################################################################

# Loading input data containing sequence abundances and subsequent input data customization
community <- read_tsv("data/mothur/raw.trim.contigs.good.unique.good.filter.unique.precluster.pick.nr_v138.wang.tax.summary") %>%
  filter(!str_detect(taxon, "^Eukaryota")) %>%
  filter(taxon!="Root")

# Removing chloroplast and mitochondrial sequences and subtracting their number from higher taxonomic levels to which
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

# Selecting the taxonomic level for plotting
select <- filter(community, taxlevel==6) %>%
  filter_at(6:ncol(.), any_vars(. >= 1)) %>%
  mutate_at(5:ncol(.), list(~. / sum(.) * 100))

# Removing the last digit from rankID for the next step
select <- select %>%
  mutate(rankID=if_else(taxon=="uncultured", str_replace(rankID, "\\.\\d+$", ""), rankID)) %>%
  mutate(rankID=if_else(taxon=="uncultured_ge", str_replace(rankID, "\\.\\d+\\.\\d+$", ""), rankID))

# Adding data to "uncultured" or "uncultured_ge" taxa describing the higher taxonomic level to which they belong
uncultured <- select(filter(community, rankID %in%
  filter(select, taxon=="uncultured" | taxon=="uncultured_ge")$rankID), rankID, taxon) %>%
  rename(rankID_uncultured=rankID, taxon_uncultured=taxon)

select <- left_join(select, uncultured, by=c("rankID"="rankID_uncultured")) %>%
  mutate(taxon=if_else(taxon=="uncultured", paste0(taxon, "_", taxon_uncultured), taxon)) %>%
  mutate(taxon=if_else(taxon=="uncultured_ge",
                       paste0(str_replace(taxon, "_ge$", ""), "_", taxon_uncultured), taxon)) %>%
  select(-taxon_uncultured)

# Transposing community data
select <- select %>%
  select(-taxlevel, -rankID, -daughterlevels, -total) %>%
  gather(key=ID, value=abundance, 2:ncol(.)) %>% 
  spread(key=names(.)[1], value="abundance")

# Loading input environmental data
env <- read_tsv("data/raw/environmental_data.csv") %>%
  mutate(ID=as.character(ID))

# Combining environmental and community data
env_select <- inner_join(env, select, by="ID")

# Deleting column containing sample labels
env_select <- env_select %>%
  select(-ID) %>%
  as.matrix()

# Calculating Spearman's rho rank correlation coefficients
rcorr <- rcorr(env_select, type="spearman")

# Extracting R-values
r_values <- rcorr$r[lower.tri(rcorr$r)]
r_values <- data.frame(t(combn(rownames(rcorr$r), 2))) %>%
  rename(V1=X1, V2=X2) %>%
  tibble(r_values=r_values)  
                         
# Extracting p-values
p_values <- rcorr$P[lower.tri(rcorr$P)]
p_values <- data.frame(t(combn(rownames(rcorr$P), 2))) %>%
  rename(V1=X1, V2=X2) %>%
  tibble(p_values=p_values)

# Calculating q-values
gobj <- qvalue(p_values$p_values)
q_values <- gobj$qvalues

# Combining R-values, p-values and q-values
r_p_q_values <- tibble(r_values, p_values=p_values$p_values, q_values)

# Selecting correlations for plotting
r_values_select <- r_p_q_values %>%
  filter((r_values > 0.80 | r_values < -0.80) &
          p_values < 0.007 &
          q_values < 0.001) %>%
  select(V1, V2, r_values)

# Creating a node list
V1 <- r_values_select %>%
  distinct(V1) %>%
  rename(label=V1)

V2 <- r_values_select %>%
  distinct(V2) %>%
  rename(label=V2)

nodes <- full_join(V1, V2, by="label") %>%
  rowid_to_column("id") %>%
  mutate(env_taxa=if_else(label %in% colnames(env), "env", "taxa"))

# Generating italic names for taxonomic groups
nodes <- nodes %>%
  mutate(names=case_when(           nodes$label=="SAR116_clade_unclassified" ~ "atop(plain('SAR116 Clade'), plain('(NR)'))",
                         str_detect(nodes$label, "uncultured") ~ paste0("atop(plain('Uncultured'), italic('", str_remove(nodes$label, "uncultured_"), "'))"),
                         str_detect(nodes$label, "unclassified") ~ paste0("atop(italic('", str_remove(nodes$label, "_unclassified"), "'), plain('(NR)'))"),
                                    nodes$label=="OCS116_clade_ge" ~ "atop(plain('OCS116'), plain('Clade'))",
                                    nodes$label=="HIMB11" ~ "plain('HIMB11')",
                                    nodes$label=="Clade_Ia" ~ "atop(plain('SAR11 Clade'), plain('(Clade Ia)'))",
                                    nodes$label=="Clade_Ib" ~ "atop(plain('SAR11 Clade'), plain('(Clade Ib)'))",
                                    nodes$label=="Clade_III_ge" ~ "atop(plain('SAR11 Clade'), plain('(Clade III)'))",
                                    nodes$label=="Marinimicrobia_(SAR406_clade)_ge" ~ "italic('Marinimicrobia')",
                                    nodes$label=="PS1_clade_ge" ~ "plain('PS1 Clade')",
                                    nodes$label=="SUP05_cluster" ~ "atop(plain('SUP05'), plain('Cluster'))",
                                    nodes$label=="Marine_Group_II_ge" ~ "atop(plain('Marine'), plain('Group II'))",
                                    nodes$label=="Candidatus_Nitrosopumilus" ~ "atop(plain('\"')*italic('Candidatus'), plain('Nitrosopumilus\"'))",
                                    nodes$label=="Candidatus_Actinomarina" ~ "atop(plain('\"')*italic('Candidatus'), plain('Actinomarina\"'))",
                                    nodes$label=="OM75_clade" ~ "atop(plain('OM75'), plain('Clade'))",
                                    nodes$label=="OM182_clade_ge" ~ "atop(plain('OM182'), plain('Clade'))",           
                                    nodes$label=="T" ~ "plain('T')",
                                    nodes$label=="TIN" ~ "plain('TIN')",
                                    nodes$label=="NO3" ~ "plain('NO'['3']^textstyle('-'))",
                                    TRUE ~ paste0("italic('", nodes$label, "')")))

# Creating an edge list
edges <- r_values_select %>%
  left_join(nodes, by=c("V1"="label")) %>%
  rename(from=id)

edges <- edges %>% 
  left_join(nodes, by=c("V2"="label")) %>% 
  rename(to=id)

edges <- edges %>%
  select(from, to, r_values) %>%
  mutate(corr=if_else(r_values > 0, "POSITIVE", "NEGATIVE")) %>%
  mutate(r_values, r_values=format(r_values, digits=2, nsmall=2))
  
# Creating the network object
network <- tbl_graph(nodes=nodes, edges=edges, directed=FALSE)

# Setting seed value for reproducibility
set.seed(19760620)

# Plotting the association network
p <- ggraph(network, layout="fr") +
  geom_edge_link(aes(label=r_values, edge_linetype=corr), color="gray50", label_size=4,
                 label_colour="black", family="Times", show.legend=FALSE,
                 check_overlap=TRUE) +
  scale_edge_linetype_manual(values=c("POSITIVE"="solid", "NEGATIVE"="dotted")) +
  geom_node_point(aes(fill=env_taxa), shape=21, size=22, color="gray50", show.legend=FALSE) +
  scale_fill_manual(values=c("env"="gray30", "taxa"="gray")) +
  geom_node_text(aes(label=names, color=env_taxa), size=4, family="Times", show.legend=FALSE,
                 parse=TRUE, check_overlap=TRUE) +
  scale_color_manual(values=c("env"="white", "taxa"="black"), labels=names) +
  theme_graph()

# Customizing node position
coord <- p$data %>%
  mutate(x=if_else(label=="Litoricola", x+0.5, x), y=if_else(label=="Litoricola", y+0.7, y)) %>%
  mutate(x=if_else(label=="Marinimicrobia_(SAR406_clade)_ge", x-0.1, x), y=if_else(label=="Marinimicrobia_(SAR406_clade)_ge", y-0.3, y)) %>%
  mutate(x=if_else(label=="Clade_Ib", x-0.6, x), y=if_else(label=="Clade_Ib", y, y)) %>%
  mutate(x=if_else(label=="Marine_Group_II_ge", x+1.2, x), y=if_else(label=="Marine_Group_II_ge", y, y)) %>%
  mutate(x=if_else(label=="Clade_Ia", x+0.9, x), y=if_else(label=="Clade_Ia", y+0.5, y)) %>%
  mutate(x=if_else(label=="SUP05_cluster", x+0.6, x), y=if_else(label=="SUP05_cluster", y-0.8, y)) %>%
  mutate(x=if_else(label=="uncultured_Ectothiorhodospiraceae", x+0.2, x), y=if_else(label=="uncultured_Ectothiorhodospiraceae", y-2.2, y)) %>%
  mutate(x=if_else(label=="OCS116_clade_ge", x, x), y=if_else(label=="OCS116_clade_ge", y+0.2, y)) %>%
  mutate(x=if_else(label=="Microbacteriaceae_unclassified", x+0.7, x), y=if_else(label=="Microbacteriaceae_unclassified", y+0.4, y)) %>%
  mutate(x=if_else(label=="OM75_clade", x+0.5, x), y=if_else(label=="OM75_clade", y-0.8, y)) %>%
  mutate(x=if_else(label=="OM182_clade_ge", x+0.8, x), y=if_else(label=="OM182_clade_ge", y-1.4, y)) %>%
  mutate(x=if_else(label=="Candidatus_Actinomarina", x, x), y=if_else(label=="Candidatus_Actinomarina", y-2.6, y)) %>%
  mutate(x=if_else(label=="uncultured_Defluviicoccales", x-0.9, x), y=if_else(label=="uncultured_Defluviicoccales", y-1.6, y)) %>%
  mutate(x=if_else(label=="NO3", x-0.9, x), y=if_else(label=="NO3", y-1.5, y)) %>%
  mutate(x=if_else(label=="TIN", x-0.6, x), y=if_else(label=="TIN", y, y)) %>%
  mutate(x=if_else(label=="Luteolibacter", x+1.8, x), y=if_else(label=="Luteolibacter", y+3.2,  y)) %>%
  mutate(x=if_else(label=="T", x-0.2, x), y=if_else(label=="T", y+1.0, y)) %>%
  mutate(x=if_else(label=="Balneola", x+0.3, x), y=if_else(label=="Balneola", y-0.2, y)) %>%
  mutate(x=if_else(label=="SAR116_clade_unclassified", x+0.4, x), y=if_else(label=="SAR116_clade_unclassified", y+0.2, y)) %>%
  mutate(x=if_else(label=="Clade_III_ge", x+0.2, x), y=if_else(label=="Clade_III_ge", y+0.3, y)) %>%
  mutate(x=if_else(label=="uncultured_Balneolaceae", x-0.3, x), y=if_else(label=="uncultured_Balneolaceae", y+1.4, y)) %>%
  mutate(x=if_else(label=="Planktomarina", x-0.8, x), y=if_else(label=="Planktomarina", y-0.6, y)) %>%
  mutate(x=if_else(label=="PS1_clade_ge", x-0.3, x), y=if_else(label=="PS1_clade_ge", y+1.8, y)) %>%
  mutate(x=if_else(label=="Candidatus_Nitrosopumilus", x+0.2, x), y=if_else(label=="Candidatus_Nitrosopumilus", y+0.3, y)) %>%
  mutate(x=if_else(label=="uncultured_Flavobacteriaceae", x-0.9, x), y=if_else(label=="uncultured_Flavobacteriaceae", y+1.3, y)) %>%
  mutate(x=if_else(label=="HIMB11", x-0.6, x), y=if_else(label=="HIMB11", y, y))
p$data <- coord

# Saving the association network
ggsave("results/figures/association_network.jpg", p, width=297, height=210, units="mm")