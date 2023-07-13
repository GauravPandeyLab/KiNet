library(tidyverse)

all_edges <- read_csv('data/ksi.csv') %>% 
  rename(from=Kinase,to=Substrate)
colors <- read_csv('data/groups.csv')
all_nodes <- read_csv('data/proteins.csv') %>% 
  rename(id=UniProt) %>% 
  left_join(colors %>% rename(color=Color),by=join_by(Group == Group))
  
pathways_df <- read_csv('data/pathway_gene_sets.csv')
pathway_categories <- pathways_df$Category %>% unique