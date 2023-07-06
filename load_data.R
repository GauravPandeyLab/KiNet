library(tidyverse)
library(igraph)
library(visNetwork)

ksi_gene <- read.csv('data/ksi_gene.csv')
ksi_gene_t <- data.frame(ksi_gene)
k_to_group <- read.csv('data/k_to_group.csv')
ks_df <- read.csv('data/ks_df.csv')
ks_edge_list <- read.csv('data/ks_edge_list.csv')
ks_group <- read.csv('data/ks_group.csv')
group_colors <- c("#1f77b4","#aec7e8","#ff7f0e","#ffbb78","#2ca02c","#98df8a","#d62728","#ff9896","#9467bd","#c5b0d5","#808080")
group_names <- c("TK","OTHER","CAMK","CMGC","AGC","STE","TKL","ATYPICAL","CK1","RGC","Non-Kinase")

get_color <- function(name){
  gp = ks_group$Group[which(ks_group$KS == name)]
  group_colors[which(group_names == gp)]
}

get_group <- function(name){
  ks_group$Group[which(ks_group$KS == name)]
}

kinase_set <- unique(ksi_gene$Kinase)
substrate_set <- unique(ksi_gene$Substrate)

ks_set <- unique(ks_group$KS)

ks_id <- 1:length(ks_group$KS)
ks_node_id_map <- c(as.character(ks_group$id))


groups <- unique(ks_group$Group)
ks_ig <- graph_from_edgelist(as.matrix(ks_edge_list[,c("From", "To")]), directed = TRUE)

unip <- read_delim('data/uniprot.txt',delim="\t")
lapply(unip['Gene Names'], function(s){strsplit(s," ")}) -> x
unip$Primary <- unlist(map(x[[1]],1))

unip2 <- read_delim('data/uniprot2.tsv',delim="\t")

pathways_df <- read_csv('data/pathway_gene_sets.csv')
pathway_categories <- pathways_df$Category %>% unique


get_protein_info <- function(name){
  rec <- unip %>% filter(Primary==name)
  entry <- rec %>% pull('Entry')
  misc_names <- rec %>% pull('Gene Names') %>% strsplit(" ") %>% unlist
  s <- unip2 %>% filter(Entry == entry) %>% pull('Protein names')
  x1 <- c(1,gregexpr("\\(",s)[[1]]+1)
  x2 <- c(x1[2]-3,gregexpr("\\)",s)[[1]]-1)
  elems <- mapply(substr,rep(s,each=length(x1)),x1,x2,USE.NAMES=FALSE)
  outs = list(fullname=elems[1],uniprot=entry)
  elems <- elems[-1]
  if(substring(elems[1],1,3)=='EC ') {
    outs$ec <- substring(elems[1],4)
    elems <- elems[-1]
  }
  outs$remaining <- c(elems, misc_names) %>% unique %>% setdiff(c(name))
  return(outs)
}