library(flexdashboard)
library(ggplot2)
library(reshape2)
library(scales)
library(tidyverse)
library(plotly)
library(FDRestimation)
library(ggrepel)
library(hrbrthemes)
library(dplyr)
library(shiny)
library(jsonlite)
library(DT)
library(patchwork)
library(network)
library(sna)
library(igraph)
library(ggraph)
library(tidygraph)
library(visNetwork)
library(stringr)

ksi_gene <- read.csv('data/ksi_gene.csv')

ksi_gene_t <- data.frame(ksi_gene)
# print("test")
colnames(ksi_gene_t) <- c('From', 'To', 'Group')

k_to_group <- data.frame(ksi_gene[,c('Kinase', 'Group')])
k_to_group <- k_to_group %>% distinct(Kinase, Group, .keep_all = TRUE)



kinase_set <- unique(ksi_gene$Kinase)
substrate_set <- unique(ksi_gene$Substrate)

ks_set <- union(kinase_set, substrate_set)

ks_df <- data.frame('KS'=ks_set)

colnames(k_to_group)[which(names(k_to_group) == 'Kinase')] <- 'KS'
ks_group <- merge(x = ks_df, y = k_to_group, by = "KS",
                  all = TRUE)

ks_group$Group[is.na(ks_group$Group)] <- 'Non-Kinase'
ks_id <- 1:length(ks_group$KS)
ks_group$id <- ks_id




ks_edge_list <- data.frame(ksi_gene_t[,c('From', 'To')])
ks_node_id_map <- c(as.character(ks_group$id))
names(ks_node_id_map) <- paste0("^", ks_group$KS, "$")
ks_edge_list$from <- as.integer(str_replace_all(string = ks_edge_list$From,
                                                pattern= ks_node_id_map))

ks_edge_list$to <- as.integer(str_replace_all(string = ks_edge_list$To,
                                              pattern= ks_node_id_map))