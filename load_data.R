library(tidyverse)


all_edges <- read_csv('data/ksi.csv') %>% 
  rename(from=Kinase,to=Substrate)
colors <- read_csv('data/groups.csv')

all_nodes <- read_csv('data/proteins.csv') %>% 
  rename(id=UniProt,group=Group) %>%
  left_join(colors %>% rename(color=Color),by=join_by(group == Group)) %>%
  select(id,everything()) %>%
  arrange(GeneName)

synonyms <- readLines('data/synonyms.txt')
#all_proteins <- all_nodes$GeneID
  
pathways_df <- read_csv('data/pathway_gene_sets.csv')

pathway_categories <- pathways_df$Category %>% unique
default_pathway_category <- "Environmental Information Processing"
pathway_choices <- pathways_df %>% split(pathways_df[['Category']])

get_pathway_genes <- function(categ,name){
  ch <- pathway_choices[[categ]] %>% filter(Term==name) %>% pull(Genes) %>% strsplit(split=" ") 
  intersect(ch[[1]],all_nodes$GeneName)
}
get_gene_name <- function(protein) {
  all_nodes %>% filter(id==protein) %>% pull(GeneName)
}
get_edges <- function(proteins) {
  E1 <- all_edges %>% filter(from %in% proteins)
  E2 <- all_edges %>% filter(to %in% proteins)
  return(rbind(E1,E2))
}
get_one_degree <- function(genes) {
  proteins <- all_nodes %>% filter(GeneName %in% genes) %>% pull(id)
  E1 <- all_edges %>% filter(from %in% proteins)
  E2 <- all_edges %>% filter(to %in% proteins)
  g=list()
  g$edges=rbind(E1,E2)
  ids <- c(E1 %>% pull(to), E2 %>% pull(from), proteins) %>% unique
  g$nodes <- all_nodes %>% filter(id %in% ids) %>% mutate(label=GeneName)
  return(g)
}
get_zero_degree <- function(genes) {
  proteins <- all_nodes %>% filter(GeneName %in% genes) %>% pull(id)
  E1 <- all_edges %>% filter(from %in% proteins) %>% filter(to %in% proteins)
  nodes <- all_nodes %>% filter(id %in% proteins) %>% mutate(label=GeneName)
  g = list(nodes=nodes,edges=E1)
  return(g)
}


get_pathway_graph <- function(genes,bool_one_degree,bool_remove_disconnected) {
  proteins <- all_nodes %>% filter(GeneName %in% genes) %>% pull(id)
  g <- if(bool_one_degree) get_one_degree(genes) else get_zero_degree(genes)
  if(bool_remove_disconnected) {
    nonself <- g$edges %>% filter(from != to)
    ids <- c(nonself %>% pull(from),nonself %>% pull(to))
    g$nodes <- g$nodes %>% filter(id %in% ids)
    g$edges <- g$edges %>% filter(from %in% ids) %>% filter(to %in% ids)
  }
  return(g)
}

get_legend_graph <- function(){
  nodes <- colors %>% mutate(id=Group) %>% rename(label=Group,color=Color) %>% select(id,everything())
  g <- list(nodes=nodes,edges=data.frame())
  return(g)
}

get_synonyms <- function(id) {
  pos <- which(synonyms==id)
  loc <- which(synonyms[pos:length(synonyms)]=="")[1]
  s <- synonyms[(pos+1):(pos+loc-2)]
  return(s)
}


render_gene_info <- function(geneName) {
  if (is.null(geneName)) {return('')}
  info <- all_nodes %>% filter(GeneName == geneName)
  L0 <- list(tags$h4(geneName))
  L1 <-  list(tags$h5('Name'),tags$b(info$FullName))
  L2 <- list(tags$h5('Synonyms'))
  for(x in get_synonyms(info$id)) {
    L2 <- c(L2,list(x,tags$br()))
  }
  unip_url <- 'https://www.uniprot.org/uniprotkb/%s/entry'
  
  L3 <- list(tags$h5('Links'))
  L3 <- c(L3,list('UniProt: ',tags$a(info$id,href= sprintf(unip_url,info$id)), tags$br()))
  if(!is.null(info$HGNC)){
    hgnc_url <- "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:%s"
    L3 <- c(L3,list('HGNC: ',tags$a(info$HGNC,href= sprintf(hgnc_url,info$HGNC)),tags$br()))
  }
  if(!is.null(info$GeneID)){
    gene_url <- "https://www.ncbi.nlm.nih.gov/gene/%s"
    L3 <- c(L3,list('NCBI Gene: ',tags$a(info$GeneID,href= sprintf(gene_url,info$GeneID)),tags$br()))
  }
  tagList(L0,L1,L2,L3)
}
