library(tidyverse)
library(igraph)


all_edges <- read_csv('data/ksi_display.csv') %>% 
  rename(from=Kinase,to=Substrate)
#widths <- data.frame(NSites =all_edges$NSites %>% unique %>% sort) %>% mutate(width=3*ntile(NSites,10))
widths <- data.frame(NSites =all_edges$NSites %>% unique %>% sort) %>% mutate(width=5)
all_edges <- merge(x=all_edges,y=widths,by="NSites",all.x=TRUE)
colors <- read_csv('data/groups.csv')

all_nodes <- read_csv('data/proteins.csv') %>% 
  rename(id=UniProt,group=Group) %>%
  left_join(colors %>% rename(color=Color),by=join_by(group == Group)) %>% 
  select(id,everything()) %>%
  arrange(GeneName) %>%
  mutate(size=20)

all_kinases <- all_nodes %>% filter(group != 'Non-kinase') %>% pull(id) %>% unique
n_kinases <- all_kinases %>% length
n_nonkinases <- setdiff(all_edges %>% pull(to) %>% unique, all_kinases) %>% length
n_interactions <- all_edges %>% nrow
n_kki <- all_edges %>% filter(to %in% all_kinases) %>% filter(from!=to)  %>% nrow
n_auto <-all_edges %>% filter(from==to)  %>% nrow

source_data <- read_csv('data/ksi_source.csv')
####### PATHWAY GENE SETS
pathways_df <- read_csv('data/pathway_gene_sets.csv')
pathway_categories <- pathways_df$Category %>% unique
default_pathway_category <- "Environmental Information Processing"
pathway_choices <- pathways_df %>% split(pathways_df[['Category']])
pathways2 <- pathways_df %>% select(Category,Term)
all_pathways <- lapply(split(pathways2, pathways2$Category), function(d) {d$Term})
get_pathway_genes <- function(pathway) {
  ch <- pathways_df %>% filter(Term==pathway) %>% pull(Genes) %>% strsplit(split=" ")
  intersect(ch[[1]],all_nodes$GeneName)
}
####### DOMAIN GENE SETS
domains_df <- read_csv('data/domain_gene_sets.csv')
domain_categories <- domains_df$Category %>% unique
default_domain_category <- "Cytoplasmic"
domain_choices <- domains_df %>% split(domains_df[['Category']])
domains2 <- domains_df %>% select(Category,Term)
all_domains <- lapply(split(domains2, domains2$Category), function(d) {d$Term})
get_domain_genes <- function(domain) {
  ch <- domains_df %>% filter(Term==domain) %>% pull(Genes) %>% strsplit(split=" ")
  intersect(ch[[1]],all_nodes$GeneName)
}

######### DISCLAIMERS
last_date = as.Date('2023-10-04',"%Y-%m-%d")
dataset_last_update = format(last_date,"%B %d, %Y")
nproteins = nrow(all_nodes)
nproteins_all = "20423"
uniprot_release = 'release 2023_04'
npathways = nrow(pathways_df)
npathways_all = "320"
kegg_name = "KEGG_2021_Human"
gseapy_release = "v1.0.6"
interpro_name = "InterPro_Domains_2019"
ndomains = nrow(domains_df)
ndomains_all = "1071"

######### DEFAULTS
get_gene_name <- function(protein) {
  all_nodes %>% filter(id==protein) %>% pull(GeneName)
}
get_edges <- function(proteins) {
  E1 <- all_edges %>% filter(from %in% proteins)
  E2 <- all_edges %>% filter(to %in% proteins)
  return(rbind(E1,E2))
}
get_one_degree <- function(genes,largeCenterNode=FALSE) {
  proteins <- all_nodes %>% filter(GeneName %in% genes) %>% pull(id)
  E1 <- all_edges %>% filter(from %in% proteins)
  E2 <- all_edges %>% filter(to %in% proteins)
  g=list()
  g$edges=rbind(E1,E2) %>% mutate(id=row_number())
  ids <- c(E1 %>% pull(to), E2 %>% pull(from), proteins) %>% unique
  g$nodes <- all_nodes %>% filter(id %in% ids) %>% mutate(label=GeneName)
  g$nodes$id <- make.unique(g$nodes$id, sep = "_")
  g$value <- 2
  if(largeCenterNode) {
    condition <- g$nodes$label %in% genes
    #g$nodes$size[condition] <- 50
    g$nodes$shape[condition] <- "circle"
    g$nodes$value[condition] <- 5
  }
  return(g)
}
get_zero_degree <- function(genes) {
  proteins <- all_nodes %>% filter(GeneName %in% genes) %>% pull(id)
  E1 <- all_edges %>% filter(from %in% proteins) %>% filter(to %in% proteins)
  nodes <- all_nodes %>% filter(id %in% proteins) %>% mutate(label=GeneName)
  g = list(nodes=nodes,edges=E1)
  return(g)
}


get_pathway_graph <- function(genes,bool_one_degree,bool_include_disconnected) {
  proteins <- all_nodes %>% filter(GeneName %in% genes) %>% pull(id)
  g <- if(bool_one_degree) get_one_degree(genes) else get_zero_degree(genes)
  g$nodes$id <- make.unique(g$nodes$id, sep = "_")
  if(bool_include_disconnected==FALSE) {
    nonself <- g$edges %>% filter(from != to)
    ids <- c(nonself %>% pull(from),nonself %>% pull(to))
    g$nodes <- g$nodes %>% filter(id %in% ids)
    g$edges <- g$edges %>% filter(from %in% ids) %>% filter(to %in% ids) 
  }
  g$edges <- g$edges %>% mutate(id=row_number())
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

render_line_break_separated <- function(elems){
  L <- list()
  for(elem in elems) {
    L <- c(L,list(elem,tags$br()))
  }
  return(L)
}

render_edge_link <- function(idx) {
  info <- all_nodes %>% filter(id==idx)
  unip_url <- 'https://www.uniprot.org/uniprotkb/%s/entry'
  L1 <- tags$a(info$GeneName,href= sprintf(unip_url,info$id))
  list(L1,tags$p(info$FullName))
}

epsdrefs <- read_csv('data/EPSDReference.csv')
render_edge_info <- function(edge) {
  if (is.null(edge)) {return('')}
  kinase <- edge %>% pull(from)
  substrate <- edge %>% pull(to)
  L1 <- tags$div(tags$h4('Kinase ',style="display: inline;"),render_edge_link(kinase))
  L2 <- tags$div(tags$h4('Substrate ',style="display: inline;"),render_edge_link(substrate))
  nsites <- edge %>% pull(NSites)
  if(nsites==1) {
    L3 <- list(tags$h4('Sites'),tags$p('1 site known'))
    } else {
      L3 <- list(tags$h4('Sites'),tags$p(nsites,'sites known'))
    }
  sites <- edge %>% pull(Sites)
  L4 <- tags$p(gsub('"','',sites))
  ppp_url <- 'https://www.phosphosite.org/uniprotAccAction?id=%s'
  L50 <- tags$h4('Sources')
  L51 <- tags$a('PhosphoSitePlus',href= sprintf(ppp_url,substrate))
  ipt_url <- 'https://research.bioinformatics.udel.edu/iptmnet/entry/%s/'
  L52 <- tags$a('iPTMNet',href= sprintf(ipt_url,substrate))
  
  epsd <- epsdrefs %>% filter(Substrate==substrate) %>% pull('EPSD')
  L53 <- ''
  if(length(epsd)>=1) {
    epsd_url <- 'https://epsd.biocuckoo.cn/View.php?id=%s'
    L53 <- tags$a('EPSD',href=sprintf(epsd_url,epsd[[1]]))
  }
  Lbr <- tags$br()
  L5 <- tags$div(L50,L51,Lbr,L52,Lbr,L53)
  
  tagList(L1,L2,L3,L4,L5)
}

get_igraph <-function(g) {
  nodes <- g$nodes
  edges <- g$edges
  n <- nodes %>% mutate(name=label) %>% select(name,label,color,group)
  e <- data.frame(from=sapply(edges %>% pull(from),get_gene_name),to=sapply(edges %>% pull(to),get_gene_name))
  h<- graph.data.frame(e,directed=TRUE,vertices=n)
  E(h)$color <- tail_of(h,E(h))$color
  return(h)
}

get_gml <- function(h) {
  ndf <- data.frame(id = seq(1,length(V(h))),label = V(h)$label,fill = V(h)$color)
  edf <- data.frame(source=head_of(h,E(h)), target=tail_of(h,E(h)), color=E(h)$color)
  nodelines <- do.call("sprintf",c('node [ id %s label "%s" graphics [ fill "%s" ] ]', ndf))
  edgelines <- do.call("sprintf",c('edge [ source %s target %s graphics [ fill "%s" targetArrow "standard"] ]',edf))
  lines <- c("graph","[","directed 1",nodelines,edgelines,"]")
  return(lines)
} 
####### VIS
vis_default <- function(g) {
  visNetwork(nodes=g$nodes,edges=g$edges,physics=T) %>%
  #visIgraphLayout(layout="layout_with_graphopt",charge=0.01,mass=100) %>% 
  #visIgraphLayout(layout="layout_in_circle") %>%
  visNodes(font=list(size=50),opacity=0.8) %>% 
  visEdges(arrows="to",smooth = list(enabled = T, type = 'dynamic'),color=list(opacity=0.5),selectionWidth=10) %>%
  visOptions(highlightNearest = list(enabled=T,degree=list(from=1,to=1),algorithm="hierarchical")) %>% 
  visExport(type="png",label="Screenshot visible region as PNG")
}

layout_choices <- list("Default"='layout_nicely',
                       "Circle"='layout_in_circle',
                       "Grid"='layout_on_grid')
#visOptions(highlightNearest = list(enabled=TRUE,labelOnly=F)) %>% 
####### EXPORT
export_choices <- list("Proteins as CSV"='nodes.csv',
                       "Interactions as CSV"='edges.csv',
                       "Interactions with sources as CSV"='edges_sources.csv',
                       "Network as GML"='network.gml',
                       "Network as GRAPHML"='network.graphml',
                       "Network as DOT"='network.dot')

export_method <- function(filename,g,file) {
  if(filename=='nodes.csv') {write.csv(g$nodes,file,row.names=F)}
  if(filename=='edges.csv') {export_edges_csv(g$edges,file)}
  if(filename=='network.gml') {export_network_gml(g,file)}
  if(filename=='network.graphml') {export_network_graphml(g,file)}
  if(filename=='network.dot') {export_network_dot(g,file)}
  if(filename=='edges_sources.csv') {export_edges_sources_csv(g,file)}
}

export_edges_csv <- function(edges,file) {
  from_names <- sapply(edges %>% pull(from), get_gene_name)
  to_names <- sapply(edges %>% pull(to), get_gene_name)
  data <- edges %>% mutate(from_label = from_names,to_label=to_names)
  write.csv(data,file,row.names=F)
}

export_network_gml <- function(g,file) {
  h <- get_igraph(g)
  lines <- get_gml(h)
  writeLines(lines,con=file,sep="\n")
}

export_network_graphml <- function(g,file) {
  h <- get_igraph(g)
  write_graph(h,file,format="graphml")
}

export_network_dot <- function(g,file) {
  h <- get_igraph(g)
  write_graph(h,file,format="dot")
}

export_edges_sources_csv <- function(g,file) {
  df <- g$edges %>% select(from,to) %>% rename(Kinase=from,Substrate=to)
  data <- merge(x=source_data,y=df)
  write.csv(data,file,row.names=F)
}