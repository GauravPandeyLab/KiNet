library(tidyverse)
library(igraph)
library(xtable)
library(rvest)
library(dplyr)
library(stringr)

all_edges <- read_csv('data/ksi_display.csv') %>% 
  rename(from=Kinase,to=Substrate)
#widths <- data.frame(NSites =all_edges$NSites %>% unique %>% sort) %>% mutate(width=3*ntile(NSites,10))
#widths <- data.frame(NSites =all_edges$NSites %>% unique %>% sort) %>% mutate(width=5)
#all_edges <- merge(x=all_edges,y=widths,by="NSites",all.x=TRUE)
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

source_data <- read_csv('data/ksi_source_evi_ref.csv')

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
# last_date = as.Date('2023-10-04',"%Y-%m-%d")
last_date = as.Date('2024-07-09',"%Y-%m-%d")
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

#### References data
human_iptmNet <- read.csv("data/iptmNet_humanOnly.csv", header=FALSE)
# epsd_references <- read.csv("data/epsd_sources.txt", sep="\t")
ppp_references <- read.csv("data/ppp_only_human_refs.csv")

### TODO: combine the source data to export ###


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
  if(largeCenterNode & nrow(g$nodes)>10) {
    condition <- g$nodes$label %in% genes
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
  L3 <- c(L3,list('UniProt: ',tags$a(info$id,href= sprintf(unip_url,info$id), 
                                     target="_blank", rel="noreferrer noopener"), tags$br()))
  if(!is.null(info$HGNC)){
    hgnc_url <- "https://www.genenames.org/data/gene-symbol-report/#!/hgnc_id/HGNC:%s"
    L3 <- c(L3,list('HGNC: ',tags$a(info$HGNC,href= sprintf(hgnc_url,info$HGNC), 
                                    target="_blank", rel="noreferrer noopener"),tags$br()))
  }
  if(!is.null(info$GeneID)){
    gene_url <- "https://www.ncbi.nlm.nih.gov/gene/%s"
    L3 <- c(L3,list('NCBI Gene: ',tags$a(info$GeneID,href= sprintf(gene_url,info$GeneID), 
                                         target="_blank", rel="noreferrer noopener"),tags$br()))
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
  L1 <- tags$a(info$GeneName,href= sprintf(unip_url,info$id), 
               target="_blank", rel="noreferrer noopener")
  list(L1,tags$p(info$FullName))
}


# Render edge information
epsdrefs <- read_csv('data/EPSDReference.csv')
render_edge_info <- function(edge) {
  if (is.null(edge)) {
    print('no edge selected?')
    return('')
  }
  kinase <- edge %>% pull(from)
  substrate <- edge %>% pull(to)
  L1 <- tags$div(tags$h4('Kinase ',style="display: inline;"),render_edge_link(kinase))
  
  L2 <- tags$div(tags$h4('Substrate ',style="display: inline;"),render_edge_link(substrate))
  
  
  # Add site & sequence table
  nsites <- edge %>% pull(NSites)
  if(nsites==1) {
    L3 <- list(tags$h4('Phosphorylation Sites'),tags$p('1 site known'))
    } else {
      L3 <- list(tags$h4('Phosphorylation Sites'),tags$p(nsites,'sites known'))
    }
  sites <- edge %>% pull(Sites)
  L4 <- tags$p(gsub('"','',sites))
  
  view_bool = FALSE
  
  edge_df <- edge %>% select(from,to) %>% rename(Kinase=from,Substrate=to)
  detail_source_data <- merge(x=source_data,y=edge_df)
  sources_KSI <- c(detail_source_data[['PrimarySource']], 
                   detail_source_data[['SecondarySource']])
  
  # Add uniprot API to parse the site (or download beforehand)
  site_str_list <- unlist(strsplit(sites, ", "))
  site_num_locs <- unlist(lapply(substring(site_str_list, 2), as.integer))
  print(substrate)
  
  ### References
  site_refs <- c() 
  iptmNet_substrate_bool <- human_iptmNet$V11 == substrate
  iptmNet_kinase_bool <- human_iptmNet$V7 == kinase
  # epsd_substrate_bool <- epsd_references$UniProt.ID == substrate
  
  
  ppp_references_kinase_bool <- ppp_references$kinase_uniprot == kinase
  ppp_references_substrate_bool <- ppp_references$substrate_uniprot == substrate
  pmid_template <- "https://pubmed.ncbi.nlm.nih.gov/"
  ### Separate reference list by sources?
  print('site list')
  print(site_str_list)
  ppp_refs_by <- ""
  source_list <- c()
  refs_list <- c()
  refs_pmid_list <- c()
  for (site_name in site_str_list) {
    site_refs_perSite <- c()
    site_refsText_perSite <- c()
    if ('iPTMNet' %in% sources_KSI) {
      site_bool <- human_iptmNet$V6 == site_name
      human_iptmNet_site_ref <- human_iptmNet[(site_bool & iptmNet_substrate_bool & iptmNet_kinase_bool), ]
      combined_ref <- paste(human_iptmNet_site_ref$V10, collapse = ",")
      pmid_name_ <- unlist(strsplit(combined_ref, ","))
      # print(combined_ref)
      site_refs_perSite <- c(site_refs_perSite, pmid_name_)
      site_refsText_perSite <- c(site_refsText_perSite, pmid_name_)
      # refs_list <- c(refs_list, pmid_name_)
    }
    
    if ('EPSD' %in% sources_KSI) {
      # source_list <- c(source_list, 'EPSD*')
      # site_num <- as.integer(substring(site_name, 2))
      # print('EPSD references')
      # print(site_num)
      # site_bool <- epsd_references$Position == site_num
      # epsd_site_ref <- epsd_references[(site_bool & epsd_substrate_bool), ]
      # 
      # epsd_combined_ref <- paste(epsd_site_ref$Reference, collapse = ";")
      # print(epsd_combined_ref)
      # pmid_name_ <- str_trim(unlist(strsplit(epsd_combined_ref, ";")))
      # pmid_name_ <- pmid_name_[pmid_name_ != ""]
      # 
      # site_refs_perSite <- c(site_refs_perSite, pmid_name_)
      # if (length(pmid_name_) > 0 ) {
      #   site_refsText_perSite <- c(site_refsText_perSite, paste0('<i>*', c(pmid_name_), '</i>'))
      # }
      # site_refs <- c(site_refs, "")
    } 
    
    if ('PhosphoSitePlus' %in% sources_KSI) {
      # read owl?
      # source_list <- c(source_list, 'PhosphoSitePlus')
      site_num_int <- as.integer(substring(site_name, 2))
      print(site_num_int)
      ppp_references_site_bool <- ppp_references$site_num == site_num_int
      ppp_site_ref <- ppp_references[(ppp_references_site_bool & ppp_references_kinase_bool & ppp_references_substrate_bool),]
      combined_ref <- paste(ppp_site_ref$xref_all, collapse = ", ")
      combined_ref <- gsub("\\[", "", combined_ref)
      combined_ref <- gsub("\\]", "", combined_ref)
      combined_ref <- gsub("\\'", "", combined_ref)
      combined_ref <- gsub("pubmed", "", combined_ref)
      # combined_ref <- gsub("[", "", combined_ref)
      # combined_ref <- gsub("", "", combined_ref)
      print(combined_ref)
      pmid_name_ <- unlist(strsplit(combined_ref, ", "))
      pmid_name_ <- pmid_name_[pmid_name_ != ""]
      pmid_name_ <- str_sort(pmid_name_)
      # 
      site_refs_perSite <- c(site_refs_perSite, pmid_name_)
      if (length(pmid_name_) > 0 ) {
        site_refsText_perSite <- c(site_refsText_perSite, pmid_name_)
      }
      ppp_refs_by <- paste(ppp_refs_by, ppp_site_ref$evidence)
      # refs_list <- c(refs_list, pmid_name_)
      # site_refs <- c(site_refs, "")
    }
    
    if (length(site_refs_perSite) > 0) {
      # pmid_name <- site_refs_perSite[!duplicated(site_refs_perSite)]
      # refs_list <- 
      pmid_url <- paste0(pmid_template, site_refs_perSite)
      
      site_ref_df <- data.frame(c(pmidName = list(site_refs_perSite), 
                                  pmidURL = list(pmid_url),
                                  pmidHTMLText = list(site_refsText_perSite)))
      
      site_ref_df <- site_ref_df[!duplicated(site_ref_df$pmidName),]
      site_ref_df <- site_ref_df[order(site_ref_df$pmidName, decreasing = TRUE),]
      
      html_with_url <- apply(site_ref_df, 1, function(x){paste0(c("<a href='", x[['pmidURL']], 
                                                                  "' target='_blank', rel='noreferrer noopener'>", 
                                                                  x[['pmidHTMLText']], '</a>' 
                                                                  # href=x[['pmidURL']], 
                                                                  # target="_blank", 
                                                                  # rel="noreferrer noopener"
      ), collapse = "")})
      html_with_url_str <- c()
      for (htmlIdx in 1:length(html_with_url)) {
        html_with_url_str <- c(html_with_url_str, as.character(html_with_url[[htmlIdx]]))
      }
      print(paste(html_with_url_str, collapse = ', '))
      site_refs <- c(site_refs, 
                     # paste(html_with_url_str, collapse = ',<br>')
                     paste(html_with_url_str, collapse = '<br>')
      )
      refs_pmid_list <- c(refs_pmid_list, site_ref_df$pmidName)
      refs_list <- c(refs_list, html_with_url_str)
      print(site_refs)
    }
    else {
      site_refs <- c(site_refs, "")
    }
  }
  
  uniprot_fasta_url <- paste(c("https://rest.uniprot.org/uniprotkb/", substrate, ".fasta"), collapse="")
  # print(uniprot_fasta_url)
  page <- read_html(uniprot_fasta_url)
  print(page)
  fasta_seq <- page %>% html_nodes("p") %>% html_text()
  seq_with_n <- unlist(strsplit(fasta_seq, "\n"))
  
  seq_in_one <- paste(seq_with_n[2:length(seq_with_n)], collapse = "")
  
  
  ### Parse Sequences data
  site_seqs <- c() 
  for ( site_num in site_num_locs) {
    site0 <- site_num-7
    if (site0 < 1) {site0 <- 1}
    seqs_HTML <- paste(c(substring(seq_in_one, site0, site_num-1),
                       "<b>", substring(seq_in_one, site_num, site_num), '</b>',
                       substring(seq_in_one, site_num+1, site_num+7)), 
                       collapse = "")
    print(seqs_HTML)
    site_seqs <- c(site_seqs, seqs_HTML)
  }
  
  
  
  print(site_seqs)
  siteTable <- data.frame(c(Site = strsplit(sites, ", "),
                            Sequence = list(site_seqs)
                            # Reference = list(site_refs)
                            )
                          )
  
  T1 <- div(id = "div-table-1",
            renderTable(xtable(siteTable),
                        sanitize.text.function=function(x){x}
                        # , options = list(dom = 't')
                        )
        )
  
  refTable <- data.frame(c())
  
  
  # 
  # actionButton('view_hide_table', 'view/hide site sequence table')
  # observeEvent(input$view_hide_table, {
  #   toggle("button-test")
  # })
  ppp_url <- 'https://www.phosphosite.org/uniprotAccAction?id=%s'
  L50 <- tags$h4('Database Sources')
  # Include only the sources exist
  

  # print(sources_KSI)
  Lbr <- tags$br()
  
  # L5C <- c(L50)
  evidence_types <- c()
  if ('iPTMNet' %in% sources_KSI) {
    # print('iptmnet')
    ipt_url <- 'https://research.bioinformatics.udel.edu/iptmnet/entry/%s/'
    L52 <- tagList(tags$strong('iPTMNet',
                          # href= sprintf(ipt_url,substrate), 
                          # target="_blank", rel="noreferrer noopener"
                          ), 
                   # tags$span(': text mining'),
                   tags$br()
                   )
    source_list <- c(source_list, 'iPTMNet')
    # source_list <- source_list
    # L5C <- c(L5C, L52, Lbr)
    evidence_types <- c(evidence_types, 'Text mining')
  } else { L52 <- '' }
  # L70 <- tags$h5('')
  
  L70 <- tags$h5('*References from all source databases. However, EPSD does not provide site/interaction-specific references. Hence, they are omitted from this table.')
  if ('EPSD' %in% sources_KSI) {
    print('epsd')
    evidence_types <- c(evidence_types, 'Unspecified experimental method')
    
    source_list <- c(source_list, 'EPSD')
    epsd <- epsdrefs %>% filter(Substrate==substrate) %>% pull('EPSD')
    L53 <- ''
    if(length(epsd)>=1) {
      epsd_url <- 'https://epsd.biocuckoo.cn/View.php?id=%s'
      L70 <- tags$h5('*References from all source databases. However, EPSD does not provide site/interaction-specific references. Hence, they are omitted from this table.')
      # L71 <- tags$h6('')
      # L7 <- tags$div(L70, L71)
      L53 <- tagList(tags$strong('EPSD*',
                            # href=sprintf(epsd_url,epsd[[1]]), 
                            # target="_blank", rel="noreferrer noopener"
                            ), 
                     # tags$span(': identified by experiment'),
                     tags$br(),
                     L70,
                     tags$br(),
                     )
    }
    # L5C <- c(L5C, L53, Lbr)
  } else { L53 <- '' }
  
  if ('PhosphoSitePlus' %in% sources_KSI) {
    ppp_evidence_string_list <- c()
    source_list <- c(source_list, 'PhosphoSitePlus')
    
    # if (grepl('identified', ppp_refs_by)){
    #   ppp_evidence_string <- c(ppp_evidence_string_list, ": identified by")
    # }
    if (grepl('antibody', ppp_refs_by)) { 
      ppp_evidence_string_list <- c(ppp_evidence_string_list, "Antibody")
      # evidence_types <- c(evidence_types, 'Antibody')
    }
    if (grepl('mass spectrometry', ppp_refs_by)) { 
      ppp_evidence_string_list <- c(ppp_evidence_string_list, "Mass spectrometry")
      # evidence_types <- c(evidence_types, 'Mass spectrometry')
    }
    if (grepl('western blot', ppp_refs_by)) { 
      ppp_evidence_string_list <- c(ppp_evidence_string_list, "Western blot")
      # evidence_types <- c(evidence_types, 'Western blot')
    }
    if (grepl('mutation analysis', ppp_refs_by)) { 
      ppp_evidence_string_list <- c(ppp_evidence_string_list, "Mutation analysis")
      # evidence_types <- c(evidence_types, 'Mutation analysis')
    }
    
    if (length(ppp_evidence_string_list) == 0){
      ppp_htmlText <- 'N/A'
      evidence_types <- c(evidence_types, 'N/A')
    } else {
      ppp_htmlText <- paste('', paste(ppp_evidence_string_list, collapse='<br>'))
      evidence_types <- c(evidence_types, ppp_htmlText)
    }
    
    
    L51 <- tagList(tags$strong('PhosphoSitePlus',
                               # href= sprintf(ppp_url,substrate), 
                               # target="_blank", rel="noreferrer noopener"
    ), 
    tags$span(ppp_htmlText),
    tags$br()
    )
    # L5C <- c(L5C, L51, Lbr)
    
  } else { L51 <- '' }
  
  L5 <- tags$div(L50,L51,L52,L53)
  # L5 <-tags$div(list(L5C))
  
  L60 <- tags$h4('Citation')
  L61 <- tags$a('')
  num_evidence_types <- length(evidence_types)
  evidence_types_string <- ""
  # if (num_evidence_types > 2) {
  #   evidence_types_string <- paste(c(paste(evidence_types[1:(num_evidence_types-1)], 
  #                                          collapse = ',<br>'), 
  #                                    evidence_types[num_evidence_types]),
  #                                  collapse = ',<br>or ')
  #   
  # } else if (num_evidence_types == 2){
  #   evidence_types_string <- paste(evidence_types,
  #                                  collapse = ',<br>or ')
  # } else if (num_evidence_types == 1){
  #   evidence_types_string <- evidence_types[1]
  # if (num_evidence_types > 2) {
  #   evidence_types_string <- paste(c(paste(evidence_types[1:(num_evidence_types-1)], 
  #                                          collapse = ',<br>'), 
  #                                    evidence_types[num_evidence_types]),
  #                                  collapse = ',<br>or ')
  #   
  # } else if (num_evidence_types == 2){
  #   evidence_types_string <- paste(evidence_types,
  #                                  collapse = ',<br>or ')
  # } else if (num_evidence_types == 1){
    # evidence_types_string <- evidence_types[1]
  #   evidence_types_string <- paste(evidence_types,
  #                                  collapse = ',<br>or ')
  # } else {
  #   evidence_types_string <- "N/A"
  # }
  
  ref_df <- data.frame(c(pmid = list(refs_pmid_list),
                         pmidHTML = list(refs_list)))
  print(ref_df)
  ref_df <- ref_df[!duplicated(ref_df$pmid),]
  
  if (nrow(ref_df) > 0 ){
    ref_df <- ref_df[order(ref_df$pmid, decreasing = FALSE),]
    ref_string <- paste(ref_df$pmidHTML, 
                             collapse = '<br>')
  } else {
    ref_string <- "N/A"
  }
  
  # ref_df$pmidHTML
  # refs_list <- unique(refs_list)
  
  # refsTable <- data.frame(c(`Source Database` = list(paste(source_list, 
  #                                                          collapse = '<br>')),
  refsTable <- data.frame(c(`Source Database` = list(source_list),
                            `Evidence Type` = list(evidence_types))
                          )
  colnames(refsTable) <- c("Source \n Database", "Evidence \n Type")
  
  
  T2 <- div(id = "div-table-1",
            renderTable(xtable(refsTable),
                        sanitize.text.function=function(x){x}
                        # , options = list(dom = 't')
            )
  )
  # L8 <- 
  L80 <- tagList(tags$h4("Supporting References (PMID)*"), 
                 HTML(ref_string))
  
  tagList(L1,L2,L3,T1,L50,T2,L80,L70)
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
  visEdges(arrows="to",smooth = list(enabled = T, type = 'dynamic'),color=list(opacity=0.5),selectionWidth=10,width=5) %>%
  visOptions(highlightNearest = list(enabled=T,degree=list(from=1,to=1),algorithm="hierarchical")) %>% 
  visExport(type="png",label="Screenshot visible region as PNG")
}

layout_choices <- list("Default"='layout_with_fr',
                       "Circle"='layout_in_circle',
                       "Grid"='layout_on_grid')
#visOptions(highlightNearest = list(enabled=TRUE,labelOnly=F)) %>% 
####### EXPORT
export_choices <- list("Interactions (with sources) as CSV"='edges_sources.csv',
                       "Proteins as CSV"='nodes.csv',
                       "Interactions as CSV"='edges.csv',
                       "Network as GML"='network.gml',
                       "Network as GRAPHML"='network.graphml',
                       "Network as DOT"='network.dot')

export_method <- function(filename,g,file) {
  if(filename=='nodes.csv') {export_nodes_csv(g$nodes,file)}
  if(filename=='edges.csv') {export_edges_csv(g$edges,file)}
  if(filename=='network.gml') {export_network_gml(g,file)}
  if(filename=='network.graphml') {export_network_graphml(g,file)}
  if(filename=='network.dot') {export_network_dot(g,file)}
  if(filename=='edges_sources.csv') {export_edges_sources_csv(g,file)}
}

export_nodes_csv <- function(nodes,file) {
  data <- nodes %>% select(c('id','GeneName','group'))
  write.csv(data,file,row.names=F)
}


# Add the gene names to export file
export_edges_csv <- function(edges,file) {
  from_names <- sapply(edges %>% pull(from), get_gene_name)
  to_names <- sapply(edges %>% pull(to), get_gene_name)
  data <- edges %>% mutate(from_label = from_names,to_label=to_names)
  data <- data %>% select(c('from','to','from_label','to_label','NSites','Sites'))
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
  # df <- g$edges %>% select(from,to) %>% rename(Kinase=from,Substrate=to)
  # print(g)
  edges <- g$edges
  from_names <- sapply(edges %>% pull(from), get_gene_name)
  to_names <- sapply(edges %>% pull(to), get_gene_name)
  data_ <- edges %>% mutate(from_label = from_names,to_label=to_names)
  data_ <- data_ %>% select(c('from','to','from_label','to_label'))
  
  df <- data_ %>% rename(Kinase=from,Substrate=to, 
                        "Kinase Name"=from_label, "Substrate Name"=to_label,
                        
                        )
  # df <- df[,c('Kinase', 'KinaseName', 'Substrate', 'SubstrateName', 'Site', 'Source Database', 'Evidence','Reference (PMID)')]
  data <- merge(x=source_data,y=df)
  data <- data %>% rename("Source Database"=PrimarySource,
                          "Reference (PMID)"=pmid_ref,
                          "Evidence"=evidence)
  print(colnames(data))
  data <- data[,c('Kinase', 'Kinase Name', 'Substrate', 'Substrate Name', 'Site', 'Source Database', 'Evidence','Reference (PMID)')]
  write.csv(data,file,row.names=F)
}