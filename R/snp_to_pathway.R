#Some code borrowed from Nick

#'@export
retrieve_go_ids <- function(all_genes){
  ensembl <- biomaRt::useMart("ensembl",dataset="hsapiens_gene_ensembl")

  goids = biomaRt::getBM(attributes = c('ensembl_gene_id', 'go_id','name_1006','hgnc_symbol'),
                         filters = 'ensembl_gene_id',
                         values = all_genes,
                         mart = ensembl)
}

#'@title Annotate rs IDs to genomic loci
#'@export
annotate_snp_loc <- function(rsid,dbname="hg19"){
  con <- RMySQL::dbConnect(RMySQL::MySQL(),
                           host="genome-mysql.cse.ucsc.edu",
                           username="genome",
                           dbname=dbname)
  variant_df <- dplyr::tbl(con,"snp150")
  anno_df <- dplyr::filter(variant_df,name %in% rsid) %>%
    dplyr::select(name,chrom,pos=chromEnd,observed,func) %>% #Pull out the SNPs we need
    dplyr::collect() %>% #grab the whole thing from the database
    dplyr::mutate(chrom=str_replace(chrom,pattern = "chr([0-9]+)",replacement = "\\1")) %>%
    dplyr::filter(!str_detect(chrom,"_")) %>%
    tidyr::separate(observed,into=c("ref","alt"),sep="/",extra="drop")

  RMySQL::dbDisconnect(con)
  functional_categories <- c("nonsense","missense","ncRNA")

  anno_df <- mutate(anno_df, has_func=map_lgl(str_split(func,","), ~any(.x %in% functional_categories)))
  return(anno_df)

}

#'@title Extract eQTLs from GTEx data
#'@param anno_df is produced by annotate_snp_loc
#'@param gtex_files file list. Should be a list of tsv files with columns chr, pos, and gene_id
#'@description You can find files in /project2/compbio/external_public_supp/GTEx/V8/GTEx_Analysis_v8_eQTL/*chrpos.tsv.gz
#'@export
get_eqtl_genes <- function(anno_df,gtex_files, tissue_names){

  eqtl_genes <- purrr::map_df(seq_along(gtex_files), function(i){
                              f <- gtex_files[i]
                              n <- tissue_names[i]

                              eqtl <- read_tsv(f)
                              x <- distinct(anno_df,chrom,pos,name) %>%
                                   mutate(chr = paste0("chr", chrom)) %>%
                                   inner_join(eqtl,by=c("chr","pos")) %>%
                                   mutate(tissue = n)
                              })
  return(eqtl_genes)
}


#'@export
get_enhancer_genes <- function(anno_df, hacer_df,
                               criteria = c("FANTOM5", "4DGenome"),
                               margin=0){
  stopifnot(all(criteria %in% c("FANTOM5", "50kb", "4DGenome")))
  stopifnot(length(criteria) >= 1)

  hacer_df <- hacer_df %>% rename(gene_FANTOM5 = `associated_gene-FANTOM5`,
                                  gene_50kb = `associated_gene-50kb`,
                                  gene_4DGenome = `associated_gene-4DGenome`)

  df <- map_df( seq(22), function(i){
    x <- filter(hacer_df, chr==paste0("chr", i)) %>%
         select(start, end) %>%
         mutate(start = start -margin, end = end + margin) %>%
         Intervals()
    snp_pos <- filter(anno_df, chrom==i) %>% with(.,pos )
    snps <- filter(anno_df, chrom==i) %>% with(.,name )
    enhancers <- filter(hacer_df, chr==paste0("chr", i)) %>% with(., Enhancer_ID)
    overlap_ix <- interval_overlap(snp_pos, x)
    len <- sapply(overlap_ix, length)
    tibble(snp = rep(snps, len),
                 Enhancer_ID = enhancers[unlist(overlap_ix)])
  })
  enhancer_info <- left_join(df, hacer_df, by="Enhancer_ID")
  #FANOM5
  enhancer_genes <- map_df(criteria, function(cri){
    col <- paste0("gene_", cri)
    x <- enhancer_info %>% rename(gene = UQ(col) ) %>%
          filter(!is.na(gene)) %>%
          select(snp, Enhancer_ID, chr, start, end, gene) %>%
          mutate(source = cri) %>%
          mutate(gene = strsplit(as.character(gene), ",")) %>%
          unnest(gene)
    rmix <- grep("NA;", x$gene)
    x[-rmix,]
  })
  return(enhancer_genes)
}

#'@export
nearest_gene <- function(anno_df,gene_df){
  snp_chrom <- distinct(anno_df,chrom)
  gene_chrom <- distinct(gene_df,chrom)

  snpgene_chrom <- inner_join(snp_chrom,gene_chrom)

  genemap <- map_df(snpgene_chrom$chrom,function(ch){
      inner_join(
      filter(anno_df,chrom==ch),
      filter(gene_df,chrom==ch),by="chrom")  %>%
      mutate(gene_dist=ifelse(pos>=gene_start & pos <=gene_end,
                                     0,
                                     .Machine$integer.max)) %>%
      mutate(gene_dist=pmin(gene_dist,pmin(abs(pos-gene_end),
                                                  abs(pos-gene_start)))) %>%
      group_by(name) %>%
      filter(gene_dist==min(gene_dist)) %>%
      distinct(gene_name,.keep_all=T) %>%  ungroup()
  })
  return(genemap)
}

#'@title Get genes with TSS within a distance of variants
#'@export
genes_within <- function(anno_df, gene_df, dist=50000){
  df <- map_df( seq(22), function(i){
    x <- filter(gene_df, chrom== i) %>%
      mutate(end = gene_start + dist, start = gene_start -dist) %>%
      select(start, end) %>%
      Intervals()
    snp_pos <- filter(anno_df, chrom==i) %>% with(.,pos )
    snps <- filter(anno_df, chrom==i) %>% with(.,name )
    genes <- filter(gene_df, chrom==i) %>% with(., ensembl_gene_stable_id)
    overlap_ix <- interval_overlap(snp_pos, x)
    len <- sapply(overlap_ix, length)
    tibble(snp = rep(snps, len),
           gene = genes[unlist(overlap_ix)])
  })
  return(df)
}


#'@export
retrieve_genes <- function(getSymbol=T){

  con <- RMySQL::dbConnect(RMySQL::MySQL(),
                           host="genome-euro-mysql.soe.ucsc.edu",
                           username="genome",
                           dbname="hg19")

  nc_genelink_df <- tbl(con,"knownToRefSeq")
  nc_genelink_df
  entrez_genelink_df <- tbl(con,"knownToLocusLink") %>% rename(ucsc_name=name,entrez_gene=value)
  kt_ensembl <- tbl(con,"knownToEnsembl")
  egene <- tbl(con,"ensGene")
  es_genelink_df <- inner_join(kt_ensembl,egene,by=c("value"="name")) %>%
    dplyr::select(ucsc_name=name,ensembl_transcript=value,ensembl_gene_stable_id=name2) %>%
    dplyr::select(-ensembl_transcript) %>%
    inner_join(entrez_genelink_df)
  kgene <- tbl(con,"knownGene")

  if(getSymbol){
    ks <- tbl(con,"kgXref")
    kgene <- inner_join(kgene,ks,by=c("name"="kgID")) %>%
      dplyr::select(ucsc_name=name,
                    chrom,
                    geneSymbol,
                    gene_start=txStart,
                    gene_end=txEnd)
  }else{
    kgene <- kgene %>%
      dplyr::select(ucsc_name=name,
                    chrom,
                    gene_start=txStart,
                    gene_end=txEnd)
  }
  gene_df <- kgene %>%
    inner_join(es_genelink_df) %>%
    dplyr::select(-ucsc_name)  %>%
    collect() %>%
    dplyr::filter(chrom %in% paste("chr",1:22,sep="")) %>%
    dplyr::mutate(chrom=stringr::str_replace(chrom,pattern = "chr([0-9]+)",replacement = "\\1")) %>%
    rename(gene_name=entrez_gene) %>% distinct(gene_name,.keep_all = T)
  RMySQL::dbDisconnect(con)
  return(gene_df)
}

