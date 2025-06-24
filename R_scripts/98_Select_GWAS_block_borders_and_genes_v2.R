
suppressMessages(library("plyr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("data.table", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("crayon", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggplot2", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("farver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("labeling", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("optparse", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("dplyr", lib.loc = "/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("withr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("backports", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("broom", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rstudioapi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cli", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tzdb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("svglite", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggeasy", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("sandwich", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("digest", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("tidyverse", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


opt = NULL

options(warn = 1)

intersect_with_ALL = function(option_list)
{
  suppressMessages(library("BiocGenerics", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("S4Vectors", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("IRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("GenomeInfoDb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("GenomicRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("Biobase", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("AnnotationDbi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("GO.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  suppressMessages(library("rtracklayer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
  library("R.oo", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/")
  library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/")
  
  
  opt_in = option_list
  opt <<- option_list
  
  cat("All options:\n")
  printList(opt)
  
  
  #### READ and transform type ----
  
  type = opt$type
  
  cat("TYPE_\n")
  cat(sprintf(as.character(type)))
  cat("\n")
  
  #### READ and transform out ----
  
  out = opt$out
  
  cat("out_\n")
  cat(sprintf(as.character(out)))
  cat("\n")
  
  #### LOOP PRINTING ----
  
  path_files<-paste(out,'Build_files','/',sep='')
  
  if (file.exists(path_files)){
    
    
  }else{
    
    dir.create(file.path(path_files))
    
  }#path_files
  

  #### READ variants_file ----
  
  variants_file<-readRDS(file=opt$variants_file)
  
  cat("variants_file_0\n")
  cat(str(variants_file))
  cat("\n")
  cat(str(unique(variants_file$VAR)))
  cat("\n")
  
  #### READ GWAS block sizes ----
  
  BLOCK_PRE<-readRDS(file=opt$BLOCK_PRE)
  
 
  
  cat("BLOCK_PRE_0\n")
  cat(str(BLOCK_PRE))
  cat("\n")
 
  #### READ PCHiC_HITS ----
  
  PCHiC_HITS<-readRDS(file=opt$PCHiC_HITS)
  
  cat("PCHiC_HITS_0\n")
  cat(str(PCHiC_HITS))
  cat("\n")
  cat(str(unique(PCHiC_HITS$VAR)))
  cat("\n")
  cat(str(unique(PCHiC_HITS$oeID)))
  cat("\n")
  
  #### READ ensembl_gtf ----
  
  ensembl_gtf = readGFF(opt$ensembl_gtf)
  
  
  cat("ensembl_gtf_0\n")
  cat(str(ensembl_gtf))
  cat("\n")
  
  
  ###### LOOP AS_IDs -----
  
  path_files<-paste(out,'Build_files','/',sep='')
  
  if (file.exists(path_files)){
    
    
  }else{
    
    dir.create(file.path(path_files))
    
  }#path_files
  
  
  AS_ID_array<-unique(variants_file$AS_ID)
  
  cat("AS_ID_array_0\n")
  cat(str(AS_ID_array))
  cat("\n")
  
  
  # AS_ID_array<-'plt__411'

  
  DEBUG<-0
  
  
  Result_FINAL_genes<-data.frame()
  Result_FINAL_boundaries<-data.frame()
 
  
  
  for(i in 1:length(AS_ID_array))
  {
    AS_ID_sel<-AS_ID_array[i]
    
    cat("------------------------------------------------------------------------------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(AS_ID_sel)))
    cat("\n")
    
    if(AS_ID_sel == 'plt__411')
    {

      DEBUG<-1

    }else{

      DEBUG<-0
    }
    
    
    variants_file_sel<-variants_file[which(variants_file$AS_ID == AS_ID_sel),]
    
    if(DEBUG == 1)
    {
      cat("variants_file_sel_0\n")
      cat(str(variants_file_sel))
      cat("\n")
    }
    
    chr_sel<-gsub("^chr","",unique(variants_file_sel$chr))
    
    if(DEBUG == 1)
    {
      cat("chr_sel_0\n")
      cat(str(chr_sel))
      cat("\n")
    }
    
    ensembl_gtf_chr_sel<-ensembl_gtf[which(ensembl_gtf$seqid == chr_sel &
                                             ensembl_gtf$type == 'gene'),]
    
    if(DEBUG == 1)
    {
      cat("ensembl_gtf_chr_sel_0\n")
      cat(str(ensembl_gtf_chr_sel))
      cat("\n")
    }
    
    gr_ensembl_genes <- GRanges(
      seqnames = as.character(ensembl_gtf_chr_sel$seqid),
      strand=ensembl_gtf_chr_sel$strand,
      name2=ensembl_gtf_chr_sel$gene_name,
      ranges=IRanges(
        start=ensembl_gtf_chr_sel$start,
        end=ensembl_gtf_chr_sel$end,
        names=ensembl_gtf_chr_sel$gene_id))
    
    if(DEBUG == 1)
    {
      cat("gr_ensembl_genes_0\n")
      cat(str(gr_ensembl_genes))
      cat("\n")
    }
    
    
    BLOCK_PRE_sel<-BLOCK_PRE[which(BLOCK_PRE$AS_ID == AS_ID_sel),]
    
    
    if(DEBUG == 1)
    {
      cat("BLOCK_PRE_sel_0\n")
      cat(str(BLOCK_PRE_sel))
      cat("\n")
    }
    
    
    temp<-data.frame()
    
    if(dim(BLOCK_PRE_sel)[1] >0)
    {
      
      gr_BLOCK_PRE <- GRanges(
        seqnames = as.character(gsub("^chr","",BLOCK_PRE_sel$chr)),
        strand="*",
        ranges=IRanges(
          start=BLOCK_PRE_sel$START,
          end=BLOCK_PRE_sel$END,
          names=BLOCK_PRE_sel$AS_ID))
      
      if(DEBUG == 1)
      {
        cat("gr_BLOCK_PRE_0\n")
        cat(str(gr_BLOCK_PRE))
        cat("\n")
      }
      
      
      #### find overlap with SNP ----
      
      m <- findOverlaps(gr_BLOCK_PRE,gr_ensembl_genes)
      
      if(DEBUG == 1)
      {
        cat("m\n")
        cat(str(m))
        cat("\n")
      }
      
      subjectHits_ensembl_genes<-subjectHits(m)
      
      if(DEBUG == 1)
      {
        cat("subjectHits_ensembl_genes\n")
        cat(str(subjectHits_ensembl_genes))
        cat("\n")
      }
      
      queryHits_BLOCK_PRE<-queryHits(m)
      
      if(DEBUG == 1)
      {
        cat("queryHits_BLOCK_PRE\n")
        cat(str(queryHits_BLOCK_PRE))
        cat("\n")
      }
      
      BLOCK_PRE_2 <- data.frame(chr=as.character(seqnames(gr_BLOCK_PRE)),
                           start=as.integer(start(gr_BLOCK_PRE)),
                           end=as.integer(end(gr_BLOCK_PRE)),
                           AS_ID=names(gr_BLOCK_PRE), stringsAsFactors = F)
      
      if(DEBUG == 1)
      {
        cat("BLOCK_PRE_2_0\n")
        cat(str(BLOCK_PRE_2))
        cat("\n")
      }
      
      BLOCK_PRE_2_hits<-BLOCK_PRE_2[queryHits_BLOCK_PRE,]
      
      if(DEBUG == 1)
      {
        cat("BLOCK_PRE_2_hits_0\n")
        cat(str(BLOCK_PRE_2_hits))
        cat("\n")
      }
      
      ensembl_genes_2 <- data.frame(chr=as.character(seqnames(gr_ensembl_genes)),
                            start=as.integer(start(gr_ensembl_genes)),
                            end=as.integer(end(gr_ensembl_genes)),
                            ensembl_gene_id=as.character(names(gr_ensembl_genes)),
                            strand=as.character(strand(gr_ensembl_genes)),
                            hgnc=as.character(gr_ensembl_genes$name2),
                            stringsAsFactors = F)
      
    
      if(DEBUG == 1)
      {
        cat("ensembl_genes_2_0\n")
        cat(str(ensembl_genes_2))
        cat("\n")
      }
      
      
      
      ensembl_genes_2_hits<-ensembl_genes_2[subjectHits_ensembl_genes,]
      
      if(dim(ensembl_genes_2_hits)[1] >0)
      {
        ensembl_genes_2_hits$AS_ID<-AS_ID_sel
        
        if(DEBUG == 1)
        {
          cat("ensembl_genes_2_hits_0\n")
          cat(str(ensembl_genes_2_hits))
          cat("\n")
          cat(str(unique(ensembl_genes_2_hits$ensembl_gene_id)))
          cat("\n")
        }
        
        temp<-rbind(ensembl_genes_2_hits,temp)
        
      }#dim(ensembl_genes_2_hits)[1] >0
        
        
        
        
        
      PCHiC_HITS_sel<-PCHiC_HITS[which(PCHiC_HITS$AS_ID == AS_ID_sel),] 
        
        
        
      if(dim(PCHiC_HITS_sel)[1] >0)
      {
        if(DEBUG == 1)
        {
          cat("PCHiC_HITS_sel_0\n")
          cat(str(PCHiC_HITS_sel))
          cat("\n")
          cat(sprintf(as.character(PCHiC_HITS_sel$ensembl_gene_id)))
          cat("\n")
        }
        
        PCHiC_HITS_sel<-unique(as.data.frame(cSplit(PCHiC_HITS_sel,sep = ';', direction = "long",
                                                    splitCols = "ensembl_gene_id"),stringsAsFactors=F))
        if(DEBUG == 1)
        {
          cat("PCHiC_HITS_sel_1\n")
          cat(str(PCHiC_HITS_sel))
          cat("\n")
          cat(sprintf(as.character(PCHiC_HITS_sel$ensembl_gene_id)))
          cat("\n")
        }
        
        FLAG_NA_bait<-dim(PCHiC_HITS_sel[is.na(PCHiC_HITS_sel$ensembl_gene_id),])[1]
        
        if(DEBUG == 1)
        {
          cat("FLAG_NA_bait:\n")
          cat(str(FLAG_NA_bait))
          cat("\n")
        
        }
        
       
        ensembl_genes_3 <- data.frame(chr=as.character(ensembl_gtf_chr_sel$seqid),
                                      start=as.integer(ensembl_gtf_chr_sel$start),
                                      end=as.integer(ensembl_gtf_chr_sel$end),
                                      ensembl_gene_id=as.character(ensembl_gtf_chr_sel$gene_id),
                                      strand=as.character(ensembl_gtf_chr_sel$strand),
                                      hgnc=as.character(ensembl_gtf_chr_sel$gene_name),
                                      stringsAsFactors = F)
        
        if(DEBUG == 1)
        {
          cat("ensembl_genes_3_0\n")
          cat(str(ensembl_genes_3))
          cat("\n")
        }
        
        ensembl_genes_3_HITS<-ensembl_genes_3[which(ensembl_genes_3$ensembl_gene_id%in%PCHiC_HITS_sel$ensembl_gene_id),]
        
        if(DEBUG == 1)
        {
          cat("ensembl_genes_3_HITS_0\n")
          cat(str(ensembl_genes_3_HITS))
          cat("\n")
        }
        
        if(dim(ensembl_genes_3_HITS)[1] >0)
        {
          ensembl_genes_3_HITS$AS_ID<-AS_ID_sel
          
          if(DEBUG == 1)
          {
            cat("ensembl_genes_3_HITS_0\n")
            cat(str(ensembl_genes_3_HITS))
            cat("\n")
            cat(str(unique(ensembl_genes_3_HITS$ensembl_gene_id)))
            cat("\n")
          }
          
          temp<-unique(rbind(ensembl_genes_3_HITS,temp))
          
        }#dim(ensembl_genes_3_HITS)[1] >0
       
      }#dim(PCHiC_HITS_sel)[1] >0 (1)
          
        
      if(DEBUG == 1)
      {
        cat("temp_0\n")
        cat(str(temp))
        cat("\n")
        cat(str(unique(temp$ensembl_gene_id)))
        cat("\n")
        cat(sprintf(as.character((temp$hgnc))))
        cat("\n")
      }
      
      ALL_POS<-c(BLOCK_PRE_sel$START,BLOCK_PRE_sel$END,temp$start,temp$end)
      
      
      ABSOLUTE_START<-min(ALL_POS)
      ABSOLUTE_END<-max(ALL_POS)
      
      if(DEBUG == 1)
      {
        cat("ABSOLUTE_START\n")
        cat(sprintf(as.character(ABSOLUTE_START)))
        cat("\n")
        
        cat("ABSOLUTE_END\n")
        cat(sprintf(as.character(ABSOLUTE_END)))
        cat("\n")
        
      }
      
      
      recomposed_block<-as.data.frame(cbind(chr_sel,ABSOLUTE_START,ABSOLUTE_END,AS_ID_sel))
      
      colnames(recomposed_block)<-c('chr','start','end','AS_ID')
      
      recomposed_block$start<-as.integer(recomposed_block$start)
      recomposed_block$end<-as.integer(recomposed_block$end)
      
      if(DEBUG == 1)
      {
        cat("recomposed_block_0\n")
        cat(str(recomposed_block))
        cat("\n")
        
        # quit(status = 1)
      }
      
      Result_FINAL_boundaries<-rbind(Result_FINAL_boundaries,recomposed_block)
      Result_FINAL_genes<-rbind(Result_FINAL_genes,temp)
     
    }#dim(BLOCK_PRE_sel)[1] >0
    
   
    
  }# i in 1:length(AS_ID_array)
  
  
 

  #### SAVE ----

 
  
  if(dim(Result_FINAL_genes)[1] >0)
  {
    Result_FINAL_genes$chr<-as.character(paste('chr',Result_FINAL_genes$chr, sep=''))
    
    Result_FINAL_genes$chr<-factor(Result_FINAL_genes$chr,
                                        levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                                                 "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                                 "chr22","chr23","chrX","chrY"), ordered=T)
    
    Result_FINAL_genes<-droplevels(Result_FINAL_genes)
    
    Result_FINAL_genes<-Result_FINAL_genes[order(Result_FINAL_genes$chr),]
    
    cat("Result_FINAL_genes_0\n")
    cat(str(Result_FINAL_genes))
    cat("\n")
    cat(str(unique(Result_FINAL_genes$AS_ID)))
    cat("\n")
    cat(str(unique(Result_FINAL_genes$ensembl_gene_id)))
    cat("\n")
    
    setwd(path_files)
    
    saveRDS(Result_FINAL_genes, file=paste('Intersected_genes',".rds",sep=''))
    
    
    
  }#dim(Result_FINAL_genes)[1] >0
  
  
  if(dim(Result_FINAL_boundaries)[1] >0)
  {
    Result_FINAL_boundaries$chr<-as.character(paste('chr',Result_FINAL_boundaries$chr, sep=''))
    
    Result_FINAL_boundaries$chr<-factor(Result_FINAL_boundaries$chr,
                       levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                                "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                "chr22","chr23","chrX","chrY"), ordered=T)
    
    Result_FINAL_boundaries<-droplevels(Result_FINAL_boundaries)
    
    Result_FINAL_boundaries<-Result_FINAL_boundaries[order(Result_FINAL_boundaries$chr),]
    
    cat("Result_FINAL_boundaries_0\n")
    cat(str(Result_FINAL_boundaries))
    cat("\n")
    cat(str(unique(Result_FINAL_boundaries$AS_ID)))
    cat("\n")
    
    setwd(path_files)
    
    saveRDS(Result_FINAL_boundaries, file=paste('Boundaries_AS_POST_genes',".rds",sep=''))
    
    
    
  }#dim(Result_FINAL_boundaries)[1] >0
  

  
}



printList = function(l, prefix = "    ") {
  list.df = data.frame(val_name = names(l), value = as.character(l))
  list_strs = apply(list.df, MARGIN = 1, FUN = function(x) { paste(x, collapse = " = ")})
  cat(paste(paste(paste0(prefix, list_strs), collapse = "\n"), "\n"))
}


#### main script ----

main = function() {
  cmd_line = commandArgs()
  cat("Command line:\n")
  cat(paste(gsub("--file=", "", cmd_line[4], fixed=T),
            paste(cmd_line[6:length(cmd_line)], collapse = " "),
            "\n\n"))
  option_list <- list(
    make_option(c("--variants_file"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--BLOCK_PRE"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ensembl_gtf"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--PCHiC_HITS"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--type"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--out"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required.")
  )
  parser = OptionParser(usage = "140__Rscript_v106.R
                        --subset type
                        --TranscriptEXP FILE.txt
                        --cadd FILE.txt
                        --ncboost FILE.txt
                        --type type
                        --out filename",
                        option_list = option_list)
  opt <<- parse_args(parser)
  
  intersect_with_ALL(opt)

  
}


###########################################################################

system.time( main() )