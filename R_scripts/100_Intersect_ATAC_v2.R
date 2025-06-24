
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
  
  #### READ ATAC_RAW ----
  
  
  VAR_vector<-c('chr18_60920854_C_T','chr3_128317978_C_T')
  
  my_DF<- data.frame(matrix(vector(), length(VAR_vector), 5,
                                                dimnames=list(c(),
                                                              c("VAR","chr","pos","ref","alt"))),  stringsAsFactors=F)
  
  my_DF$VAR<-VAR_vector
  
  my_DF$chr<-gsub("_.+$","",my_DF$VAR)
  my_DF$pos<-gsub("^[^_]+_","",my_DF$VAR)
  my_DF$pos<-as.integer(gsub("_.+$","",my_DF$pos))
  my_DF$ref<-gsub("^[^_]+_[^_]+_","",my_DF$VAR)
  my_DF$ref<-gsub("_.+$","",my_DF$ref)
  my_DF$alt<-gsub("^[^_]+_[^_]+_[^_]+_","",my_DF$VAR)
  
  gr_VARS <- GRanges(
    seqnames = as.character(gsub("chr","",my_DF$chr)),
    ranges=IRanges(
      start=as.numeric(my_DF$pos),
      end=as.numeric(my_DF$pos),
      name=my_DF$VAR))
  
  cat("gr_VARS_0\n")
  cat(str(gr_VARS))
  cat("\n")
  
  

  
  #### READ ATAC_peaks ----
  
  
  ATAC_peaks<-as.data.frame(fread(file=opt$ATAC_peaks, sep="\t", header =F), stringAsFactors=F)
  
  colnames(ATAC_peaks)<-c('chr','start','end')
  
  
  cat("ATAC_peaks_0\n")
  cat(str(ATAC_peaks))
  cat("\n")
  
  ATAC_counts<-as.data.frame(fread(file=opt$ATAC_counts, sep="\t", header =T), stringAsFactors=F)
  
  colnames(ATAC_counts)[which(colnames(ATAC_counts) == 'GMP-A')]<-'GMP_A'
  colnames(ATAC_counts)[which(colnames(ATAC_counts) == 'GMP-B')]<-'GMP_B'
  colnames(ATAC_counts)[which(colnames(ATAC_counts) == 'GMP-C')]<-'GMP_C'
  
  
  cat("ATAC_counts_0\n")
  cat(str(ATAC_counts))
  cat("\n")

  
  ATAC_DEF<-cbind(ATAC_peaks,ATAC_counts)
  
  cat("ATAC_DEF_0\n")
  cat(str(ATAC_DEF))
  cat("\n")
  
           
  
  gr_ATAC_DEF <- GRanges(
    seqnames = as.character(gsub("chr","",ATAC_DEF$chr)),
    name2=as.character(ATAC_DEF$B),
    name3=as.character(ATAC_DEF$CD4),
    name4=as.character(ATAC_DEF$CD8),
    name5=as.character(ATAC_DEF$CLP),
    name6=as.character(ATAC_DEF$CMP),
    name7=as.character(ATAC_DEF$Ery),
    name8=as.character(ATAC_DEF$GMP_A),
    name9=as.character(ATAC_DEF$GMP_B),
    name10=as.character(ATAC_DEF$GMP_C),
    name11=as.character(ATAC_DEF$HSC),
    name12=as.character(ATAC_DEF$LMPP),
    name13=as.character(ATAC_DEF$mDC),
    name14=as.character(ATAC_DEF$Mega),
    name15=as.character(ATAC_DEF$MEP),
    name16=as.character(ATAC_DEF$Mono),
    name17=as.character(ATAC_DEF$MPP),
    name18=as.character(ATAC_DEF$NK),
    name19=as.character(ATAC_DEF$pDC),
    ranges=IRanges(
      start=ATAC_DEF$start,
      end=ATAC_DEF$end))
  

  #### READ GWAS block sizes ----
  
  BLOCK_POST<-readRDS(file=opt$BLOCK_POST)
  
 
  
  cat("BLOCK_POST_0\n")
  cat(str(BLOCK_POST))
  cat("\n")
 
 
  
  
  ###### LOOP AS_IDs -----
  
  path_files<-paste(out,'Build_files','/',sep='')
  
  if (file.exists(path_files)){
    
    
  }else{
    
    dir.create(file.path(path_files))
    
  }#path_files
  
  
  AS_ID_array<-unique(BLOCK_POST$AS_ID)
  
  cat("AS_ID_array_0\n")
  cat(str(AS_ID_array))
  cat("\n")
  
  

  
  DEBUG<-0
  
  
  Result_FINAL_ATAC<-data.frame()
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
  
    
    
    BLOCK_POST_sel<-BLOCK_POST[which(BLOCK_POST$AS_ID == AS_ID_sel),]
    
    
    if(DEBUG == 1)
    {
      cat("BLOCK_POST_sel_0\n")
      cat(str(BLOCK_POST_sel))
      cat("\n")
    }
    
    chr_sel<-gsub("^chr","",unique(BLOCK_POST_sel$chr))
    
    if(DEBUG == 1)
    {
      cat("chr_sel_0\n")
      cat(str(chr_sel))
      cat("\n")
    }
    
    
    temp<-data.frame()
    
    if(dim(BLOCK_POST_sel)[1] >0)
    {
      
      gr_BLOCK_POST <- GRanges(
        seqnames = as.character(gsub("^chr","",as.character(BLOCK_POST_sel$chr))),
        strand="*",
        ranges=IRanges(
          start=BLOCK_POST_sel$start,
          end=BLOCK_POST_sel$end,
          names=BLOCK_POST_sel$AS_ID))
      
      if(DEBUG == 1)
      {
        cat("gr_BLOCK_POST_0\n")
        cat(str(gr_BLOCK_POST))
        cat("\n")
      }
      
      #### find overlap with SNP ----
      
      m <- findOverlaps(gr_BLOCK_POST,gr_ATAC_DEF)
      
      if(DEBUG == 1)
      {
        cat("m\n")
        cat(str(m))
        cat("\n")
      }
      
      subjectHits_ATAC_DEF<-subjectHits(m)
      
      if(DEBUG == 1)
      {
        cat("subjectHits_ATAC_DEF\n")
        cat(str(subjectHits_ATAC_DEF))
        cat("\n")
      }
      
      queryHits_BLOCK_POST<-queryHits(m)
      
      if(DEBUG == 1)
      {
        cat("queryHits_BLOCK_POST\n")
        cat(str(queryHits_BLOCK_POST))
        cat("\n")
      }
      
      BLOCK_POST_df <- data.frame(chr=as.character(seqnames(gr_BLOCK_POST)),
                           start=as.integer(start(gr_BLOCK_POST)),
                           end=as.integer(end(gr_BLOCK_POST)),
                           AS_ID=names(gr_BLOCK_POST), stringsAsFactors = F)
      
      if(DEBUG == 1)
      {
        cat("BLOCK_POST_df_0\n")
        cat(str(BLOCK_POST_df))
        cat("\n")
      }
      
      BLOCK_POST_df_hits<-BLOCK_POST_df[queryHits_BLOCK_POST,]
      
      if(DEBUG == 1)
      {
        cat("BLOCK_POST_df_hits_0\n")
        cat(str(BLOCK_POST_df_hits))
        cat("\n")
      }
      
      ATAC_DEF_df <- data.frame(chr=as.character(seqnames(gr_ATAC_DEF)),
                               start=as.integer(start(gr_ATAC_DEF)),
                               end=as.integer(end(gr_ATAC_DEF)),
                               B=as.numeric(gr_ATAC_DEF$name2),
                               CD4=as.numeric(gr_ATAC_DEF$name3),
                               CD8=as.numeric(gr_ATAC_DEF$name4),
                               CLP=as.numeric(gr_ATAC_DEF$name5),
                               CMP=as.numeric(gr_ATAC_DEF$name6),
                               Ery=as.numeric(gr_ATAC_DEF$name7),
                               GMP_A=as.numeric(gr_ATAC_DEF$name8),
                               GMP_B=as.numeric(gr_ATAC_DEF$name9),
                               GMP_C=as.numeric(gr_ATAC_DEF$name10),
                               HSC=as.numeric(gr_ATAC_DEF$name11),
                               LMPP=as.numeric(gr_ATAC_DEF$name12),
                               mDC=as.numeric(gr_ATAC_DEF$name13),
                               Mega=as.numeric(gr_ATAC_DEF$name14),
                               MEP=as.numeric(gr_ATAC_DEF$name15),
                               Mono=as.numeric(gr_ATAC_DEF$name16),
                               MPP=as.numeric(gr_ATAC_DEF$name17),
                               NK=as.numeric(gr_ATAC_DEF$name18),
                               pDC=as.numeric(gr_ATAC_DEF$name19),
                               stringsAsFactors = F)
      
      
      if(DEBUG == 1)
      {
        cat("ATAC_DEF_df_0\n")
        cat(str(ATAC_DEF_df))
        cat("\n")
      }
      
      
      
      ATAC_DEF_df_hits<-ATAC_DEF_df[subjectHits_ATAC_DEF,]
      
      if(dim(ATAC_DEF_df_hits)[1] >0)
      {
        ATAC_DEF_df_hits$AS_ID<-AS_ID_sel
        
        if(DEBUG == 1)
        {
          cat("ATAC_DEF_df_hits_0\n")
          cat(str(ATAC_DEF_df_hits))
          cat("\n")
        }
        
        
        Result_FINAL_ATAC<-rbind(Result_FINAL_ATAC,ATAC_DEF_df_hits)
        
        
        ALL_POS<-c(BLOCK_POST_sel$start,BLOCK_POST_sel$end,BLOCK_POST_sel$END,ATAC_DEF_df_hits$start,ATAC_DEF_df_hits$end)
        
        
       
        
      }else{
        
        ALL_POS<-c(BLOCK_POST_sel$START,BLOCK_POST_sel$END)
        
        
        
      }#dim(ATAC_DEF_df_hits)[1] >0
      
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
        
      
      }
      
      Result_FINAL_boundaries<-rbind(Result_FINAL_boundaries,recomposed_block)
      
   
      
    }#dim(BLOCK_POST_sel)[1] >0
  }# i in 1:length(AS_ID_array)
  
  
 

  #### SAVE ----

 
  
  if(dim(Result_FINAL_ATAC)[1] >0)
  {
    Result_FINAL_ATAC$chr<-as.character(paste('chr',Result_FINAL_ATAC$chr, sep=''))
    
    Result_FINAL_ATAC$chr<-factor(Result_FINAL_ATAC$chr,
                                        levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                                                 "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                                 "chr22","chr23","chrX","chrY"), ordered=T)
    
    Result_FINAL_ATAC<-droplevels(Result_FINAL_ATAC)
    
    Result_FINAL_ATAC<-Result_FINAL_ATAC[order(Result_FINAL_ATAC$chr),]
    
    cat("Result_FINAL_ATAC_0\n")
    cat(str(Result_FINAL_ATAC))
    cat("\n")
    cat(str(unique(Result_FINAL_ATAC$AS_ID)))
    cat("\n")
    cat(str(unique(Result_FINAL_ATAC$ensembl_gene_id)))
    cat("\n")
    
    setwd(path_files)
    
    saveRDS(Result_FINAL_ATAC, file=paste('Intersected_ATAC',".rds",sep=''))
    
    
    
  }#dim(Result_FINAL_ATAC)[1] >0
  
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
    
    saveRDS(Result_FINAL_boundaries, file=paste('Boundaries_AS_POST_genes_POST_ATAC',".rds",sep=''))
    
    
    
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
    make_option(c("--ATAC_peaks"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ATAC_counts"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--BLOCK_POST"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ensembl_gtf"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ensembl_gtf_sel_CORRECTED"), type="character", default=NULL, 
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