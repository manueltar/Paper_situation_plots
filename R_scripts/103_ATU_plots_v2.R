
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
suppressMessages(library("BiocGenerics", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("S4Vectors", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("IRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomeInfoDb", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GenomicRanges", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("Biobase", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("AnnotationDbi", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("GO.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("org.Hs.eg.db", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("TxDb.Hsapiens.UCSC.hg19.knownGene", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rtracklayer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("R.oo", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggtranscript", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggpubr", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


opt = NULL

options(warn = 1)

data_wrangling = function(option_list)
{

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
  
  #### READ and transform path_INTERVAL ----
  
  path_INTERVAL = opt$path_INTERVAL
  
  cat("path_INTERVAL_\n")
  cat(sprintf(as.character(path_INTERVAL)))
  cat("\n")
  
  #### READ and transform path_BP ----
  
  path_BP = opt$path_BP
  
  cat("path_BP_\n")
  cat(sprintf(as.character(path_BP)))
  cat("\n")
  
  #### Read transposed expression ----
  
  
  
  Transposed_Isoform_Expression<-readRDS(file=opt$Transposed_Isoform_Expression)
  colnames(Transposed_Isoform_Expression)[which(colnames(Transposed_Isoform_Expression) == "Sample_id")]<-"sample_id"
  
  # cat("Transposed_Isoform_Expression_0\n")
  # cat(str(Transposed_Isoform_Expression))
  # cat("\n")
  
  #### READ ensembl_gtf ----
  
  ensembl_gtf = readGFF(opt$ensembl_gtf)
  
  
  # cat("ensembl_gtf_0\n")
  # cat(str(ensembl_gtf))
  # cat("\n")
  # cat(sprintf(as.character(names(summary(as.factor(ensembl_gtf$type))))))
  # cat("\n")
  # cat(sprintf(as.character(summary(as.factor(ensembl_gtf$type)))))
  # cat("\n")
  
  
  #### READ Table_S6 ----
  
  Table_S6<-readRDS(file=opt$Table_S6)
  
  
  cat("Table_S6_0\n")
  cat(str(Table_S6))
  cat("\n")
  cat(str(unique(Table_S6$VAR)))
  cat("\n")
  
  Table_S6_subset<-Table_S6[which(Table_S6$Mechanistic_Class%in%c('ATU','DE + ATU')),]
  
  cat("Table_S6_subset_0\n")
  cat(str(Table_S6_subset))
  cat("\n")
  cat(str(unique(Table_S6_subset$VAR)))
  cat("\n")
  
  
  Table_S6_subset<-unique(as.data.frame(cSplit(Table_S6_subset,sep = ';', direction = "long",
                                              splitCols = "Whole_blood_DTU_HGNC_string"),stringsAsFactors=F))
  
  # cat("Table_S6_subset_1\n")
  # cat(str(Table_S6_subset))
  # cat("\n")
  # cat(sprintf(as.character(unique(Table_S6_subset$Whole_blood_DTU_HGNC_string))))
  # cat("\n")
  
  Table_S6_subset<-unique(as.data.frame(cSplit(Table_S6_subset,sep = ';', direction = "long",
                                               splitCols = "Neutrophil_DTU_HGNC_string"),stringsAsFactors=F))
  
  # cat("Table_S6_subset_2\n")
  # cat(str(Table_S6_subset))
  # cat("\n")
  # cat(sprintf(as.character(unique(Table_S6_subset$Neutrophil_DTU_HGNC_string))))
  # cat("\n")
 
  Table_S6_subset<-unique(as.data.frame(cSplit(Table_S6_subset,sep = ';', direction = "long",
                                               splitCols = "Monocyte_DTU_HGNC_string"),stringsAsFactors=F))
  
  # cat("Table_S6_subset_3\n")
  # cat(str(Table_S6_subset))
  # cat("\n")
  # cat(sprintf(as.character(unique(Table_S6_subset$Monocyte_DTU_HGNC_string))))
  # cat("\n")
  
  Table_S6_subset<-unique(as.data.frame(cSplit(Table_S6_subset,sep = ';', direction = "long",
                                               splitCols = "Tcell_DTU_HGNC_string"),stringsAsFactors=F))
  
  # cat("Table_S6_subset_4\n")
  # cat(str(Table_S6_subset))
  # cat("\n")
  # cat(sprintf(as.character(unique(Table_S6_subset$Tcell_DTU_HGNC_string))))
  # cat("\n")
  
  Table_S6_subset_ATU<-unique(Table_S6_subset[,c(which(colnames(Table_S6_subset) == 'VAR'),which(colnames(Table_S6_subset) == 'rs'),
                                          which(colnames(Table_S6_subset) == 'Whole_blood_DTU_HGNC_string'),
                                          which(colnames(Table_S6_subset) == 'Neutrophil_DTU_HGNC_string'),
                                          which(colnames(Table_S6_subset) == 'Monocyte_DTU_HGNC_string'),
                                          which(colnames(Table_S6_subset) == 'Tcell_DTU_HGNC_string'))])
  
  
  # cat("Table_S6_subset_ATU_0\n")
  # cat(str(Table_S6_subset_ATU))
  # cat("\n")
  
  Table_S6_subset_ATU.m<-melt(Table_S6_subset_ATU, id.vars=c('VAR','rs'), variable.name='RNASeq_source', value.name='HGNC')
  
  # cat("Table_S6_subset_ATU.m_0\n")
  # cat(str(Table_S6_subset_ATU.m))
  # cat("\n")
  # cat(sprintf(as.character(unique(Table_S6_subset_ATU.m$HGNC))))
  # cat("\n")
  
  Table_S6_subset_ATU.m$RNASeq_source<-gsub("_DTU_HGNC_string","",Table_S6_subset_ATU.m$RNASeq_source)
  
  cat(sprintf(as.character(names(summary(as.factor(Table_S6_subset_ATU.m$RNASeq_source))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_S6_subset_ATU.m$RNASeq_source)))))
  cat("\n")
  
  
  # check<-Table_S6_subset_ATU.m[which(Table_S6_subset_ATU.m$VAR == 'chr3_71355240_G_C'),]
  # 
  # cat("check_0\n")
  # cat(str(check))
  # cat("\n")
  # cat(sprintf(as.character(unique(check$HGNC))))
  # cat("\n")
  
  
  
  
  Table_S6_subset_ATU.m_NO_NA<-Table_S6_subset_ATU.m[which(Table_S6_subset_ATU.m$HGNC != "NA"),]
  
  cat("Table_S6_subset_ATU.m_NO_NA_0\n")
  cat(str(Table_S6_subset_ATU.m_NO_NA))
  cat("\n")
  cat(sprintf(as.character(unique(Table_S6_subset_ATU.m_NO_NA$HGNC))))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_S6_subset_ATU.m_NO_NA$RNASeq_source))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_S6_subset_ATU.m_NO_NA$RNASeq_source)))))
  cat("\n")
  
  Table_S6_subset_ATU.m_NO_NA<-unique(Table_S6_subset_ATU.m_NO_NA[,-which(colnames(Table_S6_subset_ATU.m_NO_NA) == 'RNASeq_source')])
  
  
 
  # check<-Table_S6_subset_ATU.m_NO_NA[which(Table_S6_subset_ATU.m_NO_NA$VAR == 'chr3_71355240_G_C'),]
  # 
  # cat("check_2\n")
  # cat(str(check))
  # cat("\n")
  # cat(sprintf(as.character(unique(check$HGNC))))
  # cat("\n")
  
  #### READ Table_S7 ----
  
  Table_S7<-readRDS(file=opt$Table_S7)
  
  cat("Table_S7_0\n")
  cat(str(Table_S7))
  cat("\n")
  cat(str(unique(Table_S7$VAR)))
  cat("\n")
  cat(str(unique(Table_S7$ensembl_gene_id)))
  cat("\n")
  
  Table_S7$Analysis[which(Table_S7$Analysis == 'DTU')]<-'ATU'
  
  Table_S7_ATU<-Table_S7[which(Table_S7$Analysis == 'ATU'),]
  
 
  cat("Table_S7_ATU_0\n")
  cat(str(Table_S7_ATU))
  cat("\n")
  cat(str(unique(Table_S7_ATU$VAR)))
  cat("\n")
  cat(str(unique(Table_S7_ATU$ensembl_gene_id)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_S7_ATU$RNASeq_source))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_S7_ATU$RNASeq_source)))))
  cat("\n")
  # 
  # check<-Table_S7_ATU[which(Table_S7_ATU$VAR == 'chr3_71355240_G_C'),]
  # 
  # cat("check_3\n")
  # cat(str(check))
  # cat("\n")
  # cat(sprintf(as.character(unique(check$HGNC))))
  # cat("\n")
  # cat(str(unique(check$transcript_id)))
  # cat("\n")
  
  
  Table_S7_ATU_subset<-merge(Table_S7_ATU,
                             Table_S6_subset_ATU.m_NO_NA,
                             by=c("VAR",'rs',"HGNC"))
  
  
  
  cat("Table_S7_ATU_subset_0\n")
  cat(str(Table_S7_ATU_subset))
  cat("\n")
  cat(str(unique(Table_S7_ATU_subset$VAR)))
  cat("\n")
  cat(str(unique(Table_S7_ATU_subset$ensembl_gene_id)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_S7_ATU_subset$RNASeq_source))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_S7_ATU_subset$RNASeq_source)))))
  cat("\n")
  cat(sprintf(as.character(unique(Table_S7_ATU_subset$HGNC))))
  cat("\n")
  
  # check<-Table_S7_ATU_subset[which(Table_S7_ATU_subset$VAR == 'chr3_71355240_G_C'),]
  # 
  # cat("check_4\n")
  # cat(str(check))
  # cat("\n")
  # cat(sprintf(as.character(unique(check$HGNC))))
  # cat("\n")
  # cat(str(unique(check$transcript_id)))
  # cat("\n")
  
  
 
  
  #### subset the ensembl gtf ----
  
  
  ensembl_gtf_sel<-unique(ensembl_gtf[which(ensembl_gtf$gene_id%in%Table_S7_ATU_subset$ensembl_gene_id &
                                       ensembl_gtf$type%in%c('gene','transcript','exon','CDS','five_prime_utr','three_prime_utr')),])
  
  # cat("ensembl_gtf_sel_0\n")
  # cat(str(ensembl_gtf_sel))
  # cat("\n")
  # cat(str(unique(ensembl_gtf_sel$gene_id)))
  # cat("\n")
  # cat(sprintf(as.character(names(summary(as.factor(ensembl_gtf_sel$transcript_biotype))))))
  # cat("\n")
  # cat(sprintf(as.character(summary(as.factor(ensembl_gtf_sel$transcript_biotype)))))
  # cat("\n")
  
  
  
  ######## LOOP TO PRINT ------
  
  
  path_graphs<-paste(out,'ATU_graphs','/',sep='')
  
  if (file.exists(path_graphs)){
    
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  
  
  VARS<-unique(Table_S7_ATU_subset$VAR)
  
  # VARS<-'chr7_101499930_G_A'
  
  # VARS<-'chr1_92925654_G_C'
  
  cat("VARS_0\n")
  cat(str(VARS))
  cat("\n")
  
  DEBUG<-1
  
  for(i in 1:length(VARS))
  {
    VAR_sel<-VARS[i]
    
    
    
    Table_S7_ATU_subset_sel<-Table_S7_ATU_subset[which(Table_S7_ATU_subset$VAR == VAR_sel),]
    
    if(DEBUG == 1)
    {
      cat("Table_S7_ATU_subset_sel_0\n")
      cat(str(Table_S7_ATU_subset_sel))
      cat("\n")
      
    }
    
    rs_sel<-unique(as.character(Table_S7_ATU_subset_sel$rs))
    
    if(DEBUG == 1)
    {
      cat("rs_sel_0\n")
      cat(str(rs_sel))
      cat("\n")
    }
    
    ENSG_array<-unique(Table_S7_ATU_subset_sel$ensembl_gene_id)
    
    if(DEBUG == 1)
    {
      cat("ENSG_array_0\n")
      cat(str(ENSG_array))
      cat("\n")
    }
    
    path_graphs<-paste(out,'ATU_graphs','/',paste(rs_sel,VAR_sel, sep='__'),'/',sep='')
    
    if (file.exists(path_graphs)){
      
      
    }else{
      
      dir.create(file.path(path_graphs))
      
    }#path_graphs
    
    for(k in 1:length(ENSG_array))
    {
      ENSG_array_sel<-ENSG_array[k]
      
      Table_S7_ATU_subset_sel_ENSG_sel<-Table_S7_ATU_subset_sel[which(Table_S7_ATU_subset_sel$ensembl_gene_id == ENSG_array_sel),]
      
      Table_S7_ATU_subset_sel_ENSG_sel$REF<-NA
      
      Table_S7_ATU_subset_sel_ENSG_sel$REF[is.na(Table_S7_ATU_subset_sel_ENSG_sel$adjusted_minus_logpval)]<-'reference'
      Table_S7_ATU_subset_sel_ENSG_sel$REF[!is.na(Table_S7_ATU_subset_sel_ENSG_sel$adjusted_minus_logpval)]<-'not_reference'
      
      Table_S7_ATU_subset_sel_ENSG_sel$REF<-factor(Table_S7_ATU_subset_sel_ENSG_sel$REF,
                                                   levels=c('not_reference','reference'),
                                                   ordered=T)
      
      
      Table_S7_ATU_subset_sel_ENSG_sel[order(Table_S7_ATU_subset_sel_ENSG_sel$RNASeq_source,Table_S7_ATU_subset_sel_ENSG_sel$REF),]
      
      order_transcript_id<-unique(Table_S7_ATU_subset_sel_ENSG_sel$transcript_id)
      
      newG<-c("HOM_REF","HET")
      
      if(DEBUG == 1)
      {
        cat("Table_S7_ATU_subset_sel_ENSG_sel_0\n")
        cat(str(Table_S7_ATU_subset_sel_ENSG_sel))
        cat("\n")
        cat(str(unique(Table_S7_ATU_subset_sel_ENSG_sel$ensembl_gene_id)))
        cat("\n")
        cat(sprintf(as.character(order_transcript_id)))
        cat("\n")
        
        cat("ORDER transcript_id\n")
        cat(sprintf(as.character(Table_S7_ATU_subset_sel_ENSG_sel$transcript_id)))
        cat("\n")
        
        
      }
      
      ensembl_gtf_sel_ENSG_sel<-ensembl_gtf_sel[which(ensembl_gtf_sel$gene_id == ENSG_array_sel &
                                                        ensembl_gtf_sel$type == 'gene'),]
      
      if(DEBUG == 1)
      {
        cat("ensembl_gtf_sel_ENSG_sel_0\n")
        cat(str(ensembl_gtf_sel_ENSG_sel))
        cat("\n")
      }
      
      HGNC_sel<-unique(ensembl_gtf_sel_ENSG_sel$gene_name)
      
      if(DEBUG == 1)
      {
        cat("HGNC_sel_0\n")
        cat(str(HGNC_sel))
        cat("\n")
      }
      
      
      chr_sel<-unique(as.character(ensembl_gtf_sel_ENSG_sel$seqid))
      
      if(DEBUG == 1)
      {
        cat("chr_sel_0\n")
        cat(str(chr_sel))
        cat("\n")
      }
      
      start_sel<-unique(ensembl_gtf_sel_ENSG_sel$start)
      
      if(DEBUG == 1)
      {
        cat("start_sel_0\n")
        cat(str(start_sel))
        cat("\n")
      }
      
      end_sel<-unique(ensembl_gtf_sel_ENSG_sel$end)
      
      if(DEBUG == 1)
      {
        cat("end_sel_0\n")
        cat(str(end_sel))
        cat("\n")
      }
      
      strand_sel<-unique(ensembl_gtf_sel_ENSG_sel$strand)
      
      if(DEBUG == 1)
      {
        cat("strand_sel_0\n")
        cat(str(strand_sel))
        cat("\n")
      }
      
      cat("---------------------------------------------------------------------------------->\t")
      cat(sprintf(as.character(i)))
      cat("\t")
      cat(sprintf(as.character(rs_sel)))
      cat("\t")
      cat(sprintf(as.character(VAR_sel)))
      cat("\t")
      cat(sprintf(as.character(k)))
      cat("\t")
      cat(sprintf(as.character(ENSG_array_sel)))
      cat("\t")
      cat(sprintf(as.character(HGNC_sel)))
      cat("\t")
      cat(sprintf(as.character(chr_sel)))
      cat("\t")
      cat(sprintf(as.character(start_sel)))
      cat("\t")
      cat(sprintf(as.character(end_sel)))
      cat("\t")
      cat(sprintf(as.character(strand_sel)))
      cat("\t")
      cat("\n")
      
      path_graphs<-paste(out,'ATU_graphs','/',paste(rs_sel,VAR_sel, sep='__'),'/',HGNC_sel,'/',sep='')
      
      if (file.exists(path_graphs)){
        
        
      }else{
        
        dir.create(file.path(path_graphs))
        
      }#path_graphs
      
      
      #### transcript plot ----
      
      transcript_id_array<-unique(Table_S7_ATU_subset_sel_ENSG_sel$transcript_id)
      
      if(DEBUG == 1)
      {
        cat("transcript_id_array_0\n")
        cat(str(transcript_id_array))
        cat("\n")
      }
      
      
      ensembl_gtf_sel_transcript_sel<-ensembl_gtf_sel[which(ensembl_gtf_sel$transcript_id %in% transcript_id_array &
                                                              ensembl_gtf_sel$type == 'transcript'),]
      
      if(DEBUG == 1)
      {
        cat("ensembl_gtf_sel_transcript_sel_0\n")
        cat(str(ensembl_gtf_sel_transcript_sel))
        cat("\n")
      }
      
      
      check<-Table_S7_ATU_subset_sel_ENSG_sel[-which(Table_S7_ATU_subset_sel_ENSG_sel$transcript_id%in%ensembl_gtf_sel_transcript_sel$transcript_id),]
      
      if(DEBUG == 1)
      {
        cat("check_0\n")
        cat(str(check))
        cat("\n")
      }
      
      if(dim(check)[1]>0)
      {
        
        ### None
        
      }#dim(check)[1]>0
      
      
      REP<-data.frame()
      for(l in 1:length(transcript_id_array))
      {
        
        transcript_id_array_sel<-transcript_id_array[l]
        
        
        # cat("---------------------------------------------------------------------------------->\t")
        # cat(sprintf(as.character(l)))
        # cat("\t")
        # cat(sprintf(as.character(transcript_id_array_sel)))
        # cat("\n")
        
        ensembl_gtf_sel_transcript_sel<-ensembl_gtf_sel[which(ensembl_gtf_sel$transcript_id == transcript_id_array_sel &
                                                                ensembl_gtf_sel$type != 'gene'),]
        # 
        # if(DEBUG == 1)
        # {
        #   cat("ensembl_gtf_sel_transcript_sel_0\n")
        #   cat(str(ensembl_gtf_sel_transcript_sel))
        #   cat("\n")
        # }
        
        ensembl_gtf_sel_transcript_sel$COORD<-l
        
        REP<-rbind(ensembl_gtf_sel_transcript_sel,REP)
        
      }#l in 1:length(transcript_id_array)
      
      
      colnames(REP)[which(colnames(REP) == 'seqid')]<-'seqnames'
      
      
      if(DEBUG == 1)
      {
        cat("REP_0\n")
        cat(str(REP))
        cat("\n")
      }
      
      REP$transcript_biotype<-factor(REP$transcript_biotype,
                                     levels=c("protein_coding","retained_intron","nonsense_mediated_decay","processed_transcript","antisense"),
                                     ordered=T)
      
      REP$transcript_id<-factor(REP$transcript_id,
                                levels=rev(order_transcript_id),
                                ordered=T)
      
      
      
      if(DEBUG == 1)
      {
        cat("REP_1\n")
        cat(str(REP))
        cat("\n")
      }
      
      
      vector.fill<-c(brewer.pal(length(levels(REP$transcript_biotype)), "Dark2"))
      
      if(DEBUG == 1)
      {
        cat("vector.fill_0\n")
        cat(str(vector.fill))
        cat("\n")
      }
      
      
      
      # extract exons
      
      REP_exons <- REP %>% dplyr::filter(type == "exon")
      
      REP_exons_rescaled <- shorten_gaps(
        REP_exons, 
        to_intron(REP_exons, "transcript_name"), 
        group_var = "transcript_name"
      )
      
      transcript_plot<-REP_exons_rescaled %>%
        dplyr::filter(type == "exon") %>%
        ggplot(aes(
          xstart = start,
          xend = end,
          y = transcript_id
        )) +
        geom_range(
          aes(fill = transcript_biotype)
        ) +
        geom_intron(
          data = REP_exons_rescaled %>% dplyr::filter(type == "intron"),
          aes(strand=strand),
          arrow = grid::arrow(ends = "last", length = grid::unit(0.25, "cm")),
          arrow.min.intron.length = 200
        )
      
      transcript_plot<-transcript_plot+
        scale_fill_manual(values=vector.fill,drop=F)+
        theme_classic()+
        theme(axis.title.y=element_blank(),
              axis.title.x=element_blank(),
              axis.text.y=element_text(angle=0,size=6, color="black", family="sans"),
              axis.text.x=element_text(angle=0,size=6, color="black", family="sans"),
              axis.line.x = element_line(size = 0.2),
              axis.ticks = element_line(size = 0.2),
              axis.line.y = element_line(size = 0.2))+
        theme(legend.title = element_blank(),
              legend.text = element_text(size=6),
              legend.key.size = unit(0.25, 'cm'), #change legend key size
              legend.key.height = unit(0.25, 'cm'), #change legend key height
              legend.key.width = unit(0.25, 'cm'), #change legend key width
              legend.position="bottom")
      
      # scale_x_continuous(name=NULL,breaks=c(start_sel,end_sel),
      #                    labels=as.character(c(start_sel,end_sel)),
      #                    limits=c(start_sel-1,end_sel+1))+
      
      setwd(path_graphs)
      
      svgname<-paste(paste('Transcripts_analysed', sep='_'),".svg",sep='')
      makesvg = TRUE
      
      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= transcript_plot,
               device="svg",
               height=4, width=4)
      }
      #### violin plots of expression  BP----
      
      
      
      BLUEPRINT_files<-paste(path_BP,VAR_sel,'/', sep='')
      filename<-paste(BLUEPRINT_files,'BP_Transcript_EXP_for_LM.rds',sep='')
      if(DEBUG == 1)
      {
        cat("BLUEPRINT_files\n")
        cat(sprintf(as.character(BLUEPRINT_files)))
        cat("\n")
      }
      
      
      
      if (file.exists(filename)){
        
        setwd(BLUEPRINT_files)
        
        EXP_BP_df = readRDS(file='BP_Transcript_EXP_for_LM.rds')
        
        if(DEBUG == 1)
        {
          cat("EXP_BP_df\n")
          cat(str(EXP_BP_df))
          cat("\n")
          cat(str(unique(EXP_BP_df$ensembl_gene_id)))
          cat("\n")
          cat(str(unique(EXP_BP_df$transcript_id)))
          cat("\n")
        }
        
        EXP_BP_df_sel<-EXP_BP_df[which(EXP_BP_df$ensembl_gene_id == ENSG_array_sel &
                                                     EXP_BP_df$transcript_id%in%order_transcript_id),]
        
        if(DEBUG == 1)
        {
          cat("EXP_BP_df_sel\n")
          cat(str(EXP_BP_df_sel))
          cat("\n")
          cat(str(unique(EXP_BP_df_sel$ensembl_gene_id)))
          cat("\n")
          cat(str(unique(EXP_BP_df_sel$transcript_id)))
          cat("\n")
          
          cat("EXP_BP_df_sel$transcript_id\n")
          cat(sprintf(as.character(unique(EXP_BP_df_sel$transcript_id))))
          cat("\n")
          
          cat("order_transcript_id\n")
          cat(sprintf(as.character(order_transcript_id)))
          cat("\n")
          
          cat("newG\n")
          cat(sprintf(as.character(newG)))
          cat("\n")
        }
        
      
        EXP_BP_df_sel$FPKM<-2^(EXP_BP_df_sel$value)
        
        EXP_BP_df_sel$transcript_id<-factor(EXP_BP_df_sel$transcript_id,
                                                  levels=order_transcript_id,
                                                  ordered=T)
        
        EXP_BP_df_sel$Cell_Type<-factor(EXP_BP_df_sel$Cell_Type, 
                                              levels=c("Monocyte","Neutrophil","Tcell"), 
                                              ordered=T)
        
        
        EXP_BP_df_sel$Genotype<-factor(EXP_BP_df_sel$Genotype,
                                             levels=newG,
                                             ordered=T)
        
        
        if(DEBUG == 1)
        {
          cat("EXP_BP_df_sel\n")
          cat(str(EXP_BP_df_sel))
          cat("\n")
          cat(sprintf(as.character(names(summary(EXP_BP_df_sel$Cell_Type)))))
          cat("\n")
          cat(sprintf(as.character(summary(EXP_BP_df_sel$Cell_Type))))
          cat("\n")
          cat(sprintf(as.character(names(summary(EXP_BP_df_sel$transcript_id)))))
          cat("\n")
          cat(sprintf(as.character(summary(EXP_BP_df_sel$transcript_id))))
          cat("\n")
          cat(sprintf(as.character(names(summary(EXP_BP_df_sel$Genotype)))))
          cat("\n")
          cat(sprintf(as.character(summary(EXP_BP_df_sel$Genotype))))
          cat("\n")
        }
        
        
        A<-round(summary(EXP_BP_df_sel$value[!is.na(EXP_BP_df_sel$value)]),2)
        
        
        # cat("summary_GENE_EXP\n")
        # cat(sprintf(as.character(names(A))))
        # cat("\n")
        # cat(sprintf(as.character(A)))
        # cat("\n")
        
        step<-abs(A[6]-A[1])/4
        
        if(step == 0)
        {
          
          step<-1
        }
        
        breaks.Rank<-sort(unique(c(A[6],seq(from= A[1], to=A[6],by=step))))
        labels.Rank<-as.character(round(breaks.Rank,0))
        
        
        if(DEBUG == 1)
        {
          cat("labels.Rank:\t")
          cat(sprintf(as.character(labels.Rank)))
          cat("\n")
          
          
          
          # quit(status = 1)
        }
        
        
        
        
        
        
        
        BLUEPRINT_violin<-ggplot(data=EXP_BP_df_sel,
                                 aes(y=Genotype, x=value, fill=Genotype)) +
          geom_violin()+
          stat_summary(fun = median, fun.min = median, fun.max = median,
                       geom = "crossbar", width = 0.5)+
          scale_x_continuous(name="log2 FPKM",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
          scale_fill_manual(values=c("#E7B800","#00AFBB","#FC4E07"),drop=F)+
          scale_y_discrete(name=NULL, drop=F)
        
        
        
        BLUEPRINT_violin<-BLUEPRINT_violin+
          facet_grid(transcript_id ~ Cell_Type, scales='free_x', space='free_x', switch="y", drop=F)+
          theme_cowplot(font_size = 4)+
          theme( strip.background = element_blank(),
                 strip.placement = "outside",
                 strip.text = element_text(size=6),
                 panel.spacing = unit(0.2, "lines"),
                 panel.background=element_rect(fill="white"),
                 panel.border=element_rect(colour="white",size=0,5),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())+
          theme_classic()+
          theme(axis.title.y=element_blank(),
                axis.title.x=element_text(size=8, color="black", family="sans"),
                axis.text.y=element_blank(),
                axis.text.x=element_text(angle=0,size=6, color="black", family="sans"),
                axis.line.x = element_line(size = 0.2),
                axis.ticks.x = element_blank(),
                axis.ticks.y = element_line(size = 0.2),
                axis.line.y = element_line(size = 0.2))+
          theme(legend.title = element_text(size=6),
                legend.text = element_text(size=6),
                legend.key.size = unit(0.25, 'cm'), #change legend key size
                legend.key.height = unit(0.25, 'cm'), #change legend key height
                legend.key.width = unit(0.25, 'cm'), #change legend key width
                legend.position="bottom")+
          ggeasy::easy_center_title()
        
        
        setwd(path_graphs)
        
        svgname<-paste(paste('BluePrint_violin', sep='_'),".svg",sep='')
        makesvg = TRUE
        
        if (makesvg == TRUE)
        {
          ggsave(svgname, plot= BLUEPRINT_violin,
                 device="svg",
                 height=4, width=6)
        }
        
        
        array_transcript_id<-unique(EXP_BP_df_sel$transcript_id)
        
        
        if(DEBUG == 1)
        {
          cat("array_transcript_id_0\n")
          cat(str(array_transcript_id))
          cat("\n")
        }
        
        for(iteration_array_transcript_id in 1:length(array_transcript_id)){
          
          array_transcript_id_sel<-array_transcript_id[iteration_array_transcript_id]
          
          if(DEBUG == 1)
          {
            cat(sprintf(as.character(iteration_array_transcript_id)))
            cat("\t")
            cat(sprintf(as.character(array_transcript_id_sel)))
            cat("\n")
          }
          
          EXP_BP_df_sel_transcript_id_sel<-droplevels(EXP_BP_df_sel[which(EXP_BP_df_sel$transcript_id == array_transcript_id_sel),])
          
          if(DEBUG == 1)
          {
            cat("EXP_BP_df_sel_transcript_id_sel_0\n")
            cat(str(EXP_BP_df_sel_transcript_id_sel))
            cat("\n")
          }
          
          
          array_Cell_Type<-unique(EXP_BP_df_sel_transcript_id_sel$Cell_Type)
          
          if(DEBUG == 1)
          {
            cat("array_Cell_Type_0\n")
            cat(str(array_Cell_Type))
            cat("\n")
          }
          
          for(iteration_array_Cell_Type in 1:length(array_Cell_Type)){
            
            array_Cell_Type_sel<-array_Cell_Type[iteration_array_Cell_Type]
            
            if(DEBUG == 1)
            {
              cat(sprintf(as.character(iteration_array_Cell_Type)))
              cat("\t")
              cat(sprintf(as.character(array_Cell_Type_sel)))
              cat("\n")
            }
            
            EXP_BP_df_sel_transcript_id_sel_Cell_Type_sel<-droplevels(EXP_BP_df_sel_transcript_id_sel[which(EXP_BP_df_sel_transcript_id_sel$Cell_Type == array_Cell_Type_sel),])
            
            if(DEBUG == 1)
            {
              cat("EXP_BP_df_sel_transcript_id_sel_Cell_Type_sel_0\n")
              cat(str(EXP_BP_df_sel_transcript_id_sel_Cell_Type_sel))
              cat("\n")
            }
            
            A<-round(summary(EXP_BP_df_sel_transcript_id_sel_Cell_Type_sel$value[!is.na(EXP_BP_df_sel_transcript_id_sel_Cell_Type_sel$value)]),2)
            
            
            # cat("summary_GENE_EXP\n")
            # cat(sprintf(as.character(names(A))))
            # cat("\n")
            # cat(sprintf(as.character(A)))
            # cat("\n")
            
            step<-abs(A[6]-A[1])/4
            
            if(step == 0)
            {
              
              step<-1
            }
            
            breaks.Rank<-sort(unique(c(0,A[6],seq(from= A[1], to=A[6],by=step))))
            labels.Rank<-as.character(round(breaks.Rank,0))
            
            
            if(DEBUG == 1)
            {
              cat("labels.Rank:\t")
              cat(sprintf(as.character(labels.Rank)))
              cat("\n")
              
              
              
              # quit(status = 1)
            }
            
            
            
            BP_per_transcript_violin<-ggplot(data=EXP_BP_df_sel_transcript_id_sel_Cell_Type_sel,
                                             aes(y=value, x=Genotype, fill=Genotype)) +
              geom_violin()+
              stat_summary(fun = median, fun.min = median, fun.max = median,
                           geom = "crossbar", width = 0.5)+
              scale_y_continuous(name="log2 FPKM",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
              scale_fill_manual(values=c("#E7B800","#00AFBB","#FC4E07"),drop=F)+
              scale_x_discrete(name=NULL, drop=F)
            
            
            
            BP_per_transcript_violin<-BP_per_transcript_violin+
              facet_grid(transcript_id ~ Cell_Type, scales='free_x', space='free_x', switch="y", drop=F)+
              theme_cowplot(font_size = 4)+
              theme( strip.background = element_blank(),
                     strip.placement = "outside",
                     strip.text = element_text(size=6),
                     panel.spacing = unit(0.2, "lines"),
                     panel.background=element_rect(fill="white"),
                     panel.border=element_rect(colour="white",size=0,5),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())+
              theme_classic()+
              theme(axis.title.y=element_text(size=8, color="black", family="sans"),
                    axis.title.x=element_blank(),
                    axis.text.y=element_text(angle=0,size=6, color="black", family="sans"),
                    axis.text.x=element_blank(),
                    axis.line.x = element_line(size = 0.2),
                    axis.ticks.x = element_blank(),
                    axis.ticks.y = element_line(size = 0.2),
                    axis.line.y = element_line(size = 0.2))+
              theme(legend.title = element_text(size=6),
                    legend.text = element_text(size=6),
                    legend.key.size = unit(0.25, 'cm'), #change legend key size
                    legend.key.height = unit(0.25, 'cm'), #change legend key height
                    legend.key.width = unit(0.25, 'cm'), #change legend key width
                    legend.position="hidden")+
              ggeasy::easy_center_title()
            
            
            setwd(path_graphs)
            
            svgname<-paste(paste('BP_per_transcript_violin',array_transcript_id_sel,array_Cell_Type_sel, sep='_'),".svg",sep='')
            makesvg = TRUE
            
            if (makesvg == TRUE)
            {
              ggsave(svgname, plot= BP_per_transcript_violin,
                     device="svg",
                     height=2, width=2)
            }
            
          }#iteration_array_Cell_Type in 1:length(array_Cell_Type)
        }#iteration_array_transcript_id in 1:length(array_transcript_id
        
        
       
      }else{
        
        BLUEPRINT_violin<-NA
        
        }#file.exists(BLUEPRINT_files)
      #### violin plots of expression  BLOOD----
      
     
      
      INTERVAL_files<-paste(path_INTERVAL,VAR_sel,'/', sep='')
      
      if(DEBUG == 1)
      {
        cat("INTERVAL_files\n")
        cat(sprintf(as.character(INTERVAL_files)))
        cat("\n")
      }
      
      if (file.exists(INTERVAL_files)){
        
        setwd(INTERVAL_files)
        
        filename<-paste("INTERVAL_covariates_and_PEER_factors_",VAR_sel,".rds", sep='')
        
        if (file.exists(filename)){
          
          INTERVAL_covariates_and_PEER_factors_sel<-readRDS(file=filename)
          
          if(DEBUG == 1)
          {
            cat("INTERVAL_covariates_and_PEER_factors_sel\n")
            cat(str(INTERVAL_covariates_and_PEER_factors_sel))
            cat("\n")
          }
          
          Transposed_Isoform_Expression_sel<-Transposed_Isoform_Expression[,c(which(colnames(Transposed_Isoform_Expression) == "sample_id"),
                                                                              which(colnames(Transposed_Isoform_Expression) %in% order_transcript_id))]
          if(DEBUG == 1)
          {
            cat("Transposed_Isoform_Expression_sel_1\n")
            cat(str(Transposed_Isoform_Expression_sel))
            cat("\n")
          }
          
          Transposed_Isoform_Expression_sel.m<-melt(Transposed_Isoform_Expression_sel, id.vars=c("sample_id"),value.name = "log2TPM",variable.name = "transcript_id")
          Transposed_Isoform_Expression_sel.m$TPM<--1+2^(Transposed_Isoform_Expression_sel.m$log2TPM)
          
          
          if(DEBUG == 1)
          {
            cat("Transposed_Isoform_Expression_sel.m_0\n")
            cat(str(Transposed_Isoform_Expression_sel.m))
            cat("\n")
            
            cat("sample_id\n")
            cat(str(unique(Transposed_Isoform_Expression_sel.m$sample_id)))
            cat("\n")
            
            cat("transcript_id\n")
            cat(str(unique(Transposed_Isoform_Expression_sel.m$transcript_id)))
            cat("\n")
            
            # ########################################
            # quit(status = 1)
          }
          
          ### calculate summatory TPM per sample
          
          Transposed_Isoform_Expression_sel.m.dt<-data.table(Transposed_Isoform_Expression_sel.m,
                                                             key=c("sample_id"))
          
          if(DEBUG == 1)
          {
            cat("Transposed_Isoform_Expression_sel.m.dt_0\n")
            cat(str(Transposed_Isoform_Expression_sel.m.dt))
            cat("\n")
          }
          
          
          Summary_table_GENE_EXP<-as.data.frame(Transposed_Isoform_Expression_sel.m.dt[, .(sum_GENE_EXP=sum(TPM)),
                                                                                       by=key(Transposed_Isoform_Expression_sel.m.dt)],stringsAsFactors=F)
          
          if(DEBUG == 1)
          {
            
            cat("Summary_table_GENE_EXP_0\n")
            cat(str(Summary_table_GENE_EXP))
            cat("\n")
            
            # quit(status = 1)
          }
          
          
          Transposed_Isoform_Expression_sel.m<-merge(Transposed_Isoform_Expression_sel.m,
                                                     Summary_table_GENE_EXP,
                                                     by="sample_id")
          
          
          if(DEBUG == 1)
          {
            cat("Transposed_Isoform_Expression_sel.m_1\n")
            cat(str(Transposed_Isoform_Expression_sel.m))
            cat("\n")
            
          }
          
          Transposed_Isoform_Expression_sel.m.dt<-data.table(Transposed_Isoform_Expression_sel.m,
                                                             key=c("sample_id","transcript_id"))
          
          Ratio_df<-as.data.frame(Transposed_Isoform_Expression_sel.m.dt[, .(log2TPM=log2TPM,
                                                                             TPM=TPM,
                                                                             sum_GENE_EXP=sum_GENE_EXP,
                                                                             Ratio=TPM/sum_GENE_EXP),by=key(Transposed_Isoform_Expression_sel.m.dt)],stringsAsFactors=F)
          
          if(DEBUG == 1)
          {
            
            cat("Ratio_df_0\n")
            cat(str(Ratio_df))
            cat("\n")
            #
            # cat("distrib_ratios\n")
            # cat(sprintf(as.character(names(summary(Ratio_df$Ratio)))))
            # cat("\n")
            # cat(sprintf(as.character(summary(Ratio_df$Ratio))))
            # cat("\n")
            #
            #
            # quit(status = 1)
          }
          
          
          
          
          
          #### Merge with covariates matrix & keep HET ----
          
          
          
          Ratio_df<-merge(Ratio_df,
                          INTERVAL_covariates_and_PEER_factors_sel,
                          by=c("sample_id"),
                          all=T)
          
          
          
          
          if(DEBUG == 1)
          {
            
            cat("Ratio_df_3\n")
            cat(str(Ratio_df))
            cat("\n")
            
            # quit(status = 1)
            
          }
          
          Ratio_df_HET<-droplevels(Ratio_df[which(Ratio_df$Genotype != "HOM"),])
          
          Ratio_df_HET$transcript_id<-factor(Ratio_df_HET$transcript_id,
                                             levels=order_transcript_id,
                                             ordered=T)
          
          Ratio_df_HET$Cell_Type<-'Whole blood'
          
          
          Ratio_df_HET$Genotype<-factor(Ratio_df_HET$Genotype,
                                        levels=newG,
                                        ordered=T)
          
          
          
          if(DEBUG == 1)
          {
            
            cat("Ratio_df_HET_0\n")
            cat(str(Ratio_df_HET))
            cat("\n")
            
            
            # quit(status = 1)
            
          }
          
          ##### Violin visual TPM ----
          
          A<-round(summary(Ratio_df_HET$TPM[!is.na(Ratio_df_HET$TPM)]),2)
          
          if(DEBUG == 1)
          {
            cat("summary_TPM\n")
            cat(sprintf(as.character(names(A))))
            cat("\n")
            cat(sprintf(as.character(A)))
            cat("\n")
          }
          
          step<-abs(A[6]-A[1])/4
          
          if(step == 0)
          {
            
            step<-1
          }
          
          breaks.Rank<-sort(unique(c(A[6],seq(from= A[1], to=A[6],by=step))))
          labels.Rank<-as.character(round(breaks.Rank,1))
          
          
          
          if(DEBUG == 1)
          {
            cat("labels.Rank:\t")
            cat(sprintf(as.character(labels.Rank)))
            cat("\n")
            
            
            # quit(status = 1)
          }
          
          
          
          
          Blood_violin<-ggplot(data=Ratio_df_HET,
                               aes(y=Genotype, x=TPM, fill=Genotype)) +
            geom_violin()+
            stat_summary(fun = median, fun.min = median, fun.max = median,
                         geom = "crossbar", width = 0.5)+
            scale_x_continuous(name="TPM",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
            scale_fill_manual(values=c("#E7B800","#00AFBB","#FC4E07"),drop=F)+
            scale_y_discrete(name=NULL, drop=F)
          
          
          
          Blood_violin<-Blood_violin+
            facet_grid(transcript_id ~ Cell_Type, scales='free_x', space='free_x', switch="y", drop=F)+
            theme_cowplot(font_size = 4)+
            theme( strip.background = element_blank(),
                   strip.placement = "outside",
                   strip.text = element_text(size=6),
                   panel.spacing = unit(0.2, "lines"),
                   panel.background=element_rect(fill="white"),
                   panel.border=element_rect(colour="white",size=0,5),
                   panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank())+
            theme_classic()+
            theme(axis.title.y=element_blank(),
                  axis.title.x=element_text(size=8, color="black", family="sans"),
                  axis.text.y=element_blank(),
                  axis.text.x=element_text(angle=0,size=6, color="black", family="sans"),
                  axis.line.x = element_line(size = 0.2),
                  axis.ticks.x = element_line(size = 0.2),
                  axis.ticks.y = element_blank(),
                  axis.line.y = element_line(size = 0.2))+
            theme(legend.title = element_text(size=6),
                  legend.text = element_text(size=6),
                  legend.key.size = unit(0.25, 'cm'), #change legend key size
                  legend.key.height = unit(0.25, 'cm'), #change legend key height
                  legend.key.width = unit(0.25, 'cm'), #change legend key width
                  legend.position="bottom")+
            ggeasy::easy_center_title()
          
          
          setwd(path_graphs)
          
          svgname<-paste(paste('Whole_Blood_violin', sep='_'),".svg",sep='')
          makesvg = TRUE
          
          if (makesvg == TRUE)
          {
            ggsave(svgname, plot= Blood_violin,
                   device="svg",
                   height=4, width=2)
          }
          
          array_transcript_id<-unique(Ratio_df_HET$transcript_id)
          
          
          if(DEBUG == 1)
          {
            cat("array_transcript_id_0\n")
            cat(str(array_transcript_id))
            cat("\n")
          }
          
          for(iteration_array_transcript_id in 1:length(array_transcript_id)){
            
            array_transcript_id_sel<-array_transcript_id[iteration_array_transcript_id]
            
            if(DEBUG == 1)
            {
              cat(sprintf(as.character(iteration_array_transcript_id)))
              cat("\t")
              cat(sprintf(as.character(array_transcript_id_sel)))
              cat("\n")
            }
            
            Ratio_df_HET_transcript_id_sel<-droplevels(Ratio_df_HET[which(Ratio_df_HET$transcript_id == array_transcript_id_sel),])
            
            if(DEBUG == 1)
            {
              cat("Ratio_df_HET_transcript_id_sel_0\n")
              cat(str(Ratio_df_HET_transcript_id_sel))
              cat("\n")
            }
            
            
            A<-round(summary(Ratio_df_HET_transcript_id_sel$TPM[!is.na(Ratio_df_HET_transcript_id_sel$TPM)]),2)
            
            
            # cat("summary_GENE_EXP\n")
            # cat(sprintf(as.character(names(A))))
            # cat("\n")
            # cat(sprintf(as.character(A)))
            # cat("\n")
            
            step<-abs(A[6]-A[1])/4
            
            if(step == 0)
            {
              
              step<-1
            }
            
            breaks.Rank<-sort(unique(c(A[6],seq(from= A[1], to=A[6],by=step))))
            labels.Rank<-as.character(round(breaks.Rank,0))
            
            
            if(DEBUG == 1)
            {
              cat("labels.Rank:\t")
              cat(sprintf(as.character(labels.Rank)))
              cat("\n")
              
              
              
              # quit(status = 1)
            }
            
            
            
            
            
            
            
            WB_per_transcript_violin<-ggplot(data=Ratio_df_HET_transcript_id_sel,
                                             aes(y=TPM, x=Genotype, fill=Genotype)) +
              geom_violin()+
              stat_summary(fun = median, fun.min = median, fun.max = median,
                           geom = "crossbar", width = 0.5)+
              scale_y_continuous(name="TPM",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
              scale_fill_manual(values=c("#E7B800","#00AFBB","#FC4E07"),drop=F)+
              scale_x_discrete(name=NULL, drop=F)
            
            
            
            WB_per_transcript_violin<-WB_per_transcript_violin+
              facet_grid(transcript_id ~ Cell_Type, scales='free_x', space='free_x', switch="y", drop=F)+
              theme_cowplot(font_size = 4)+
              theme( strip.background = element_blank(),
                     strip.placement = "outside",
                     strip.text = element_text(size=6),
                     panel.spacing = unit(0.2, "lines"),
                     panel.background=element_rect(fill="white"),
                     panel.border=element_rect(colour="white",size=0,5),
                     panel.grid.major = element_blank(),
                     panel.grid.minor = element_blank())+
              theme_classic()+
              theme(axis.title.y=element_text(size=8, color="black", family="sans"),
                    axis.title.x=element_blank(),
                    axis.text.y=element_text(angle=0,size=6, color="black", family="sans"),
                    axis.text.x=element_blank(),
                    axis.line.x = element_line(size = 0.2),
                    axis.ticks.x = element_blank(),
                    axis.ticks.y = element_line(size = 0.2),
                    axis.line.y = element_line(size = 0.2))+
              theme(legend.title = element_text(size=6),
                    legend.text = element_text(size=6),
                    legend.key.size = unit(0.25, 'cm'), #change legend key size
                    legend.key.height = unit(0.25, 'cm'), #change legend key height
                    legend.key.width = unit(0.25, 'cm'), #change legend key width
                    legend.position="hidden")+
              ggeasy::easy_center_title()
            
            
            setwd(path_graphs)
            
            svgname<-paste(paste('WB_per_transcript_violin',array_transcript_id_sel, sep='_'),".svg",sep='')
            makesvg = TRUE
            
            if (makesvg == TRUE)
            {
              ggsave(svgname, plot= WB_per_transcript_violin,
                     device="svg",
                     height=2, width=2)
            }
            
          }#iteration_array_transcript_id in 1:length(array_transcript_id
          
          
          
          
        }#file.exists(filename)
      }else{
        
        Blood_violin<-NA
        
        }#file.exists(INTERVAL_files)
      
      graph_DEF<-plot_grid(transcript_plot,Blood_violin,BLUEPRINT_violin,
                           nrow = 1,
                           ncol = 3,
                           rel_widths=c(1,0.5,1))
      
      
      
      
      setwd(path_graphs)
      
      svgname<-paste(paste("GRAPH_DEF",'transcripts', sep='_'),".svg",sep='')
      makesvg = TRUE
      
      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= graph_DEF,
               device="svg",
               height=3, 
               width=6)
        
      }
     
    }#k in 1:length(ENSG_array)
  }#i in 1:length(VARS)
  
  
  
  # #######################################################
  # quit(status = 1)
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
    make_option(c("--Table_S7"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--path_INTERVAL"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--path_BP"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Transposed_Isoform_Expression"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Table_S6"), type="numeric", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ensembl_gtf"), type="numeric", default=NULL, 
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
  
  data_wrangling(opt)
  
  
}


###########################################################################

system.time( main() )
  