
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

intersect_with_PCHiC = function(option_list)
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
  
 
  
  #### READ GWAS block sizes ----
  
  GWAS_blocks<-as.data.frame(fread(file=opt$GWAS_blocks, sep=',', header=T), stringsAsFactors=F)
  
  
  
  cat("GWAS_blocks_0\n")
  cat(str(GWAS_blocks))
  cat("\n")
  
  #### Read annotation files----
  
  PCHiC_original_file<-as.data.frame(fread(file=opt$PCHiC_original_file, sep="\t", header=T), stringsAsFactors=F)
  
  cat("PCHiC_original_file_0\n")
  cat(str(PCHiC_original_file))
  cat("\n")
 
  
  #### READ variants_file ----
  
  variants_file<-readRDS(file=opt$variants_file)
  
  cat("variants_file_0\n")
  cat(str(variants_file))
  cat("\n")
  cat(str(unique(variants_file$VAR)))
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
  
  
  # AS_ID_array<-'rbc__248'
  
  DEBUG<-0
  
  
  Result_FINAL<-data.frame()
  Result_FINAL_PCHiC<-data.frame()
  
  
  for(i in 1:length(AS_ID_array))
  {
    AS_ID_sel<-AS_ID_array[i]
    
    cat("------------------------------------------------------------------------>\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(AS_ID_sel)))
    cat("\n")
    
    # if(AS_ID_sel == 'chr1_158613314_G_A')
    # {
    # 
    #   DEBUG<-1
    # 
    # }else{
    # 
    #   DEBUG<-0
    # }
    
    
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
    
    
    phenotype_sel<-unique(variants_file_sel$phenotype)
    
    if(DEBUG == 1)
    {
      cat("phenotype_sel_0\n")
      cat(str(phenotype_sel))
      cat("\n")
    }
    
    
    
    block_sel<-unique(variants_file_sel$block_no)
    
    if(DEBUG == 1)
    {
      cat("block_sel_0\n")
      cat(str(block_sel))
      cat("\n")
    }
    
    
    GWAS_blocks_sel<-GWAS_blocks[which(GWAS_blocks$pheno == phenotype_sel &
                                         GWAS_blocks$block == block_sel &
                                         GWAS_blocks$chr == chr_sel),]
    
    GWAS_blocks_sel$AS_ID<-AS_ID_sel
    
    if(DEBUG == 1)
    {
      cat("GWAS_blocks_sel_0\n")
      cat(str(GWAS_blocks_sel))
      cat("\n")
    }
    
    if(dim(GWAS_blocks_sel)[1] >0)
    {
      
      gr_BLOCK <- GRanges(
        seqnames = as.character(gsub("^chr","",as.character(GWAS_blocks_sel$chr))),
        strand="*",
        ranges=IRanges(
          start=GWAS_blocks_sel$beg,
          end=GWAS_blocks_sel$end,
          names=GWAS_blocks_sel$AS_ID))
      
      if(DEBUG == 1)
      {
        cat("gr_BLOCK_0\n")
        cat(str(gr_BLOCK))
        cat("\n")
      }
      
      
      PCHiC_original_file_chr_sel<-PCHiC_original_file[which(PCHiC_original_file$baitChr == chr_sel &
                                                               PCHiC_original_file$oeChr == chr_sel),]
      
      if(DEBUG == 1)
      {
        cat("PCHiC_original_file_chr_sel_0\n")
        cat(str(PCHiC_original_file_chr_sel))
        cat("\n")
      }
      
    
      
      gr_oe <- GRanges(
        seqnames = as.character(PCHiC_original_file_chr_sel$oeChr),
        name20=as.character(PCHiC_original_file_chr_sel$baitID),
        ranges=IRanges(
          start=PCHiC_original_file_chr_sel$oeStart,
          end=PCHiC_original_file_chr_sel$oeEnd,
          name=PCHiC_original_file_chr_sel$oeID))
      
      
      gr_bait <- GRanges(
        seqnames = as.character(PCHiC_original_file_chr_sel$baitChr),
        name2=as.character(PCHiC_original_file_chr_sel$Mon),
        name3=as.character(PCHiC_original_file_chr_sel$Mac0),
        name4=as.character(PCHiC_original_file_chr_sel$Mac1),
        name5=as.character(PCHiC_original_file_chr_sel$Mac2),
        name6=as.character(PCHiC_original_file_chr_sel$Neu),
        name7=as.character(PCHiC_original_file_chr_sel$MK),
        name8=as.character(PCHiC_original_file_chr_sel$EP),
        name9=as.character(PCHiC_original_file_chr_sel$Ery),
        name10=as.character(PCHiC_original_file_chr_sel$FoeT),
        name11=as.character(PCHiC_original_file_chr_sel$nCD4),
        name12=as.character(PCHiC_original_file_chr_sel$tCD4),
        name13=as.character(PCHiC_original_file_chr_sel$aCD4),
        name14=as.character(PCHiC_original_file_chr_sel$naCD4),
        name15=as.character(PCHiC_original_file_chr_sel$nCD8),
        name16=as.character(PCHiC_original_file_chr_sel$tCD8),
        name17=as.character(PCHiC_original_file_chr_sel$nB),
        name18=as.character(PCHiC_original_file_chr_sel$tB),
        name19=as.character(PCHiC_original_file_chr_sel$baitName),
        name20=as.character(PCHiC_original_file_chr_sel$oeID),
        ranges=IRanges(
          start=PCHiC_original_file_chr_sel$baitStart,
          end=PCHiC_original_file_chr_sel$baitEnd,
          name=PCHiC_original_file_chr_sel$baitID))
      
      #### find overlap with gr_oe ----
      
      m <- findOverlaps(gr_BLOCK,gr_oe)
      
      if(DEBUG == 1)
      {
        cat("m\n")
        cat(str(m))
        cat("\n")
      }
      
      subjectHits_oe<-subjectHits(m)
      
      if(DEBUG == 1)
      {
        cat("subjectHits_oe\n")
        cat(str(subjectHits_oe))
        cat("\n")
      }
      
      queryHits_BLOCK<-queryHits(m)
      
      if(DEBUG == 1)
      {
        cat("queryHits_BLOCK\n")
        cat(str(queryHits_BLOCK))
        cat("\n")
      }
      
      BLOCK_df <- data.frame(chr=as.character(seqnames(gr_BLOCK)),
                             start=as.integer(start(gr_BLOCK)),
                             end=as.integer(end(gr_BLOCK)),
                             AS_ID=names(gr_BLOCK), stringsAsFactors = F)
     
      
      BLOCK_df_hits<-unique(BLOCK_df[queryHits_BLOCK,])
      
      if(DEBUG == 1)
      {
        cat("BLOCK_df_hits_0\n")
        cat(str(BLOCK_df_hits))
        cat("\n")
      }
      
      oe_df <- data.frame(oeChr=as.character(seqnames(gr_oe)),
                            oeStart=as.integer(start(gr_oe)),
                            oeEnd=as.integer(end(gr_oe)),
                            oeID=names(gr_oe),
                            baitID=as.character(gr_oe$name20),
                            stringsAsFactors = F)
      
      
    
      
      oe_df_hits<-unique(oe_df[subjectHits_oe,])
      
      if(dim(oe_df_hits)[1] >0)
      {
        oe_df_hits$AS_ID<-AS_ID_sel
        
        if(DEBUG == 1)
        {
          cat("oe_df_hits_0\n")
          cat(str(oe_df_hits))
          cat("\n")
          cat(str(unique(oe_df_hits$oeID)))
          cat("\n")
          cat("Min POS\n")
          cat(sprintf(as.character(min(oe_df_hits$oeStart))))
          cat("\n")
          cat("Max POS\n")
          cat(sprintf(as.character(max(oe_df_hits$oeEnd))))
          cat("\n")
        }
        
        #### find overlap with gr_bait ----
        
        m <- findOverlaps(gr_BLOCK,gr_bait)
        
        if(DEBUG == 1)
        {
          cat("m\n")
          cat(str(m))
          cat("\n")
        }
        
        subjectHits_bait<-subjectHits(m)
        
        if(DEBUG == 1)
        {
          cat("subjectHits_bait\n")
          cat(str(subjectHits_bait))
          cat("\n")
        }
        
        queryHits_BLOCK<-queryHits(m)
        
        if(DEBUG == 1)
        {
          cat("queryHits_BLOCK\n")
          cat(str(queryHits_BLOCK))
          cat("\n")
        }
        
        BLOCK_df <- data.frame(chr=as.character(seqnames(gr_BLOCK)),
                               start=as.integer(start(gr_BLOCK)),
                               end=as.integer(end(gr_BLOCK)),
                               AS_ID=names(gr_BLOCK), stringsAsFactors = F)
        
        # if(DEBUG == 1)
        # {
        #   cat("BLOCK_df_0\n")
        #   cat(str(BLOCK_df))
        #   cat("\n")
        # }
        
        BLOCK_df_hits<-unique(BLOCK_df[queryHits_BLOCK,])
        
        if(DEBUG == 1)
        {
          cat("BLOCK_df_hits_0\n")
          cat(str(BLOCK_df_hits))
          cat("\n")
        }
        
        bait_df <- data.frame(baitChr=as.character(seqnames(gr_bait)),
                              baitStart=as.integer(start(gr_bait)),
                              baitEnd=as.integer(end(gr_bait)),
                              baitID=names(gr_bait),
                              baitName=as.character(gr_bait$name19),
                              oeID=as.character(gr_bait$name20),
                              Mon=as.numeric(gr_bait$name2),
                              Mac0=as.numeric(gr_bait$name3),
                              Mac1=as.numeric(gr_bait$name4),
                              Mac2=as.numeric(gr_bait$name5),
                              Neu=as.numeric(gr_bait$name6),
                              MK=as.numeric(gr_bait$name7),
                              EP=as.numeric(gr_bait$name8),
                              Ery=as.numeric(gr_bait$name9),
                              FoeT=as.numeric(gr_bait$name10),
                              nCD4=as.numeric(gr_bait$name11),
                              tCD4=as.numeric(gr_bait$name12),
                              aCD4=as.numeric(gr_bait$name13),
                              naCD4=as.numeric(gr_bait$name14),
                              nCD8=as.numeric(gr_bait$name15),
                              tCD8=as.numeric(gr_bait$name16),
                              nB=as.numeric(gr_bait$name17),
                              tB=as.numeric(gr_bait$name18),
                              stringsAsFactors = F)
        
        
        if(DEBUG == 1)
        {
          cat("bait_df_0\n")
          cat(str(bait_df))
          cat("\n")
          cat(str(unique(bait_df$baitID)))
          cat("\n")
        }
        
        
        
        bait_df_hits<-bait_df[subjectHits_bait,]
        
        if(dim(bait_df_hits)[1] >0)
        {
          bait_df_hits$AS_ID<-AS_ID_sel
          
          if(DEBUG == 1)
          {
            cat("bait_df_hits_0\n")
            cat(str(bait_df_hits))
            cat("\n")
            cat(str(unique(bait_df_hits$baitID)))
            cat("\n")
          }
          
          
          bait_df_hits<-merge(oe_df_hits,
                              bait_df_hits,
                              by=c("oeID","baitID","AS_ID"))
          
          if(DEBUG == 1)
          {
            cat("bait_df_hits_1\n")
            cat(str(bait_df_hits))
            cat("\n")
            cat(str(unique(bait_df_hits$baitID)))
            cat("\n")
            cat("bait_min\n")
            cat(sprintf(as.character(min(bait_df_hits$baitStart))))
            cat("\n")
            cat("bait_max\n")
            cat(sprintf(as.character(max(bait_df_hits$baitEnd))))
            cat("\n")
            cat("oe_min\n")
            cat(sprintf(as.character(min(bait_df_hits$oeStart))))
            cat("\n")
            cat("oe_max\n")
            cat(sprintf(as.character(max(bait_df_hits$oeEnd))))
            cat("\n")
          }
          
          Result_FINAL_PCHiC<-rbind(Result_FINAL_PCHiC,bait_df_hits)
          
          
          TOTAL_pos<-sort(unique(c(BLOCK_df$start,BLOCK_df$end,bait_df_hits$baitStart,bait_df_hits$baitEnd,bait_df_hits$oeStart,bait_df_hits$oeEnd)))
          
          MIN_POS<-min(TOTAL_pos)
          MAX_POS<-max(TOTAL_pos)
          
          if(DEBUG == 1)
          {
            cat("TOTAL_pos_0\n")
            cat(str(TOTAL_pos))
            cat("\n")
            cat(sprintf(as.character(MIN_POS)))
            cat("\n")
            cat(sprintf(as.character(MAX_POS)))
            cat("\n")
          }
          
        
         
          
          
        }else{
          
          
          
          TOTAL_pos<-sort(unique(c(BLOCK_df$start,BLOCK_df$end,oe_df_hits$oeStart,oe_df_hits$oeEnd)))
        
          
          MIN_POS<-min(TOTAL_pos)
          MAX_POS<-max(TOTAL_pos)
          
          if(DEBUG == 1)
          {
            cat("TOTAL_pos_0\n")
            cat(str(TOTAL_pos))
            cat("\n")
            cat(sprintf(as.character(MIN_POS)))
            cat("\n")
            cat(sprintf(as.character(MAX_POS)))
            cat("\n")
          }
          
          
        }#dim(bait_df_hits)[1] >0
      
        
        
        
        
      }else{
        
        TOTAL_pos<-sort(unique(c(BLOCK_df$start,BLOCK_df$end)))
        
        
        MIN_POS<-min(TOTAL_pos)
        MAX_POS<-max(TOTAL_pos)
        
        }#dim(oe_df_hits)[1] >0
    
      temp<-as.data.frame(cbind(AS_ID_sel,paste('chr',chr_sel, sep=''),MIN_POS,MAX_POS), stringsAsFactors=F)
      colnames(temp)<-c('AS_ID','chr','START','END')
      
      temp$START<-as.integer(temp$START)
      temp$END<-as.integer(temp$END)
      
      if(DEBUG == 1)
      {
        cat("temp_0\n")
        cat(str(temp))
        cat("\n")
        
      }
      
    Result_FINAL<-rbind(temp,Result_FINAL)
      
     
    }# dim(GWAS_blocks_sel)[1] >0
  }#i in 1:length(AS_ID_array)
  
  
  #### SAVE ----
  
  if(dim(Result_FINAL_PCHiC)[1] >0)
  {
    
   
    cat("Result_FINAL_PCHiC_0\n")
    cat(str(Result_FINAL_PCHiC))
    cat("\n")
    cat(str(unique(Result_FINAL_PCHiC$AS_ID)))
    cat("\n")
    
    setwd(path_files)
    
    saveRDS(Result_FINAL_PCHiC, file=paste('PCHiC_AS_pre_ENSEMBL',".rds",sep=''))
    
    
    
  }#dim(Result_FINAL_PCHiC)[1] >0
  
  
  if(dim(Result_FINAL)[1] >0)
  {
    
    Result_FINAL$chr<-factor(Result_FINAL$chr,
                             levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                                      "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                      "chr22","chr23","chrX","chrY"), ordered=T)
    
    Result_FINAL<-droplevels(Result_FINAL)
    
    Result_FINAL<-Result_FINAL[order(Result_FINAL$chr),]
    
    cat("Result_FINAL_0\n")
    cat(str(Result_FINAL))
    cat("\n")
    cat(str(unique(Result_FINAL$AS_ID)))
    cat("\n")
    
    setwd(path_files)
    
    saveRDS(Result_FINAL, file=paste('Boundaries_AS_PRE_genes',".rds",sep=''))
    
    
    
  }#dim(Result_FINAL)[1] >0
}

ENSEMBL_gene_id_to_PCHiC = function(option_list)
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
  
  
  #### READ ensembl_gtf ----
  
  ensembl_gtf = readGFF(opt$ensembl_gtf)
  
  
  cat("ensembl_gtf_0\n")
  cat(str(ensembl_gtf))
  cat("\n")
  
  
 
  #### READ genes_without_symbol_CORRECTED ----
  
  
  genes_without_symbol_CORRECTED<-as.data.frame(fread(file=opt$genes_without_symbol_CORRECTED, sep="\t", header =T), stringAsFactors=F)
  
  colnames(genes_without_symbol_CORRECTED)[which(colnames(genes_without_symbol_CORRECTED) == 'Alias')]<-'baitName'
  colnames(genes_without_symbol_CORRECTED)[which(colnames(genes_without_symbol_CORRECTED) == 'SYMBOL')]<-'baitSYMBOL'
  
  
  cat("genes_without_symbol_CORRECTED_0\n")
  cat(str(genes_without_symbol_CORRECTED))
  cat("\n")
  cat(str(unique(genes_without_symbol_CORRECTED$baitName)))
  cat("\n")
  cat(str(unique(genes_without_symbol_CORRECTED$baitSYMBOL)))
  cat("\n")
  
  genes_without_symbol_CORRECTED_NO_NA<-genes_without_symbol_CORRECTED[!is.na(genes_without_symbol_CORRECTED$baitSYMBOL),]
  
  cat("genes_without_symbol_CORRECTED_NO_NA_0\n")
  cat(str(genes_without_symbol_CORRECTED_NO_NA))
  cat("\n")
  cat(str(unique(genes_without_symbol_CORRECTED_NO_NA$baitName)))
  cat("\n")
  cat(str(unique(genes_without_symbol_CORRECTED_NO_NA$baitSYMBOL)))
  cat("\n")
  
  
  #### READ ensembl_gtf_sel_CORRECTED ----
  
  
  ensembl_gtf_sel_CORRECTED<-as.data.frame(fread(file=opt$ensembl_gtf_sel_CORRECTED, sep="\t", header =T), stringAsFactors=F)
  
  
  cat("ensembl_gtf_sel_CORRECTED_0\n")
  cat(str(ensembl_gtf_sel_CORRECTED))
  cat("\n")
  cat(str(unique(ensembl_gtf_sel_CORRECTED$gene_name)))
  cat("\n")
  cat(str(unique(ensembl_gtf_sel_CORRECTED$ensembl_gene_id)))
  cat("\n")
  
  
  
  
  ###### LOOP AS_IDs -----
  
  path_files<-paste(out,'Build_files','/',sep='')
  
  if (file.exists(path_files)){
    
    
  }else{
    
    dir.create(file.path(path_files))
    
  }#path_files
  
  setwd(path_files)
  

  Result_FINAL_PCHiC<-readRDS(file=paste('PCHiC_AS_pre_ENSEMBL',".rds",sep=''))
  
  cat("Result_FINAL_PCHiC_0\n")
  cat(str(Result_FINAL_PCHiC))
  cat("\n")
  cat(str(unique(Result_FINAL_PCHiC$AS_ID)))
  cat("\n")
  
  ##### csplit long ----
  
  DEBUG <- 1
  
  Result_FINAL_PCHiC<-unique(as.data.frame(cSplit(Result_FINAL_PCHiC,sep = ';', direction = "long",
                                              splitCols = "baitName"),stringsAsFactors=F))
  if(DEBUG == 1)
  {
    cat("Result_FINAL_PCHiC_1\n")
    cat(str(Result_FINAL_PCHiC))
    cat("\n")
   
  }
  
  ##### FLAG NA ----
  
  
 
  FLAG_NA_bait<-dim(Result_FINAL_PCHiC[is.na(Result_FINAL_PCHiC$baitName),])[1]
  
  if(DEBUG == 1)
  {
    cat("FLAG_NA_bait:\n")
    cat(str(FLAG_NA_bait))
    cat("\n")
    
  }
  
  
  if(FLAG_NA_bait > 0)
  {
    
    Result_FINAL_PCHiC<-Result_FINAL_PCHiC[!is.na(Result_FINAL_PCHiC$baitName),]
    
    if(DEBUG == 1)
    {
      cat("Result_FINAL_PCHiC_2_Hello_world\n")
      cat(str(Result_FINAL_PCHiC))
      cat("\n")
     
    }
  }else{
    
    if(DEBUG == 1)
    {
      cat("Result_FINAL_PCHiC_2\n")
      cat(str(Result_FINAL_PCHiC))
      cat("\n")
      
    }
    
  }#FLAG_NA_bait > 0
  
  
  ##### merges ----
  
  Result_FINAL_PCHiC<-merge(Result_FINAL_PCHiC,
                            ensembl_gtf_sel_CORRECTED,
                            by='baitName',
                            all.x=T)
  
  
  if(DEBUG == 1)
  {
    cat("Result_FINAL_PCHiC_POST_MERGE_3\n")
    cat(str(Result_FINAL_PCHiC))
    cat("\n")
    
  }
  
  Result_FINAL_PCHiC<-merge(Result_FINAL_PCHiC,
                            genes_without_symbol_CORRECTED_NO_NA,
                            by='baitName',
                            all.x=T)
  
  
  Result_FINAL_PCHiC$baitSYMBOL[is.na(Result_FINAL_PCHiC$baitSYMBOL)]<-Result_FINAL_PCHiC$baitName[is.na(Result_FINAL_PCHiC$baitSYMBOL)]
  
  if(DEBUG == 1)
  {
    cat("Result_FINAL_PCHiC_POST_MERGE_4\n")
    cat(str(Result_FINAL_PCHiC))
    cat("\n")
    
  }
  
  ##### MapIDs ----
  
  Result_FINAL_PCHiC$ensembl_gene_id[is.na(Result_FINAL_PCHiC$ensembl_gene_id)] <- mapIds(org.Hs.eg.db, keys=Result_FINAL_PCHiC$baitSYMBOL[is.na(Result_FINAL_PCHiC$ensembl_gene_id)], keytype="SYMBOL",
                                                                                          column="ENSEMBL", multiVals="first")
  
  if(DEBUG == 1)
  {
    cat("Result_FINAL_PCHiC_4\n")
    cat(str(Result_FINAL_PCHiC))
    cat("\n")
    cat(str(unique(Result_FINAL_PCHiC$AS_ID)))
    cat("\n")
  }
  
  Result_FINAL_PCHiC_NA<-Result_FINAL_PCHiC[is.na(Result_FINAL_PCHiC$ensembl_gene_id),]
  
  if(dim(Result_FINAL_PCHiC_NA)[1] >0)
  {
    
    cat("----------------------------------------------------------------------------------->Result_FINAL_PCHiC_NA_bait_name\n")
    cat(sprintf(as.character(unique(Result_FINAL_PCHiC_NA$baitName))))
    cat("\n")

  }
  
  Result_FINAL_PCHiC_NO_NA<-Result_FINAL_PCHiC[!is.na(Result_FINAL_PCHiC$ensembl_gene_id),]
  
  if(DEBUG == 1)
  {
    cat("Result_FINAL_PCHiC_NO_NA_0\n")
    cat(str(Result_FINAL_PCHiC_NO_NA))
    cat("\n")
    cat(str(unique(Result_FINAL_PCHiC_NO_NA$AS_ID)))
    cat("\n")
  }
  
  
  ##### Redo the collapse ----
  
  
  keys_to_keep<-colnames(Result_FINAL_PCHiC_NO_NA)[-which(colnames(Result_FINAL_PCHiC_NO_NA)%in%c("baitName","baitSYMBOL","ensembl_gene_id"))]
  
  if(DEBUG == 1)
  {
    cat("keys_to_keep_4\n")
    cat(str(keys_to_keep))
    cat("\n")
  }
  
  Result_FINAL_PCHiC_NO_NA.dt<-data.table(Result_FINAL_PCHiC_NO_NA, key=keys_to_keep)
  
  if(DEBUG == 1)
  {
    cat("Result_FINAL_PCHiC_NO_NA.dt_4\n")
    cat(str(Result_FINAL_PCHiC_NO_NA.dt))
    cat("\n")
  }
  
  Result_FINAL_PCHiC_NO_NA_collapse<-as.data.frame(Result_FINAL_PCHiC_NO_NA.dt[,.(baitName=paste(baitName, collapse=";"),
                                                                                 baitSYMBOL=paste(baitSYMBOL, collapse=";"),
                                                                                 ensembl_gene_id=paste(ensembl_gene_id, collapse=";")), by=key(Result_FINAL_PCHiC_NO_NA.dt)], stringsAsFactors=F)
  
  if(DEBUG == 1)
  {
    cat("Result_FINAL_PCHiC_NO_NA_collapse_0\n")
    cat(str(Result_FINAL_PCHiC_NO_NA_collapse))
    cat("\n")
    cat(str(unique(Result_FINAL_PCHiC_NO_NA_collapse$AS_ID)))
    cat("\n")
  }
  
  setwd(path_files)
  
  saveRDS(Result_FINAL_PCHiC_NO_NA_collapse, file=paste('PCHiC_AS_POST_ENSEMBL',".rds",sep=''))
  
  
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
    make_option(c("--GWAS_blocks"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--PCHiC_original_file"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ensembl_gtf"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--genes_without_symbol_CORRECTED"), type="character", default=NULL, 
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
  
  intersect_with_PCHiC(opt)
  ENSEMBL_gene_id_to_PCHiC(opt)

  
}


###########################################################################

system.time( main() )