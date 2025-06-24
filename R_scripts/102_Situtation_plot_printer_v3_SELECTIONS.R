
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
suppressMessages(library("RColorBrewer", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("cowplot", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggrepel", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggarchery", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("ggnewscale", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("splitstackshape", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("viridis", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))

opt = NULL

options(warn = 1)

Printer_function = function(option_list)
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
  
  #### Read rescue_haplotypes----
  
  rescue_haplotypes<-as.data.frame(fread(file=opt$rescue_haplotypes, sep="\t", header=T), stringsAsFactors=F)
  
 
  
  cat("rescue_haplotypes_0\n")
  cat(str(rescue_haplotypes))
  cat("\n")
  
  
  #### Read annotation files----
  
  DEBUG<-0
  
  dataframe_to_print_SIT_plots<-as.data.frame(fread(file=opt$dataframe_to_print_SIT_plots, sep="\t", header=T), stringsAsFactors=F)
  
  colnames(dataframe_to_print_SIT_plots)[which(colnames(dataframe_to_print_SIT_plots) == 'AS')]<-'AS_ID'
  
  dataframe_to_print_SIT_plots$beta_span<-gsub("\"","",dataframe_to_print_SIT_plots$beta_span)
  dataframe_to_print_SIT_plots$beta_span<-gsub("\'","",dataframe_to_print_SIT_plots$beta_span)
  
  
  cat("dataframe_to_print_SIT_plots_0\n")
  cat(str(dataframe_to_print_SIT_plots))
  cat("\n")
  cat(str(unique(dataframe_to_print_SIT_plots$AS_ID)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(dataframe_to_print_SIT_plots$RNA_source))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(dataframe_to_print_SIT_plots$RNA_source)))))
  cat("\n")
  
  ##### csplit long ENSG ----
  
  
  dataframe_subset_ENSG<-unique(dataframe_to_print_SIT_plots[,c(which(colnames(dataframe_to_print_SIT_plots) == 'AS_ID'),
                                                                 which(colnames(dataframe_to_print_SIT_plots) == 'Genes'))])
  if(DEBUG == 1)
  {
    cat("dataframe_subset_ENSG_0\n")
    cat(str(dataframe_subset_ENSG))
    cat("\n")
    cat(str(unique(dataframe_subset_ENSG$AS_ID)))
    cat("\n")
  }
  
  dataframe_subset_ENSG_NO_NA<-dataframe_subset_ENSG[!is.na(dataframe_subset_ENSG$Genes),]
  
  cat("dataframe_subset_ENSG_NO_NA_0\n")
  cat(str(dataframe_subset_ENSG_NO_NA))
  cat("\n")
  cat(str(unique(dataframe_subset_ENSG_NO_NA$AS_ID)))
  cat("\n")
  
  dataframe_subset_ENSG_NO_NA<-unique(as.data.frame(cSplit(dataframe_subset_ENSG_NO_NA,sep = ';', direction = "long",
                                                            splitCols = "Genes"),stringsAsFactors=F))
  
  colnames(dataframe_subset_ENSG_NO_NA)[which(colnames(dataframe_subset_ENSG_NO_NA) == 'Genes')]<-'ensembl_gene_id'
  
  
 
  cat("dataframe_subset_ENSG_NO_NA_1\n")
  cat(str(dataframe_subset_ENSG_NO_NA))
  cat("\n")
  cat(str(unique(dataframe_subset_ENSG_NO_NA$AS_ID)))
  cat("\n")
  
  
 
  
  ##### csplit long PCHiC ----
  
 
  dataframe_subset_PCHiC<-unique(dataframe_to_print_SIT_plots[,c(which(colnames(dataframe_to_print_SIT_plots) == 'AS_ID'),
                                                                 which(colnames(dataframe_to_print_SIT_plots) == 'oeID'),
                                                                 which(colnames(dataframe_to_print_SIT_plots) == 'PCHiC_cells'))])
  if(DEBUG == 1)
  {
    cat("dataframe_subset_PCHiC_0\n")
    cat(str(dataframe_subset_PCHiC))
    cat("\n")
    cat(str(unique(dataframe_subset_PCHiC$AS_ID)))
    cat("\n")
  }
  
  dataframe_subset_PCHiC_NO_NA<-dataframe_subset_PCHiC[!is.na(dataframe_subset_PCHiC$oeID),]
  
  cat("dataframe_subset_PCHiC_NO_NA_0\n")
  cat(str(dataframe_subset_PCHiC_NO_NA))
  cat("\n")
  cat(str(unique(dataframe_subset_PCHiC_NO_NA$AS_ID)))
  cat("\n")
  
  
  
 
  dataframe_subset_PCHiC_NO_NA<-unique(as.data.frame(cSplit(dataframe_subset_PCHiC_NO_NA,sep = ';', direction = "long",
                                                  splitCols = "oeID"),stringsAsFactors=F))
 
  if(DEBUG == 1)
  {
    cat("dataframe_subset_PCHiC_NO_NA_1\n")
    cat(str(dataframe_subset_PCHiC_NO_NA))
    cat("\n")
    cat(str(unique(dataframe_subset_PCHiC_NO_NA$AS_ID)))
    cat("\n")
  }
    
  dataframe_subset_PCHiC_NO_NA<-unique(as.data.frame(cSplit(dataframe_subset_PCHiC_NO_NA,sep = ';', direction = "long",
                                                            splitCols = "PCHiC_cells"),stringsAsFactors=F))
  
  colnames(dataframe_subset_PCHiC_NO_NA)[which(colnames(dataframe_subset_PCHiC_NO_NA) == 'PCHiC_cells')]<-'Cell_Type'
  
  
  cat("dataframe_subset_PCHiC_NO_NA_2\n")
  cat(str(dataframe_subset_PCHiC_NO_NA))
  cat("\n")
  cat(str(unique(dataframe_subset_PCHiC_NO_NA$AS_ID)))
  cat("\n")
  
  
  
  ##### csplit long ATAC ----
  
  dataframe_subset_ATAC<-unique(dataframe_to_print_SIT_plots[,c(which(colnames(dataframe_to_print_SIT_plots) == 'AS_ID'),
                                                                 which(colnames(dataframe_to_print_SIT_plots) == 'ATAC_cells'))])
  if(DEBUG == 1)
  {
    cat("dataframe_subset_ATAC_0\n")
    cat(str(dataframe_subset_ATAC))
    cat("\n")
    cat(str(unique(dataframe_subset_ATAC$AS_ID)))
    cat("\n")
  }
  
  dataframe_subset_ATAC_NO_NA<-dataframe_subset_ATAC[!is.na(dataframe_subset_ATAC$ATAC_cells),]
  
  cat("dataframe_subset_ATAC_NO_NA_0\n")
  cat(str(dataframe_subset_ATAC_NO_NA))
  cat("\n")
  cat(str(unique(dataframe_subset_ATAC_NO_NA$AS_ID)))
  cat("\n")
  
  
  
  
  dataframe_subset_ATAC_NO_NA<-unique(as.data.frame(cSplit(dataframe_subset_ATAC_NO_NA,sep = ';', direction = "long",
                                                            splitCols = "ATAC_cells"),stringsAsFactors=F))
  
  colnames(dataframe_subset_ATAC_NO_NA)[which(colnames(dataframe_subset_ATAC_NO_NA) == 'ATAC_cells')]<-'Cell_Type'
  
  
  cat("dataframe_subset_ATAC_NO_NA_1\n")
  cat(str(dataframe_subset_ATAC_NO_NA))
  cat("\n")
  cat(str(unique(dataframe_subset_ATAC_NO_NA$AS_ID)))
  cat("\n")
  

  
  #### READ ATAC_tracks ----
  
  ATAC_tracks<-readRDS(file=opt$ATAC_tracks)
  
  cat("ATAC_tracks_0\n")
  cat(str(ATAC_tracks))
  cat("\n")
  cat(str(unique(ATAC_tracks$AS_ID)))
  cat("\n")
  
  ATAC_tracks_subset<-ATAC_tracks[which(ATAC_tracks$AS_ID%in%dataframe_subset_ATAC_NO_NA$AS_ID),]
  
  cat("ATAC_tracks_subset_0\n")
  cat(str(ATAC_tracks_subset))
  cat("\n")
  cat(str(unique(ATAC_tracks_subset$AS_ID)))
  cat("\n")
  
  ATAC_tracks_subset.m<-melt(ATAC_tracks_subset, 
                          id.vars=c('chr','start','end','AS_ID'), variable.name="Cell_Type", value.name="counts")
  
  if(DEBUG == 1)
  {
    cat("ATAC_tracks_subset.m_0\n")
    cat(str(ATAC_tracks_subset.m))
    cat("\n")
    cat(str(unique(ATAC_tracks_subset.m$AS_ID)))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(ATAC_tracks_subset.m$Cell_Type))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(ATAC_tracks_subset.m$Cell_Type)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(ATAC_tracks_subset.m$counts)))))
    cat("\n")
    cat(sprintf(as.character(summary(ATAC_tracks_subset.m$counts))))
    cat("\n")
  }
  
  ATAC_tracks_subset.m<-merge(dataframe_subset_ATAC_NO_NA,
                           ATAC_tracks_subset.m,
                           by=c("AS_ID","Cell_Type"))
  
  cat("ATAC_tracks_subset.m_1\n")
  cat(str(ATAC_tracks_subset.m))
  cat("\n")
  cat(str(unique(ATAC_tracks_subset.m$AS_ID)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(ATAC_tracks_subset.m$Cell_Type))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(ATAC_tracks_subset.m$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(ATAC_tracks_subset.m$counts)))))
  cat("\n")
  cat(sprintf(as.character(summary(ATAC_tracks_subset.m$counts))))
  cat("\n")
    
 
  
  
  #### READ PCHiC_HITS ----
  
  PCHiC_HITS<-readRDS(file=opt$PCHiC_HITS)
  
  cat("PCHiC_HITS_0\n")
  cat(str(PCHiC_HITS))
  cat("\n")
  cat(str(unique(PCHiC_HITS$AS_ID)))
  cat("\n")
  
  
  PCHiC_HITS_subset<-PCHiC_HITS[which(PCHiC_HITS$AS_ID%in%dataframe_subset_PCHiC_NO_NA$AS_ID),]
  
  cat("PCHiC_HITS_subset_0\n")
  cat(str(PCHiC_HITS_subset))
  cat("\n")
  cat(str(unique(PCHiC_HITS_subset$AS_ID)))
  cat("\n")
  
  indx.scale.CS<-c(which(colnames(PCHiC_HITS_subset) == 'AS_ID'),
                   which(colnames(PCHiC_HITS_subset) == 'oeStart'),which(colnames(PCHiC_HITS_subset) == 'oeEnd'),which(colnames(PCHiC_HITS_subset) == 'oeID'),
              which(colnames(PCHiC_HITS_subset) == 'baitChr'),which(colnames(PCHiC_HITS_subset) == 'baitStart'),which(colnames(PCHiC_HITS_subset) == 'baitEnd'),which(colnames(PCHiC_HITS_subset) == 'baitID'),
              which(colnames(PCHiC_HITS_subset) == 'baitSYMBOL'),
              which(colnames(PCHiC_HITS_subset) == 'Mon'),which(colnames(PCHiC_HITS_subset) == 'Mac0'),which(colnames(PCHiC_HITS_subset) == 'Mac1'),
              which(colnames(PCHiC_HITS_subset) == 'Mac2'),which(colnames(PCHiC_HITS_subset) == 'Neu'),which(colnames(PCHiC_HITS_subset) == 'MK'),
              which(colnames(PCHiC_HITS_subset) == 'EP'),which(colnames(PCHiC_HITS_subset) == 'Ery'),which(colnames(PCHiC_HITS_subset) == 'FoeT'),
              which(colnames(PCHiC_HITS_subset) == 'nCD4'),which(colnames(PCHiC_HITS_subset) == 'tCD4'),which(colnames(PCHiC_HITS_subset) == 'aCD4'),
              which(colnames(PCHiC_HITS_subset) == 'naCD4'),which(colnames(PCHiC_HITS_subset) == 'nCD8'),which(colnames(PCHiC_HITS_subset) == 'tCD8'),
              which(colnames(PCHiC_HITS_subset) == 'nB'),which(colnames(PCHiC_HITS_subset) == 'tB'))
  
  
  
  scale_df<-unique(PCHiC_HITS_subset[,indx.scale.CS])
  

  if(DEBUG == 1)
  {
 
    cat("scale_df_0\n")
    cat(str(scale_df))
    cat("\n")
    cat(str(unique(scale_df$AS_ID)))
    cat("\n")
  }

  
  scale_df.m<-melt(scale_df, id.vars=c("AS_ID","oeStart","oeEnd","oeID","baitChr","baitStart","baitEnd","baitID","baitSYMBOL"),
                 variable.name="Cell_Type", value.name="ChicagoScore")
  
  
  if(DEBUG == 1)
  {
    cat("scale_df.m_0\n")
    cat(str(scale_df.m))
    cat("\n")
    cat(str(unique(scale_df.m$AS_ID)))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(scale_df.m$Cell_Type))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(scale_df.m$Cell_Type)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(scale_df.m$ChicagoScore)))))
    cat("\n")
    cat(sprintf(as.character(summary(scale_df.m$ChicagoScore))))
    cat("\n")
  }
  
  scale_df.m_subset<-merge(dataframe_subset_PCHiC_NO_NA,
                           scale_df.m,
                           by=c("AS_ID","Cell_Type","oeID"))
  
  cat("scale_df.m_subset_0\n")
  cat(str(scale_df.m_subset))
  cat("\n")
  cat(str(unique(scale_df.m_subset$AS_ID)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(scale_df.m_subset$Cell_Type))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(scale_df.m_subset$Cell_Type)))))
  cat("\n")
  cat(sprintf(as.character(names(summary(scale_df.m_subset$ChicagoScore)))))
  cat("\n")
  cat(sprintf(as.character(summary(scale_df.m_subset$ChicagoScore))))
  cat("\n")
  
  breaks_CS<-unique(sort(unique(c(0,5,max(scale_df.m_subset$ChicagoScore[!is.na(scale_df.m_subset$ChicagoScore)])))))
  labels_CS<-as.character(round(breaks_CS,0))
  
 
  cat("breaks_CS\n")
  cat(sprintf(as.character(breaks_CS)))
  cat("\n")
  cat("labels_CS\n")
  cat(sprintf(as.character(labels_CS)))
  cat("\n")
 
  
  
  check<-scale_df.m[-which(scale_df.m$AS_ID%in%scale_df.m_subset$AS_ID),]
  
  if(DEBUG == 1)
  {
    cat("check_0\n")
    cat(str(check))
    cat("\n")
    cat(str(unique(check$AS_ID)))
    cat("\n")
    cat(sprintf(as.character(names(summary(as.factor(check$Cell_Type))))))
    cat("\n")
    cat(sprintf(as.character(summary(as.factor(check$Cell_Type)))))
    cat("\n")
    cat(sprintf(as.character(names(summary(check$ChicagoScore)))))
    cat("\n")
    cat(sprintf(as.character(summary(check$ChicagoScore))))
    cat("\n")
  }
  # setwd(out)
  # 
  # write.table(check,file='test.tsv', sep="\t", quote=F, row.names = F)
  # 
  # 
  # ####################################################################
  # quit(status = 1)
  
  #### vector_colors ----
  
  vector_colors_ALL_but_index<-brewer.pal(9, "Blues")[c(9,7,5,3)]
  vector_colors_index<-brewer.pal(9, "PuRd")[c(5,7)]
  
  vector_colors<-c(vector_colors_ALL_but_index,vector_colors_index)
  
  cat("vector_colors_0\n")
  cat(str(vector_colors))
  cat("\n")
 
 
  #### READ variants_file ----
  
  variants_file<-readRDS(file=opt$variants_file)
  
  cat("variants_file_0\n")
  cat(str(variants_file))
  cat("\n")
  cat(str(unique(variants_file$VAR)))
  cat("\n")
  
 
 
  
  #### READ GWAS block sizes ----
  
  BLOCK_POST<-readRDS(file=opt$BLOCK_POST)
  
 
  
  cat("BLOCK_POST_0\n")
  cat(str(BLOCK_POST))
  cat("\n")
  
 
  
  #### READ gene_tracks ----
  
  gene_tracks<-readRDS(file=opt$gene_tracks)
  
  cat("gene_tracks_0\n")
  cat(str(gene_tracks))
  cat("\n")
  
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
  
  ###### LOOP AS_IDs -----
  
  
  
  path_graphs<-paste(out,'Graphs_selected','/',sep='')
  
  if (file.exists(path_graphs)){
    
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  
  
  VARS<-unique(dataframe_to_print_SIT_plots$VAR)
  
  # VARS<-'chr7_101499930_G_A'
  
  # VARS<-'chr2_219020958_C_T'
  
  cat("VARS_0\n")
  cat(str(VARS))
  cat("\n")
  
  
  # 'chr18_60920854_C_T'
  # VARS<-'chr8_41589736_T_G'
  # VARS<-'chr18_60920854_C_T'
  
  # VARS<-'chr3_184091102_T_G'
  
  # VARS<-'chr3_46354444_C_T'
  
  DEBUG<-0
  
  for(i in 1:length(VARS))
  {
    VAR_sel<-VARS[i]
    
    
   
    dataframe_to_print_SIT_plots_sel<-dataframe_to_print_SIT_plots[which(dataframe_to_print_SIT_plots$VAR == VAR_sel),]
    
    if(DEBUG == 1)
    {
      cat("dataframe_to_print_SIT_plots_sel_0\n")
      cat(str(dataframe_to_print_SIT_plots_sel))
      cat("\n")
    }
   
    
    AS_ID_array<-unique(dataframe_to_print_SIT_plots_sel$AS_ID)
    
    # AS_ID_array<-'mpv__423'
    
    if(DEBUG == 1)
    {
      cat("AS_ID_array_0\n")
      cat(str(AS_ID_array))
      cat("\n")
    }
    
    rs_sel<-unique(dataframe_to_print_SIT_plots_sel$rs)
    
    if(DEBUG == 1)
    {
      cat("rs_sel_0\n")
      cat(str(rs_sel))
      cat("\n")
    }
    
   
   
    path_graphs<-paste(out,'Graphs_selected','/',paste(rs_sel,VAR_sel, sep='__'),'/',sep='')

    if (file.exists(path_graphs)){


    }else{

      dir.create(file.path(path_graphs))

    }#path_graphs
    
   
    
    # AS_ID_array<-'ret__395'
    # AS_ID_array<-'rbc__248'
    # neut__229
    #VAR_sel == 'chr2_219020958_C_T' & AS_ID_sel == 'baso__92
    
    for(k in 1:length(AS_ID_array))
    {
      AS_ID_sel<-AS_ID_array[k]
      
      if(VAR_sel == 'chr1_92925654_G_C' & AS_ID_sel == 'mono__27')
      {

        DEBUG<-1

      }else{

        DEBUG<-0
      }
      
      dataframe_to_print_SIT_plots_sel_AS_sel<-dataframe_to_print_SIT_plots[which(dataframe_to_print_SIT_plots$VAR == VAR_sel &
                                                                                    dataframe_to_print_SIT_plots$AS_ID == AS_ID_sel),]
      
      
      if(DEBUG == 1)
      {
        cat("dataframe_to_print_SIT_plots_sel_AS_sel_0\n")
        cat(str(dataframe_to_print_SIT_plots_sel_AS_sel))
        cat("\n")
      }
      
      source_sel<-unique(dataframe_to_print_SIT_plots_sel_AS_sel$RNA_source)
      
      if(DEBUG == 1)
      {
        cat("source_sel_0\n")
        cat(str(source_sel))
        cat("\n")
      }
     
      
      BLOCK_POST_sel<-BLOCK_POST[which(BLOCK_POST$AS_ID == AS_ID_sel),]
      
      if(DEBUG == 1)
      {
        cat("BLOCK_POST_sel_0\n")
        cat(str(BLOCK_POST_sel))
        cat("\n")
      }
      
      chr_sel<-unique(BLOCK_POST_sel$chr)
      
      if(DEBUG == 1)
      {
        cat("chr_sel_0\n")
        cat(str(chr_sel))
        cat("\n")
      }
      
      start_sel<-unique(BLOCK_POST_sel$start)
      
      if(DEBUG == 1)
      {
        cat("start_sel_0\n")
        cat(str(start_sel))
        cat("\n")
      }
      
      end_sel<-unique(BLOCK_POST_sel$end)
      
      if(DEBUG == 1)
      {
        cat("end_sel_0\n")
        cat(str(end_sel))
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
      cat(sprintf(as.character(AS_ID_sel)))
      cat("\t")
      cat(sprintf(as.character(chr_sel)))
      cat("\t")
      cat(sprintf(as.character(start_sel)))
      cat("\t")
      cat(sprintf(as.character(end_sel)))
      cat("\t")
      cat(sprintf(as.character(source_sel)))
      cat("\n")
      
      path_graphs<-paste(out,'Graphs_selected','/',paste(rs_sel,VAR_sel, sep='__'),'/',AS_ID_sel,'/',sep='')
      
      if (file.exists(path_graphs)){
        
        
      }else{
        
        dir.create(file.path(path_graphs))
        
      }#path_graphs
      
      
      variants_file_sel<-unique(variants_file[which(variants_file$AS_ID == AS_ID_sel),])
      
      variants_file_sel$is_cond_ind<-factor(variants_file_sel$is_cond_ind,
                                            levels=c("0","1"),
                                            ordered=T)
      
      if(DEBUG == 1)
      {
        cat("variants_file_sel_0\n")
        cat(str(variants_file_sel))
        cat("\n")
        cat(str(unique(variants_file_sel$VAR)))
        cat("\n")
      }
      
      #### PLOT BETA ----
      
     
      
      phenotype_sel<-unique(variants_file_sel$phenotype)
      
      if(DEBUG == 1)
      {
        cat("phenotype_sel_0\n")
        cat(str(phenotype_sel))
        cat("\n")
      }
      
      phenotype_DEF_sel<-unique(variants_file_sel$phenotype_DEF)
      
      if(DEBUG == 1)
      {
        cat("phenotype_DEF_sel_0\n")
        cat(str(phenotype_DEF_sel))
        cat("\n")
      }
      
      
      indx_min<-which(colnames(variants_file_sel) == 'min')
      indx_max<-which(colnames(variants_file_sel) == 'max')
      
      
      
      vector_SE<-c(variants_file_sel[,indx_min][!is.na(variants_file_sel[,indx_min])],
                   variants_file_sel[,indx_max][!is.na(variants_file_sel[,indx_max])])
      
      
      summary_vector_SE<-summary(vector_SE)
      
      
      
      if(DEBUG == 1)
      {
        cat("vector_SE\n")
        cat(sprintf(as.character(vector_SE)))
        cat("\n")
        cat("summary_vector_SE\n")
        cat(sprintf(as.character(names(summary_vector_SE))))
        cat("\n")
        cat(sprintf(as.character(summary_vector_SE)))
        cat("\n")
      }
      
      max_value<-summary_vector_SE[6]
      min_value<-summary_vector_SE[1]
      
   
      beta_span<-sort(as.numeric(unlist(strsplit(dataframe_to_print_SIT_plots_sel_AS_sel$beta_span, split=";"))))
      
      if(DEBUG == 1)
      {
        cat("beta_span\n")
        cat(sprintf(as.character(beta_span)))
        cat("\n")
      }
         
   
      #breaks_Beta<-unique(sort(unique(c(0,min_value,max_value,seq(min_value,max_value, by=step)))))
      #breaks_Beta<-unique(sort(unique(c(0,-0.1,0.1,min_value,max_value))))
     # breaks_Beta<-unique(sort(unique(c(0,-0.1,0.1,-0.2,0.2,0.3,-0.3,min_value,max_value))))
      breaks_Beta<-unique(sort(unique(beta_span)))
      labels_Beta<-as.character(round(breaks_Beta,2))
      
      limit_min<-min(min_value,beta_span)
      limit_max<-max(max_value,beta_span)
      
      if(DEBUG == 1)
      {
      
        cat("breaks_Beta\n")
        cat(sprintf(as.character(breaks_Beta)))
        cat("\n")
        cat("labels_Beta\n")
        cat(sprintf(as.character(labels_Beta)))
        cat("\n")
      }
      

      # geom_segment(data=variants_file_sel,
      #              aes(y=min,
      #                  yend=max,
      #                  x=pos37,
      #                  xend=pos37),
      #              color="black", size=2)+
        
      
      finemap_beta_dot_plot<-ggplot(data=variants_file_sel,
                      aes(x=pos37))+
        geom_segment(aes(y=breaks_Beta[1],
                     yend=finemap_beta,
                     x=pos37,
                     xend=pos37), size=0.5,color="gray",linetype="dashed")+
        scale_x_continuous(name=NULL,breaks=c(start_sel,end_sel),labels=as.character(c(start_sel,end_sel)),
                           limits=c(start_sel-1,end_sel+1))+
        geom_rect(data=variants_file_sel,
                  aes(ymin=min,
                      ymax=max,
                      xmin=pos37-10,
                      xmax=pos37+10),
                  color="black")+
        geom_point(data=variants_file_sel,
                   aes(x=pos37,
                       y=finemap_beta,
                       fill=Variant_classification,
                       shape=is_cond_ind),
                   size=1.5, stroke=0.5)+
        scale_shape_manual(values=c(21,23), drop=F)+
        scale_fill_manual(values=vector_colors, drop=F)+
        scale_y_continuous(name=paste("beta (s.d.)",phenotype_DEF_sel,sep=' '),breaks=breaks_Beta,labels=labels_Beta,
                           limits=c(limit_min,limit_max))+
        theme_classic()+
        theme(plot.title=element_text(size=6, color="black", family="sans"),
              axis.title.y=element_text(size=6, color="black", family="sans"),
              axis.title.x=element_text(size=6, color="black", family="sans"),
              axis.text.y=element_text(angle=0,size=6, color="black", family="sans"),
              axis.text.x=element_text(angle=0,size=6, color="black", family="sans"),
              axis.line.x = element_line(size = 0.2),
              axis.ticks.x = element_line(size = 0.2),
              axis.line.y = element_line(size = 0.2),
              axis.ticks.y = element_line(size = 0.2))+
        theme(legend.title = element_blank(),
              legend.text = element_text(size=6),
              legend.position="hidden")+
        geom_hline(yintercept=0, color="black",size=0.25)+
        ggeasy::easy_center_title()
      
      finemap_beta_dot_plot<-finemap_beta_dot_plot+
        geom_text_repel(data=variants_file_sel,
                      aes(x=pos37,
                          y=finemap_beta,
                          label=rs),
                      box.padding = 1,
                      color='black',
                      max.overlaps = Inf,
                      show.legend = FALSE,
                      size=2)
      
      if(DEBUG == 1)
      {
        cat("Effect size plot\n")
        
      }
      
      
      finemap_beta_dot_plot2<-ggplot(data=variants_file_sel,
                                     aes(x=pos37))+
        scale_x_continuous(name=NULL,breaks=c(start_sel,end_sel),labels=as.character(c(start_sel,end_sel)),
                           limits=c(start_sel-1,end_sel+1))+
        geom_point(data=variants_file_sel,
                   aes(x=pos37,
                       y=finemap_beta,
                       color=Variant_classification,
                       fill=Variant_classification,
                       shape=is_cond_ind),
                   size=1.5, stroke=0.5)+
        scale_shape_manual(values=c(21,23), drop=F)+
        scale_color_manual(values=vector_colors, drop=F)+
        scale_fill_manual(values=vector_colors, drop=F)+
        scale_y_continuous(name=paste("beta (s.d.)",phenotype_DEF_sel,sep=' '),breaks=breaks_Beta,labels=labels_Beta,
                           limits=c(limit_min,limit_max))+
        geom_hline(yintercept=0, color="black",size=0.25)+
        geom_text_repel(data=variants_file_sel,
                        aes(x=pos37,
                            y=finemap_beta,
                            label=rs),
                        box.padding = 1,
                        color='black',
                        max.overlaps = Inf,
                        show.legend = FALSE,
                        size=2)+
        theme_classic()+
        theme(plot.title=element_text(size=6, color="black", family="sans"),
              axis.title.y=element_text(size=6, color="black", family="sans"),
              axis.title.x=element_text(size=6, color="black", family="sans"),
              axis.text.y=element_text(angle=0,size=6, color="black", family="sans"),
              axis.text.x=element_text(angle=0,size=6, color="black", family="sans"),
              axis.line.x = element_line(size = 0.2),
              axis.ticks.x = element_line(size = 0.2),
              axis.line.y = element_line(size = 0.2),
              axis.ticks.y = element_line(size = 0.2))+
        theme(legend.title = element_text(size=6),
              legend.text = element_text(size=6),
              legend.key.size = unit(0.5, 'cm'), #change legend key size
              legend.key.height = unit(0.5, 'cm'), #change legend key height
              legend.key.width = unit(0.5, 'cm'), #change legend key width
              legend.position="right")+
        ggeasy::easy_center_title()

      if(DEBUG == 1)
      {
        cat("Effect size plot 2\n")
        
      }
      
      
  
      setwd(path_graphs)
      
      svgname<-paste('Finemap_Beta_COORD_order_2',".svg",sep='')
      makesvg = TRUE
      
      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= finemap_beta_dot_plot2,
               device="svg",
               height=5, width=5)
      }
      
      if(DEBUG == 1)
      {
        cat("Finemap_Beta_COORD_order plot 2\n")
        
      }
     
      
      
      #### gene tracks plot----
      
     
      
      gene_tracks_sel<-gene_tracks[which(gene_tracks$AS_ID == AS_ID_sel),]
      
      if(DEBUG == 1)
      {
        cat("gene_tracks_sel_0\n")
        cat(str(gene_tracks_sel))
        cat("\n")
      }
      
      if(dim(gene_tracks_sel)[1] >0)
      {
        finemap_beta_dot_plot<-finemap_beta_dot_plot+
          theme(axis.text.x=element_blank(),
                axis.line.x = element_blank(),
                axis.ticks.x = element_blank())
        
        
        
        Table_S7_DE_sel<-Table_S7[which(Table_S7$VAR == VAR_sel &
                                       Table_S7$Analysis == 'DE'),]
        
        
        if(dim(Table_S7_DE_sel)[1] >0)
        {
          
        }else{
          rescue_haplotypes_sel<-rescue_haplotypes[which(rescue_haplotypes$Proxy_VAR == VAR_sel &
                                                           rescue_haplotypes$ensembl_gene_id%in%gene_tracks_sel$ensembl_gene_id &
                                                           rescue_haplotypes$comparison == 'HOM_REF__HET|HET'),]
          
          if(DEBUG == 1)
          {
            cat("rescue_haplotypes_sel_0\n")
            cat(str(rescue_haplotypes_sel))
            cat("\n")
            cat(str(unique(rescue_haplotypes_sel$ensembl_gene_id)))
            cat("\n")
            cat(str(unique(rescue_haplotypes_sel$VAR)))
            cat("\n")
            cat(str(unique(rescue_haplotypes_sel$Proxy_VAR)))
            cat("\n")
            cat(str(unique(rescue_haplotypes_sel$comparison)))
            cat("\n")
          }
          
          my_DF<- data.frame(matrix(vector(), dim(rescue_haplotypes_sel)[1], 11,
                                    dimnames=list(c(),
                                                  colnames(Table_S7))),  stringsAsFactors=F)
          
          
          if(DEBUG == 1)
          {
            cat("my_DF_0\n")
            cat(str(my_DF))
            cat("\n")
            
          }
          
          my_DF$ensembl_gene_id<-rescue_haplotypes_sel$ensembl_gene_id
          my_DF$RNASeq_source<-'Whole blood'
          my_DF$Analysis<-'DE'
          my_DF$HGNC<-rescue_haplotypes_sel$HGNC
          my_DF$Beta<-rescue_haplotypes_sel$Block_PCHiC_Beta
          my_DF$Beta_Z_score<-rescue_haplotypes_sel$Block_PCHiC_Beta_Z_score
          my_DF$adjusted_pval<-rescue_haplotypes_sel$Block_PCHiC_minuslogpvalue
          my_DF$adjusted_minus_logpval<-rescue_haplotypes_sel$Block_PCHiC_minuslogpvalue
          
          
          if(DEBUG == 1)
          {
            cat("my_DF_1\n")
            cat(str(my_DF))
            cat("\n")
            
          }
            
          Table_S7_DE_sel<-my_DF
          
       
        }#dim(Table_S7_DE_sel)[1] >0
        
        if(DEBUG == 1)
        {
          cat("Table_S7_DE_sel_0\n")
          cat(str(Table_S7_DE_sel))
          cat("\n")
        }
        
        Table_S7_DE_sel<-Table_S7_DE_sel[which(Table_S7_DE_sel$RNASeq_source == source_sel),]
        
        if(DEBUG == 1)
        {
          cat("Table_S7_DE_sel_0.5\n")
          cat(str(Table_S7_DE_sel))
          cat("\n")
        }
        
       
        
        Table_S7_DE_sel$Beta_Z_score<-round(Table_S7_DE_sel$Beta_Z_score,2)
        
        A<-round(summary(Table_S7_DE_sel$Beta_Z_score[!is.na(Table_S7_DE_sel$Beta_Z_score)]),4)
        
        max_Beta_Z_score<-max(A)
        min_Beta_Z_score<-min(A)
        
        abs_max_Beta_Z_score<-abs(max_Beta_Z_score)
        abs_min_Beta_Z_score<-abs(min_Beta_Z_score)
        
        if(DEBUG == 1)
        {
          cat("max_Beta_Z_score_0\n")
          cat(str(max_Beta_Z_score))
          cat("\n")
          
          cat("min_Beta_Z_score_0\n")
          cat(str(min_Beta_Z_score))
          cat("\n")
          
          cat("abs_max_Beta_Z_score_0\n")
          cat(str(abs_max_Beta_Z_score))
          cat("\n")
          
          cat("abs_min_Beta_Z_score_0\n")
          cat(str(abs_min_Beta_Z_score))
          cat("\n")
        }
        
        if(abs_max_Beta_Z_score > abs_min_Beta_Z_score)
        {
          step<-abs(max_Beta_Z_score - -1*max_Beta_Z_score)/4
          breaks.EXP<-sort(unique(c(0,seq(from= -1*abs_max_Beta_Z_score, to=abs_max_Beta_Z_score,by=step))))
          breaks.EXP<-sort(unique(c(0,-1*abs_max_Beta_Z_score, abs_max_Beta_Z_score)))
          labels.EXP<-as.character(round(breaks.EXP,2))
          
        }else{
          
          step<-abs(min_Beta_Z_score - -1*min_Beta_Z_score)/4
          breaks.EXP<-sort(unique(c(0,seq(from= -1*abs_min_Beta_Z_score, to=abs_min_Beta_Z_score,by=step))))
          breaks.EXP<-sort(unique(c(0,-1*abs_min_Beta_Z_score,abs_min_Beta_Z_score,by=step)))
          labels.EXP<-as.character(round(breaks.EXP,2))
          
        }#abs_max_Beta_Z_score > abs_min_Beta_Z_score
       
        
        
        if(DEBUG == 1)
        {
          cat("labels.EXP:\t")
          cat(sprintf(as.character(labels.EXP)))
          cat("\n")
          
          
          
          # quit(status = 1)
        }
        
        
     
        
        
        
        indx.negative.strand<-which(gene_tracks_sel$strand == '-')
        
        if(DEBUG == 1)
        {
          cat("indx.negative.strand_0\n")
          cat(str(indx.negative.strand))
          cat("\n")
        }
        
    
        
        indx.positive.strand<-which(gene_tracks_sel$strand == '+')
        
        if(DEBUG == 1)
        {
          cat("indx.positive.strand_0\n")
          cat(str(indx.positive.strand))
          cat("\n")
        }
        
    
        
        gene_tracks_sel$COORD<-NA
        
        gene_tracks_sel$COORD[indx.positive.strand]<-1
        gene_tracks_sel$COORD[indx.negative.strand]<-0
        
        if(DEBUG == 1)
        {
          cat("gene_tracks_sel_1\n")
          cat(str(gene_tracks_sel))
          cat("\n")
        }
        
        gene_tracks_sel$definitive_start<-NA
        gene_tracks_sel$definitive_end<-NA
        
        gene_tracks_sel$definitive_start[indx.positive.strand]<-gene_tracks_sel$start[indx.positive.strand]
        gene_tracks_sel$definitive_end[indx.positive.strand]<-gene_tracks_sel$end[indx.positive.strand]
        
        
        gene_tracks_sel$definitive_start[indx.negative.strand]<-gene_tracks_sel$end[indx.negative.strand]
        gene_tracks_sel$definitive_end[indx.negative.strand]<-gene_tracks_sel$start[indx.negative.strand]
        
        gene_tracks_sel$ymin<-gene_tracks_sel$COORD
        gene_tracks_sel$ymax<-gene_tracks_sel$COORD+0.2
        
        if(DEBUG == 1)
        {
          cat("gene_tracks_sel_2\n")
          cat(str(gene_tracks_sel))
          cat("\n")
          
          cat("ymin_0\n")
          cat(sprintf(as.character(names(summary(gene_tracks_sel$ymin)))))
          cat("\n")
          cat(sprintf(as.character(summary(gene_tracks_sel$ymin))))
          cat("\n")
          
          cat("ymax_0\n")
          cat(sprintf(as.character(names(summary(gene_tracks_sel$ymax)))))
          cat("\n")
          cat(sprintf(as.character(summary(gene_tracks_sel$ymax))))
          cat("\n")
        }
        
        
        
        
        
        

        random_vec <- round(runif(n=dim(gene_tracks_sel)[1], min=-0.25, max=0.25),2)
        
        gene_tracks_sel$ymin<-gene_tracks_sel$ymin+random_vec
        gene_tracks_sel$ymax<-gene_tracks_sel$ymax+random_vec
        
        if(DEBUG == 1)
        {
          # cat("jitter_pos_0\n")
          # cat(str(jitter_pos))
          # cat("\n")
          
          cat("random_vec_0\n")
          cat(str(random_vec))
          cat("\n")
          
          cat("ymin_1\n")
          cat(sprintf(as.character(names(summary(gene_tracks_sel$ymin)))))
          cat("\n")
          cat(sprintf(as.character(summary(gene_tracks_sel$ymin))))
          cat("\n")
          
          cat("ymax_1\n")
          cat(sprintf(as.character(names(summary(gene_tracks_sel$ymax)))))
          cat("\n")
          cat(sprintf(as.character(summary(gene_tracks_sel$ymax))))
          cat("\n")
         
          
        }
        
        Table_S7_DE_sel<-merge(Table_S7_DE_sel,
                               gene_tracks_sel,
                               by='ensembl_gene_id',
                               all.y=T)
        
        Table_S7_DE_sel$SIG<-NA
        
        Table_S7_DE_sel$SIG[which(Table_S7_DE_sel$adjusted_minus_logpval >= 1.3)]<-'YES'
        Table_S7_DE_sel$SIG[which(Table_S7_DE_sel$adjusted_minus_logpval < 1.3)]<-'NO'
        
        Table_S7_DE_sel$SIG<-factor(Table_S7_DE_sel$SIG,
                                    levels=c('YES','NO'),
                                    ordered=T)
        
        
        
        
        
        if(DEBUG == 1)
        {
          cat("Table_S7_DE_sel_1\n")
          cat(str(Table_S7_DE_sel))
          cat("\n")
        }
        
        # geom_point(data=variants_file_sel,
        #            aes(x=pos37,
        #                y=1.5,
        #                fill=Variant_classification,
        #                shape=is_cond_ind),
        #            size=1.5, stroke=0.5)+
        #   scale_shape_manual(values=c(21,23), drop=F)+
        #   scale_fill_manual(values=vector_colors, drop=F)+
        #   new_scale("fill")+
        
        # geom_rect(data=gene_tracks_sel,
        #           aes(ymin=ymin,
        #               ymax=ymax,
        #               xmin=definitive_start,
        #               xmax=definitive_end), fill="gray", color='white',size=0.2)+
        
        
        graph_genes<-ggplot()+
          scale_x_continuous(name=NULL,breaks=c(start_sel,end_sel),
                             labels=as.character(c(start_sel,end_sel)),
                             limits=c(start_sel-1,end_sel+1))+
          scale_y_continuous(name=NULL, 
                             breaks=seq(0,1,by=1),
                             labels=as.character(c('gDNA -','gDNA +')), 
                             limits=c(-0.5,1.5))+
          geom_segment(data=variants_file_sel,
                       aes(y=-0.5,
                           yend=1.5,
                           x=pos37,
                           xend=pos37), size=0.5,color="gray",linetype="dashed")+
          geom_rect(data=Table_S7_DE_sel,
                    aes(ymin=ymin,
                        ymax=ymax,
                        xmin=definitive_start,
                        xmax=definitive_end,
                        fill=Beta_Z_score), color=NA,size=0.2)+
          scale_fill_gradient2(
            midpoint = 0,
            low = "blue",
            mid = "white",
            high = "red",
            breaks=breaks.EXP,labels=labels.EXP,
            limits=c(breaks.EXP[1],breaks.EXP[length(breaks.EXP)]),
            name=paste('Beta',"Z-score",sep="\n"),na.value = "gray")+
          new_scale("fill")+
          geom_rect(data=Table_S7_DE_sel[is.na(Table_S7_DE_sel$SIG),],
                    aes(ymin=ymin,
                        ymax=ymax,
                        xmin=definitive_start,
                        xmax=definitive_end), fill= NA, color='black',size=0.05)+
          geom_rect(data=Table_S7_DE_sel[which(Table_S7_DE_sel$SIG == 'NO'),],
                    aes(ymin=ymin,
                        ymax=ymax,
                        xmin=definitive_start,
                        xmax=definitive_end), fill= NA, color='black',size=0.05)+
          geom_rect(data=Table_S7_DE_sel[which(Table_S7_DE_sel$SIG == 'YES'),],
                    aes(ymin=ymin,
                        ymax=ymax,
                        xmin=definitive_start,
                        xmax=definitive_end), fill= NA, color='black',size=0.25)+
          geom_hline(yintercept=0.5, color="black",size=0.25)+
          theme_classic()+
          theme(plot.title=element_text(size=6, color="black", family="sans"),
                axis.title.y=element_blank(),
                axis.title.x=element_text(size=6, color="black", family="sans"),
                axis.text.y=element_text(angle=0,size=6, color="black", family="sans"),
                axis.text.x=element_text(angle=0,size=6, color="black", family="sans"),
                axis.line.x = element_line(size = 0.2),
                axis.ticks = element_line(size = 0.2),
                axis.line.y = element_line(size = 0.2))+
          theme(legend.title = element_text(size=6),
                legend.text = element_text(size=6),
                legend.key.size = unit(0.25, 'cm'), #change legend key size
                legend.key.height = unit(0.25, 'cm'), #change legend key height
                legend.key.width = unit(0.25, 'cm'), #change legend key width
                legend.position="bottom")+
          ggeasy::easy_center_title()
        
        
        
       
        
        dataframe_subset_ENSG_NO_NA_sel<-dataframe_subset_ENSG_NO_NA[which(dataframe_subset_ENSG_NO_NA$AS_ID == AS_ID_sel),]
        
        if(DEBUG == 1)
        {
          cat("dataframe_subset_ENSG_NO_NA_sel_2\n")
          cat(str(dataframe_subset_ENSG_NO_NA_sel))
          cat("\n")
        }
        
        
        indx.selected_genes<-which(gene_tracks_sel$ensembl_gene_id%in%dataframe_subset_ENSG_NO_NA_sel$ensembl_gene_id)
        
        if(DEBUG == 1)
        {
          cat("indx.selected_genes_0\n")
          cat(str(indx.selected_genes))
          cat("\n")
        }
        
        graph_genes2<-graph_genes+
          geom_text_repel(data=gene_tracks_sel,
                          aes(x=definitive_start,
                              y=ymax,
                              label=hgnc),
                          # Repel away from the left edge, not from the right.
                          xlim = c(NA, Inf),
                          # Do not repel from top or bottom edges.
                          ylim = c(-Inf, Inf),
                          family="sans",
                          color='black',
                          fontface='bold.italic',
                          size=2)
        
        graph_genes<-graph_genes+
          geom_text_repel(data=gene_tracks_sel[indx.selected_genes,],
                          aes(x=definitive_start,
                              y=ymax,
                              label=hgnc),
                          # Repel away from the left edge, not from the right.
                          xlim = c(NA, Inf),
                          # Do not repel from top or bottom edges.
                          ylim = c(-Inf, Inf),
                          family="sans",
                          color='black',
                          fontface='bold.italic',
                          size=2)
        
        setwd(path_graphs)
        
        svgname<-paste(paste("Genes",'COORD_order', sep='_'),".svg",sep='')
        makesvg = TRUE
        
        if (makesvg == TRUE)
        {
          ggsave(svgname, plot= graph_genes2,
                 device="svg",
                 height=dataframe_to_print_SIT_plots_sel_AS_sel$Height, width=dataframe_to_print_SIT_plots_sel_AS_sel$Width)
        }
        
        if(DEBUG == 1)
        {
          cat("genes plot\n")
          
        }
        
        
     
        
      }else{
        
        graph_genes<-NA
        
      }#dim(gene_tracks_sel)[1] >0
      
      
     
      
      # quit(status = 1)
     
     
    
      #### ATAC PLOT ----
      
      
      
      
      ATAC_tracks_subset.m_sel<-ATAC_tracks_subset.m[which(ATAC_tracks_subset.m$AS_ID == AS_ID_sel),]
      
      if(DEBUG == 1)
      {
        cat("ATAC_tracks_subset.m_sel_0\n")
        cat(str(ATAC_tracks_subset.m_sel))
        cat("\n")
      }
      
      if(dim(ATAC_tracks_subset.m_sel)[1] >0)
      {
        finemap_beta_dot_plot<-finemap_beta_dot_plot+
          theme(axis.text.x=element_blank(),
                axis.line.x = element_blank(),
                axis.ticks.x = element_blank())
        
        graph_genes<-graph_genes+
          theme(axis.text.x=element_blank(),
                axis.line.x = element_blank(),
                axis.ticks.x = element_blank())

        dataframe_subset_ATAC_NO_NA_sel<-dataframe_subset_ATAC_NO_NA[which(dataframe_subset_ATAC_NO_NA$AS_ID == AS_ID_sel),]
        
        if(DEBUG == 1)
        {
          cat("dataframe_subset_ATAC_NO_NA_sel_0\n")
          cat(str(dataframe_subset_ATAC_NO_NA_sel))
          cat("\n")
        }
        
        
       
      
        
        ATAC_tracks_subset.m_sel$Cell_Type<-factor(ATAC_tracks_subset.m_sel$Cell_Type,
                                                   levels = unique(dataframe_subset_ATAC_NO_NA_sel$Cell_Type),
                                                   ordered=T)
        
        
        if(DEBUG == 1)
        {
          cat("ATAC_tracks_subset.m_sel_2\n")
          cat(str(ATAC_tracks_subset.m_sel))
          cat("\n")
          cat(sprintf(as.character(names(summary(ATAC_tracks_subset.m_sel$Cell_Type)))))
          cat("\n")
          cat(sprintf(as.character(summary(ATAC_tracks_subset.m_sel$Cell_Type))))
          cat("\n")
          
        }
        
        
        
        ATAC_tracks_subset.m_sel$counts_by_10<-round(ATAC_tracks_subset.m_sel$counts,0)/10
        
        if(DEBUG == 1)
        {
          cat("ATAC_tracks_subset.m_sel_3\n")
          cat(str(ATAC_tracks_subset.m_sel))
          cat("\n")
          cat(sprintf(as.character(names(summary(ATAC_tracks_subset.m_sel$Cell_Type)))))
          cat("\n")
          cat(sprintf(as.character(summary(ATAC_tracks_subset.m_sel$Cell_Type))))
          cat("\n")
          
        }
        
        
        breaks_counts<-unique(sort(unique(c(0,50,max(ATAC_tracks_subset.m_sel$counts_by_10)))))
        labels_counts<-as.character(round(breaks_counts,1))
        
        if(DEBUG == 1)
        {
          cat("step_counts\n")
          cat(sprintf(as.character(step)))
          cat("\n")
          cat("breaks_counts\n")
          cat(sprintf(as.character(breaks_counts)))
          cat("\n")
          cat("labels_counts\n")
          cat(sprintf(as.character(labels_counts)))
          cat("\n")
        }
         
        ATAC_dot_plot<-ggplot()+
          geom_vline(xintercept=variants_file_sel$pos37, color="gray", linetype='dashed',size=0.5)+
          geom_rect(data=ATAC_tracks_subset.m_sel,
                    aes(xmin=start,
                        xmax=end,
                        ymin=0,
                        ymax=counts_by_10,
                        color=counts_by_10),
                        size=0.3)+
          scale_y_continuous(name="ATAC-seq",breaks=NULL,
                             limits=c(breaks_counts[1],breaks_counts[length(breaks_counts)]+1))+
          scale_x_continuous(name=NULL,breaks=c(start_sel,end_sel),labels=as.character(c(start_sel,end_sel)),
                             limits=c(start_sel-1,end_sel+1))+
          scale_color_viridis(
            alpha = 1,
            begin = 0,
            end = 1,
            direction = 1,
            discrete = FALSE,
            option = "D",
            aesthetics = "color",
            name=paste('ATAC counts',"x10",sep="\n"),
            na.value = NA,
            breaks=breaks_counts,
            labels=labels_counts,
            limits=c(breaks_counts[1],
                     breaks_counts[length(breaks_counts)])
          )
          
          # scale_color_gradient2(
          #   midpoint = 50,
          #   low = "black",
          #   mid = "yellow",
          #   high = "red",
          #   breaks=breaks_counts,labels=labels_counts,
          #   limits=c(breaks_counts[1],breaks_counts[length(breaks_counts)]),name=paste('ATAC reads',"x10",sep="\n"),na.value = "white")
        
        # scale_y_continuous(name="ATAC reads x10",breaks=breaks_counts,labels=labels_counts,
        #                    limits=c(breaks_counts[1],breaks_counts[length(breaks_counts)]+1))+
        
       
        
        
        if(DEBUG == 1)
        {
          cat("ATAC plot I \n")
          
        }
        
        # setwd(path_graphs)
        # 
        # svgname<-paste(paste("test", sep='_'),".svg",sep='')
        # makesvg = TRUE
        # 
        # if (makesvg == TRUE)
        # {
        #   ggsave(svgname, plot= ATAC_dot_plot,
        #          device="svg",
        #          height=10, width=4)
        # }
        
        
        ATAC_dot_plot<-ATAC_dot_plot+
          facet_grid(Cell_Type ~ ., scales='free_x', space='free_x', switch="y") +
          theme_cowplot(font_size = 8)+
          theme( strip.background = element_blank(),
                 strip.placement = "outside",
                 strip.text = element_text(size=6),
                 panel.spacing = unit(0.2, "lines"),
                 panel.background=element_rect(fill="white"),
                 panel.border=element_rect(colour="white",size=0,5),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())+
          theme_classic()+
          theme(plot.title=element_text(size=6, color="black", family="sans"),
                axis.title.y=element_text(size=6, color="black", family="sans"),
                axis.title.x=element_text(size=6, color="black", family="sans"),
                axis.text.y=element_text(angle=0,size=6, color="black", family="sans"),
                axis.text.x=element_text(angle=0,size=6, color="black", family="sans"),
                axis.line.x = element_line(size = 0.2),
                axis.ticks = element_line(size = 0.2),
                axis.line.y = element_line(size = 0.2))+
          theme(legend.title = element_text(size=6),
                legend.text = element_text(size=6),
                legend.key.size = unit(0.25, 'cm'), #change legend key size
                legend.key.height = unit(0.25, 'cm'), #change legend key height
                legend.key.width = unit(0.25, 'cm'), #change legend key width
                legend.position="bottom")+
          ggeasy::easy_center_title()
        
        
        
        
        if(DEBUG == 1)
        {
          cat("ATAC plot II \n")
          
        }
        
        setwd(out)

        write.table(ATAC_tracks_subset.m_sel,file='test.tsv', sep="\t",quote=F,row.names = F)

        setwd(path_graphs)

        svgname<-paste(paste("ATAC",'COORD_order', sep='_'),".svg",sep='')
        makesvg = TRUE

        if (makesvg == TRUE)
        {
          ggsave(svgname, plot= ATAC_dot_plot,
                 device="svg",
                 height=dataframe_to_print_SIT_plots_sel_AS_sel$Height, width=dataframe_to_print_SIT_plots_sel_AS_sel$Width)
        }
        
      }else{
        
        ATAC_dot_plot<-NA
      }#dim(ATAC_tracks_subset.m_sel)[1] >0
      
      
      
     
      
      #### PCHiC plot ----
      
     
      
      
      
      scale_df.m_subset_sel<-scale_df.m_subset[which(scale_df.m_subset$AS_ID == AS_ID_sel),]
      
      if(DEBUG == 1)
      {
        cat("scale_df.m_subset_sel_0\n")
        cat(str(scale_df.m_subset_sel))
        cat("\n")
      }
      
      if(dim(scale_df.m_subset_sel)[1] > 0)
      {
        
        if(dim(ATAC_tracks_subset.m_sel)[1] >0)
        {
          ATAC_dot_plot<-ATAC_dot_plot+
            theme(axis.text.x=element_blank(),
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank())
          
        }else{
          graph_genes<-graph_genes+
            theme(axis.text.x=element_blank(),
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank())  
          finemap_beta_dot_plot<-finemap_beta_dot_plot+
            theme(axis.text.x=element_blank(),
                  axis.line.x = element_blank(),
                  axis.ticks.x = element_blank())
          
          }#dim(ATAC_tracks_subset.m_sel)[1] >0
        
        
        dataframe_subset_PCHiC_NO_NA_sel<-dataframe_subset_PCHiC_NO_NA[which(dataframe_subset_PCHiC_NO_NA$AS_ID == AS_ID_sel),]
        
        if(DEBUG == 1)
        {
          cat("dataframe_subset_PCHiC_NO_NA_sel_0\n")
          cat(str(dataframe_subset_PCHiC_NO_NA_sel))
          cat("\n")
        }
        
        #### arc ----
        
        
        
        indx.arc<-c(which(colnames(scale_df.m_subset_sel) == 'AS_ID'),
                    which(colnames(scale_df.m_subset_sel) == 'oeStart'),which(colnames(scale_df.m_subset_sel) == 'oeEnd'),which(colnames(scale_df.m_subset_sel) == 'oeID'),
                    which(colnames(scale_df.m_subset_sel) == 'baitChr'),which(colnames(scale_df.m_subset_sel) == 'baitStart'),which(colnames(scale_df.m_subset_sel) == 'baitEnd'),which(colnames(scale_df.m_subset_sel) == 'baitID'),
                    which(colnames(scale_df.m_subset_sel) == 'baitSYMBOL'),
                    which(colnames(scale_df.m_subset_sel) == 'Cell_Type'),which(colnames(scale_df.m_subset_sel) == 'ChicagoScore'))
        
        
        
        arc_df<-unique(scale_df.m_subset_sel[,indx.arc])
        arc_df$chr<-paste('chr',arc_df$baitChr, sep='')
        arc_df$oe_middle_point<-as.integer(arc_df$oeStart + ((arc_df$oeEnd - arc_df$oeStart)/2))
        
        
        if(DEBUG == 1)
        {
          cat("arc_df_0\n")
          cat(str(arc_df))
          cat("\n")
        }
      
        
        arc_df$Cell_Type<-factor(arc_df$Cell_Type,
                                          levels = unique(dataframe_subset_PCHiC_NO_NA_sel$Cell_Type),
                                          ordered=T)
        
        if(DEBUG == 1)
        {
          cat("arc_df_1\n")
          cat(str(arc_df))
          cat("\n")
          cat(sprintf(as.character(names(summary(arc_df$Cell_Type)))))
          cat("\n")
          cat(sprintf(as.character(summary(arc_df$Cell_Type))))
          cat("\n")
        }
        
        arc_df$start<-NA
        
        arc_df$start<-arc_df$oe_middle_point
        
        # arc_df$start[which(arc_df$oe_middle_point >= arc_df$baitStart)]<-arc_df$baitStart
        # arc_df$start[which(arc_df$oe_middle_point < arc_df$baitStart)]<-arc_df$oe_middle_point
        
        arc_df$end<-NA
        
        arc_df$end<-arc_df$baitStart
        
        # arc_df$end[which(arc_df$oe_middle_point >= arc_df$baitStart)]<-arc_df$oe_middle_point
        # arc_df$end[which(arc_df$oe_middle_point < arc_df$baitStart)]<-arc_df$baitStart
        # arc_df$name<-paste(arc_df$baitSYMBOL,arc_df$oeID, sep='__')
        arc_df$name<-arc_df$baitSYMBOL
        
        
        if(DEBUG == 1)
        {
          cat("arc_df_2\n")
          cat(str(arc_df))
          cat("\n")
          cat(sprintf(as.character(names(summary(arc_df$Cell_Type)))))
          cat("\n")
          cat(sprintf(as.character(summary(arc_df$Cell_Type))))
          cat("\n")
          
        }
        
        indx.int<-c(which(colnames(arc_df) == 'chr'),which(colnames(arc_df) == 'start'),which(colnames(arc_df) == 'end'),
                    which(colnames(arc_df) == 'name'),which(colnames(arc_df) == 'Cell_Type'),
                    which(colnames(arc_df) == 'ChicagoScore'))
        
        arc_df_FINAL<-unique(arc_df[,indx.int])
        
        arc_df_FINAL$Feature<-'arc'
        arc_df_FINAL$y_start<-2
        arc_df_FINAL$y_end<-1
        
        if(DEBUG == 1)
        {
          cat("arc_df_FINAL_0\n")
          cat(str(arc_df_FINAL))
          cat("\n")
        }
        
        #### oe  and bait ----
        
        indx.oe<-c(which(colnames(scale_df.m_subset_sel) == 'oeChr'),which(colnames(scale_df.m_subset_sel) == 'oeStart'),which(colnames(scale_df.m_subset_sel) == 'oeEnd'),which(colnames(scale_df.m_subset_sel) == 'oeID'))
        
        oe_df<-unique(scale_df.m_subset_sel[,indx.oe])
        oe_df$chr<-paste('chr',oe_df$oeChr, sep='')
        oe_df$start<-oe_df$oeStart
        oe_df$end<-oe_df$oeEnd
        oe_df$name<-paste('oeID',as.character(oe_df$oeID),sep=' ')
        
        if(DEBUG == 1)
        {
          cat("oe_df_0\n")
          cat(str(oe_df))
          cat("\n")
          
        }
        
        indx.int<-c(which(colnames(oe_df) == 'chr'),which(colnames(oe_df) == 'start'),which(colnames(oe_df) == 'end'),
                    which(colnames(oe_df) == 'name'))
        
        oe_df_subset<-unique(oe_df[,indx.int])
        
        if(DEBUG == 1)
        {
          cat("oe_df_subset_0\n")
          cat(str(oe_df_subset))
          cat("\n")
          cat(sprintf(as.character(oe_df_subset$name)))
          cat("\n")
        }
        
        oe_df_subset_expanded<-do.call("rbind", replicate(length(levels(arc_df_FINAL$Cell_Type)), oe_df_subset, simplify = FALSE))
        oe_df_subset_expanded$Cell_Type<-rep(levels(arc_df_FINAL$Cell_Type), dim(oe_df_subset)[1])
        
        oe_df_subset_expanded$Cell_Type<-factor(oe_df_subset_expanded$Cell_Type,
                                                levels = levels(arc_df_FINAL$Cell_Type),
                                                ordered=T)
        
        oe_df_subset_expanded$ChicagoScore<-NA
        oe_df_subset_expanded$Feature<-'oe'
        oe_df_subset_expanded$y_start<-2
        oe_df_subset_expanded$y_end<-2
        
        if(DEBUG == 1)
        {
          cat("oe_df_subset_expanded_0\n")
          cat(str(oe_df_subset_expanded))
          cat("\n")
          cat(sprintf(as.character(names(summary(oe_df_subset_expanded$Cell_Type)))))
          cat("\n")
          cat(sprintf(as.character(summary(oe_df_subset_expanded$Cell_Type))))
          cat("\n")
          
        }
        
        
        indx.bait<-c(which(colnames(scale_df.m_subset_sel) == 'baitChr'),which(colnames(scale_df.m_subset_sel) == 'baitStart'),
                     which(colnames(scale_df.m_subset_sel) == 'baitEnd'),which(colnames(scale_df.m_subset_sel) == 'baitID'),which(colnames(scale_df.m_subset_sel) == 'baitSYMBOL'))
        
        bait_df<-unique(scale_df.m_subset_sel[,indx.bait])
        bait_df$chr<-paste('chr',bait_df$baitChr, sep='')
        bait_df$start<-bait_df$baitStart
        bait_df$end<-bait_df$baitEnd
        # bait_df$name<-bait_df$baitSYMBOL #paste(bait_df$baitID,bait_df$baitSYMBOL, sep= '__')
        bait_df$name<-paste('bait:',bait_df$baitSYMBOL,sep="\n")
        
        
        if(DEBUG == 1)
        {
          cat("bait_df_\n")
          cat(str(bait_df))
          cat("\n")
          
        }
        
        indx.int<-c(which(colnames(bait_df) == 'chr'),which(colnames(bait_df) == 'start'),which(colnames(bait_df) == 'end'),
                    which(colnames(bait_df) == 'name'))
        
        bait_df_subset<-unique(bait_df[,indx.int])
        
        if(DEBUG == 1)
        {
          cat("bait_df_subset_0\n")
          cat(str(bait_df_subset))
          cat("\n")
          
        }
        
        
        bait_df_subset_expanded<-do.call("rbind", replicate(length(levels(arc_df_FINAL$Cell_Type)), bait_df_subset, simplify = FALSE))
        bait_df_subset_expanded$Cell_Type<-rep(levels(arc_df_FINAL$Cell_Type), dim(bait_df_subset)[1])
        
        bait_df_subset_expanded$Cell_Type<-factor(bait_df_subset_expanded$Cell_Type,
                                                  levels = levels(arc_df_FINAL$Cell_Type),
                                                  ordered=T)
        
        bait_df_subset_expanded$ChicagoScore<-NA
        bait_df_subset_expanded$Feature<-'oe'
        bait_df_subset_expanded$y_start<-1
        bait_df_subset_expanded$y_end<-1
        
        if(DEBUG == 1)
        {
          cat("bait_df_subset_expanded_0\n")
          cat(str(bait_df_subset_expanded))
          cat("\n")
          cat(sprintf(as.character(names(summary(bait_df_subset_expanded$Cell_Type)))))
          cat("\n")
          cat(sprintf(as.character(summary(bait_df_subset_expanded$Cell_Type))))
          cat("\n")
          
        }
        
        
        REP<-rbind(oe_df_subset_expanded,bait_df_subset_expanded,arc_df_FINAL)
        
        if(DEBUG == 1)
        {
          cat("REP_0\n")
          cat(str(REP))
          cat("\n")
          cat(sprintf(as.character(names(summary(REP$Cell_Type)))))
          cat("\n")
          cat(sprintf(as.character(summary(REP$Cell_Type))))
          cat("\n")
          
        }
        
        #### PCHiC graph ----
        
     
        
        if(DEBUG == 1)
        {
          cat("Graph_Part_START:\n")
          
        }
        
        graph_PCHiC<-ggplot()+
          geom_curve(data=REP[which(REP$Feature == 'arc'),],
                     aes(x=start,
                         xend=end,
                         y=y_start,
                         yend=y_end,
                         color=ChicagoScore),
                     curvature = -0.025, size=1)+
          geom_segment(data=REP[which(REP$Feature == 'oe'),],
                       aes(x=start,
                           xend=end,
                           y=y_start,
                           yend=y_end),
                       color="black", size=2)+
          geom_segment(data=REP[which(REP$Feature == 'Bait'),],
                       aes(x=start,
                           xend=end,
                           y=y_start,
                           yend=y_end),
                       color="black", size=2)+
          scale_x_continuous(name=NULL,breaks=c(start_sel,end_sel),
                             labels=as.character(c(start_sel,end_sel)),
                             limits=c(start_sel-1,end_sel+1))+
          scale_y_continuous(name='PCHi-C',
                             breaks=NULL,
                             limits=c(0.75,2.25))+
          scale_color_gradient2(
            low = "black",
            mid = "yellow",
            high = "red",
            midpoint = 5,
            breaks=breaks_CS,labels=labels_CS,
            limits=c(breaks_CS[1],breaks_CS[length(breaks_CS)]),name=paste('Chicago',"score",sep="\n"),na.value = "white")
          
        
        if(DEBUG == 1)
        {
          cat("PCHiC plot I \n")
          
        }
        
        graph_PCHiC<-graph_PCHiC+
          facet_grid(Cell_Type ~ ., scales='free_x', space='free_x', switch="y") +
          theme_cowplot(font_size = 8)+
          theme( strip.background = element_blank(),
                 strip.placement = "outside",
                 strip.text = element_text(size=6),
                 panel.spacing = unit(0.2, "lines"),
                 panel.background=element_rect(fill="white"),
                 panel.border=element_rect(colour="white",size=0.5),
                 panel.grid.major = element_blank(),
                 panel.grid.minor = element_blank())+
          theme_classic()+
          theme(plot.title=element_text(size=6, color="black", family="sans"),
                axis.title.y=element_text(size=6, color="black", family="sans"),
                axis.title.x=element_text(size=6, color="black", family="sans"),
                axis.text.y=element_text(angle=0,size=6, color="black", family="sans"),
                axis.text.x=element_text(angle=0,size=6, color="black", family="sans"),
                axis.line.x = element_line(size = 0.2),
                axis.ticks = element_line(size = 0.2),
                axis.line.y = element_line(size = 0.2))+
          theme(legend.title = element_text(size=6),
                legend.text = element_text(size=6),
                legend.key.size = unit(0.25, 'cm'), #change legend key size
                legend.key.height = unit(0.25, 'cm'), #change legend key height
                legend.key.width = unit(0.25, 'cm'), #change legend key width
                legend.position="bottom")+
          geom_text_repel(data=bait_df_subset,
                          aes(x=end+1,
                              y=1,
                              label=name),
                          box.padding = 0.25,
                          color='black',
                          max.overlaps = Inf,
                          show.legend = FALSE,
                          size=2)+
          geom_text_repel(data=oe_df_subset,
                          aes(x=start-1,
                              y=2,
                              label=name),
                          box.padding = 0.25,
                          color='black',
                          max.overlaps = Inf,
                          show.legend = FALSE,
                          size=2)+
          geom_vline(xintercept=variants_file_sel$pos37, color="gray", linetype='dashed',size=0.5)+
          ggeasy::easy_center_title()
        
        
        
        
        if(DEBUG == 1)
        {
          cat("PCHiC plot II \n")
          
        }
        
        
        
        
        

        setwd(path_graphs)

        svgname<-paste(paste("PCHiC",'COORD_order', sep='_'),".svg",sep='')
        makesvg = TRUE

        if (makesvg == TRUE)
        {
          ggsave(svgname, plot= graph_PCHiC,
                 device="svg",
                 height=dataframe_to_print_SIT_plots_sel_AS_sel$Height, width=dataframe_to_print_SIT_plots_sel_AS_sel$Width)
        }
        
       
        
      }else{
        
        graph_PCHiC<-NA
        
      }#dim(scale_df.m_subset_sel)[1] > 0
      
      
      
      
      #### FINAL GRAPH----

      setwd(path_graphs)
      
      svgname<-paste(paste("Finemap_Beta",'COORD_order', sep='_'),".svg",sep='')
      makesvg = TRUE
      
      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= finemap_beta_dot_plot,
               device="svg",
               height=dataframe_to_print_SIT_plots_sel_AS_sel$Height, width=dataframe_to_print_SIT_plots_sel_AS_sel$Width)
      }
      
      if(DEBUG == 1)
      {
        cat("dataframe_to_print_SIT_plots_sel_AS_sel_REMEMBER\n")
        cat(str(dataframe_to_print_SIT_plots_sel_AS_sel))
        cat("\n")
        
      }
      
      
      graph_DEF<-plot_grid(finemap_beta_dot_plot,graph_genes,ATAC_dot_plot,graph_PCHiC,
                           nrow = 4,
                           ncol = 1,
                           rel_heights=c(dataframe_to_print_SIT_plots_sel_AS_sel$weight_A,
                                         dataframe_to_print_SIT_plots_sel_AS_sel$weight_B,
                                         dataframe_to_print_SIT_plots_sel_AS_sel$weight_C,
                                         dataframe_to_print_SIT_plots_sel_AS_sel$weight_D))

      


      setwd(path_graphs)

      svgname<-paste(paste("GRAPH_DEF",'COORD_order', sep='_'),".svg",sep='')
      makesvg = TRUE

      if (makesvg == TRUE)
      {
        ggsave(svgname, plot= graph_DEF,
               device="svg",
               height=dataframe_to_print_SIT_plots_sel_AS_sel$Height, width=dataframe_to_print_SIT_plots_sel_AS_sel$Width)
        
      }
      
      # if( DEBUG == 1)
      # {
      #   quit(status = 1)
      # }
      
    }# k in 1:length(AS_ID_array)
  }#i in 1:length(VARS)
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
    make_option(c("--dataframe_to_print_SIT_plots"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--BLOCK_POST"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--rescue_haplotypes"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Table_S7"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ATAC_tracks"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--gene_tracks"), type="character", default=NULL, 
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
  
  Printer_function(opt)

  
}


###########################################################################

system.time( main() )