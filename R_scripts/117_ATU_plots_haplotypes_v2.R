
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
  
  #### READ ATU_haplotype_variants ----
  
  ATU_haplotype_variants<-unlist(strsplit(opt$ATU_haplotype_variants, split=','))
  
  
  cat("ATU_haplotype_variants_0\n")
  cat(str(ATU_haplotype_variants))
  cat("\n")
  
  #### READ ENSG_selected ----
  
  ENSG_selected<-unlist(strsplit(opt$ENSG_selected, split=','))
  
  
  cat("ENSG_selected_0\n")
  cat(str(ENSG_selected))
  cat("\n")
  
  
  #### READ transcripts_selected ----
  
  transcripts_selected<-unlist(strsplit(opt$transcripts_selected, split=','))
  
  
  cat("transcripts_selected_0\n")
  cat(str(transcripts_selected))
  cat("\n")
  
  
 
  
  #### READ and transform path_INTERVAL ----
  
  path_INTERVAL = opt$path_INTERVAL
  
  cat("path_INTERVAL_\n")
  cat(sprintf(as.character(path_INTERVAL)))
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
  
  
  
  
  
 
  
  
  #### subset the ensembl gtf ----
  
  
  ensembl_gtf_sel<-unique(ensembl_gtf[which(ensembl_gtf$gene_id%in%ENSG_selected &
                                       ensembl_gtf$type%in%c('gene','transcript','exon','CDS','five_prime_utr','three_prime_utr')),])
  
  rm(ensembl_gtf)
  
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
  
  
  path_graphs<-paste(out,'ATU_graphs_Haplotypes','/',sep='')
  
  if (file.exists(path_graphs)){
    
    
  }else{
    
    dir.create(file.path(path_graphs))
    
  }#path_graphs
  
  
  
  
  
 
  
  DEBUG<-1
  
  for(i in 1:length(ATU_haplotype_variants))
  {
    ATU_haplotype_variants_sel<-ATU_haplotype_variants[i]
    filename<-paste(path_INTERVAL,ATU_haplotype_variants_sel, sep='')
    
    #chr2_219020958_C_T/Haplotypes/chr2_219020958_C_T__chr2_219214529_A_T/INTERVAL_covariates_and_PEER_factors_Haplotype_chr2_219020958_C_T__chr2_219214529_A_T.rds
    
    
    haplotype_sel<-gsub("^[^/]+/[^/]+/","",ATU_haplotype_variants_sel)
    haplotype_sel<-gsub("/.+$","",haplotype_sel)
    
    cat("---------------------------------------------------------------------------------->\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(ATU_haplotype_variants_sel)))
    cat("\t")
    cat(sprintf(as.character(haplotype_sel)))
    cat("\n")
    
    INTERVAL_covariates_and_PEER_factors_sel<-readRDS(file=filename)
    
    if(DEBUG == 1)
    {
      cat("INTERVAL_covariates_and_PEER_factors_sel_0\n")
      cat(str(INTERVAL_covariates_and_PEER_factors_sel))
      cat("\n")
      cat(sprintf(as.character(names(summary(INTERVAL_covariates_and_PEER_factors_sel$Haplotype)))))
      cat("\n")
      cat(sprintf(as.character(summary(INTERVAL_covariates_and_PEER_factors_sel$Haplotype))))
      cat("\n")
      
    }
 
    
    path_graphs<-paste(out,'ATU_graphs_Haplotypes','/',haplotype_sel,'/',sep='')
    
    if (file.exists(path_graphs)){
      
      
    }else{
      
      dir.create(file.path(path_graphs))
      
    }#path_graphs
    
    #### ENSG ----
    
    ENSG_array_sel<-ENSG_selected[i]
    
    cat(sprintf(as.character(ENSG_array_sel)))
    cat("\n")
    
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
    
    path_graphs<-paste(out,'ATU_graphs_Haplotypes','/',haplotype_sel,'/',HGNC_sel,'/',sep='')
    
    if (file.exists(path_graphs)){
      
      
    }else{
      
      dir.create(file.path(path_graphs))
      
    }#path_graphs
    
    
    transcripts_array_sel<-transcripts_selected[i]
    
    cat(sprintf(as.character(transcripts_array_sel)))
    cat("\n")
    
    
    transcript_set<-rev(unlist(strsplit(transcripts_array_sel, split="__")))
    
    if(DEBUG == 1)
    {
      cat("transcript_set_0\n")
      cat(str(transcript_set))
      cat("\n")
    }
    
    #### transcript plot ----
    
  
    
    ensembl_gtf_sel_transcript_sel<-ensembl_gtf_sel[which(ensembl_gtf_sel$transcript_id %in% transcript_set &
                                                            ensembl_gtf_sel$type == 'transcript'),]
    
    if(DEBUG == 1)
    {
      cat("ensembl_gtf_sel_transcript_sel_0\n")
      cat(str(ensembl_gtf_sel_transcript_sel))
      cat("\n")
    }
    
   
    
    
    REP<-data.frame()
    for(l in 1:length(transcript_set))
    {
      
      transcript_set_sel<-transcript_set[l]
      
      
      # cat("---------------------------------------------------------------------------------->\t")
      # cat(sprintf(as.character(l)))
      # cat("\t")
      # cat(sprintf(as.character(transcript_set_sel)))
      # cat("\n")
      
      ensembl_gtf_sel_transcript_sel<-ensembl_gtf_sel[which(ensembl_gtf_sel$transcript_id == transcript_set_sel &
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
      
    }#l in 1:length(transcript_set)
    
    
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
                              levels=rev(transcript_set),
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
            legend.position="bottom")+
      guides(fill=guide_legend(nrow=2,byrow=TRUE))
    
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
    
    
    #### Blood violin plot ----
    
   
    Transposed_Isoform_Expression_sel<-Transposed_Isoform_Expression[,c(which(colnames(Transposed_Isoform_Expression) == "sample_id"),
                                                                        which(colnames(Transposed_Isoform_Expression) %in% transcript_set))]
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
    
    Ratio_df_HET<-droplevels(Ratio_df[which(Ratio_df$Haplotype%in%c('HOM_REF','HET|HOM_REF','HOM_REF|HET','HET|HET')),])
    
    Ratio_df_HET$transcript_id<-factor(Ratio_df_HET$transcript_id,
                                       levels=transcript_set,
                                       ordered=T)
    
    Ratio_df_HET$Cell_Type<-'Whole blood'
    
    
    Ratio_df_HET$Haplotype<-factor(Ratio_df_HET$Haplotype,
                                  levels=c('HOM_REF','HET|HOM_REF','HOM_REF|HET','HET|HET'),
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
    
    
    vector.haplotypes<-c(brewer.pal(4, "Spectral"))
    
    
    
    Blood_violin<-ggplot(data=Ratio_df_HET,
                         aes(x=Haplotype, y=TPM, fill=Haplotype)) +
      geom_violin()+
      stat_summary(fun = median, fun.min = median, fun.max = median,
                   geom = "crossbar", width = 0.5)+
      scale_y_continuous(name="TPM",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
      scale_fill_manual(values=vector.haplotypes,drop=F)+
      scale_x_discrete(name=NULL, drop=F)
    
    
    
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
      theme(axis.title.y=element_text(size=8, color="black", family="sans"),
            axis.title.x=element_blank(),
            axis.text.y=element_text(angle=0,size=6, color="black", family="sans"),
            axis.text.x=element_blank(),
            axis.line.x = element_line(size = 0.2),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(size = 0.2),
            axis.line.y = element_line(size = 0.2))+
      theme(legend.title = element_text(size=6, color="black", family="sans"),
            legend.text = element_text(size=6, color="black", family="sans"),
            legend.key.size = unit(0.4, 'cm'), #change legend key size
            legend.key.height = unit(0.4, 'cm'), #change legend key height
            legend.key.width = unit(0.4, 'cm'), #change legend key width
            legend.position="bottom")+
      guides(fill=guide_legend(nrow=2,byrow=TRUE))+
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
    
    
    
  
  
  graph_DEF<-plot_grid(transcript_plot,Blood_violin,
                       nrow = 1,
                       ncol = 2,
                       rel_widths=c(1,1))
  
  
  
  
  setwd(path_graphs)
  
  svgname<-paste(paste("GRAPH_DEF",'transcripts', sep='_'),".svg",sep='')
  makesvg = TRUE
  
  if (makesvg == TRUE)
  {
    ggsave(svgname, plot= graph_DEF,
           device="svg",
           height=3, 
           width=3)
    
  } 
  
 
  
  
  array_transcript_ids<-levels(Ratio_df_HET$transcript_id)
  
  if(DEBUG == 1)
  {
    
    cat("array_transcript_ids_0\n")
    cat(str(array_transcript_ids))
    cat("\n")
  }
  
  
  for(iteration_transcript_ids in 1:length(array_transcript_ids)){
    
    transcript_id_sel<-array_transcript_ids[iteration_transcript_ids]
    
    cat("--------------------------------------------->\t")
    cat(sprintf(as.character(iteration_transcript_ids)))
    cat("\t")
    cat(sprintf(as.character(transcript_id_sel)))
    cat("\n")
    
    Ratio_df_HET_transcript_id_sel<-droplevels(Ratio_df_HET[which(Ratio_df_HET$transcript_id == transcript_id_sel),])
    
    if(DEBUG == 1)
    {
      
      cat("Ratio_df_HET_transcript_id_sel_0\n")
      cat(str(Ratio_df_HET_transcript_id_sel))
      cat("\n")
    }
    
    
    ##### Violin visual TPM  sel ----
    
    A<-round(summary(Ratio_df_HET_transcript_id_sel$TPM[!is.na(Ratio_df_HET_transcript_id_sel$TPM)]),2)
    
    if(DEBUG == 1)
    {
      cat("summary_TPM\n")
      cat(sprintf(as.character(names(A))))
      cat("\n")
      cat(sprintf(as.character(A)))
      cat("\n")
    }
    
    step<-abs(A[6]-A[1])/3
    
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
    
    
    vector.haplotypes<-c(brewer.pal(4, "Spectral"))
    
    
    
    Blood_violin<-ggplot(data=Ratio_df_HET_transcript_id_sel,
                         aes(x=Haplotype, y=TPM, fill=Haplotype)) +
      geom_violin()+
      stat_summary(fun = median, fun.min = median, fun.max = median,
                   geom = "crossbar", width = 0.5)+
      scale_y_continuous(name="TPM",breaks=breaks.Rank,labels=labels.Rank, limits=c(breaks.Rank[1],breaks.Rank[length(breaks.Rank)]))+
      scale_fill_manual(values=vector.haplotypes,drop=F)+
      scale_x_discrete(name=NULL, drop=F)
    
    
    
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
      theme(axis.title.y=element_text(size=8, color="black", family="sans"),
            axis.title.x=element_blank(),
            axis.text.y=element_text(angle=0,size=6, color="black", family="sans"),
            axis.text.x=element_blank(),
            axis.line.x = element_line(size = 0.2),
            axis.ticks.x = element_blank(),
            axis.ticks.y = element_line(size = 0.2),
            axis.line.y = element_line(size = 0.2))+
      theme(legend.title = element_text(size=6, color="black", family="sans"),
            legend.text = element_text(size=6, color="black", family="sans"),
            legend.key.size = unit(0.4, 'cm'), #change legend key size
            legend.key.height = unit(0.4, 'cm'), #change legend key height
            legend.key.width = unit(0.4, 'cm'), #change legend key width
            legend.position="bottom")+
      guides(fill=guide_legend(nrow=2,byrow=TRUE))+
      ggeasy::easy_center_title()
    
    
    setwd(path_graphs)
    
    svgname<-paste(paste('Whole_Blood_violin',transcript_id_sel, sep='_'),".svg",sep='')
    makesvg = TRUE
    
    if (makesvg == TRUE)
    {
      ggsave(svgname, plot= Blood_violin,
             device="svg",
             height=3, width=3)
    }
    
    
  }#iteration_transcript_ids in 1:length(array_transcript_ids)
}#i in 1:length(ATU_haplotype_variants

  
  
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
    make_option(c("--ATU_haplotype_variants"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ENSG_selected"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--transcripts_selected"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--path_INTERVAL"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--Transposed_Isoform_Expression"), type="character", default=NULL, 
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
  