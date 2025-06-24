
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
suppressMessages(library("liftOver", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))
suppressMessages(library("rCNV", lib.loc="/home/manuel.tardaguila/R/x86_64-pc-linux-gnu-library/4.1/"))


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
  
  
  
  #### READ Table_S1 ----
  
  
  Table_S1<-as.data.frame(readRDS(file=opt$Table_S1) , stringsAsFactors=F)
 

  cat("Table_S1_0\n")
  cat(str(Table_S1))
  cat("\n")
  cat(str(unique(Table_S1$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_S1$Variant_classification))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_S1$Variant_classification)))))
  cat("\n")
  
  Table_S1_subset<-Table_S1[which(Table_S1$Variant_classification == 'index_variants'),]
  
  cat("Table_S1_subset_0\n")
  cat(str(Table_S1_subset))
  cat("\n")
  cat(str(unique(Table_S1_subset$VAR)))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(Table_S1_subset$Variant_classification))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(Table_S1_subset$Variant_classification)))))
  cat("\n")
 
  #### READ ALL_db ----
  
  ALL_db<-as.data.frame(fread(file=opt$ALL_db, sep="\t", header =T), stringAsFactors=F)
  
  ALL_db<-unique(ALL_db[,-which(colnames(ALL_db) == 'maf_origin')])
  
  ALL_db$AS_ID<-paste(ALL_db$phenotype,ALL_db$block_no, sep="__")
  
  # Reverse the alignment
  
  ALL_db$finemap_z<--1*ALL_db$finemap_z
  ALL_db$finemap_beta<--1*ALL_db$finemap_beta
  
  #### Nicole's new phenotype labels ----
  
  ALL_db$phenotype_DEF<-revalue(ALL_db$phenotype,
                                       c("plt"="PLT#",
                                         "mpv"="MPV",
                                         "pdw"="PDW",
                                         "pct"="PCT",
                                         "H-IPF"="H-IPF",
                                         "IPF_perc"="IPF%",
                                         "IPF"="IPF#",
                                         "P-LCR"="P-LCR",
                                         "rbc"="RBC#",
                                         "mcv"="MCV",
                                         "hct"="HCT",
                                         "mch"="MCH",
                                         "mchc"="MCHC",
                                         "hgb"="HGB",
                                         "rdw_cv"='RDW',
                                         "MacrorR"='MacroR',
                                         "MicroR"='MicroR',
                                         "RBC-He"='RBC-He',
                                         "Delta-He"='Delta-He',
                                         "Hyer-He"='Hyper-He',
                                         "RDW-SD"='RDW-SD',
                                         "RPI"='RPI',
                                         "HFR"='HFR',
                                         "MFR"='MFR',
                                         "LFR"='LFR',
                                         "IG"='IG#',
                                         "IG_perc"='IG%',
                                         "ret"='RET#',
                                         "ret_p"='RET%',
                                         "irf"='IRF',
                                         "hlr"='HLSR#',
                                         "hlr_p"='HLSR%',
                                         "mrv"='MRV',
                                         "mscv"='MSCV',
                                         "RET-FSC"='RET-FSC',
                                         "RET-RBC-FSC"='RET-FSC',
                                         "IRF-FSC"='IRF-FSC',
                                         "RET-He"='RET-He',
                                         "RET-UPP"='RET-UPP',
                                         "lymph"='LYMPH#',
                                         "lymph_p"='LYMPH%',
                                         "LY-FSC"="LY-FSC",
                                         "LY-FSC-DW"="LY-FSC-DW",
                                         "LY-SSC"="LY-SSC",
                                         "LY-SSC-DW"="LY-SSC-DW",
                                         "LY-SFL"="LY-SFL",
                                         "LY-SFL-DW"="LY-SFL-DW",
                                         "baso"='BASO#',
                                         "baso_p"='BASO%',
                                         "eo"='EO#',
                                         "eo_p"='EO%',
                                         "mono"='MONO#',
                                         "mono_p"='MONO%',
                                         "MO-FSC"="MO-FSC",
                                         "MO-FSC-DW"="MO-FSC-DW",
                                         "MO-SSC"="MO-SSC",
                                         "MO-SSC-DW"="MO-SSC-DW",
                                         "MO-SFL"="MO-SFL",
                                         "MO-SFL-DW"="MO-SFL-DW",
                                         "neut"='NEUT#',
                                         "neut_p"='NEUT%',
                                         "NE-FSC"="NE-FSC",
                                         "NE-FSC-DW"="NE-FSC-DW",
                                         "NE-SSC"="NE-SSC",
                                         "NE-SSC-DW"="NE-SSC-DW",
                                         "NE-SFL"="NE-SFL",
                                         "NE-SFL-DW"="NE-SFL-DW",
                                         "wbc"='WBC#'))
  
  
  ALL_db$phenotype_DEF<-factor(as.character(ALL_db$phenotype_DEF),
                                      levels=c("PLT#","MPV","PDW","PCT","H-IPF","IPF%","IPF#","P-LCR",
                                               "RBC#","MCV","HCT","MCH","MCHC","HGB",'RDW','MacroR','MicroR','RBC-He','Delta-He','Hyper-He','RDW-SD','RPI','HFR','MFR','LFR','IG#','IG%',
                                               'RET#','RET%','IRF','HLSR#','HLSR%','MRV','MSCV',"RET-FSC","RET-RBC-FSC","IRF-FSC","RET-He","RET-UPP",
                                               'LYMPH#','LYMPH%',"LY-FSC","LY-FSC-DW","LY-SSC","LY-SSC-DW","LY-SFL","LY-SFL-DW",
                                               'BASO#','BASO%','EO#','EO%',
                                               'MONO#','MONO%',"MO-FSC","MO-FSC-DW","MO-SSC","MO-SSC-DW","MO-SFL","MO-SFL-DW",
                                               'NEUT#','NEUT%',"NE-FSC","NE-FSC-DW","NE-SSC","NE-SSC-DW","NE-SFL","NE-SFL-DW",
                                               'WBC#'), ordered = T)
  
  ALL_db<-droplevels(ALL_db)
  
  
  cat("ALL_db_0\n")
  cat(str(ALL_db))
  cat("\n")
  cat(str(unique(ALL_db$VAR)))
  cat("\n")
  
  
  
  ALL_db$chr<-gsub("_.+$","",ALL_db$VAR)
  ALL_db$pos<-gsub("^[^_]+_","",ALL_db$VAR)
  ALL_db$pos<-as.integer(gsub("_.+$","",ALL_db$pos))
  ALL_db$ref<-gsub("^[^_]+_[^_]+_","",ALL_db$VAR)
  ALL_db$ref<-gsub("_.+$","",ALL_db$ref)
  ALL_db$alt<-gsub("^[^_]+_[^_]+_[^_]+_","",ALL_db$VAR)
  
  cat("ALL_db_1\n")
  cat(str(ALL_db))
  cat("\n")
  cat(sprintf(as.character(names(summary(as.factor(ALL_db$chr))))))
  cat("\n")
  cat(sprintf(as.character(summary(as.factor(ALL_db$chr)))))
  cat("\n")
  
  ALL_db$chr<-factor(ALL_db$chr,
                       levels=c("chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11",
                                "chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21",
                                "chr22","chr23","chrX","chrY"), ordered=T)
  
  
  ALL_db<-ALL_db[order(ALL_db$chr,ALL_db$pos),]
  
  cat("ALL_db_2\n")
  cat(str(ALL_db))
  cat("\n")
  cat(str(unique(ALL_db$chr)))
  cat("\n")
  cat(str(unique(ALL_db$pos)))
  cat("\n")
  cat(str(unique(ALL_db$ref)))
  cat("\n")
  cat(str(unique(ALL_db$alt)))
  cat("\n")
  
  #### Get only the AS of index variants ----
  
  
  
  ALL_db_subset<-unique(ALL_db[which(ALL_db$VAR%in%Table_S1_subset$VAR |
                                       ALL_db$VAR == 'chr9_135864513_C_T'),])
  
  cat("ALL_db_subset_0\n")
  cat(str(ALL_db_subset))
  cat("\n")
  cat(str(unique(ALL_db_subset$VAR)))
  cat("\n")
  
  VAR_vector<-unique(ALL_db_subset$VAR)
  
  cat("VAR_vector_0\n")
  cat(str(VAR_vector))
  cat("\n")
  
 
 
  
  #### LOOP PRINTING ----
  
  path_files<-paste(out,'Build_files','/',sep='')
  
  if (file.exists(path_files)){
    
    
  }else{
    
    dir.create(file.path(path_files))
    
  }#path_files
  
  
  
  DEBUG<-0
  
  Result_FINAL<-data.frame()
  
  for(i in 1:length(VAR_vector))
  {
    VAR_sel<-VAR_vector[i]
    
    cat("------------------------------------------------------------------------>\t")
    cat(sprintf(as.character(i)))
    cat("\t")
    cat(sprintf(as.character(VAR_sel)))
    cat("\n")
    
    # if(VAR_sel == 'chr1_158613314_G_A')
    # {
    # 
    #   DEBUG<-1
    # 
    # }else{
    # 
    #   DEBUG<-0
    # }
    
    
    ALL_db_subset_sel<-ALL_db_subset[which(ALL_db_subset$VAR == VAR_sel),]
    
    if(DEBUG == 1)
    {
      cat("ALL_db_subset_sel_0\n")
      cat(str(ALL_db_subset_sel))
      cat("\n")
    }
    
    
    rsid_sel<-unique(ALL_db_subset_sel$rs)
    
    if(DEBUG == 1)
    {
      cat("rsid_sel_0\n")
      cat(str(rsid_sel))
      cat("\n")
    }
    
  
    
    AS_ID_array_sel<-unique(ALL_db_subset_sel$AS_ID)
    
    if(DEBUG == 1)
    {
      cat("AS_ID_array_sel_0\n")
      cat(str(AS_ID_array_sel))
      cat("\n")
    }
    
    Result_AS<-data.frame()
    
    for(k in 1:length(AS_ID_array_sel))
    {
      AS_ID_sel<-AS_ID_array_sel[k]
      
      cat("--------------------->\t")
      cat(sprintf(as.character(k)))
      cat("\t")
      cat(sprintf(as.character(AS_ID_sel)))
      cat("\n")
      
     
      
      ALL_db_restricted_sel<-ALL_db[which(ALL_db$AS_ID == AS_ID_sel),]
      
      if(DEBUG == 1)
      {
        cat("ALL_db_restricted_sel_0\n")
        cat(str(ALL_db_restricted_sel))
        cat("\n")
        cat(str(unique(ALL_db_restricted_sel$VAR)))
        cat("\n")
      }
      
      phenotype_sel<-unique(ALL_db_restricted_sel$phenotype)
      
      if(DEBUG == 1)
      {
        cat("phenotype_sel_0\n")
        cat(str(phenotype_sel))
        cat("\n")
      }
      
      phenotype_DEF_sel<-unique(ALL_db_restricted_sel$phenotype_DEF)
      
      if(DEBUG == 1)
      {
        cat("phenotype_DEF_sel_0\n")
        cat(str(phenotype_DEF_sel))
        cat("\n")
      }
      
      
      if(VAR_sel == 'chr2_74920648_G_A')
      {
        ALL_db_restricted_sel_cond_ind<-ALL_db_restricted_sel[which(ALL_db_restricted_sel$rs%in%c('rs183347927','rs568602257')),]
        
        # ALL_db_restricted_sel_cond_ind$factor_finemap_z<-factor()
        
        if(DEBUG == 1)
        {
          cat("ALL_db_restricted_sel_cond_ind_0\n")
          cat(str(ALL_db_restricted_sel_cond_ind))
          cat("\n")
          cat(str(unique(ALL_db_restricted_sel_cond_ind$VAR)))
          cat("\n")
          cat(str(unique(ALL_db_restricted_sel_cond_ind$rs)))
          cat("\n")
        }
        
      }else{
        
        if(VAR_sel == 'chr2_219020958_C_T')
        {
          ALL_db_restricted_sel_cond_ind<-ALL_db_restricted_sel[which(ALL_db_restricted_sel$is_cond_ind == 1 |
                                                                        ALL_db_restricted_sel$rs == 'rs192280464'),]
          
          if(DEBUG == 1)
          {
            cat("ALL_db_restricted_sel_cond_ind_0\n")
            cat(str(ALL_db_restricted_sel_cond_ind))
            cat("\n")
            cat(str(unique(ALL_db_restricted_sel_cond_ind$VAR)))
            cat("\n")
          }
        }else{
          
          if(VAR_sel == 'chr9_135920196_C_T')
          {
            ALL_db_restricted_sel_cond_ind<-ALL_db_restricted_sel[which(ALL_db_restricted_sel$is_cond_ind == 1 |
                                                                          ALL_db_restricted_sel$phenotype == 'mpv'),]
            
 
            
            if(DEBUG == 1)
            {
              cat("ALL_db_restricted_sel_cond_ind_0\n")
              cat(str(ALL_db_restricted_sel_cond_ind))
              cat("\n")
              cat(str(unique(ALL_db_restricted_sel_cond_ind$VAR)))
              cat("\n")
            } 
            
          }else{
            
            if(VAR_sel == 'chr2_24091099_C_T')
            {
              ALL_db_restricted_sel_cond_ind<-ALL_db_restricted_sel[which(ALL_db_restricted_sel$is_cond_ind == 1 |
                                                                            ALL_db_restricted_sel$rs == 'rs138903557'),]
              
              if(DEBUG == 1)
              {
                cat("ALL_db_restricted_sel_cond_ind_0\n")
                cat(str(ALL_db_restricted_sel_cond_ind))
                cat("\n")
                cat(str(unique(ALL_db_restricted_sel_cond_ind$VAR)))
                cat("\n")
              }
              
            }else{
              
              if(VAR_sel == 'chr22_50949811_T_C')
              {
                ALL_db_restricted_sel_cond_ind<-ALL_db_restricted_sel[which(ALL_db_restricted_sel$is_cond_ind == 1 |
                                                                              ALL_db_restricted_sel$rs == 'rs74624637'),]
                
                if(DEBUG == 1)
                {
                  cat("ALL_db_restricted_sel_cond_ind_0\n")
                  cat(str(ALL_db_restricted_sel_cond_ind))
                  cat("\n")
                  cat(str(unique(ALL_db_restricted_sel_cond_ind$VAR)))
                  cat("\n")
                }
                
              }else{
                
                ALL_db_restricted_sel_cond_ind<-ALL_db_restricted_sel[which(ALL_db_restricted_sel$is_cond_ind == 1),]
                
                if(DEBUG == 1)
                {
                  cat("ALL_db_restricted_sel_cond_ind_0\n")
                  cat(str(ALL_db_restricted_sel_cond_ind))
                  cat("\n")
                  cat(str(unique(ALL_db_restricted_sel_cond_ind$VAR)))
                  cat("\n")
                }  
              }#VAR_sel == 'chr22_50949811_T_C'
            }#VAR_sel == 'chr2_24091099_C_T'
          }#VAR_sel == 'chr9_135920196_C_T'
        }#VAR_sel == 'chr2_219020958_C_T'
      }#if(VAR_sel == 'chr2_74920648_G_A')
      
      
      
      if(dim(ALL_db_restricted_sel_cond_ind)[1] >0)
      {
        ALL_db_restricted_sel_cond_ind<-merge(ALL_db_restricted_sel_cond_ind,
                                              Table_S1,
                                              by=c("VAR","rs"),
                                              all.x=T)
        
        if(DEBUG == 1)
        {
          cat("ALL_db_restricted_sel_cond_ind_1\n")
          cat(str(ALL_db_restricted_sel_cond_ind))
          cat("\n")
          cat(str(unique(ALL_db_restricted_sel_cond_ind$VAR)))
          cat("\n")
        }
        
        
        #### calculate the segment that is going to be the SE -----
        
        ALL_db_restricted_sel_cond_ind$min<-ALL_db_restricted_sel_cond_ind$finemap_beta - ALL_db_restricted_sel_cond_ind$finemap_se
        ALL_db_restricted_sel_cond_ind$max<-ALL_db_restricted_sel_cond_ind$finemap_beta + ALL_db_restricted_sel_cond_ind$finemap_se
        
        ALL_db_restricted_sel_cond_ind<-ALL_db_restricted_sel_cond_ind[order(ALL_db_restricted_sel_cond_ind$pos37),]
        rs_vector<-ALL_db_restricted_sel_cond_ind$rs
        
        
        if(DEBUG == 1)
        {
          cat("ALL_db_restricted_sel_cond_ind_2\n")
          cat(str(ALL_db_restricted_sel_cond_ind))
          cat("\n")
          cat(str(unique(ALL_db_restricted_sel_cond_ind$VAR)))
          cat("\n")
          cat("rs_vector\n")
          cat(sprintf(as.character(rs_vector)))
          cat("\n")
        }
        
        ALL_db_restricted_sel_cond_ind$DUMMY<-factor(ALL_db_restricted_sel_cond_ind$rs,
                                                     levels=rs_vector,
                                                     ordered=T)
        
        if(DEBUG == 1)
        {
          cat("ALL_db_restricted_sel_cond_ind_3\n")
          cat(str(ALL_db_restricted_sel_cond_ind))
          cat("\n")
          cat("DUMMY\n")
          cat(sprintf(as.character(names(summary(ALL_db_restricted_sel_cond_ind$DUMMY)))))
          cat("\n")
          cat(sprintf(as.character(summary(ALL_db_restricted_sel_cond_ind$DUMMY))))
          cat("\n")
        }
        
       Result_AS<-rbind(Result_AS,ALL_db_restricted_sel_cond_ind) 
        
       # quit(status = 1)
      }# dim(ALL_db_restricted_sel_cond_ind)[1] >0
    }#k in 1:length(AS_ID_array_sel)
    
    
    Result_FINAL<-rbind(Result_AS,Result_FINAL) 
    
    if(DEBUG == 1)
    {
      cat("Result_FINAL_0\n")
      cat(str(Result_FINAL))
      cat("\n")
    }
    
  }# i in 1:length(VAR_vector)
  
  
  if(dim(Result_FINAL)[1] >0)
  {
    cat("Result_FINAL_1\n")
    cat(str(Result_FINAL))
    cat("\n")
    cat(str(unique(Result_FINAL$VAR)))
    cat("\n")
    
    setwd(path_files)
    
    saveRDS(Result_FINAL,file='AS_selector.rds')
  }
 
  
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
    make_option(c("--Table_S1"), type="character", default=NULL, 
                metavar="type", 
                help="Path to tab-separated input file listing regions to analyze. Required."),
    make_option(c("--ALL_db"), type="character", default=NULL, 
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