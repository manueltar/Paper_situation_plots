#!/bin/bash
 
set -e

eval "$(conda shell.bash hook)"

Rscripts_path=$(echo "/home/manuel.tardaguila/Scripts/R/")

output_dir=$1

mkdir -p $output_dir

Log_files=$(echo "$output_dir""/""Log_files/")


mkdir -p $Log_files


#### AS_selector_and_beta_printer ####

module load R/4.1.0

type=$(echo "AS_selector_and_beta_printer")
outfile_AS_selector_and_beta_printer=$(echo "$Log_files""outfile_1_""$type"".log")
touch $outfile_AS_selector_and_beta_printer
echo -n "" > $outfile_AS_selector_and_beta_printer
name_AS_selector_and_beta_printer=$(echo "$type""_job")


Rscript_AS_selector_and_beta_printer=$(echo "$Rscripts_path""95_AS_selector_and_beta_se_printer.R")

Table_S1=$(echo "/group/soranzo/manuel.tardaguila/Paper_bits/FIX_TABLES/Provisional_Tables/Table_S1_Provisional.rds")
ALL_db=$(echo "/group/soranzo/manuel.tardaguila/Rebuild_sit_plots/Reference_files/ALL_db.tsv")



myjobid_AS_selector_and_beta_printer=$(sbatch --job-name=$name_AS_selector_and_beta_printer --output=$outfile_AS_selector_and_beta_printer --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=3 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_AS_selector_and_beta_printer  --ALL_db $ALL_db --Table_S1 $Table_S1 --type $type --out $output_dir")
myjobid_seff_AS_selector_and_beta_printer=$(sbatch --dependency=afterany:$myjobid_AS_selector_and_beta_printer --open-mode=append --output=$outfile_AS_selector_and_beta_printer --job-name=$(echo "seff""_""$name_AS_selector_and_beta_printer") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_AS_selector_and_beta_printer >> $outfile_AS_selector_and_beta_printer")

# #### Intersect_VAR_PCHiC ####

module load R/4.1.0


type=$(echo "Intersect_VAR_PCHiC")
outfile_Intersect_VAR_PCHiC=$(echo "$Log_files""outfile_2_""$type"".log")
touch $outfile_Intersect_VAR_PCHiC
echo -n "" > $outfile_Intersect_VAR_PCHiC
name_Intersect_VAR_PCHiC=$(echo "$type""_job")


Rscript_Intersect_VAR_PCHiC=$(echo "$Rscripts_path""96_Intersect_original_PCHiC_ALL_v2.R")


GWAS_blocks=$(echo "/group/soranzo/manuel.tardaguila/Rebuild_sit_plots/Reference_files/Vuckovic_Blocks.csv.gz")
variants_file=$(echo "$output_dir""/""Build_files""/""AS_selector.rds")
PCHiC_original_file=$(echo "/group/soranzo/manuel.tardaguila/Rebuild_sit_plots/Reference_files/PCHiC_peak_matrix_cutoff5.tsv")
ensembl_gtf=$(echo "/group/soranzo/manuel.tardaguila/Rebuild_sit_plots/Reference_files/Homo_sapiens.GRCh37.87.gtf.gz")
genes_without_symbol_CORRECTED=$(echo "/group/soranzo/manuel.tardaguila/Rebuild_sit_plots/Reference_files/genes_without_symbol_CORRECTED.tsv")
ensembl_gtf_sel_CORRECTED=$(echo "/group/soranzo/manuel.tardaguila/Rebuild_sit_plots/Reference_files/ensembl_gtf_sel_CORRECTED.tsv")


#  --dependency=afterany:$myjobid_AS_selector_and_beta_printer

myjobid_Intersect_VAR_PCHiC=$(sbatch --dependency=afterany:$myjobid_AS_selector_and_beta_printer --job-name=$name_Intersect_VAR_PCHiC --output=$outfile_Intersect_VAR_PCHiC --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_Intersect_VAR_PCHiC --PCHiC_original_file $PCHiC_original_file --variants_file $variants_file --GWAS_blocks $GWAS_blocks --ensembl_gtf $ensembl_gtf --genes_without_symbol_CORRECTED $genes_without_symbol_CORRECTED --ensembl_gtf_sel_CORRECTED $ensembl_gtf_sel_CORRECTED --type $type --out $output_dir")
myjobid_seff_Intersect_VAR_PCHiC=$(sbatch --dependency=afterany:$myjobid_Intersect_VAR_PCHiC --open-mode=append --output=$outfile_Intersect_VAR_PCHiC --job-name=$(echo "seff""_""$name_Intersect_VAR_PCHiC") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Intersect_VAR_PCHiC >> $outfile_Intersect_VAR_PCHiC")





#### Intersect_Genes ####

module load R/4.1.0


type=$(echo "Intersect_Genes")
outfile_Intersect_Genes=$(echo "$Log_files""outfile_3_""$type"".log")
touch $outfile_Intersect_Genes
echo -n "" > $outfile_Intersect_Genes
name_Intersect_Genes=$(echo "$type""_job")


Rscript_Intersect_Genes=$(echo "$Rscripts_path""98_Select_GWAS_block_borders_and_genes_v2.R")



variants_file=$(echo "$output_dir""/""Build_files""/""AS_selector.rds")
PCHiC_HITS=$(echo "$output_dir""/""Build_files""/""PCHiC_AS_POST_ENSEMBL.rds")
BLOCK_PRE=$(echo "$output_dir""/""Build_files""/""Boundaries_AS_PRE_genes.rds")
ensembl_gtf=$(echo "/group/soranzo/manuel.tardaguila/Rebuild_sit_plots/Reference_files/Homo_sapiens.GRCh37.87.gtf.gz")



myjobid_Intersect_Genes=$(sbatch --dependency=afterany:$myjobid_Intersect_VAR_PCHiC  --job-name=$name_Intersect_Genes --output=$outfile_Intersect_Genes --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=3 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_Intersect_Genes --ensembl_gtf $ensembl_gtf --PCHiC_HITS $PCHiC_HITS --BLOCK_PRE $BLOCK_PRE --variants_file $variants_file --type $type --out $output_dir")
myjobid_seff_Intersect_Genes=$(sbatch --dependency=afterany:$myjobid_Intersect_Genes --open-mode=append --output=$outfile_Intersect_Genes --job-name=$(echo "seff""_""$name_Intersect_Genes") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Intersect_Genes >> $outfile_Intersect_Genes")



#### Intersect_ATAC ####

module load R/4.1.0


type=$(echo "Intersect_ATAC")
outfile_Intersect_ATAC=$(echo "$Log_files""outfile_4_""$type"".log")
touch $outfile_Intersect_ATAC
echo -n "" > $outfile_Intersect_ATAC
name_Intersect_ATAC=$(echo "$type""_job")


Rscript_Intersect_ATAC=$(echo "$Rscripts_path""100_Intersect_ATAC_v2.R")


ATAC_peaks=$(echo "/group/soranzo/manuel.tardaguila/Rebuild_sit_plots/Reference_files/29August2017_EJCsamples_allReads_500bp.bed")
ATAC_counts=$(echo "/group/soranzo/manuel.tardaguila/Rebuild_sit_plots/Reference_files/29August2017_EJCsamples_allReads_500bp.counts.txt")
BLOCK_POST=$(echo "$output_dir""/""Build_files""/""Boundaries_AS_POST_genes.rds")


myjobid_Intersect_ATAC=$(sbatch  --dependency=afterany:$myjobid_Intersect_Genes --job-name=$name_Intersect_ATAC --output=$outfile_Intersect_ATAC --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=3 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_Intersect_ATAC --ATAC_peaks $ATAC_peaks --ATAC_counts $ATAC_counts --BLOCK_POST $BLOCK_POST  --type $type --out $output_dir")
myjobid_seff_Intersect_ATAC=$(sbatch --dependency=afterany:$myjobid_Intersect_ATAC --open-mode=append --output=$outfile_Intersect_ATAC --job-name=$(echo "seff""_""$name_Intersect_ATAC") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Intersect_ATAC >> $outfile_Intersect_ATAC")





# #### Printing_selected_plots ####

# module load R/4.1.0


# type=$(echo "Printing_selected_plots")
# outfile_Printing_selected_plots=$(echo "$Log_files""outfile_5_""$type"".log")
# touch $outfile_Printing_selected_plots
# echo -n "" > $outfile_Printing_selected_plots
# name_Printing_selected_plots=$(echo "$type""_job")


# Rscript_Printing_selected_plots=$(echo "$Rscripts_path""102_Situtation_plot_printer_v3_SELECTIONS.R")


# Table_S7=$(echo "/group/soranzo/manuel.tardaguila/Paper_bits/FIX_TABLES/Provisional_Tables/Table_S7_Provisional.rds")
# variants_file=$(echo "$output_dir""/""Build_files""/""AS_selector.rds")
# PCHiC_HITS=$(echo "$output_dir""/""Build_files""/""PCHiC_AS_POST_ENSEMBL.rds")
# BLOCK_POST=$(echo "$output_dir""/""Build_files""/""Boundaries_AS_POST_genes_POST_ATAC.rds")
# ATAC_tracks=$(echo "$output_dir""/""Build_files""/""Intersected_ATAC.rds")
# gene_tracks=$(echo "$output_dir""/""Build_files""/""Intersected_genes.rds")
# dataframe_to_print_SIT_plots=$(echo "/group/soranzo/manuel.tardaguila/Rebuild_sit_plots/dataframe_to_print_SIT_plots_3.txt")
# rescue_haplotypes=$(echo "/group/soranzo/manuel.tardaguila/INTERVAL_results/chr9_135920196_C_T_ALL_BY_ALL_DE_LM_FPKM_results_Haplotypes.tsv")

# #--dependency=afterany:$myjobid_Intersect_ATAC

# myjobid_Printing_selected_plots=$(sbatch --job-name=$name_Printing_selected_plots --output=$outfile_Printing_selected_plots --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=1024M --parsable --wrap="Rscript $Rscript_Printing_selected_plots --ATAC_tracks $ATAC_tracks --PCHiC_HITS $PCHiC_HITS --BLOCK_POST $BLOCK_POST --variants_file $variants_file --gene_tracks $gene_tracks --Table_S7 $Table_S7 --dataframe_to_print_SIT_plots $dataframe_to_print_SIT_plots --rescue_haplotypes $rescue_haplotypes --type $type --out $output_dir")

# myjobid_seff_Printing_selected_plots=$(sbatch --dependency=afterany:$myjobid_Printing_selected_plots --open-mode=append --output=$outfile_Printing_selected_plots --job-name=$(echo "seff""_""$name_Printing_selected_plots") --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_Printing_selected_plots >> $outfile_Printing_selected_plots")




# ####################################################################################################################################################

# #### ATU_plots ############################# ATU_plots


# type=$(echo "ATU_plots")
# outfile_ATU_plots=$(echo "$Log_files""outfile_6_""$type"".log")
# touch $outfile_ATU_plots
# echo -n "" > $outfile_ATU_plots
# name_ATU_plots=$(echo "$type""_job")
# seff_name=$(echo "seff""_""$type")

# Rscript_ATU_plots=$(echo "$Rscripts_path""103_ATU_plots_v2.R")

# Table_S7=$(echo "/group/soranzo/manuel.tardaguila/Paper_bits/FIX_TABLES/Provisional_Tables/Table_S7_Provisional.rds")
# Table_S6=$(echo "/group/soranzo/manuel.tardaguila/Paper_bits/FIX_TABLES/Provisional_Tables/Table_S6_Provisional.rds")
# ensembl_gtf=$(echo "/group/soranzo/manuel.tardaguila/Rebuild_sit_plots/Reference_files/Homo_sapiens.GRCh37.87.gtf.gz")
# Transposed_Isoform_Expression=$(echo "/group/soranzo/manuel.tardaguila/INTERVAL_results/Transposed_Isoform_Expression_df.rds")
# path_INTERVAL=$(echo "/group/soranzo/manuel.tardaguila/INTERVAL_results/")
# path_BP=$(echo "/group/soranzo/manuel.tardaguila/BP_results/")


# myjobid_ATU_plots=$(sbatch --job-name=$name_ATU_plots --output=$outfile_ATU_plots --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_ATU_plots --Table_S7 $Table_S7 --Table_S6 $Table_S6 --Transposed_Isoform_Expression $Transposed_Isoform_Expression --ensembl_gtf $ensembl_gtf --type $type --out $output_dir --path_INTERVAL $path_INTERVAL --path_BP $path_BP")
# myjobid_seff_ATU_plots=$(sbatch --dependency=afterany:$myjobid_ATU_plots --open-mode=append --output=$outfile_ATU_plots --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_ATU_plots >> $outfile_ATU_plots")


# ####################################################################################################################################################

# #### ATU_plots_haplotypes ############################# ATU_plots_haplotypes


# type=$(echo "8_ATU_plots_haplotypes")
# outfile_ATU_plots_haplotypes=$(echo "$Log_files""outfile_7_""$type"".log")
# touch $outfile_ATU_plots_haplotypes
# echo -n "" > $outfile_ATU_plots_haplotypes
# name_ATU_plots_haplotypes=$(echo "$type""_job")
# seff_name=$(echo "seff""_""$type")

# Rscript_ATU_plots_haplotypes=$(echo "$Rscripts_path""117_ATU_plots_haplotypes_v2.R")


# ATU_haplotype_variants=$(echo 'chr2_219020958_C_T/Haplotypes/chr2_219020958_C_T__chr2_219214529_A_T/INTERVAL_covariates_and_PEER_factors_Haplotype_chr2_219020958_C_T__chr2_219214529_A_T.rds,chr9_135920196_C_T/Haplotypes/chr9_135920196_C_T__chr9_135864513_C_T/INTERVAL_covariates_and_PEER_factors_Haplotype_chr9_135920196_C_T__chr9_135864513_C_T.rds,chr22_50949811_T_C/Haplotypes/chr22_50949811_T_C__chr22_50967912_C_T/INTERVAL_covariates_and_PEER_factors_Haplotype_chr22_50949811_T_C__chr22_50967912_C_T.rds')
# ENSG_selected=$(echo 'ENSG00000018280,ENSG00000165702,ENSG00000025708')
# transcripts_selected=$(echo 'ENST00000490872__ENST00000233202__ENST00000354352,ENST00000372123__ENST00000372122__ENST00000339463,ENST00000395678__ENST00000487577__ENST00000252029')


# ensembl_gtf=$(echo "/group/soranzo/manuel.tardaguila/Rebuild_sit_plots/Reference_files/Homo_sapiens.GRCh37.87.gtf.gz")
# Transposed_Isoform_Expression=$(echo "/group/soranzo/manuel.tardaguila/INTERVAL_results/Transposed_Isoform_Expression_df.rds")
# path_INTERVAL=$(echo "/group/soranzo/manuel.tardaguila/INTERVAL_results/")


# myjobid_ATU_plots_haplotypes=$(sbatch --job-name=$name_ATU_plots_haplotypes --output=$outfile_ATU_plots_haplotypes --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=4 --mem-per-cpu=4096M --parsable --wrap="Rscript $Rscript_ATU_plots_haplotypes --ENSG_selected $ENSG_selected --ATU_haplotype_variants $ATU_haplotype_variants --Transposed_Isoform_Expression $Transposed_Isoform_Expression --ensembl_gtf $ensembl_gtf --type $type --transcripts_selected $transcripts_selected --out $output_dir --path_INTERVAL $path_INTERVAL")
# myjobid_seff_ATU_plots_haplotypes=$(sbatch --dependency=afterany:$myjobid_ATU_plots_haplotypes --open-mode=append --output=$outfile_ATU_plots_haplotypes --job-name=$seff_name --partition=cpuq --time=24:00:00 --nodes=1 --ntasks-per-node=1 --mem-per-cpu=128M --parsable --wrap="seff $myjobid_ATU_plots_haplotypes >> $outfile_ATU_plots_haplotypes")
