library(tidyverse)


load_oncotated <- function(onc_file){
  fr <- read_tsv(onc_file,comment="#",col_types = cols(.default='c'))
  fr %>% dplyr::filter(Tumor_Sample_Barcode=="TUMOR")
}
onc_files <- list.files('./data/tumor_variants',pattern="*.oncotated.txt",full.names=T)
fr <- lapply(onc_files,load_oncotated)
fr <- do.call(rbind,fr) 
fr <- fr %>% separate(allelic_depth,into=c('ref_depth','alt_depth'),sep=',',remove=F) %>%
  mutate_at(c('tumor_f','ref_depth','alt_depth','t_lod_fstar__INFO__'),funs(as.numeric))
fr <- fr %>% rename("cgc_syndrome"=`CGC_Cancer Syndrome`,
                    "cgc_somatic"=`CGC_Cancer Somatic Mut`,
                    "mutsig_results"=`MUTSIG_Published_Results`,
                    "exac_path"=`ExAC_clinvar_pathogenic`,
                    "cosmic_overlapping"=`COSMIC_overlapping_mutations`)



# spike at 12.5%, I suspect this is the cellularity.



# everthing but RNA
all_vars_filter <- fr %>% dplyr::filter(!stringr::str_detect(Variant_Classification,'RNA')) %>%
  mutate(all_filters=ifelse(germline_risk=="PASS" & t_lod_fstar=="PASS" & homologous_mapping_event=="PASS" &
                              alt_allele_in_normal=="PASS" & triallelic_site=="PASS" & clustered_events=="PASS" &
                              multi_event_alt_allele_in_normal=="PASS" & str_contraction=="PASS","PASS","FAIL")) %>%
  mutate(tlod_only=ifelse(germline_risk=="PASS" &  t_lod_fstar=="FAIL" & homologous_mapping_event=="PASS" &
                            alt_allele_in_normal=="PASS" & triallelic_site=="PASS" & clustered_events=="PASS" &
                            multi_event_alt_allele_in_normal=="PASS" & str_contraction=="PASS","PASS","FAIL")) %>%
  dplyr::select(Hugo_Symbol,Variant_Classification,Chromosome,Start_position,Codon_Change,Protein_Change,all_filters,tlod_only,tumor_f,
                alt_depth,ref_depth,germline_risk,t_lod_fstar__INFO__,
                cgc_syndrome,cgc_somatic,mutsig_results,cosmic_overlapping)
save(all_vars_filter,file="./all_vars_filter.RData")
# Coding variant only
# Oncotator breaks out the FILTER field in MuTect2 VCF output into separate columns for whether a sample PASS or FAIL each filter.
filter_fr <- fr %>% dplyr::filter(stringr::str_detect(Variant_Classification,'Missense|Nonsense|Ins|Del')) %>%
  mutate(all_filters=ifelse(germline_risk=="PASS" & t_lod_fstar=="PASS" & homologous_mapping_event=="PASS" &
                              alt_allele_in_normal=="PASS" & triallelic_site=="PASS" & clustered_events=="PASS" &
                              multi_event_alt_allele_in_normal=="PASS" & str_contraction=="PASS","PASS","FAIL")) %>%
  mutate(tlod_only=ifelse(germline_risk=="PASS" &  t_lod_fstar=="FAIL" & homologous_mapping_event=="PASS" &
                              alt_allele_in_normal=="PASS" & triallelic_site=="PASS" & clustered_events=="PASS" &
                              multi_event_alt_allele_in_normal=="PASS" & str_contraction=="PASS","PASS","FAIL")) %>%
  dplyr::select(Hugo_Symbol,Variant_Classification,Chromosome,Start_position,Codon_Change,Protein_Change,all_filters,tlod_only,tumor_f,
                alt_depth,ref_depth,germline_risk,t_lod_fstar__INFO__,
                cgc_syndrome,cgc_somatic,mutsig_results,cosmic_overlapping)
save(filter_fr,file="./filter_fr.RData")