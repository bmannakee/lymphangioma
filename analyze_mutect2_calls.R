load("all_vars_filter.RData")
load("filter_fr.RData")
library(tidyverse)
library(cowplot)


# Plot the variant allele frequency for passed variants
pl1 <- ggplot(all_vars_filter %>% dplyr::filter(all_filters=="PASS"),aes(tumor_f)) +
  geom_histogram(binwidth = .01,color='black',fill="white") + 
  scale_x_continuous(name="Variant Allele Frequency",breaks=c(0,0.05,0.1,0.15,0.2,0.25,0.3,0.5),limits=c(0,.5)) +
  labs(title="Tumor Content",
       subtitle="All variants called by MuTect including non-coding",
       y="Count") +
  theme(axis.text.x=element_text(angle=45)) + 
  theme_classic(base_size=15)

#filter_fr %>% dplyr::filter(tlod_only=="PASS" & alt_depth > 2 & cgc_somatic=="yes") %>% print(n=nrow(.))
#filter_fr %>% dplyr::filter((all_filters=="PASS" | tlod_only=="PASS") & cgc_somatic=="yes" & alt_depth > 2) %>% print(n=nrow(.))
filter_fr %>% dplyr::filter((all_filters=="PASS" | tlod_only=="PASS") & 
                              (cgc_somatic=="yes" | !is.na(cosmic_overlapping) | !is.na(mutsig_results)) & 
                              alt_depth > 1 &
                              tumor_f > 0.06) %>%
  replace_na(list(cgc_somatic="no",
                  cgc_syndrome="",
                  mutsig_results="")) %>%
  rename("Gene_Name"=Hugo_Symbol,
         "VAF"=tumor_f,
         "Evidence_level"=t_lod_fstar__INFO__) -> final_fr
