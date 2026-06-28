install.packages("cbioportalR")
library(cbioportalR)
library(dplyr)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))
library(purrr)
library(ggplot2)
library(readxl)
BiocManager::install("maftools")
library(maftools)
library(ComplexHeatmap)
install.packages("grid")
library(grid)
library(tidyr)

# query
set_cbioportal_db("public")

all_studies <- available_studies()

get_mutations_by_study(
  study_id = selected_studies,
  molecular_profile_id = NULL,
  add_hugo = TRUE,
  base_url = NULL
)

# 321 exome studies
selected_studies <- c('all_stjude_2015','acyc_mgh_2016','ccrcc_dfci_2019','chol_jhu_2013','chol_nccs_2013','chol_nus_2012','cllsll_icgc_2011','coad_caseccc_2015',
                      'ctcl_columbia_2015','hnsc_tcga_pub','ihch_ismms_2015','metastatic_solid_tumors_mich_2017','pediatric_dkfz_2017','brain_cptac_gdc','breast_cptac_gdc',
                      'coad_cptac_gdc','luad_cptac_gdc','lusc_cptac_gdc','ohnca_cptac_gdc','ovary_cptac_gdc','pancreas_cptac_gdc','kirp_tcga','blca_tcga_pub',
                      'acbc_mskcc_2015','blca_tcga_pan_can_atlas_2018','brca_mbcproject_2022','ccrcc_irc_2014','ccrcc_utokyo_2013','coadread_genentech','cellline_nci60',
                      'cll_iuopa_2015','coadread_dfci_2016','cll_broad_2015','brca_tcga_pan_can_atlas_2018','cesc_tcga_pan_can_atlas_2018','chol_tcga_pan_can_atlas_2018',
                      'ccle_broad_2019','coad_cptac_2019','coadread_cass_2020','cll_broad_2022','coadread_tcga','coadread_tcga_pub','desm_broad_2015','dlbc_broad_2012',
                      'cscc_hgsc_bcm_2014','coadread_tcga_pan_can_atlas_2018','dlbc_tcga_pan_can_atlas_2018','difg_glass_2019','cscc_ucsf_2021','cscc_ranson_2022',
                      'difg_glass','esca_tcga','es_dfarber_broad_2014','esca_tcga_pan_can_atlas_2018','gbm_tcga','gbm_tcga_pan_can_atlas_2018','kirc_tcga','kirc_tcga_pub',
                      'laml_tcga','laml_tcga_pub','meso_tcga','luad_tcga_pub','mbl_sickkids_2016','mixed_pipseq_2017','mds_mskcc_2020','paad_tcga','nhl_bcgsc_2011',
                      'nhl_bcgsc_2013','pcpg_tcga','prad_broad_2013','prad_tcga_pub','pcpg_tcga_pub','pptc_2019','pog570_bcgsc_2020','stad_tcga_pub','thca_tcga_pub',
                      'hnsc_tcga_pan_can_atlas_2018','kich_tcga_pan_can_atlas_2018','kirc_tcga_pan_can_atlas_2018','kirp_tcga_pan_can_atlas_2018',
                      'laml_tcga_pan_can_atlas_2018','lihc_tcga_pan_can_atlas_2018','luad_tcga_pan_can_atlas_2018','lusc_tcga_pan_can_atlas_2018',
                      'meso_tcga_pan_can_atlas_2018','ov_tcga_pan_can_atlas_2018','paad_tcga_pan_can_atlas_2018','pcpg_tcga_pan_can_atlas_2018',
                      'prad_tcga_pan_can_atlas_2018','sarc_tcga_pan_can_atlas_2018','skcm_tcga_pan_can_atlas_2018','stad_tcga_pan_can_atlas_2018',
                      'tgct_tcga_pan_can_atlas_2018','thca_tcga_pan_can_atlas_2018','thym_tcga_pan_can_atlas_2018','ucec_tcga_pan_can_atlas_2018',
                      'ucs_tcga_pan_can_atlas_2018','uvm_tcga_pan_can_atlas_2018','coad_silu_2022','acc_tcga_pan_can_atlas_2018','pancan_mappyacts_2022',
                      'thyroid_gatci_2024','acc_tcga_gdc','acc_tcga','blca_tcga','ampca_bcm_2016','blca_dfarber_mskcc_2014','blca_bgi','all_stjude_2013','acyc_mskcc_2013',
                      'acyc_jhu_2016','acyc_sanger_2013','all_stjude_2016','angs_project_painter_2018','all_phase2_target_2018_pub','aml_target_2018_pub','bfn_duke_nus_2015',
                      'blca_cornell_2016','aml_ohsu_2018','angs_painter_2020','aml_ohsu_2022','alal_target_gdc','aml_target_gdc','bll_target_gdc','nbl_target_gdc',
                      'os_target_gdc','brca_broad','brca_bccrc','brca_igr_2015','brca_mbcproject_wagle_2017','blca_tcga_pub_2017','brca_jup_msk_2020','brain_cptac_2020',
                      'brca_cptac_2020','brca_hta9_htan_2022','brca_dfci_2020','brca_tcga','brca_sanger','brca_tcga_pub2015','brca_tcga_pub','brca_smc_2018',
                      'brca_pareja_msk_2020','pancan_pcawg_2020','chl_sccc_2023','dlbc_tcga','esca_broad','escc_icgc','es_iocurie_2014','gbc_shanghai_2014','egc_tmucih_2015',
                      'dlbcl_duke_2017','dlbcl_dfci_2018','gbm_columbia_2019','egc_trap_ccr_msk_2023','gbm_tcga_pub2013','hnsc_broad','hcc_inserm_fr_2015',
                      'gbm_mayo_pdx_sarkaria_2019','hccihch_pku_2019','gbm_cptac_2021','hcc_meric_2021','hcc_clca_2024','kich_tcga','lgg_tcga','hnsc_tcga','kich_tcga_pub',
                      'hnsc_jhu','hnsc_mdanderson_2013','ihch_smmu_2014','lihc_tcga','lgg_ucsf_2014','lgggbm_tcga_pub','lihc_amc_prv','lihc_riken','luad_mskcc_2015',
                      'luad_broad','liad_inserm_fr_2014','lcll_broad_2013','luad_cptac_2020','luad_tcga','lusc_tcga','lusc_tcga_pub','mbl_broad_2012','mbl_icgc',
                      'mbl_pcgp','luad_oncosg_2020','lung_smc_2016','mbl_dkfz_2017','lusc_cptac_2021','lung_nci_2022','mm_broad','mcl_idibips_2013','mds_tokyo_2011',
                      'mel_tsam_liang_2017','mel_ucla_2016','mixed_allen_2018','mel_dfci_2019','mbn_sfu_2023','npc_nusingapore','nepc_wcm_2016','nbl_ucologne_2015',
                      'nbl_broad_2013','mrt_bcgsc_2016','nbl_target_2018_pub','mpn_cimr_2013','nsclc_mskcc_2015','nsclc_mskcc_2018','stad_oncosg_2018',
                      'mng_utoronto_2021','mpcproject_broad_2021','ov_tcga','mpnst_mskcc','nbl_amc_2012','nccrcc_genentech_2014','ov_tcga_pub','paac_jhu_2014',
                      'paad_icgc','paad_utsw_2015','nsclc_tcga_broad_2016','paad_qcmg_uq_2016','pact_jhu_2011','nsclc_tracerx_2017','paad_cptac_2021','nst_nfosi_ntap',
                      'panet_jhu_2011','pcnsl_mayo_2015','prad_broad','crc_hta11_htan_2021','panet_shanghai_2013','plmeso_nyu_2015','prad_cpcg_2017','panet_arcnet_2017',
                      'past_dkfz_heidelberg_2013','prad_eururol_2017','prad_tcga','prad_fhcrc','prad_mich','prad_mskcc_2014','prad_su2c_2015','prad_p1000',
                      'prad_su2c_2019','prostate_dkfz_2018','prad_mskcc_cheny1_organoids_2014','prostate_pcbm_swiss_2019','sarc_tcga','sclc_clcgp','sclc_jhu','skcm_broad',
                      'rms_nih_2014','sarc_tcga_pub','sclc_cancercell_gardner_2017','rt_target_2018_pub','sclc_ucologne_2015','skcm_tcga','stad_tcga','skcm_broad_dfarber',
                      'skcm_yale','stad_pfizer_uhongkong','skcm_broad_brafresist_2012','skcm_mskcc_2014','skcm_tcga_pub_2015','skcm_dfci_2015','tgct_tcga','thym_tcga',
                      'thca_tcga','stad_uhongkong','stad_utokyo','tet_nci_2014','stes_tcga_pub','stmyec_wcm_2022','ucec_tcga','ucs_tcga','ucs_jhu_2014','ucec_tcga_pub',
                      'um_qimr_2016','uccc_nih_2017','ucec_cptac_2020','uvm_tcga','wt_target_2018_pub','vsc_cuk_2018','utuc_cornell_baylor_mdacc_2019','utuc_igbmc_2021',
                      'rcc_cptac_gdc','lgg_tcga_pan_can_atlas_2018','uec_cptac_gdc','wt_target_gdc','brca_aurora_2023','pancan_pdmr_2025','schw_ctf_synodos_2025',
                      'blca_tcga_gdc','brca_tcga_gdc','cesc_tcga_gdc','chol_tcga_gdc','chrcc_tcga_gdc','ccrcc_tcga_gdc','aml_tcga_gdc','dlbclnos_tcga_gdc','esca_tcga_gdc',
                      'difg_tcga_gdc','coad_tcga_gdc','gbm_tcga_gdc','hnsc_tcga_gdc','hcc_tcga_gdc','luad_tcga_gdc','lusc_tcga_gdc','hgsoc_tcga_gdc','plmeso_tcga_gdc',
                      'paad_tcga_gdc','mnet_tcga_gdc','prad_tcga_gdc','nsgct_tcga_gdc','prcc_tcga_gdc','read_tcga_gdc','soft_tissue_tcga_gdc','skcm_tcga_gdc',
                      'stad_tcga_gdc','thpa_tcga_gdc','thym_tcga_gdc','ucec_tcga_gdc','ucs_tcga_gdc','um_tcga_gdc'
)

# ARHG genes (52)
arhg_pattern <- "^ARHG"
arhg_genes <- c(
  "ARHGEF1", "ARHGEF2", "ARHGEF3", "ARHGEF4", "ARHGEF5",
  "ARHGEF6", "ARHGEF7", "ARHGEF9", "ARHGEF10", "ARHGEF11",
  "ARHGEF12", "ARHGEF15", "ARHGEF16", "ARHGEF17", "ARHGEF18",
  "ARHGEF19", "ARHGEF25", "ARHGEF26", "ARHGEF28",
  "ARHGAP1", "ARHGAP4", "ARHGAP5", "ARHGAP8", "ARHGAP9",
  "ARHGAP10", "ARHGAP11A", "ARHGAP11B", "ARHGAP12", "ARHGAP15",
  "ARHGAP17", "ARHGAP18", "ARHGAP19", "ARHGAP20", "ARHGAP21",
  "ARHGAP22", "ARHGAP23", "ARHGAP24", "ARHGAP25", "ARHGAP26",
  "ARHGAP27", "ARHGAP28", "ARHGAP29", "ARHGAP30", "ARHGAP31",
  "ARHGAP32", "ARHGAP33", "ARHGAP35", "ARHGAP36", "ARHGAP39",
  "ARHGAP42", "ARHGAP44", "ARHGAP45"
)

get_arhg_mutations <- function(study_id, genes) {
  profiles <- get_molecular_profiles(cbio, studyId = study_id)
  mut_profile <- profiles %>% filter(grepl("mutations", molecularProfileId)) %>% pull(molecularProfileId)
  
  if (length(mut_profile) == 0) return(tibble())
  
  safe_query <- purrr::safely(get_mutation_data)
  result <- safe_query(cbio,
                       molecularProfileId = mut_profile,
                       sampleListId = past0(study, "_all"),
                       genes = genes)
  if (!is.null(result$result)) {
    result$result %>% mutate(studyId = study_id)
  } else {
    tibble()
  }
}

get_arhg_mutations <- function(study_id) {
  df <- tryCatch({
    get_mutations_by_study(study_id = study_id)
  }, error = function(e) {
    message("Error for study ", study_id, ": ", e$message)
    return(NULL)
  })
  
  if (!is.null(df)) {
    df %>%
      filter(grepl(arhg_pattern, hugoGeneSymbol, ignore.case = TRUE)) %>%
      mutate(studyId = study_id)
  } else {
    tibble()
  }
}
 

all_results <- map_dfr(selected_studies, get_arhg_mutations)

write.csv(all_results, file.path(BASE_DIR, "cBioPortal/Datasets/arhg_mutations_across_studies.csv"), row.names = FALSE)

all_results <- read_excel(file.path(BASE_DIR, "cBioPortal/Datasets/arhg_mutations_across_studies.xlsx"))

# bargraph for mutations across cancer types
mutations_by_cancer_type <- all_results %>%
  group_by(CancerType) %>%
  summarise(n_mutations = n())

ggplot(mutations_by_cancer_type, aes(x = reorder(CancerType, n_mutations), y = n_mutations)) + 
  geom_col(fill = "steelblue") +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "ARHG Mutations by Cancer Type",
       x = "Cancer Type",
       y = "Mutation Count")

# bargraph for myeloid

# Subset to myeloid samples
myeloid_data <- all_results %>% 
  filter(CancerType == "myeloid")

# Count mutations per gene
mutations_by_gene <- myeloid_data %>%
  group_by(hugoGeneSymbol) %>%
  summarise(n_mutations = n()) %>%
  arrange(desc(n_mutations))

# Optionally keep only top 20 genes
top_genes1 <- mutations_by_gene %>% slice_max(n_mutations, n = 20)

# Make bar graph
ggplot(top_genes1, aes(x = reorder(hugoGeneSymbol, n_mutations), y = n_mutations)) +
  geom_col(fill = "steelblue") +
  coord_flip() +  # horizontal bars
  labs(title = "Top ARHG Gene Mutations in Myeloid",
       x = "Gene",
       y = "Mutation Count") +
  theme_minimal()


# by study
# mutations_by_study <- all_results %>%
#   group_by(studyId) %>%
#   summarise(n_mutations = n())
# 
# ggplot(mutations_by_study, aes(x = reorder(studyId, n_mutations), y = n_mutations)) + 
#   geom_col(fill = "steelblue") +
#   coord_flip() +
#   labs(title = "NUmber of ARHG Mutations by Study",
#        x = "Study",
#        y = "Mutation Count")

# top 20 most mutated ARHG
top_genes <- all_results %>%
  group_by(hugoGeneSymbol) %>%
  summarise(n_mutations = n()) %>%
  arrange(desc(n_mutations)) %>%
  slice_max(n_mutations, n = 20)

ggplot(top_genes, aes(x = reorder(hugoGeneSymbol, n_mutations), y = n_mutations)) +
  geom_col(fill = "tomato") +
  coord_flip() +
  scale_y_continuous(expand = c(0, 0)) +
  labs(title = "Top 20 mutated ARHG genes",
       x = "Gene",
       y = "Mutation count")

# checking for duplicates
duplicated_patients <- all_results$patientId[duplicated(all_results$patientId)]
length(unique(duplicated_patients))

# oncoprint
# maf_ready <- all_results %>%
#   transmute(
#     Hugo_Symbol = hugoGeneSymbol,   
#     Tumor_Sample_Barcode = sampleId,  
#     Variant_Classification = mutationType, 
#     Chromosome = chr,           
#     Start_Position = startPosition,    
#     End_Position = endPosition,        
#     Reference_Allele = referenceAllele, 
#     Tumor_Seq_Allele2 = variantAllele,
#     Variant_Type = variantType
#   )
# 
# maf <- read.maf(maf_ready)
# oncoplot(maf, top = 10)

keep <- all_results %>%
  distinct(studyId, sampleId, hugoGeneSymbol) %>%
  count(studyId, sampleId, name = "mut_n") %>%
  group_by(CancerType) %>%
  slice_max(mut_n, n = 200, with_ties = FALSE) %>%
  ungroup()

sub <- all_results %>%
  semi_join(keep, by = c("studyId","sampleId")) %>%
  mutate(present = 1) %>%
  distinct(sampleId, gene, present) %>%
  pivot_wider(names_from = sampleId, values_from = present, values_fill = 0)

mat <- as.matrix(sub[,-1])
rownames(mat) <- sub$gene

oncoPrint(
  mat > 0,
  alter_fun = list(background = function(x, y, w, h) grid::grid.rect(x,y,w,h)),
  col = c("TRUE" = "firebrick"),
  column_title = "ARHG mutations – BRCA (top 200 samples)",
  use_raster = TRUE 
)


# Myeloid oncoprint
myeloid_samples <- all_results$sampleId[all_results$CancerType == "myeloid"]
length(unique(myeloid_samples))


myeloid_data <- all_results %>%
  filter(CancerType == "myeloid") %>%
  select(sampleId, hugoGeneSymbol, variantType) %>%
  summarise(n_mutations = n())

gene_counts <- myeloid_data %>%
  filter(!is.na(variantType) & variantType != "") %>%
  group_by(hugoGeneSymbol) %>%
  summarise(mutation_count = n(), .groups = 'drop') %>%
  arrange(desc(mutation_count)) %>%
  slice_head(n = 20)

top_genes_myeloid <- gene_counts$hugoGeneSymbol

myeloid_data_filtered <- myeloid_data %>%
  filter(hugoGeneSymbol %in% top_genes_myeloid)

# Collapse into matrix format
myeloid_wide <- myeloid_data_filtered %>%
  distinct(sampleId, hugoGeneSymbol, variantType) %>%
  pivot_wider(names_from = hugoGeneSymbol,
              values_from = variantType,
              values_fill = list(variantType = ""))

# Convert to matrix
mat <- as.matrix(myeloid_wide[,-1])
rownames(mat) <- myeloid_wide$sampleId

# Define colors for mutation types
mutation_colors <- c(
  "SNP" = "tomato",
  "DEL" = "purple",
  "INS" = "orange",
  " " = "white"
)

# Define alteration functions
alter_fun <- list(
  background = function(x, y, w, h) {
    grid.rect(x, y, w, h, gp = gpar(fill = "white", col = NA))
  },
  SNP = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = mutation_colors["SNP"], col = NA))
  },
  DEL = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = mutation_colors["DEL"], col = NA))
  },
  INS = function(x, y, w, h) {
    grid.rect(x, y, w*0.9, h*0.9, gp = gpar(fill = mutation_colors["INS"], col = NA))
  }
)

# Generate Oncoprint
oncoPrint(
  mat,
  alter_fun = alter_fun,
  col = mutation_colors,
  column_title = "ARHG Mutations in Myeloid Samples",
  remove_empty_columns = TRUE,
  remove_empty_rows = TRUE,
  alter_fun_is_vectorized = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE
)

