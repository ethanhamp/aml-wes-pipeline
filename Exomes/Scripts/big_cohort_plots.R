library(dplyr)
library(maftools)
library(purrr)
library(stringr)
library(tidyverse)
library(readxl)

BASE_DIR <- Sys.getenv("ARHG_BASE_DIR", unset = normalizePath("."))
ARTIFACT_GENE_FILE <- file.path(BASE_DIR, "General/artifact_genes.txt")

dev.off()

big_cohort <- read_excel(file.path(BASE_DIR, "Exomes/Datasets/big_cohort_wes.xlsx"))

combined <- big_cohort[big_cohort$`AD Total` >=20, ]

# excluded_samples and hypermutators loaded from config/ — see README

artifact_genes <- if (file.exists(ARTIFACT_GENE_FILE)) {
  lines <- readLines(ARTIFACT_GENE_FILE, warn = FALSE)
  lines[nzchar(lines) & !startsWith(lines, "#")]
} else { character(0) }

combined <- combined[!combined$Gene %in% artifact_genes, ]

# all_samples_clean <- all_samples_clean[all_samples_clean$`Alt Percentage` >0.02, ]

variant_key <- "MANE HGVS"

clean <- combined %>%
  mutate(
    Timepoint = tolower(trimws(as.character(`Sample Type`))),     # "diagnosis"/"relapse"
    Group     = as.character(File),
    HGVS      = str_squish(as.character(`MANE HGVS`)),
    # optional: normalize transcript version so NM_XXXX.6:c.xxx == NM_XXXX:c.xxx
    HGVS      = str_replace(HGVS, "(NM_\\d+)\\.\\d+(:.*)", "\\1\\2"),
    VAF_raw   = as.character(`Alt Percentage`),
    VAF       = parse_number(VAF_raw),
    VAF       = ifelse(str_detect(VAF_raw, "%"), VAF/100, VAF)
  ) %>%
  filter(Timepoint %in% c("germline","tumor"),
         !is.na(HGVS) & HGVS != "")

variant_patient_counts <- clean %>%
  distinct(`Sample ID`, HGVS) %>%    # don’t double-count the same patient across timepoints
  count(HGVS, name = "n_patients")


clean <- clean %>%
  inner_join(variant_patient_counts, by = "HGVS") %>%
  filter(n_patients >= 2)

infer_variant_type <- function(ref, alt){
  ref <- toupper(trimws(ref)); alt <- toupper(trimws(alt))
  clen <- function(x){
    x <- gsub("-", "", x, fixed = TRUE)
    x <- ifelse(grepl("^<.*>$", x), "", x)  # treat symbolic like <DEL> as 0
    nchar(x)
  }
  lr <- clen(ref); la <- clen(alt)
  out <- ifelse(is.na(lr) | is.na(la), NA_character_,
                ifelse(lr == 0 & la > 0, "INS",
                       ifelse(lr > 0 & la == 0, "DEL",
                              ifelse(lr == la & lr == 1, "SNP",
                                     ifelse(lr == la & lr == 2, "DNP",
                                            ifelse(lr == la & lr == 3, "TNP",
                                                   ifelse(lr == la & lr >= 4, "ONP",
                                                          ifelse(lr != la, ifelse(lr > la, "DEL", "INS"), NA_character_)
                                                   )))))))
  out
}

# 2) Map your column names (edit these if needed)
ref_col <- intersect(names(clean), c("Reference_Allele","REF","Ref","reference","ref"))[1]
alt_col <- intersect(names(clean), c("Tumor_Seq_Allele2","ALT","Alt","allele","alt"))[1]
stopifnot(!is.na(ref_col), !is.na(alt_col))

# 3) If ALT can have multiple entries (e.g., "A,C"), split to one row per allele,
#    then compute Variant_Type and (optionally) re-aggregate later if needed.
clean <- clean %>%
  mutate(!!alt_col := as.character(.data[[alt_col]])) %>%
  separate_rows(!!sym(alt_col), sep = ",") %>%
  mutate(
    Variant_Type = infer_variant_type(.data[[ref_col]], .data[[alt_col]])
  )

# 1) Canonical MAF labels we want to end up with
maf_levels <- c(
  "Frame_Shift_Del","Frame_Shift_Ins",
  "Nonsense_Mutation","Nonstop_Mutation",
  "Splice_Site","Translation_Start_Site",
  "In_Frame_Del","In_Frame_Ins",
  "Missense_Mutation","Silent",
  "5'UTR","3'UTR","5'Flank","3'Flank",
  "Intron","IGR","RNA"
)

# 0) Nuke any older defs so they don’t linger
if (exists("canonicalize_variant_class")) rm(canonicalize_variant_class)
if (exists("map_token_to_maf")) rm(map_token_to_maf)

# 1) Severity order (left = most severe)
severity_rank <- c(
  "Frame_Shift_Del","Frame_Shift_Ins",
  "Nonsense_Mutation","Nonstop_Mutation",
  "Splice_Site","Translation_Start_Site",
  "In_Frame_Del","In_Frame_Ins",
  "Missense_Mutation","Silent",
  "5'UTR","3'UTR","5'Flank","3'Flank",
  "Intron","IGR","RNA","Unknown"
)
sev_index <- function(x) match(x, severity_rank, nomatch = length(severity_rank))

# 2) Token -> canonical MAF label (uses Variant_Type to decide FS Ins/Del)
map_token_to_maf <- function(tok, variant_type) {
  t  <- str_trim(tolower(ifelse(is.na(tok), "", tok)))
  vt <- toupper(ifelse(is.na(variant_type), "", variant_type))
  is_ins <- vt == "INS"
  is_del <- vt == "DEL"
  
  dplyr::case_when(
    t == "" ~ NA_character_,
    
    # nonsense / nonstop / start
    t %in% c("frameshift&stopgain","stopgain&nonframeshift insertion",
             "stopgain;nonsynonymous SNV", "nonsynonymous SNV;stopgain", "stopgain") ~ "Nonsense_Mutation",
    t %in% c("stoploss","stop_lost","stop_lost_variant","nonstop_mutation",
             "Nonstop_Mutation;Nonsense_Mutation") ~ "Nonstop_Mutation",
    t %in% c("startloss;nonsynonymous SNV","initiator_codon;nonsynonymous SNV","initiator_codon",
             "startloss") ~ "Translation_Start_Site",
    
    # splice
    str_detect(t, "splice_acceptor|splice_donor") ~ "Splice_Site",
    str_detect(t, "splice_region") ~ "Splice_Site",  # change to Splice_Region if you prefer
    
    # frameshift
    str_detect(t, "frameshift") & is_del ~ "Frame_Shift_Del",
    str_detect(t, "frameshift") & is_ins ~ "Frame_Shift_Ins",
    str_detect(t, "frameshift") ~ "Frame_Shift_Ins",  # fallback if VT unknown
    
    # in-frame indels
    str_detect(t, "nonframeshift deletion|inframe_deletion") ~ "In_Frame_Del",
    str_detect(t, "nonframeshift insertion|inframe_insertion|dup|Frame_Shift_Ins&Translation_Start_Site|
               Frame_Shift_Ins&Nonstop_Mutation|In_Frame_Ins;Frame_Shift_Ins|frameshift insertion") ~ "In_Frame_Ins",
    
    # missense / synonymous
    t %in% c("missense","missense_variant","nonsynonymous snv","nonsynonymous_snv",
             "nonsynonymous SNV", "Missense_Mutation;Translation_Start_Site",
             "Missense_Mutation;initiator_codor", "Missense_Mutation;Nonstop_Mutation",
             "Missense_Mutation&Nonsense_Mutation") ~ "Missense_Mutation",
    t %in% c("synonymous snv","synonymous_snv","synonymous_variant") ~ "Silent",
    
    # UTR / flank / intron / intergenic
    str_detect(t, "3.?prime_?utr|^utr3$") ~ "3'UTR",
    str_detect(t, "5.?prime_?utr|^utr5$") ~ "5'UTR",
    str_detect(t, "^intron") ~ "Intron",
    str_detect(t, "^intergenic|\\big r\\b|intergenic_region") ~ "IGR",
    str_detect(t, "^upstream") ~ "5'Flank",
    str_detect(t, "^downstream") ~ "3'Flank",
    
    # fusions → keep out of non-syn bucket
    str_detect(t, "gene[_-]?fusion|fusion") ~ "RNA",
    
    TRUE ~ "Unknown"
  )
}

# 3) Collapse multi-annotations per row: pick most severe
canonicalize_variant_class <- function(raw_vec, variant_type_vec) {
  n <- length(raw_vec)
  out <- character(n)
  
  for (i in seq_len(n)) {
    raw_i <- ifelse(is.na(raw_vec[i]), "", raw_vec[i])
    vt_i  <- ifelse(is.na(variant_type_vec[i]), NA_character_, variant_type_vec[i])
    
    parts <- unlist(strsplit(raw_i, "[;&,]"))
    parts <- trimws(parts)
    parts <- parts[parts != ""]
    
    if (length(parts) == 0) { out[i] <- NA_character_; next }
    
    labs <- vapply(parts, function(tok) map_token_to_maf(tok, vt_i), character(1))
    labs <- unique(labs[!is.na(labs)])
    if (!length(labs)) { out[i] <- NA_character_; next }
    
    out[i] <- labs[which.min(sev_index(labs))]
  }
  
  out
}

maf_df <- clean %>%
  rename(Hugo_Symbol = Gene,
         Chromosome = Chr,
         Start_Position = Start,
         End_Position = Stop,
         Variant_Classification = Effect,
         Variant_Type = Variant_Type,
         Tumor_Sample_Barcode = Group,
         Reference_Allele = Ref,
         Tumor_Seq_Allele2 = Alt)

# 4) APPLY to your data frame (replace DF with your object)
# Requires DF has columns: Variant_Classification (raw strings) and Variant_Type (SNP/DNP/TNP/ONP/INS/DEL)
DF <- maf_df %>%
  mutate(
    Variant_Classification = canonicalize_variant_class(Variant_Classification, Variant_Type),
    
    # Optional backfill if (and only if) in THIS file blank means splice:
    Variant_Classification = if_else(
      is.na(Variant_Classification) | Variant_Classification == "Unknown",
      "Splice_Site",
      Variant_Classification
    )
  )
DF$Hugo_Symbol <- sub("^SLC.*", "SLC", DF$Hugo_Symbol)
DF$Hugo_Symbol <- sub("^ARHG.*", "ARHG", DF$Hugo_Symbol)
DF$Hugo_Symbol <- sub("^DOCK.*", "DOCK", DF$Hugo_Symbol)
DF$Hugo_Symbol <- sub("^ZNF.*", "ZNF", DF$Hugo_Symbol)
DF$Hugo_Symbol <- sub("^PCDH.*", "PCDH", DF$Hugo_Symbol)

diagnosis_df <- DF[!is.na(DF$Timepoint) & DF$Timepoint == "tumor" & DF$`Alt Percentage`>0.02, , drop = FALSE]

# QA
table(DF$Variant_Classification, useNA = "ifany")

maf <- read.maf(maf = DF)
diagnosis_maf <- read.maf(maf = diagnosis_df)

vc_cols = c(
  'red3',
  'red3',
  'red3',
  'grey0',
  'grey0',
  'steelblue',
  'grey0',
  'grey50',
  'red'
)

names(vc_cols) = c(
  'Frame_Shift_Ins',
  'Splice_Site',
  'Nonsense_Mutation',
  'Missense_Mutation',
  'In_Frame_Ins',
  'Multi_Hit',
  'In_Frame_Del', #TKD
  'Frame_Shift_Del', #ITD
  'Translation_Start_Site'
)

# -------------------------------------------------------------
getGeneSummary(maf)
plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)

oncoplot(maf = maf, top = 30)
plotVaf(maf = maf, vafCol = 'Alt Percentage')

somaticInteractions(maf = maf, top = 25, pvalue = c(0.05, 0.1))

maf.sig = oncodrive(maf = maf, AACol = 'MANE HGVS', minMut = 5, pvalMethod = 'zscore')

plotOncodrive(res = maf.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)

# diagnosis samples
plotmafSummary(maf = diagnosis_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
oncoplot(maf = diagnosis_maf, top = 30, colors = vc_cols, borderCol = 'grey80',
         bgCol = 'white', titleText = "Diagnosis samples (n=154)")
plotVaf(maf = diagnosis_maf, vafCol = 'Alt Percentage')

somaticInteractions(maf = diagnosis_maf, top = 25, pvalue = c(0.05, 0.1), fontSize = 0.5, leftMar = 6)

diag.maf.sig = oncodrive(maf = diagnosis_maf, AACol = 'MANE HGVS', minMut = 5, pvalMethod = 'zscore')

plotOncodrive(res = diag.maf.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)

# relapse samples
plotmafSummary(maf = relapse_maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
oncoplot(maf = relapse_maf, top = 30, titleText = "Relapse samples (n=29)")
plotVaf(maf = maf, vafCol = 'Alt Percentage')

somaticInteractions(maf = relapse_maf, top = 25, pvalue = c(0.05, 0.1), fontSize = 0.5)

rel.maf.sig = oncodrive(maf = relapse_maf, AACol = 'MANE HGVS', minMut = 5, pvalMethod = 'zscore')

plotOncodrive(res = rel.maf.sig, fdrCutOff = 0.1, useFraction = TRUE, labelSize = 0.5)

genes = c("NPM1", "FLT3", "TET2", "ARHG", "NRAS", "PHIP", "DNMT3A", "ASXL1", "HOMEZ", "IDH2",
          "SRSF2", "WT1", "TP53", "SLC", "TUBB8B", "ZNF", "DEAF1", "PCDH")
coOncoplot(m1 = diagnosis_maf, m2 = relapse_maf, m1Name = 'Diagnosis AML', m2Name = 'Relapse AML', genes = genes, removeNonMutated = TRUE)

coBarplot(m1 = diagnosis_maf, m2 = relapse_maf, m1Name = "Diagnosis", m2Name = "Relapse", 
          genes = genes, orderBy = "m2", legendTxtSize = 0.8)
summarise(VAF = max(VAF, na.rm = TRUE), .groups = "drop")