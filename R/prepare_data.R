# Build the small, derived CSVs the site displays, from the manuscript's supplementary tables.
#
# Input : data-raw/SupplementaryTables.xlsx  (not committed; see .gitignore)
# Output: data/site/*.csv                    (committed; the only data the site reads)
#
#   Rscript R/prepare_data.R
#
# Only values shown on the site are kept. Not exported: iPSYCH hazard ratios (ST1 HR_*),
# demographics beyond N, per-SNP tables (ST6), regional PRS estimates (ST12), medication
# young and pediatric sensitivity groups, and any per-disorder or per-CNV regional map.

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(tidyr)
})

xlsx <- "data-raw/SupplementaryTables.xlsx"
stopifnot(file.exists(xlsx))
out_dir <- "data/site"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
rd <- function(sheet) read_excel(xlsx, sheet = sheet, .name_repair = "minimal")

# ---- 1. Global effect sizes: structural variants (ST1) and diagnoses (ST2) -------------------
st1 <- rd("ST1")
names(st1)[names(st1) == "Sample size (N)"] <- "N"

# ST1 lists each variant on one row, followed by continuation rows (CNV = NA) that hold the
# UK Biobank carriers of the same variant. Carriers are summed per variant; control rows dropped.
st1 <- st1 |>
  mutate(variant = CNV) |>
  fill(variant, .direction = "down") |>
  filter(variant != "Controls")

carriers <- st1 |>
  mutate(n_num = suppressWarnings(as.numeric(N))) |>
  group_by(variant) |>
  summarise(n_carriers = sum(n_num, na.rm = TRUE), .groups = "drop")

variants <- st1 |>
  filter(!is.na(CNV), !is.na(es_TotalSA)) |>
  left_join(carriers, by = "variant") |>
  transmute(
    group = "variant",
    name = CNV,
    subtype = case_when(
      Type_CNV_SCA == "SCA" ~ "Sex chromosome aneuploidy",
      Type_Del_Dup == "Del" ~ "Deletion",
      Type_Del_Dup == "Dup" ~ "Duplication",
      TRUE ~ NA_character_
    ),
    es_SA = es_TotalSA,
    es_CT = es_MeanCT,
    source = if_else(Flag_published_computed == "computed",
                     "This study (individual-level MRI)", "Published meta-analysis"),
    pmid = as.character(PubMedID_data_source),
    n_label = if_else(Flag_published_computed == "computed",
                      paste0(n_carriers, " carriers"), "Published estimate")
  )

st2 <- rd("ST2")
# Figure 2 shows the 9 main diagnostic groups (8 disorders; clinical high risk split by later
# psychosis). Young and pediatric groups were sensitivity analyses and are not shown.
fig2_groups <- c("SCZ", "CHRPS+", "CHRPS-", "OCD", "BD", "MDD", "CD", "ASD", "ADHD")
dx_names <- c(SCZ = "Schizophrenia", "CHRPS+" = "Clinical high risk, later psychosis",
              "CHRPS-" = "Clinical high risk, no later psychosis", OCD = "Obsessive-compulsive disorder",
              BD = "Bipolar disorder", MDD = "Major depressive disorder", CD = "Conduct disorder",
              ASD = "Autism", ADHD = "ADHD")
diagnoses <- st2 |>
  filter(Disorder %in% fig2_groups) |>
  transmute(
    group = "diagnosis",
    name = unname(dx_names[Disorder]),
    subtype = "ENIGMA case-control group",
    es_SA = es_TotalSA,
    es_CT = es_MeanCT,
    source = "ENIGMA summary statistics",
    pmid = as.character(PubMedID_data_source),
    n_label = paste0(format(`Cases sample size (N)`, big.mark = ",", trim = TRUE), " cases, ",
                     format(`Controls sample size (N)`, big.mark = ",", trim = TRUE), " controls"),
    n_cases = `Cases sample size (N)`
  )

global <- bind_rows(variants, diagnoses)
stopifnot(sum(global$group == "variant") == 18, sum(global$group == "diagnosis") == 9)
write.csv(global, file.path(out_dir, "global_effects.csv"), row.names = FALSE)

# Medication subgroups as in Figure 2D-E: adult groups only (9 medicated, 5 unmedicated)
st3 <- rd("ST3")
med_names <- c("SCZ-m12" = "Schizophrenia, 1st and 2nd generation antipsychotics",
               "SCZ-m1" = "Schizophrenia, 1st generation antipsychotics",
               "SCZ-m2" = "Schizophrenia, 2nd generation antipsychotics",
               "BD-m2" = "Bipolar disorder, antiepileptics", "BD-m1" = "Bipolar disorder, lithium",
               ASD = "Autism", OCD = "Obsessive-compulsive disorder",
               MDD = "Major depressive disorder", ADHD = "ADHD", SCZ = "Schizophrenia")
medication <- st3 |>
  filter(Participant_age_group == "Adult") |>
  transmute(
    group = NPD_med_groups,
    name = unname(med_names[Disorder]),
    es_SA = es_TotalSA, es_CT = es_MeanCT,
    n_label = paste0(format(`Cases sample size (N)`, big.mark = ",", trim = TRUE), " cases, ",
                     format(`Controls sample size (N)`, big.mark = ",", trim = TRUE), " controls")
  )
stopifnot(sum(medication$group == "medicated") == 9, sum(medication$group == "unmedicated") == 5,
          !anyNA(medication$name))
write.csv(medication, file.path(out_dir, "medication_effects.csv"), row.names = FALSE)

# ---- 2. Regional latent maps (ST8, PC1 columns) and the cortical gradient -------------------
st8 <- rd("ST8")
maps <- st8 |>
  transmute(
    label = paste0("lh_", Desikan_roi),
    gradient = cortical_gradient,
    PC1_SA_CNVs, PC1_SA_SNPs, PC1_SA_NPDs, PC1_CT_CNVs, PC1_CT_NPDs
  )
stopifnot(nrow(maps) == 34, !anyDuplicated(maps$label))
write.csv(maps, file.path(out_dir, "regional_maps.csv"), row.names = FALSE)

# ---- 3. Map vs gradient correlations (ST9, PC1 rows) ----------------------------------------
st9 <- rd("ST9")
grad <- st9 |>
  filter(consensus_measure == "PC1") |>
  transmute(map, r = cor, p_spin = pspin_one_sided)
stopifnot(nrow(grad) == 5)
write.csv(grad, file.path(out_dir, "gradient_correlations.csv"), row.names = FALSE)

# Heritability vs gradient correlations, as shown in Figure 5 (values from the supplementary
# figure on regional twin and SNP heritability; regional SA adjusted for total SA).
heritability <- data.frame(
  map = c("CT_twin", "SA_twin", "CT_snph2", "SA_snph2"),
  r = c(-0.75, -0.69, -0.49, -0.47),
  significant = TRUE
)
write.csv(heritability, file.path(out_dir, "heritability_gradient.csv"), row.names = FALSE)

# ---- 4. Opposing effects: sign split and variance test (ST13), PRS simulation (ST14) ----------
is_primary <- function(version) grepl("(primary)", version, fixed = TRUE)
st13 <- rd("ST13") |>
  filter(is_primary(`SNP list version`)) |>
  transmute(disorder = Disorder, n_snps = `N SNPs`,
            variance_ratio = `Variance Ratio`, p_variance = `P (variance)`,
            pct_smaller = `% risk alleles with smaller cortical estimate`)
write.csv(st13, file.path(out_dir, "variance_test.csv"), row.names = FALSE)

st14 <- rd("ST14")
names(st14)[grepl("^Median R", names(st14))] <- "Median R2" # avoid locale issues with the superscript
st14 <- st14 |>
  filter(is_primary(`SNP list version`)) |>
  transmute(disorder = Disorder, scenario = Scenario, n_snps = `N SNPs`,
            r2_median = `Median R2`, r2_lo = `95% CI Lower`, r2_hi = `95% CI Upper`,
            n_sim = `N Simulations`)
write.csv(st14, file.path(out_dir, "prs_simulation.csv"), row.names = FALSE)

message("Wrote ", length(list.files(out_dir)), " files to ", out_dir)
