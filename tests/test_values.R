# Every value the site shows must equal the supplementary tables (tolerance 1e-8).
#   Rscript tests/test_values.R
suppressPackageStartupMessages(library(readxl))
xlsx <- "data-raw/SupplementaryTables.xlsx"
if (!file.exists(xlsx)) { message("SKIP: ", xlsx, " not present"); quit(status = 0) }
tol <- 1e-8
fail <- 0
check <- function(ok, what) { if (!isTRUE(ok)) { message("FAIL: ", what); fail <<- fail + 1 } }
rd <- function(s) as.data.frame(read_excel(xlsx, sheet = s, .name_repair = "minimal"))
csv <- function(n) read.csv(file.path("data/site", paste0(n, ".csv")), stringsAsFactors = FALSE)

g <- csv("global_effects")
st1 <- rd("ST1"); st1 <- st1[!is.na(st1$CNV) & !is.na(st1$es_TotalSA), ]
st2 <- rd("ST2")
for (k in seq_len(nrow(st1))) {
  r <- g[g$group == "variant" & g$name == st1$CNV[k], ]
  check(nrow(r) == 1 && abs(r$es_SA - st1$es_TotalSA[k]) < tol && abs(r$es_CT - st1$es_MeanCT[k]) < tol,
        paste("ST1", st1$CNV[k]))
}
fig2 <- c(SCZ = "Schizophrenia", "CHRPS+" = "Clinical high risk, later psychosis",
          "CHRPS-" = "Clinical high risk, no later psychosis", OCD = "Obsessive-compulsive disorder",
          BD = "Bipolar disorder", MDD = "Major depressive disorder", CD = "Conduct disorder",
          ASD = "Autism", ADHD = "ADHD")
check(sum(g$group == "diagnosis") == 9, "9 diagnostic groups as in Figure 2")
for (k in names(fig2)) {
  r <- g[g$group == "diagnosis" & g$name == fig2[[k]], ]; o <- st2[st2$Disorder == k, ]
  check(nrow(r) == 1 && abs(r$es_SA - o$es_TotalSA) < tol && abs(r$es_CT - o$es_MeanCT) < tol, paste("ST2", k))
}
md <- csv("medication_effects"); st3 <- rd("ST3"); st3 <- st3[st3$Participant_age_group == "Adult", ]
check(nrow(md) == nrow(st3) && sum(md$group == "medicated") == 9, "ST3 adult medication groups")
check(all(abs(sort(md$es_CT) - sort(st3$es_MeanCT)) < tol) && all(abs(sort(md$es_SA) - sort(st3$es_TotalSA)) < tol), "ST3 values")
# headline numbers as printed in Figure 2
fc <- function(a, b) median(abs(a)) / median(abs(b))
check(round(fc(g$es_SA[g$group == "variant"], g$es_SA[g$group == "diagnosis"])) == 11, "11-fold SA")
check(round(fc(g$es_CT[g$group == "variant"], g$es_CT[g$group == "diagnosis"])) == 4, "4-fold CT")
check(round(fc(md$es_CT[md$group == "medicated"], md$es_CT[md$group == "unmedicated"])) == 8, "8-fold medication")

m <- csv("regional_maps"); st8 <- rd("ST8")
for (col in c("PC1_SA_CNVs", "PC1_SA_SNPs", "PC1_SA_NPDs", "PC1_CT_CNVs", "PC1_CT_NPDs")) {
  check(max(abs(m[[col]] - st8[[col]][match(m$label, paste0("lh_", st8$Desikan_roi))])) < tol, paste("ST8", col))
}
check(max(abs(m$gradient - st8$cortical_gradient[match(m$label, paste0("lh_", st8$Desikan_roi))])) < tol, "ST8 gradient")

gc <- csv("gradient_correlations"); st9 <- rd("ST9"); st9 <- st9[st9$consensus_measure == "PC1", ]
check(all(abs(gc$r - st9$cor[match(gc$map, st9$map)]) < tol), "ST9 r")
check(all(abs(gc$p_spin - st9$pspin_one_sided[match(gc$map, st9$map)]) < tol), "ST9 p")

primary <- function(x) x[grepl("(primary)", x$`SNP list version`, fixed = TRUE), ]
vt <- csv("variance_test"); st13 <- primary(rd("ST13"))
check(nrow(vt) == nrow(st13) && all(vt$disorder == st13$Disorder) && all(vt$n_snps == st13$`N SNPs`), "ST13 primary rows")
check(all(abs(vt$variance_ratio - st13$`Variance Ratio`) < tol), "ST13 ratio")
check(all(abs(vt$p_variance - st13$`P (variance)`) < tol), "ST13 p")
check(all(abs(vt$pct_smaller - st13$`% risk alleles with smaller cortical estimate`) < tol), "ST13 sign split")

ps <- csv("prs_simulation"); st14 <- primary(rd("ST14")); names(st14)[grepl("^Median R", names(st14))] <- "med"
key <- paste(st14$Disorder, st14$Scenario, st14$`N SNPs`)
check(nrow(ps) == nrow(st14), "ST14 primary rows")
check(all(abs(ps$r2_median - st14$med[match(paste(ps$disorder, ps$scenario, ps$n_snps), key)]) < tol), "ST14 median")

# hygiene: nothing third-party in the site data
all_cols <- unlist(lapply(list.files("data/site", full.names = TRUE), function(f) names(read.csv(f))))
check(!any(grepl("HR_|iPSYCH|BETA|SNP_", all_cols)), "no third-party columns")

if (fail > 0) { message(fail, " check(s) failed"); quit(status = 1) }
message("All value checks passed")
