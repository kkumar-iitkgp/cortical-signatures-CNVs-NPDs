# cortical-signatures-CNVs-NPDs

Live site: https://kkumar-iitkgp.github.io/cortical-signatures-CNVs-NPDs/

A small static site that shows the main results of
*Copy number variants reveal divergent genetic and diagnostic cortical signatures across psychiatric disorders*
(Kumar, Jacquemont et al., Communications Medicine) for a general and scientific audience.

It shows only results reported in the paper and its supplementary tables.

## Layout

| Path | What it is |
|---|---|
| `index.qmd` | Findings page: three findings, one takeaway each |
| `about.qmd` | Methods in brief, data sources, limitations, glossary, citation |
| `config/site.yml` | Paper facts in one place: status, DOI, figure references, sample sizes |
| `config/metrics.yml` | Every definition shown on the site (tooltips, "How to read this", glossary) |
| `data/site/*.csv` | The only data the site reads, derived from the supplementary tables |
| `data/brain/dk_lh_inflated.json` | Desikan left-hemisphere mesh from ggseg3d (pinned version) |
| `R/prepare_data.R` | Supplementary tables to `data/site/` |
| `R/extract_mesh.R` | Re-creates the brain mesh file (run once) |
| `R/helpers.R` | Text helpers and the four figures |
| `tests/test_values.R` | Checks every displayed value against the supplementary tables |

## Build

Requires R (plotly, htmlwidgets, htmltools, jsonlite, yaml, knitr, rmarkdown; readxl, dplyr and
tidyr for data preparation) and Quarto 1.6 or later.

```sh
# put the workbook at data-raw/SupplementaryTables.xlsx (not committed)
Rscript R/prepare_data.R      # tables -> data/site/*.csv
Rscript tests/test_values.R   # every displayed value checked against the tables
quarto preview                # local preview
quarto render                 # writes the static site to docs/
```

The static site is the `docs/` folder (plain HTML, CSS and JavaScript). GitHub Pages serves it
directly: Settings > Pages > Deploy from a branch > `main`, folder `/docs`.

## Schizophrenia variant list

The paper's primary schizophrenia list is now the European-ancestry lead variants (128 in the
cortical GWAS); the original 259-variant list is the sensitivity analysis. `R/prepare_data.R`
reads ST10 and ST11 from the workbook; when the v3 workbook carries both schizophrenia rows,
keep the primary (128) row there and the site text updates from the CSV.

## Licence

Content CC BY 4.0; code MIT. See `LICENSE`.
