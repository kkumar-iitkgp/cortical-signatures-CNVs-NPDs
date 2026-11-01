# Extract the Desikan-Killiany left-hemisphere inflated surface mesh used for the 3D maps.
#
# Source: ggseg3d (Mowinckel & Vidal-Pineiro 2020, doi:10.1177/2515245920928009; MIT licence),
# dataset `dk_3d`, pinned to the last commit before ggseg3d moved from plotly to Three.js.
# Run once; the output (data/brain/dk_lh_inflated.json) is committed, so the site does not
# depend on any particular ggseg3d version at render time.
#
#   Rscript R/extract_mesh.R

library(jsonlite)

commit <- "23cc79bdb4689dcaba3ba3bfa201c3739229c663"
url <- sprintf("https://raw.githubusercontent.com/ggseg/ggseg3d/%s/data/dk_3d.rda", commit)
tmp <- tempfile(fileext = ".rda")
download.file(url, tmp, mode = "wb", quiet = TRUE)
load(tmp) # creates dk_3d

hemi <- dk_3d$ggseg_3d[[which(dk_3d$surf == "inflated" & dk_3d$hemi == "left")]]
hemi <- hemi[!is.na(hemi$label), ]

x <- y <- z <- numeric(0)
i <- j <- k <- integer(0)
vertex_region <- character(0)
regions <- data.frame(label = character(0), name = character(0))
offset <- 0L

for (r in seq_len(nrow(hemi))) {
  m <- hemi$mesh[[r]]
  lab <- hemi$label[r]
  nv <- nrow(m$vertices)
  x <- c(x, m$vertices$x)
  y <- c(y, m$vertices$y)
  z <- c(z, m$vertices$z)
  # ggseg3d faces are 1-based; plotly mesh3d wants 0-based
  i <- c(i, as.integer(m$faces$i) - 1L + offset)
  j <- c(j, as.integer(m$faces$j) - 1L + offset)
  k <- c(k, as.integer(m$faces$k) - 1L + offset)
  vertex_region <- c(vertex_region, rep(lab, nv))
  regions <- rbind(regions, data.frame(
    label = lab,
    name = ifelse(is.na(hemi$region[r]), "medial wall", hemi$region[r])
  ))
  offset <- offset + nv
}

out <- list(
  source = sprintf("ggseg3d dk_3d (inflated, left), commit %s", commit),
  regions = regions,
  x = round(x, 3), y = round(y, 3), z = round(z, 3),
  i = i, j = j, k = k,
  vertex_region = vertex_region
)

dir.create("data/brain", showWarnings = FALSE, recursive = TRUE)
write_json(out, "data/brain/dk_lh_inflated.json", auto_unbox = TRUE, digits = NA)
message("Wrote data/brain/dk_lh_inflated.json: ", length(x), " vertices, ",
        length(i), " faces, ", nrow(regions), " regions")
