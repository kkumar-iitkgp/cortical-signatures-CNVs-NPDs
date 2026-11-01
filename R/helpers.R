# Helpers shared by all pages: configuration, text, and the four figures.
# Figures carry no in-plot text beyond axis labels; explanation lives around them.

suppressPackageStartupMessages({
  library(yaml)
  library(jsonlite)
  library(plotly)
  library(htmlwidgets)
  library(htmltools)
})

# ---- configuration --------------------------------------------------------------------------
site    <- yaml::read_yaml("config/site.yml", fileEncoding = "UTF-8")
metrics <- yaml::read_yaml("config/metrics.yml", fileEncoding = "UTF-8")
names(metrics) <- vapply(metrics, `[[`, "", "id")

read_site_csv <- function(name) read.csv(file.path("data/site", paste0(name, ".csv")),
                                         stringsAsFactors = FALSE)

# ---- palette (shared with styles.css) -------------------------------------------------------
pal <- list(
  # as in the paper's Figure 2 (ColorBrewer Dark2)
  variant     = "#E7298A",  # rare genetic variants (CNVs)
  diagnosis   = "#1B9E77",  # psychiatric diagnoses
  medicated   = "#D95F02",
  unmedicated = "#7570B3",
  common      = "#B8860B",  # common variants (SNPs)
  heritability = "#66A61E",
  # as in Figures 3 to 5: PC1 maps purple to orange; gradient sensorimotor (blue) to association (yellow)
  pc_low  = "#5E3C99", pc_mid = "#F4F2F0", pc_high = "#E66101",
  sm      = "#3E78B2", assoc = "#E9B824", grad_mid = "#E9E6DC",
  nodata  = "#B9BDC4",
  axis    = "#6B7480",
  grid    = "rgba(128,128,128,0.18)"
)
font_family <- "'Public Sans', system-ui, -apple-system, 'Segoe UI', sans-serif"

# ---- text helpers ---------------------------------------------------------------------------
term <- function(id, field = "label") metrics[[id]][[field]]

# Markdown for the collapsed "How to read this" box: a few lines specific to the panel,
# followed by the definitions of the terms it uses, pulled from metrics.yml.
how_to_read <- function(lines, ids = character(0)) {
  defs <- vapply(ids, function(id) {
    m <- metrics[[id]]
    txt <- paste0("**", m$label, ".** ", trimws(m$long))
    if (!is.null(m$caveat)) txt <- paste0(txt, " ", trimws(m$caveat))
    txt
  }, "")
  body <- paste(c(lines, defs), collapse = "\n\n")
  paste0('::: {.callout-note collapse="true" appearance="minimal" title="How to read this"}\n',
         body, "\n:::\n")
}

provenance <- function(ref, samples) {
  paste0('<p class="provenance">From ', ref, '. Samples: ', samples, '.</p>')
}

fmt_p <- function(p) if (p < 0.001) "< 0.001" else paste0("= ", formatC(p, format = "fg", digits = 2))
fmt_r <- function(r) sub("-", "\u2212", sprintf("%.2f", r))

# ---- shared plotly layout -------------------------------------------------------------------
base_layout <- function(p, ...) {
  p |>
    layout(
      font = list(family = font_family, size = 13, color = pal$axis),
      paper_bgcolor = "rgba(0,0,0,0)", plot_bgcolor = "rgba(0,0,0,0)",
      hoverlabel = list(font = list(family = font_family, size = 13)),
      showlegend = FALSE,
      ...
    ) |>
    config(displayModeBar = FALSE, responsive = TRUE)
}

# ---- Figure 1: global effect sizes (as in Figure 2 of the paper) -----------------------------
# Generic two-group dot strip for thickness (top) and surface area (bottom).
# scale = "abs": absolute Cohen's d on a log axis with medians (the paper's view);
# scale = "signed": signed Cohen's d on a linear axis, zero line, no medians.
strip_plot <- function(d, groups, labels, scale = c("abs", "signed")) {
  scale <- match.arg(scale)
  panel <- function(measure, measure_lab, show_x) {
    dd <- d
    dd$signed <- dd[[paste0("es_", measure)]]
    dd$x <- if (scale == "abs") abs(dd$signed) else dd$signed
    dd$ypos <- match(dd$group, rev(groups)) - 1
    dd <- dd[order(dd$group, dd$x), ]
    dd$ypos <- dd$ypos + ave(dd$x, dd$group, FUN = function(v) ((seq_along(v) %% 5) - 2) * 0.075)
    dd$hover <- paste0("<b>", dd$name, "</b><br>", measure_lab, ": d = ",
                       sub("-", "\u2212", formatC(dd$signed, format = "f", digits = 2)),
                       "<br>", dd$n_label, "<extra></extra>")
    p <- plot_ly()
    if (scale == "signed") {
      p <- add_segments(p, x = 0, xend = 0, y = -0.6, yend = length(groups) - 0.4,
                        line = list(color = pal$axis, width = 1), hoverinfo = "skip")
    }
    for (g in groups) {
      gd <- dd[dd$group == g, ]
      p <- add_markers(p, data = gd, x = ~x, y = ~ypos, hovertemplate = ~hover,
                       marker = list(size = 10, color = pal[[g]], opacity = 0.9,
                                     line = list(width = 1, color = "rgba(255,255,255,0.85)")))
      if (scale == "abs") {
        y0 <- match(g, rev(groups)) - 1; m <- median(gd$x)
        p <- add_segments(p, x = m, xend = m, y = y0 - 0.32, yend = y0 + 0.32,
                          line = list(color = pal[[g]], width = 3, dash = "dot"), hoverinfo = "text",
                          text = paste0("Median |d| = ", sprintf("%.2f", m)))
      }
    }
    xax <- if (scale == "abs") {
      list(type = "log", range = c(-2.6, 0.55), tickvals = c(0.003, 0.01, 0.03, 0.1, 0.3, 1, 3),
           ticktext = c("0.003", "0.01", "0.03", "0.1", "0.3", "1", "3"),
           title = if (show_x) list(text = "Effect size |d|, log scale") else "")
    } else {
      list(range = c(-2.85, 1.0), tickvals = c(-2.5, -2, -1.5, -1, -0.5, 0, 0.5),
           ticktext = c("\u22122.5", "\u22122", "\u22121.5", "\u22121", "\u22120.5", "0", "0.5"),
           title = if (show_x) list(text = "Effect size d (smaller \u2190 0 \u2192 larger than controls)") else "")
    }
    p |> layout(
      xaxis = c(xax, list(showgrid = TRUE, gridcolor = pal$grid, zeroline = FALSE, fixedrange = TRUE)),
      yaxis = list(tickvals = seq_along(groups) - 1, ticktext = unname(labels[rev(groups)]),
                   range = c(-0.55, length(groups) - 0.45), showgrid = FALSE, zeroline = FALSE,
                   title = list(text = measure_lab, standoff = 12), fixedrange = TRUE)
    )
  }
  subplot(panel("CT", term("thickness"), FALSE), panel("SA", term("surface_area"), TRUE),
          nrows = 2, shareX = TRUE, titleY = TRUE, margin = 0.06) |>
    base_layout(height = 430, margin = list(l = 10, r = 10, t = 10, b = 50)) |>
    slim()
}

# drop the duplicate copy of the data plotly keeps for re-building
slim <- function(p) { p <- plotly_build(p); p$x$attrs <- NULL; p }

# Two views of the same data, switched by buttons (signed first, then absolute).
global_figure <- function(global, id = "global") {
  labs <- c(variant = "Genetic variants", diagnosis = "Diagnoses")
  grp <- c("variant", "diagnosis")
  btn <- function(view, label, pressed) {
    tags$button(type = "button", class = "seg", `data-view` = view,
                `aria-pressed` = tolower(pressed), label)
  }
  tagList(
    div(class = "seg-group view-toggle", role = "group", `aria-label` = "Effect size view",
        `data-target` = id,
        btn("signed", "Direction and size", FALSE),
        btn("abs", "Size only", TRUE)),
    div(id = paste0(id, "-abs"), class = "view-panel", strip_plot(global, grp, labs, "abs")),
    div(id = paste0(id, "-signed"), class = "view-panel", hidden = NA,
        strip_plot(global, grp, labs, "signed"))
  )
}

view_toggle_script <- function() {
  tags$script(HTML("
document.addEventListener('click', function(e) {
  var b = e.target.closest('.view-toggle button'); if (!b) return;
  var g = b.parentElement, id = g.dataset.target;
  g.querySelectorAll('button').forEach(function(x) { x.setAttribute('aria-pressed', x === b); });
  ['abs', 'signed'].forEach(function(v) {
    var el = document.getElementById(id + '-' + v); if (!el) return;
    el.hidden = (v !== b.dataset.view);
    if (!el.hidden) el.querySelectorAll('.js-plotly-plot').forEach(function(p) { Plotly.Plots.resize(p); });
  });
});"))
}

plot_medication <- function(med) {
  strip_plot(med, c("medicated", "unmedicated"),
             c(medicated = "Medicated", unmedicated = "Unmedicated"), "abs")
}

# ---- Figure 2: rotatable brain with map switcher ----------------------------------------------
# One mesh, values swapped in the browser; controls are plain HTML buttons (see brain_controls()).
brain_maps_spec <- function(maps, grad) {
  r_of <- function(key) grad$r[grad$map == key]
  p_of <- function(key) grad$p_spin[grad$map == key]
  pc <- list(list(0, pal$pc_low), list(0.5, pal$pc_mid), list(1, pal$pc_high))
  list(
    gradient = list(col = "gradient", flip = FALSE, ends = c("Sensorimotor", "Association"),
                    cs = list(list(0, pal$sm), list(0.5, pal$grad_mid), list(1, pal$assoc)),
                    caption = "Reference: the sensorimotor-to-association gradient."),
    SA_CNVs = list(col = "PC1_SA_CNVs", r = r_of("SA_CNVs"), p = p_of("SA_CNVs"), cs = pc),
    SA_SNPs = list(col = "PC1_SA_SNPs", r = r_of("SA_SNPs"), p = p_of("SA_SNPs"), cs = pc),
    SA_NPDs = list(col = "PC1_SA_NPDs", r = r_of("SA_NPDs"), p = p_of("SA_NPDs"), cs = pc),
    CT_CNVs = list(col = "PC1_CT_CNVs", r = r_of("CT_CNVs"), p = p_of("CT_CNVs"), cs = pc),
    CT_NPDs = list(col = "PC1_CT_NPDs", r = r_of("CT_NPDs"), p = p_of("CT_NPDs"), cs = pc)
  )
}

plot_brain <- function(maps, grad, mesh_file = "data/brain/dk_lh_inflated.json",
                       el_id = "brain") {
  mesh <- jsonlite::fromJSON(mesh_file)
  region_name <- setNames(mesh$regions$name, mesh$regions$label)
  stopifnot(all(maps$label %in% names(region_name)))

  face_region <- mesh$vertex_region[mesh$i + 1]
  face_data <- face_region %in% maps$label

  # a mesh restricted to some faces, with its vertices re-indexed from 0
  submesh <- function(faces) {
    used <- sort(unique(c(mesh$i[faces], mesh$j[faces], mesh$k[faces]))) + 1L
    idx <- rep(NA_integer_, length(mesh$x)); idx[used] <- seq_along(used) - 1L
    list(x = round(mesh$x[used], 1), y = round(mesh$y[used], 1), z = round(mesh$z[used], 1),
         i = idx[mesh$i[faces] + 1], j = idx[mesh$j[faces] + 1], k = idx[mesh$k[faces] + 1],
         region = mesh$vertex_region[used])
  }
  dm <- submesh(face_data)
  gm <- submesh(!face_data)

  # values are sent per region (34 numbers per map); the browser expands them to vertices
  spec <- brain_maps_spec(maps, grad)
  region_values <- lapply(spec, function(s) {
    v <- maps[[s$col]]
    if (isTRUE(s$flip)) v <- -v
    round(v, 4)
  })
  vidx <- match(dm$region, maps$label) - 1L
  init <- round(region_values$SA_CNVs[vidx + 1], 2)
  shade <- list(ambient = 0.9, diffuse = 0.35, specular = 0.02, roughness = 1, fresnel = 0.05)

  p <- plot_ly() |>
    add_trace(
      type = "mesh3d", x = dm$x, y = dm$y, z = dm$z, i = dm$i, j = dm$j, k = dm$k,
      intensity = init, intensitymode = "vertex",
      colorscale = list(list(0, pal$pc_low), list(0.5, pal$pc_mid), list(1, pal$pc_high)),
      cmin = -max(abs(init)), cmax = max(abs(init)),
      hovertemplate = "%{text}<extra></extra>", flatshading = FALSE, lighting = shade,
      colorbar = list(orientation = "h", len = 0.45, thickness = 10, outlinewidth = 0, x = 0.5, xanchor = "center", y = 0.02, yanchor = "top",
                      tickvals = c(-0.92, 0.92) * max(abs(init)), ticktext = c("Lower", "Higher"),
                      tickfont = list(size = 12))
    ) |>
    add_trace(
      type = "mesh3d", x = gm$x, y = gm$y, z = gm$z, i = gm$i, j = gm$j, k = gm$k,
      hoverinfo = "skip", flatshading = FALSE, lighting = shade
    ) |>
    base_layout(
      height = 420, margin = list(l = 0, r = 0, t = 0, b = 30),
      scene = list(
        xaxis = list(visible = FALSE), yaxis = list(visible = FALSE), zaxis = list(visible = FALSE),
        aspectmode = "data", bgcolor = "rgba(0,0,0,0)",
        camera = list(eye = list(x = -1.85, y = 0, z = 0.1), up = list(x = 0, y = 0, z = 1))
      )
    )
  p <- plotly_build(p)
  p$x$attrs <- NULL # drop the duplicate copy of the data plotly keeps for re-building
  nd <- p$x$data[[2]]
  nd[c("intensity", "colorscale", "colorbar", "showscale", "marker", "line")] <- NULL
  nd$color <- pal$nodata
  p$x$data[[2]] <- nd
  p$elementId <- el_id

  js <- "
function(el, x, data) {
  var state = { measure: 'SA', layer: 'CNVs', side: 'outer' };
  var root = document.getElementById(data.controls);
  if (!root) return;
  var caption = root.querySelector('[data-caption]');
  var names = data.names.map(function(n) { return n.charAt(0).toUpperCase() + n.slice(1); });
  Plotly.restyle(el, { text: [data.vidx.map(function(i) { return names[i]; })] }, [0]);
  function key() { return state.layer === 'gradient' ? 'gradient' : state.measure + '_' + state.layer; }
  function fmt(v) { return (v < 0 ? '\\u2212' : '') + Math.abs(v).toFixed(2); }
  function fmtp(p) { return p < 0.001 ? '< 0.001' : '= ' + Number(p.toPrecision(2)); }
  function render() {
    if (state.measure === 'CT' && state.layer === 'SNPs') state.layer = 'CNVs';
    var k = key(), s = data.spec[k], rv = data.values[k];
    var v = data.vidx.map(function(i) { return rv[i]; });
    var m = Math.max.apply(null, rv.map(Math.abs));
    var ends = s.ends || ['Lower', 'Higher'];
    Plotly.restyle(el, {
      intensity: [v], cmin: [-m], cmax: [m], colorscale: [s.cs],
      'colorbar.tickvals': [[-m * 0.92, m * 0.92]], 'colorbar.ticktext': [ends]
    }, [0]);
    root.querySelectorAll('button[data-layer]').forEach(function(b) {
      b.setAttribute('aria-pressed', b.dataset.layer === state.layer);
      b.disabled = (state.measure === 'CT' && b.dataset.layer === 'SNPs');
    });
    root.querySelectorAll('button[data-measure]').forEach(function(b) {
      b.setAttribute('aria-pressed', b.dataset.measure === state.measure);
      b.disabled = state.layer === 'gradient';
    });
    root.querySelectorAll('button[data-side]').forEach(function(b) {
      b.setAttribute('aria-pressed', b.dataset.side === state.side);
    });
    if (k === 'gradient') {
      caption.textContent = s.caption;
    } else {
      var where = s.r < 0 ? 'sensorimotor' : 'association';
      var stats = ' (r = ' + fmt(s.r) + ' with the gradient, spin p ' + fmtp(s.p) + ').';
      caption.textContent = data.labels[state.layer] + ', ' + data.labels[state.measure] +
        (s.p < 0.05 ? ': higher scores in ' + where + ' regions' :
                      ': a tilt toward ' + where + ' regions that did not reach significance') + stats;
    }
  }
  root.addEventListener('click', function(e) {
    var b = e.target.closest('button'); if (!b || b.disabled) return;
    if (b.dataset.layer) state.layer = b.dataset.layer;
    if (b.dataset.measure) state.measure = b.dataset.measure;
    if (b.dataset.side) {
      state.side = b.dataset.side;
      var ex = state.side === 'outer' ? -1.85 : 1.85;
      Plotly.relayout(el, { 'scene.camera.eye': { x: ex, y: 0, z: 0.1 } });
    }
    render();
  });
  render();
}"
  htmlwidgets::onRender(p, js, data = list(
    controls = paste0(el_id, "-controls"),
    spec = spec, values = region_values, vidx = vidx,
    names = unname(region_name[maps$label]),
    labels = list(CNVs = "Rare variants", SNPs = "Common variants", NPDs = "Diagnoses",
                  SA = "surface area", CT = "thickness")
  ))
}

brain_controls <- function(el_id = "brain") {
  btn <- function(attr, val, label, pressed = FALSE) {
    args <- list(type = "button", class = "seg", `aria-pressed` = tolower(pressed), label)
    args[[paste0("data-", attr)]] <- val
    do.call(tags$button, args)
  }
  div(
    id = paste0(el_id, "-controls"), class = "brain-controls",
    div(class = "seg-group", role = "group", `aria-label` = "Map",
        btn("layer", "CNVs", "Rare variants", TRUE),
        btn("layer", "SNPs", "Common variants"),
        btn("layer", "NPDs", "Diagnoses"),
        btn("layer", "gradient", "Gradient")),
    div(class = "seg-row",
        div(class = "seg-group", role = "group", `aria-label` = "Measure",
            btn("measure", "SA", "Surface area", TRUE),
            btn("measure", "CT", "Thickness")),
        div(class = "seg-group", role = "group", `aria-label` = "View",
            btn("side", "outer", "Outer side", TRUE),
            btn("side", "inner", "Inner side"))),
    p(class = "brain-caption", `data-caption` = NA, `aria-live` = "polite")
  )
}

# ---- Figure 3: shared patterns and heritability against the gradient (as in Figure 5) ---------
plot_gradient_bars <- function(grad, herit) {
  rows <- data.frame(
    key = c("CNVs", "SNPs", "NPDs", "twin", "snph2"),
    lab = c("Rare variants", "Common variants", "Diagnoses", "Twin heritability", "SNP heritability"),
    col = c(pal$variant, pal$common, pal$diagnosis, pal$heritability, pal$heritability),
    sym = c("circle", "circle", "circle", "diamond", "diamond-tall"),
    stringsAsFactors = FALSE
  )
  rows$y <- rev(seq_len(nrow(rows)))
  panel <- function(measure, measure_lab, show_x) {
    d <- rows
    d$r <- NA_real_; d$sig <- TRUE; d$note <- ""
    for (k in c("CNVs", "SNPs", "NPDs")) {
      g <- grad[grad$map == paste0(measure, "_", k), ]
      if (nrow(g)) { d$r[d$key == k] <- g$r; d$sig[d$key == k] <- g$p_spin < 0.05
                     d$note[d$key == k] <- paste0("spin p ", fmt_p(g$p_spin)) }
    }
    for (k in c("twin", "snph2")) {
      h <- herit[herit$map == paste0(measure, "_", k), ]
      d$r[d$key == k] <- h$r; d$note[d$key == k] <- "spin p < 0.05"
    }
    # background: sensorimotor (blue) to association (yellow), as in Figure 5
    xs <- seq(-1, 1, length.out = 21)
    ramp <- colorRampPalette(c(pal$sm, "#FFFFFF", pal$assoc))(20)
    shapes <- lapply(1:20, function(n) list(type = "rect", xref = "x", yref = "paper",
                                            x0 = xs[n], x1 = xs[n + 1], y0 = 0, y1 = 1,
                                            fillcolor = ramp[n], opacity = 0.22, line = list(width = 0),
                                            layer = "below"))
    p <- plot_ly()
    for (n in which(!is.na(d$r))) {
      p <- p |>
        add_segments(x = 0, xend = d$r[n], y = d$y[n], yend = d$y[n], hoverinfo = "skip",
                     line = list(color = d$col[n], width = 2.5)) |>
        add_markers(x = d$r[n], y = d$y[n],
                    hovertemplate = paste0("<b>", d$lab[n], ", ", tolower(measure_lab), "</b><br>r = ",
                                           fmt_r(d$r[n]), "<br>", d$note[n], "<extra></extra>"),
                    marker = list(size = 12, symbol = if (d$sig[n]) d$sym[n] else paste0(d$sym[n], "-open"),
                                  color = d$col[n], line = list(width = 2, color = d$col[n])))
    }
    labs <- d$lab
    if (measure == "CT") labs[d$key == "SNPs"] <- "Common variants (area only)"
    p |> layout(
      shapes = shapes,
      xaxis = list(range = c(-1, 1), zeroline = TRUE, zerolinecolor = pal$axis, zerolinewidth = 1,
                   gridcolor = pal$grid, tickvals = c(-1, -0.5, 0, 0.5, 1),
                   ticktext = c("\u22121", "\u22120.5", "0", "0.5", "1"), fixedrange = TRUE,
                   title = if (show_x) list(text = "\u2190 sensorimotor    correlation with gradient (r)    association \u2192") else ""),
      yaxis = list(tickvals = d$y, ticktext = labs, range = c(0.4, nrow(d) + 0.6), showgrid = FALSE,
                   zeroline = FALSE, fixedrange = TRUE, automargin = TRUE,
                   title = list(text = measure_lab, standoff = 10))
    )
  }
  subplot(panel("CT", term("thickness"), FALSE), panel("SA", term("surface_area"), TRUE),
          nrows = 2, shareX = TRUE, titleY = TRUE, margin = 0.05) |>
    base_layout(height = 520, margin = list(l = 10, r = 15, t = 10, b = 55)) |>
    slim()
}

# ---- Figure 4: simulated risk-score association ----------------------------------------------
plot_prs <- function(sim) {
  dis_lab <- c(BD = "Bipolar disorder", SCZ = "Schizophrenia")
  sc_lab <- c(Empirical = "observed directions", Concordant = "all one direction")
  sc_col <- c(Empirical = pal$axis, Concordant = pal$common)
  sim <- sim[order(sim$disorder, sim$scenario != "Empirical"), ]
  sim$ypos <- c(4, 3, 1.6, 0.6)[seq_len(nrow(sim))]
  sim$lab <- paste0(dis_lab[sim$disorder], "<br>", sc_lab[sim$scenario])
  sim$hover <- paste0("<b>", dis_lab[sim$disorder], "</b>, ", sc_lab[sim$scenario],
                      "<br>Median R\u00b2 = ", formatC(sim$r2_median, format = "e", digits = 1),
                      "<br>95% interval ", formatC(sim$r2_lo, format = "e", digits = 1), " to ",
                      formatC(sim$r2_hi, format = "e", digits = 1),
                      "<br>", sim$n_snps, " variants; ", format(sim$n_sim, big.mark = ","),
                      " simulations<extra></extra>")
  p <- plot_ly()
  for (n in seq_len(nrow(sim))) {
    d <- sim[n, ]; col <- sc_col[[d$scenario]]
    p <- p |>
      add_segments(x = d$r2_lo, xend = d$r2_hi, y = d$ypos, yend = d$ypos, hoverinfo = "skip",
                   line = list(color = col, width = 2.5)) |>
      add_markers(x = d$r2_median, y = d$ypos, hovertemplate = d$hover,
                  marker = list(size = 14, color = col))
  }
  ticks <- 10^(-7:-2)  # log axis
  sup <- c("\u207b\u2077", "\u207b\u2076", "\u207b\u2075", "\u207b\u2074", "\u207b\u00b3", "\u207b\u00b2")
  p |> base_layout(
    height = 320, margin = list(l = 10, r = 15, t = 10, b = 55),
    xaxis = list(type = "log", tickvals = ticks, ticktext = paste0("10", sup), range = c(-7.5, -1.6),
                 gridcolor = pal$grid, zeroline = FALSE, fixedrange = TRUE,
                 title = list(text = "R\u00b2, log scale")),
    yaxis = list(tickvals = sim$ypos, ticktext = sim$lab, range = c(0, 4.6), showgrid = FALSE,
                 zeroline = FALSE, fixedrange = TRUE, automargin = TRUE)
  )
}

# numbers used in the text, computed from the same CSVs as the figures
fold_change <- function(global, measure = "SA") {
  v <- abs(global[[paste0("es_", measure)]])
  median(v[global$group == "variant"]) / median(v[global$group == "diagnosis"])
}
