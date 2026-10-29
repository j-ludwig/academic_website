# =============================================================================
# looking_at_data_figures.R  --  data and figures for the chapter "Looking at Multivariate Data" of the ISQS 6350 course notes.
# Put this file and _common.R (shared helpers) in one folder, set it as the working
# directory, and run source("looking_at_data_figures.R"). Every figure is a function, e.g. fig_...().
# Packages: tidyverse, plotly, patchwork, mvtnorm, MASS, GGally, ggrepel, palmerpenguins, lattice
# =============================================================================

source("_common.R", local = TRUE)
suppressPackageStartupMessages({
  library(GGally)
  library(ggrepel)
  library(palmerpenguins)
})

f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)

# ---- Data ----------------------------------------------------------------------
# ---- Data -------------------------------------------------------------------
states <- as.data.frame(state.x77) |>
  rownames_to_column("state") |>
  rename(LifeExp = `Life Exp`, HSGrad = `HS Grad`) |>
  mutate(abb = state.abb, region = state.region)
vars6 <- c("Income", "Illiteracy", "LifeExp", "Murder", "HSGrad", "Frost")
labs6 <- c("Income", "Illiteracy", "Life exp.", "Murder", "HS grad", "Frost")
X6 <- as.matrix(states[, vars6]); rownames(X6) <- states$state
n_st <- nrow(X6)
peng <- na.omit(penguins)
pal_species <- c(Adelie = "#1b9e77", Chinstrap = "#d95f02", Gentoo = "#7570b3")
pal_region  <- c("Northeast" = "#1b9e77", "South" = "#d95f02",
                 "North Central" = "#7570b3", "West" = "#e7298a")

# ---- Numbers quoted in the text --------------------------------------------------
# ---- Numbers quoted in the text ------------------------------------------------
v2   <- c("Illiteracy", "Murder")
xb2  <- colMeans(X6[, v2]); S2 <- cov(X6[, v2])
# the unusual-combination example: illiteracy and high-school graduation, Hawaii
vU   <- c("Illiteracy", "HSGrad")
xbU  <- colMeans(X6[, vU]); SU <- cov(X6[, vU])
z_hi <- scale(X6[, vU])["Hawaii", ]
dmU  <- sqrt(mahalanobis(X6[, vU], xbU, SU))
r_ih <- cor(states$Illiteracy, states$HSGrad)
d2_st <- mahalanobis(X6, colMeans(X6), cov(X6))
r_all <- cor(peng$bill_length_mm, peng$bill_depth_mm)
r_within <- peng |> group_by(species) |> summarise(r = cor(bill_length_mm, bill_depth_mm)) |> deframe()
s_im <- cov(states$Illiteracy, states$Murder)
r_im <- cor(states$Illiteracy, states$Murder)

# =============================================================================
# The data matrix
# =============================================================================

# The data matrix with its actual numbers: the first rows, a gap, Texas and the last row.
# Returns an HTML table (Texas row blue, Murder column red, as in fig_heat()).
data_matrix_table <- function(show = c("Alabama", "Alaska", "Arizona", "Arkansas", "California",
                                       "Texas", "Wyoming"), font = 18) {
  fmt <- list(Income = function(v) format(v, big.mark = ","), Illiteracy = function(v) sprintf("%.1f", v),
              LifeExp = function(v) sprintf("%.1f", v), Murder = function(v) sprintf("%.1f", v),
              HSGrad = function(v) sprintf("%.1f", v), Frost = function(v) as.character(v))
  units <- c("$ per capita", "%", "years", "per 100k", "%", "days")
  td <- function(x, style = "") sprintf('<td style="padding:3px 12px;text-align:right;%s">%s</td>', style, x)
  th <- function(x, style = "") sprintf('<th style="padding:3px 12px;text-align:right;%s">%s</th>', style, x)
  red <- "background:#fde3e3;"; blue <- "background:#dbe9f6;"
  head <- paste0("<tr>", th("", "text-align:left;"),
                 paste0(mapply(function(l, u, j) th(sprintf('%s<br><span style="font-weight:normal;font-size:80%%;color:#666">%s</span>', l, u),
                                                     if (vars6[j] == "Murder") red else ""),
                               labs6, units, seq_along(vars6)), collapse = ""), "</tr>")
  row <- function(s) {
    t <- which(states$state == s); hl <- if (s == "Texas") blue else ""
    cells <- vapply(seq_along(vars6), function(j) td(fmt[[vars6[j]]](X6[t, j]),
                     paste0(hl, if (vars6[j] == "Murder") red else "")), "")
    paste0("<tr>", sprintf('<td style="padding:3px 12px;text-align:left;font-weight:bold;%s">%d %s</td>', hl, t, s),
           paste0(cells, collapse = ""), "</tr>")
  }
  dots <- paste0("<tr>", td("&#8942;", "text-align:left;"), paste0(rep(td("&#8942;"), 6), collapse = ""), "</tr>")
  first <- show[1:5]; rest <- show[-(1:5)]
  body <- c(vapply(first, row, ""), dots, row(rest[1]), dots, row(rest[2]))
  htmltools::HTML(sprintf('<table style="border-collapse:collapse;font-size:%dpx;margin:0 auto;font-variant-numeric:tabular-nums;">%s%s</table>',
                          font, head, paste0(body, collapse = "")))
}
fig_heat <- function() {
  Z6 <- scale(X6)
  heat_df <- as_tibble(Z6) |>
    mutate(abb = factor(states$abb, levels = states$abb)) |>
    pivot_longer(-abb, names_to = "variable", values_to = "z") |>
    mutate(variable = factor(variable, levels = rev(vars6), labels = rev(labs6)))
  row_t <- which(states$abb == "TX"); col_j <- 7 - which(vars6 == "Murder")
  ggplot(heat_df, aes(abb, variable, fill = z)) +
    geom_tile(colour = "white", linewidth = 0.2) +
    annotate("rect", ymin = 0.5, ymax = 6.5, xmin = row_t - 0.5, xmax = row_t + 0.5,
             fill = NA, colour = col_a1, linewidth = 1.1) +
    annotate("rect", ymin = col_j - 0.5, ymax = col_j + 0.5, xmin = 0.5, xmax = n_st + 0.5,
             fill = NA, colour = col_a2, linewidth = 1.1) +
    scale_fill_gradient2(low = "#b35806", mid = "white", high = "#542788",
                         limits = c(-3, 3), oob = scales::squish, name = "standardised\nvalue") +
    labs(x = NULL, y = NULL) +
    theme(axis.text.x = element_text(size = 8, angle = 90, vjust = 0.5), panel.grid = element_blank())
}
fig_dims <- function() {
  hl <- c("TX", "NV", "IA", "AL", "VT")
  tx <- states |> dplyr::filter(abb %in% hl)
  p_1d <- ggplot(states, aes(Murder, 0)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    geom_point(colour = "#6f7a8c", alpha = 0.6, size = 2.2) +
    geom_point(data = tx, colour = "black", size = 2.4) +
    geom_text_repel(data = tx, aes(label = abb), nudge_y = 0.3, size = 3.6, seed = 2,
                    direction = "x", min.segment.length = 0) +
    scale_y_continuous(limits = c(-0.5, 0.6), breaks = NULL) +
    labs(title = "p = 1: a line", subtitle = "Each state is a number", x = "murders per 100,000", y = NULL)
  p_2d <- ggplot(states, aes(Illiteracy, Murder)) +
    geom_point(colour = "#6f7a8c", alpha = 0.6, size = 2.2) +
    geom_point(data = tx, colour = "black", size = 2.4) +
    geom_text_repel(data = tx, aes(label = abb), size = 3.6, seed = 2, min.segment.length = 0) +
    labs(title = "p = 2: a plane", subtitle = "Each state is a point", x = "illiteracy (%)", y = "murders per 100,000")
  p_1d + p_2d + fig_dims_3d_static() + plot_layout(widths = c(1, 1, 1.15))
}
# A static 3D view drawn with ggplot: scale each variable to [0, 1], turn the cube, project onto the page.
fig_dims_3d_static <- function(azimuth = -35, elevation = 22) {
  hl <- c("TX", "NV", "IA", "AL", "VT")
  rng <- list(x = c(0.5, 3), y = c(1, 16), z = c(67.5, 74))          # illiteracy, murder, life expectancy
  unit <- function(v, r) (v - r[1]) / diff(r)
  a <- azimuth * pi / 180; e <- elevation * pi / 180
  proj <- function(x, y, z) {                                         # x, y, z in [0, 1]
    x <- x - 0.5; y <- y - 0.5; z <- z - 0.5
    u <- x * cos(a) - y * sin(a)                                      # turn around the vertical axis
    d <- x * sin(a) + y * cos(a)
    tibble(px = u, py = z * cos(e) + d * sin(e))                      # tilt and drop depth
  }
  P <- with(states, proj(unit(Illiteracy, rng$x), unit(Murder, rng$y), unit(LifeExp, rng$z))) |>
    mutate(abb = states$abb, h = abb %in% hl)
  Fl <- with(states, proj(unit(Illiteracy, rng$x), unit(Murder, rng$y), 0)) |> mutate(abb = states$abb)
  # floor and the two back walls (x = 0 and y = 1), like the corner of a room
  edges <- list(c(0,0,0, 1,0,0), c(1,0,0, 1,1,0), c(1,1,0, 0,1,0), c(0,1,0, 0,0,0),     # floor
                c(0,0,0, 0,0,1), c(0,1,0, 0,1,1), c(0,0,1, 0,1,1),                       # wall x = 0
                c(1,1,0, 1,1,1), c(0,1,1, 1,1,1))                                        # wall y = 1
  grid_f <- c(lapply(c(0.25, 0.5, 0.75), function(g) c(g,0,0, g,1,0)),
              lapply(c(0.25, 0.5, 0.75), function(g) c(0,g,0, 1,g,0)))
  seg <- function(L) purrr::map_dfr(L, function(v) {
    a0 <- proj(v[1], v[2], v[3]); a1 <- proj(v[4], v[5], v[6])
    tibble(x = a0$px, y = a0$py, xend = a1$px, yend = a1$py) })
  E <- seg(edges); G <- seg(grid_f)
  ax <- dplyr::bind_rows(proj(0.5, -0.2, 0) |> mutate(lab = "illiteracy"),
                         proj(1.22, 0.5, 0) |> mutate(lab = "murder"),
                         proj(-0.08, 1, 1.12) |> mutate(lab = "life exp."))
  hlP <- P |> dplyr::filter(h); hlF <- Fl |> dplyr::filter(abb %in% hl)
  ggplot() +
    geom_segment(data = G, aes(x, y, xend = xend, yend = yend), colour = "grey90") +
    geom_segment(data = E, aes(x, y, xend = xend, yend = yend), colour = "grey65") +
    geom_segment(data = tibble(x = hlP$px, y = hlP$py, xend = hlF$px, yend = hlF$py),
                 aes(x, y, xend = xend, yend = yend), colour = "grey45", linetype = 3) +
    geom_point(data = hlF, aes(px, py), colour = "grey60", size = 1.2) +
    geom_point(data = dplyr::filter(P, !h), aes(px, py), colour = "#6f7a8c", alpha = 0.6, size = 2.2) +
    geom_point(data = hlP, aes(px, py), colour = "black", size = 2.4) +
    geom_text_repel(data = hlP, aes(px, py, label = abb), size = 3.6, seed = 3, min.segment.length = 0.2,
                    box.padding = 0.35) +
    geom_text(data = ax, aes(px, py, label = lab), size = 3.6, colour = "grey30") +
    coord_equal(clip = "off") +
    labs(title = "p = 3: a point in space", subtitle = "Dotted: drop to the floor (the plane)", x = NULL, y = NULL) +
    theme(axis.text = element_blank(), panel.grid = element_blank())
}

# =============================================================================
# Types of variables
# =============================================================================

fig_types <- function() {
  p_nom <- ggplot(peng, aes(island)) +
    geom_bar(fill = "#6f7a8c") +
    labs(title = "Nominal", subtitle = "Island: labels,\nno order", x = NULL, y = "count")
  cut_med <- levels(ggplot2::diamonds$cut)[median(as.integer(ggplot2::diamonds$cut))]
  p_ord <- ggplot(ggplot2::diamonds, aes(cut, fill = cut == cut_med)) +
    geom_bar() +
    scale_fill_manual(values = c(`TRUE` = col_a1, `FALSE` = "#6f7a8c"), guide = "none") +
    scale_y_continuous(labels = scales::label_comma()) +
    labs(title = "Ordinal", subtitle = paste0("Diamond cut: ordered,\nmedian '", cut_med, "'"), x = NULL, y = NULL) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1))
  p_int <- ggplot(airquality, aes(Temp)) +
    geom_histogram(bins = 20, fill = "#6f7a8c", colour = "white") +
    labs(title = "Interval", subtitle = "Temperature (°F): zero\nis arbitrary", x = "°F (New York, 1973)", y = NULL)
  p_rat <- ggplot(peng, aes(body_mass_g)) +
    geom_histogram(bins = 25, fill = "#6f7a8c", colour = "white") +
    geom_vline(xintercept = mean(peng$body_mass_g), colour = col_a2, linewidth = 1) +
    labs(title = "Ratio", subtitle = "Body mass: true zero,\nmean (red) and ratios", x = "grams", y = NULL)
  p_nom + p_ord + p_int + p_rat + plot_layout(nrow = 1)
}

# The same scatterplot idea applied to different types of variables
# Euclidean distance = the length of a ruler laid between two points (Texas and California, zoomed in)
fig_euclid <- function() {
  a <- X6["California", c("Illiteracy", "Murder")]; b <- X6["Texas", c("Illiteracy", "Murder")]
  d <- sqrt(sum((a - b)^2))
  two <- tibble(state = c("California", "Texas"), x = c(a[1], b[1]), y = c(a[2], b[2]))
  ggplot() +
    annotate("segment", x = a[1], xend = b[1], y = a[2], yend = a[2], linetype = 2, colour = "grey40") +
    annotate("segment", x = b[1], xend = b[1], y = a[2], yend = b[2], linetype = 2, colour = "grey40") +
    annotate("segment", x = a[1], xend = b[1], y = a[2], yend = b[2], colour = col_a2, linewidth = 1.4) +
    annotate("text", x = (a[1] + b[1]) / 2, y = a[2] - 0.18, label = sprintf("%.1f (illiteracy)", b[1] - a[1]),
             colour = "grey30", size = 4.5) +
    annotate("text", x = b[1] + 0.08, y = (a[2] + b[2]) / 2, label = sprintf("%.1f (murder)", b[2] - a[2]),
             colour = "grey30", size = 4.5, hjust = 0) +
    annotate("text", x = (a[1] + b[1]) / 2 - 0.12, y = (a[2] + b[2]) / 2 + 0.2,
             label = sprintf("ruler: %.2f", d), colour = col_a2, size = 4.8, hjust = 1, fontface = "bold") +
    geom_point(data = two, aes(x, y), colour = "black", size = 3.2) +
    geom_text(data = two, aes(x, y, label = state), nudge_y = 0.3, size = 4.5) +
    coord_equal(xlim = c(0.3, 3.3), ylim = c(9.6, 12.8)) +
    labs(x = "illiteracy (%)", y = "murders per 100,000",
         title = "A ruler between two points",
         subtitle = "Zoomed in, both axes on the same scale")
}

# What each type of variable allows: summaries and plots (check = works, cross = meaningless, ~ = with care)
type_matrix <- local({
  cols <- c("count, mode", "median, ranks", "differences", "mean, variance", "covariance, correlation",
            "ratios (twice as much)", "bar chart of counts", "histogram, boxplot", "scatterplot axis",
            "colour hue, panels", "point size, colour gradient")
  group <- rep(c("Summaries", "Plots"), c(6, 5))
  rows <- c("Nominal (species)", "Ordinal (diamond cut)", "Interval (temperature °F)", "Ratio (body mass)")
  v <- rbind(c("y", "n", "n", "n", "n", "n",   "y", "n", "n", "y", "n"),
             c("y", "y", "n", "~", "~", "n",   "y", "~", "~", "y", "~"),
             c("y", "y", "y", "y", "y", "n",   "~", "y", "y", "~", "y"),
             c("y", "y", "y", "y", "y", "y",   "~", "y", "y", "~", "y"))
  expand_grid(r = seq_along(rows), k = seq_along(cols)) |>
    mutate(type = factor(rows[r], levels = rev(rows)), what = factor(cols[k], levels = cols),
           group = factor(group[k], levels = c("Summaries", "Plots")), v = v[cbind(r, k)])
})
fig_type_matrix <- function(base = 13) {
  sym <- c(y = "\u2713", n = "\u2717", `~` = "~")
  ggplot(type_matrix, aes(what, type)) +
    geom_tile(aes(fill = v), colour = "white", linewidth = 1.2) +
    geom_text(aes(label = sym[v], colour = v), size = base / 2.2, fontface = "bold") +
    scale_fill_manual(values = c(y = "#e3f1e6", n = "#f6e3e1", `~` = "#fbf1dc"), guide = "none") +
    scale_colour_manual(values = c(y = "#1b7a3a", n = "#b03a2e", `~` = "#a8730f"), guide = "none") +
    scale_x_discrete(position = "top") +
    facet_grid(~ group, scales = "free_x", space = "free_x") +
    labs(x = NULL, y = NULL, caption = "\u2713 works    ~ with care (few categories or values, binning)    \u2717 meaningless") +
    theme_minimal(base_size = base) +
    theme(axis.text.x.top = element_text(angle = 90, hjust = 0, vjust = 0.5), panel.grid = element_blank(),
          strip.placement = "outside", strip.text = element_text(face = "bold", size = base),
          axis.text.y = element_text(size = base), plot.caption = element_text(hjust = 0, colour = "grey30"))
}

fig_type_plots <- function() {
  cnt <- peng |> count(island, species)
  bad1 <- ggplot(peng, aes(island, species)) +
    geom_point(colour = "#6f7a8c", size = 3) +
    labs(title = "Scatterplot: two categorical",
         subtitle = sprintf("%d penguins, but only %d dots", nrow(peng), nrow(cnt)), x = "island", y = "species")
  good1 <- ggplot(cnt, aes(island, species)) +
    geom_tile(aes(fill = n), colour = "white") +
    geom_text(aes(label = n), size = 5) +
    scale_fill_gradient(low = "#eef4fa", high = col_a1, guide = "none") +
    labs(title = "Count table (heatmap)", subtitle = "Each combination shows its count", x = "island", y = "species")
  bad2 <- ggplot(peng, aes(species, body_mass_g)) +
    geom_point(colour = "#6f7a8c", size = 2, alpha = 0.6) +
    labs(title = "Scatterplot: categorical vs numeric", subtitle = "Points pile up on three lines",
         x = "species", y = "body mass (g)")
  good2 <- ggplot(peng, aes(species, body_mass_g)) +
    geom_boxplot(width = 0.5, outlier.shape = NA, fill = "grey95") +
    geom_jitter(width = 0.15, colour = "#6f7a8c", alpha = 0.5, size = 1.2) +
    labs(title = "Boxplot with jittered points", subtitle = "Centre, spread and every penguin",
         x = "species", y = "body mass (g)")
  (bad1 + good1) / (bad2 + good2)
}

# A static 3D scatterplot drawn with ggplot (variables scaled to [0, 1], cube turned and projected)
static_3d <- function(x, y, z, colour = NULL, azimuth = -35, elevation = 22,
                      labs3 = c("x", "y", "z"), title = NULL, subtitle = NULL) {
  u <- function(v) (v - min(v)) / diff(range(v))
  a <- azimuth * pi / 180; e <- elevation * pi / 180
  proj <- function(x, y, z) {
    x <- x - 0.5; y <- y - 0.5; z <- z - 0.5
    tibble(px = x * cos(a) - y * sin(a), py = z * cos(e) + (x * sin(a) + y * cos(a)) * sin(e))
  }
  P <- proj(u(x), u(y), u(z)) |> mutate(col = if (is.null(colour)) "all" else as.character(colour))
  edges <- list(c(0,0,0, 1,0,0), c(1,0,0, 1,1,0), c(1,1,0, 0,1,0), c(0,1,0, 0,0,0),
                c(0,0,0, 0,0,1), c(0,1,0, 0,1,1), c(0,0,1, 0,1,1), c(1,1,0, 1,1,1), c(0,1,1, 1,1,1))
  E <- purrr::map_dfr(edges, function(v) { a0 <- proj(v[1], v[2], v[3]); a1 <- proj(v[4], v[5], v[6])
    tibble(x = a0$px, y = a0$py, xend = a1$px, yend = a1$py) })
  corners <- list(c(0, 0), c(0, 1), c(1, 1))                  # vertical edges of the box
  top <- corners[[which.min(sapply(corners, function(cc) proj(cc[1], cc[2], 1)$px))]]
  ax <- dplyr::bind_rows(proj(0.5, -0.25, 0), proj(1.25, 0.5, 0), proj(top[1], top[2], 1.1)) |>
    mutate(lab = labs3, py = py - c(0.07, 0.07, 0))
  # drop the label of an axis that points (almost) straight at the viewer
  shown <- c(abs(cos(a)), abs(sin(a)), 1) > 0.2
  ax <- ax[shown, ]
  p <- ggplot() +
    geom_segment(data = E, aes(x, y, xend = xend, yend = yend), colour = "grey65") +
    geom_point(data = P, aes(px, py, colour = col), alpha = 0.6, size = 1.6) +
    geom_text(data = ax, aes(px, py, label = lab), size = 3.6, colour = "grey30") +
    coord_equal(clip = "off") +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme(axis.text = element_blank(), panel.grid = element_blank())
  if (is.null(colour)) p + scale_colour_manual(values = c(all = "#6f7a8c"), guide = "none")
  else p + scale_colour_manual(values = pal_species, name = NULL)
}

fig_3d_views <- function() {
  v <- function(az, el, t, s) static_3d(peng$bill_length_mm, peng$bill_depth_mm, peng$flipper_length_mm,
                                        azimuth = az, elevation = el, title = t, subtitle = s,
                                        labs3 = c("bill length", "bill depth", "flipper length"))
  v(-90, 0, "View 1", "Two groups?") + v(-20, 25, "View 2", "One blob?") + v(40, 60, "View 3", "Three groups?")
}

# The same three penguin variables in 2D: two axes plus a colour gradient for the third (numeric) variable
fig_3d_in_2d <- function() {
  ggplot(peng, aes(bill_length_mm, bill_depth_mm, colour = flipper_length_mm)) +
    geom_point(size = 2.2, alpha = 0.85) +
    scale_colour_viridis_c(option = "plasma", end = 0.92, name = "flipper\nlength (mm)") +
    labs(x = "bill length (mm)", y = "bill depth (mm)",
         title = "Two axes + a colour gradient", subtitle = "Three groups appear")
}

# Mahalanobis distances of all 50 states from the centre (six variables), largest first
fig_maha_rank <- function(k = 12, dist = sqrt(d2_st), xlab = "Mahalanobis distance (six variables)") {
  d <- tibble(state = names(dist), d = dist) |> slice_max(d, n = k) |>
    mutate(state = reorder(state, d))
  ggplot(d, aes(d, state)) +
    geom_col(fill = col_a1, width = 0.7) +
    geom_text(aes(label = sprintf("%.2f", d)), hjust = -0.15, size = 3.8) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.12))) +
    labs(x = xlab, y = NULL)
}

# =============================================================================
# Many variables and a third variable in 2D
# =============================================================================

fig_pairs <- function(base = 12) {
  ggpairs(states, columns = vars6,
          lower = list(continuous = wrap("points", alpha = 0.6, size = 1.1, colour = "#6f7a8c")),
          diag = list(continuous = wrap("densityDiag", colour = "#2c3e50")),
          upper = list(continuous = wrap("cor", size = 4, colour = "black", stars = FALSE)),
          columnLabels = labs6) +
    theme_minimal(base_size = base) +
    theme(panel.border = element_rect(fill = NA, colour = "grey85"), axis.text = element_text(size = 7),
          panel.grid.minor = element_blank())
}
fig_3d <- function(height = NULL) {
  p <- plot_ly(peng, x = ~bill_length_mm, y = ~bill_depth_mm, z = ~flipper_length_mm,
               color = ~species, colors = pal_species,
               type = "scatter3d", mode = "markers", marker = list(size = 3, opacity = 0.8),
               text = ~paste(species, island, sex, sep = ", "), hoverinfo = "text") |>
    layout(scene = list(xaxis = list(title = "bill length"), yaxis = list(title = "bill depth"),
                        zaxis = list(title = "flipper length"),
                        camera = list(eye = list(x = 1.6, y = -1.6, z = 0.8))),
           legend = list(x = 0.02, y = 0.95), margin = list(l = 0, r = 0, t = 10, b = 0))
  if (!is.null(height)) p$height <- height
  p
}
fig_simpson <- function() {
  p_grey <- ggplot(peng, aes(bill_length_mm, bill_depth_mm)) +
    geom_point(colour = "#6f7a8c", alpha = 0.6, size = 1.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.8) +
    labs(title = "Ignoring species", subtitle = sprintf("r = %.2f: longer bills are shallower?", r_all),
         x = "bill length (mm)", y = "bill depth (mm)")
  p_col <- ggplot(peng, aes(bill_length_mm, bill_depth_mm, colour = species)) +
    geom_point(alpha = 0.7, size = 1.8) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 0.8) +
    scale_colour_manual(values = pal_species) +
    labs(title = "Colour = species",
         subtitle = sprintf("Within each species r = %.2f to %.2f", min(r_within), max(r_within)),
         x = "bill length (mm)", y = "bill depth (mm)", colour = NULL) +
    theme(legend.position = c(0.84, 0.9), legend.background = element_rect(fill = "white", colour = NA))
  p_grey + p_col
}
# Bubble plot: a numeric third (and fourth) variable as point size (and fill colour)
fig_bubble <- function() {
  big <- states |> slice_max(Population, n = 6)
  ggplot(states, aes(Illiteracy, Murder)) +
    geom_point(aes(size = Population, fill = Income), shape = 21, colour = "grey30", alpha = 0.8) +
    geom_text_repel(data = big, aes(label = state), size = 3.6, seed = 3, box.padding = 0.8,
                    min.segment.length = 0) +
    scale_size_area(max_size = 16, breaks = c(1000, 5000, 20000),
                    labels = scales::label_comma(suffix = "k"), name = "population") +
    scale_fill_gradient2(low = "#b35806", mid = "#f7f7f7", high = "#542788",
                         midpoint = median(states$Income), name = "income ($)") +
    labs(x = "illiteracy (%)", y = "murders per 100,000")
}
# Facet plot, one grouping variable: facet_wrap() lays the panels out in a row (wrapping if needed)
fig_facet_wrap <- function() {
  ggplot(peng, aes(flipper_length_mm, body_mass_g, colour = species)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.6) +
    facet_wrap(~ species) +
    scale_colour_manual(values = pal_species, guide = "none") +
    labs(x = "flipper length (mm)", y = "body mass (g)") +
    theme(panel.border = element_rect(fill = NA, colour = "grey85"))
}
fig_trellis <- function() {
  ggplot(peng, aes(flipper_length_mm, body_mass_g, colour = species)) +
    geom_point(alpha = 0.7, size = 1.5) +
    geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "black", linewidth = 0.6) +
    facet_grid(sex ~ species) +
    scale_colour_manual(values = pal_species, guide = "none") +
    labs(x = "flipper length (mm)", y = "body mass (g)") +
    theme(panel.border = element_rect(fill = NA, colour = "grey85"))
}
# Trellis plot with the lattice package: the "|" in the formula reads "given"
fig_lattice <- function() {
  lattice::xyplot(body_mass_g ~ flipper_length_mm | species * sex, data = peng,
                  type = c("p", "r"), col = "#6f7a8c", pch = 16, cex = 0.6,
                  xlab = "flipper length (mm)", ylab = "body mass (g)")
}

# =============================================================================
# Summaries of the cloud
# =============================================================================

fig_quadrants <- function() {
  cen <- states |> mutate(ci = Illiteracy - xb2[1], cm = Murder - xb2[2],
                          prod = ci * cm, sign = ifelse(prod >= 0, "positive", "negative"))
  p_rect <- ggplot(cen) +
    geom_rect(aes(xmin = 0, xmax = ci, ymin = 0, ymax = cm, fill = sign), alpha = 0.12) +
    geom_hline(yintercept = 0, colour = "grey40") + geom_vline(xintercept = 0, colour = "grey40") +
    geom_point(aes(ci, cm, colour = sign), size = 2.2) +
    annotate("text", x = c(1.5, -0.9, -0.9, 1.5), y = c(8, 8, -6.5, -6.5), label = c("+", "-", "+", "-"),
             size = 9, colour = c(col_a1, col_a2, col_a1, col_a2)) +
    scale_fill_manual(values = c(positive = col_a1, negative = col_a2), guide = "none") +
    scale_colour_manual(values = c(positive = col_a1, negative = col_a2), guide = "none") +
    labs(title = "Centred data", subtitle = "Rectangle = product, sign by quadrant",
         x = "illiteracy - mean", y = "murder - mean")
  p_prod <- ggplot(cen, aes(reorder(abb, prod), prod, fill = sign)) +
    geom_col() +
    geom_hline(yintercept = s_im, linetype = 2) +
    annotate("label", x = 8, y = s_im + 0.2, label = sprintf("covariance = %.2f", s_im), size = 3.8, hjust = 0) +
    scale_fill_manual(values = c(positive = col_a1, negative = col_a2), guide = "none") +
    coord_flip() +
    labs(title = "The 50 products", subtitle = "Sum / (n - 1) = covariance", x = "states", y = "product") +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank())
  p_rect + p_prod + plot_layout(widths = c(1.2, 1))
}
fig_cor_heat <- function() {
  Rdf <- as_tibble(cor(X6), rownames = "row") |>
    pivot_longer(-row, names_to = "col", values_to = "r") |>
    mutate(row = factor(row, levels = rev(vars6), labels = rev(labs6)),
           col = factor(col, levels = vars6, labels = labs6))
  ggplot(Rdf, aes(col, row, fill = r)) +
    geom_tile(colour = "white") +
    geom_text(aes(label = sprintf("%.2f", r), colour = abs(r) > 0.6), size = 4.3) +
    scale_fill_gradient2(low = col_a2, mid = "white", high = col_a1, limits = c(-1, 1)) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    coord_equal() +
    labs(x = NULL, y = NULL) +
    theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 30, hjust = 1))
}
fig_cor_gallery <- function(nrow = 2) {
  make_cloud <- function(r) { M <- app_W %*% t(sym_root(sigma_of(1, 1, r))); tibble(x = M[, 1], y = M[, 2]) }
  w1 <- app_W[, 1]
  clouds <- list(
    "strong negative" = make_cloud(-0.9), "none" = make_cloud(0),
    "moderate positive" = make_cloud(0.5), "strong positive" = make_cloud(0.9),
    "a curve" = tibble(x = w1, y = w1^2 - 1 + 0.2 * app_W[, 2]),
    "one outlier" = tibble(x = c(app_W[1:199, 1] * 0.3, 4), y = c(app_W[1:199, 2] * 0.3, 4)))
  gal <- imap_dfr(clouds, ~ mutate(.x, panel = .y))
  lab_r <- imap_dfr(clouds, ~ tibble(panel = .y, r = cor(.x$x, .x$y)))
  gal <- gal |> left_join(lab_r, by = "panel") |>
    mutate(panel = sprintf("%s: r = %.2f", panel, r), panel = factor(panel, levels = unique(panel)))
  ells <- gal |> group_by(panel) |>
    group_modify(~ { e <- ellipse_xy(sym_root(cov(cbind(.x$x, .x$y))), 1)
                     tibble(x = e[1, ] + mean(.x$x), y = e[2, ] + mean(.x$y)) })
  ggplot(gal, aes(x, y)) +
    geom_point(colour = "#6f7a8c", alpha = 0.45, size = 0.9) +
    geom_path(data = ells, colour = "black", linewidth = 0.8) +
    geom_point(data = dplyr::filter(gal, x == 4 & y == 4), colour = col_a2, size = 2.4) +
    facet_wrap(~ panel, nrow = nrow) +
    coord_equal(xlim = c(-3.5, 4.3), ylim = c(-3.5, 4.3)) +
    labs(x = NULL, y = NULL) +
    theme(axis.text = element_blank(), panel.border = element_rect(fill = NA, colour = "grey85"),
          strip.text = element_text(size = 9))
}
fig_cor_slider <- function() {
  r_grid <- round(seq(-0.95, 0.95, by = 0.05), 2)
  rl <- function(r) sprintf("%.2f", r)
  frames_pts <- map_dfr(r_grid, function(r) {
    M <- app_W %*% t(sym_root(sigma_of(1, 1, r)))
    tibble(id = seq_len(nrow(M)), x = M[, 1], y = M[, 2], prod = M[, 1] * M[, 2], k = rl(r))
  })
  frames_ell <- map_dfr(r_grid, function(r) {
    e <- ellipse_xy(sym_root(sigma_of(1, 1, r)), 1); tibble(x = e[1, ], y = e[2, ], k = rl(r))
  })
  frames_txt <- map_dfr(r_grid, function(r) {
    pr <- frames_pts$prod[frames_pts$k == rl(r)]
    tibble(x = -3.4, y = 3.5, k = rl(r),
           text = sprintf("r = %.2f<br>positive products: %d<br>negative products: %d", r, sum(pr > 0), sum(pr < 0)))
  })
  cs <- list(c(0, "#d62728"), c(0.5, "#dddddd"), c(1, "#1f77b4"))
  p_sl <- plot_ly() |>
    add_trace(x = c(-4, 4), y = c(0, 0), type = "scatter", mode = "lines", line = list(color = "#999", width = 1),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(x = c(0, 0), y = c(-4, 4), type = "scatter", mode = "lines", line = list(color = "#999", width = 1),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(data = frames_pts, x = ~x, y = ~y, frame = ~k, ids = ~id, type = "scatter", mode = "markers",
              marker = list(size = 6, color = ~prod, colorscale = cs, cmin = -1.5, cmax = 1.5, opacity = 0.8,
                            line = list(width = 0.3, color = "#555")),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(data = frames_ell, x = ~x, y = ~y, frame = ~k, type = "scatter", mode = "lines",
              line = list(color = "black", width = 2.5), hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(data = frames_txt, x = ~x, y = ~y, frame = ~k, text = ~text, type = "scatter", mode = "text",
              textposition = "bottom right", textfont = list(size = 13, color = "#2c3e50"),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    layout(xaxis = list(range = c(-3.6, 3.6), title = "z<sub>1</sub>", zeroline = FALSE),
           yaxis = list(range = c(-3.6, 3.6), title = "z<sub>2</sub>", scaleanchor = "x", zeroline = FALSE),
           margin = list(l = 60, r = 20, t = 50, b = 130), showlegend = FALSE) |>
    animation_opts(frame = 120, transition = 0, redraw = TRUE) |>
    animation_slider(currentvalue = list(prefix = "r = ")) |>
    animation_button(x = 0, xanchor = "left", y = -0.3, yanchor = "top", label = "Play")
  jump <- lapply(c(-0.9, 0, 0.9), function(r) list(
    label = paste("r =", rl(r)), method = "animate",
    args = list(list(rl(r)), list(mode = "immediate", frame = list(duration = 0, redraw = TRUE),
                                  transition = list(duration = 0)))))
  finalize_animation(p_sl, jump, start_frame = "0.50", prefix = "r = ")
}

# =============================================================================
# Unusual combinations
# =============================================================================

fig_unusual <- function() {
  st2 <- states |> mutate(dm = dmU, hi = state == "Hawaii")
  z_long <- as_tibble(scale(X6[, vU])) |> mutate(state = states$state) |>
    pivot_longer(-state, names_to = "variable", values_to = "z") |>
    mutate(variable = factor(variable, levels = vU, labels = c("Illiteracy", "HS grad")))
  p_box <- ggplot(z_long, aes(variable, z)) +
    geom_boxplot(width = 0.45, outlier.shape = 1, fill = "grey95") +
    geom_point(data = dplyr::filter(z_long, state == "Hawaii"), colour = col_a2, size = 3.5) +
    annotate("text", x = 1.5, y = 1.55, label = "Hawaii", colour = col_a2, size = 4.5) +
    coord_cartesian(ylim = c(-2.5, 3)) +
    labs(title = "One at a time", subtitle = "Hawaii: not an outlier", x = NULL, y = "standardised value")
  ell_r <- map_dfr(1:3, function(r) { e <- ellipse_xy(sym_root(SU), r)
    tibble(x = e[1, ] + xbU[1], y = e[2, ] + xbU[2], r = factor(r)) })
  top <- dplyr::slice_max(st2, dm, n = 4)
  p_sc <- ggplot(st2, aes(Illiteracy, HSGrad)) +
    geom_path(data = ell_r, aes(x, y, group = r, linetype = r, colour = r), linewidth = 0.8) +
    geom_point(aes(colour = hi), size = 2.4, alpha = 0.8) +
    geom_label_repel(data = top, aes(label = sprintf("%s (%.2f)", state, dm)), size = 3.7, seed = 1,
                     fill = alpha("white", 0.85), label.size = 0, label.padding = 0.12,
                     min.segment.length = 0, box.padding = 0.6) +
    scale_colour_manual(values = c(`FALSE` = "#6f7a8c", `TRUE` = col_a2, "1" = "black", "2" = "grey45", "3" = "grey65"),
                        guide = "none") +
    scale_linetype_manual(values = c("1" = 1, "2" = 2, "3" = 3), guide = "none") +
    labs(title = "Together", subtitle = "Ellipses: Mahalanobis distance 1, 2, 3",
         x = "illiteracy (%)", y = "high-school graduates (%)")
  p_box + p_sc + plot_layout(widths = c(0.7, 1.3))
}
