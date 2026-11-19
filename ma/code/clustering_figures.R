# =============================================================================
# clustering_figures.R  --  data and figures for the chapter "Cluster Analysis" of the ISQS 6350 course notes.
# Put this file and _common.R (shared helpers) in one folder, set it as the working directory, and run
# source("clustering_figures.R"). Every figure is a function, e.g. fig_...().
# Packages: tidyverse, plotly, patchwork, MASS, mvtnorm, ggrepel, palmerpenguins, mclust
# =============================================================================

suppressPackageStartupMessages(library(mclust))   # loaded first, so that the map() of purrr (loaded by _common.R)
source("_common.R", local = TRUE)                  # comes before the map() of mclust
suppressPackageStartupMessages({
  library(ggrepel)
  library(palmerpenguins)
})

f0 <- function(x) formatC(x, format = "f", digits = 0, big.mark = ",")
f1 <- function(x) formatC(x, format = "f", digits = 1, big.mark = ",")
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
pct <- function(x, d = 0) paste0(formatC(100 * x, format = "f", digits = d), "%")
mns <- function(s) sub("^-", "\u2212", s)                # a proper minus sign for labels in figures
axes0 <- list(geom_hline(yintercept = 0, colour = "grey88"), geom_vline(xintercept = 0, colour = "grey88"))


# ---- Data -------------------------------------------------------------------
peng <- na.omit(penguins)                                 # 333 penguins with all four measurements and the sex
pal_species <- c(Adelie = "#1b9e77", Chinstrap = "#d95f02", Gentoo = "#7570b3")
pal_cl <- c("1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3", "4" = "#e7298a")  # clusters, as in the Introduction
bills <- as.matrix(peng[, c("bill_length_mm", "bill_depth_mm")])                 # bill length and depth, in mm
meas4 <- as.matrix(peng[, c("bill_length_mm", "bill_depth_mm", "flipper_length_mm", "body_mass_g")])
meas4_z <- scale(meas4)                                   # the four measurements, each standardised
peng_ag <- filter(peng, species != "Chinstrap") |> mutate(species = droplevels(species))   # Adelie and Gentoo
flip <- peng_ag$flipper_length_mm                         # the example with one variable: two species, two groups
lab_len <- "bill length (mm)"; lab_dep <- "bill depth (mm)"
bill_lims <- coord_equal(xlim = c(31, 61), ylim = c(13, 22), expand = FALSE)   # one scale on both axes: round looks round

# The US states (Chapter 2): six indicators
states <- as.data.frame(state.x77) |>
  rownames_to_column("state") |>
  rename(LifeExp = `Life Exp`, HSGrad = `HS Grad`) |>
  mutate(abb = state.abb)
vars6 <- c("Income", "Illiteracy", "LifeExp", "Murder", "HSGrad", "Frost")
X6 <- as.matrix(states[, vars6]); rownames(X6) <- states$abb
Z6 <- scale(X6)                                           # standardised, as for PCA with the correlation matrix

# Cluster labels are arbitrary. To compare clusters with known groups (such as the species), renumber the
# clusters one to one so that as many observations as possible carry the number of their group.
perms <- function(v) if (length(v) <= 1) list(v) else
  do.call(c, lapply(seq_along(v), function(i) lapply(perms(v[-i]), function(p) c(v[i], p))))
match_labels <- function(cl, ref) {
  ref <- as.integer(factor(ref)); k <- max(c(cl, ref))
  tab <- table(factor(cl, 1:k), factor(ref, 1:k))
  all_p <- perms(1:k)
  p <- all_p[[which.max(sapply(all_p, function(p) sum(tab[cbind(1:k, p)])))]]
  p[cl]                                                   # cluster l becomes number p[l]
}
agree <- function(cl, ref) mean(match_labels(cl, ref) == as.integer(factor(ref)))   # share in the matching group
cross_tab <- function(cl, ref) table(cluster = match_labels(cl, ref), group = ref)

# Ellipse holding a share `level` of a normal distribution with mean mu and covariance S (2 variables)
ell_df <- function(mu, S, level = 0.8) {
  e <- ellipse_xy(sym_root(S), sqrt(qchisq(level, 2)))
  tibble(x = mu[1] + e[1, ], y = mu[2] + e[2, ])
}


# =============================================================================
# 1  Groups behind the data
# =============================================================================

# The figure of the Introduction: a continuous factor, and clusters on and off the factor's line
fig_factor_to_cluster <- function() {
  set.seed(11)
  lam_b <- c(1.2, 0.8); psi_b <- 0.22
  circle_at <- function(mu, r = sqrt(psi_b) * sqrt(qchisq(0.8, 2)))
    tibble(x1 = mu[1] + r * cos(seq(0, 2 * pi, length.out = 120)),
           x2 = mu[2] + r * sin(seq(0, 2 * pi, length.out = 120)))
  f_draw <- rnorm(240)
  cont_pts <- tibble(x1 = lam_b[1] * f_draw + rnorm(240, sd = sqrt(psi_b)),
                     x2 = lam_b[2] * f_draw + rnorm(240, sd = sqrt(psi_b)), f = f_draw)
  cont_circ <- purrr::map_dfr(c(-1.5, -0.75, 0, 0.75, 1.5), function(f)
    circle_at(lam_b * f) |> mutate(g = f))
  on_centers  <- tibble(g = c("1", "2", "3"), x1 = lam_b[1] * c(-1.5, 0, 1.5), x2 = lam_b[2] * c(-1.5, 0, 1.5))
  off_centers <- tibble(g = c("1", "2", "3"), x1 = c(-1.8, 1.5, 1.3), x2 = c(0.3, 1.5, -1.6))
  draw_clusters <- function(centers, n = 80) purrr::map_dfr(seq_len(nrow(centers)), function(i)
    tibble(g = centers$g[i],
           x1 = centers$x1[i] + rnorm(n, sd = sqrt(psi_b)),
           x2 = centers$x2[i] + rnorm(n, sd = sqrt(psi_b))))
  on_pts  <- draw_clusters(on_centers)
  off_pts <- draw_clusters(off_centers)
  circles_for <- function(centers) purrr::map_dfr(seq_len(nrow(centers)), function(i)
    circle_at(c(centers$x1[i], centers$x2[i])) |> mutate(g = centers$g[i]))
  factor_line <- geom_line(data = tibble(x1 = lam_b[1] * c(-2.4, 2.4), x2 = lam_b[2] * c(-2.4, 2.4)),
                           aes(x1, x2), colour = col_f, linetype = 2)
  pal3 <- c("1" = "#1b9e77", "2" = "#d95f02", "3" = "#7570b3")
  lims <- coord_equal(xlim = c(-3.2, 3.2), ylim = c(-2.6, 2.6))
  p1 <- ggplot() + factor_line +
    geom_point(data = cont_pts, aes(x1, x2, colour = f), size = 1.3, alpha = 0.8) +
    geom_path(data = cont_circ, aes(x1, x2, group = g), colour = "grey25", linewidth = 0.6) +
    scale_colour_gradient2(low = "#b35806", mid = "#d4d4d4", high = "#542788", guide = "none") +
    lims + labs(title = "Continuous factor", subtitle = "Centers λf slide\nalong the line", x = expression(italic(x)[1]), y = expression(italic(x)[2]))
  p2 <- ggplot() + factor_line +
    geom_point(data = on_pts, aes(x1, x2, colour = g), size = 1.3, alpha = 0.8) +
    geom_path(data = circles_for(on_centers), aes(x1, x2, group = g), colour = "grey25", linewidth = 0.8) +
    scale_colour_manual(values = pal3, guide = "none") +
    lims + labs(title = "Clusters on the line", subtitle = "The center jumps\namong three points", x = expression(italic(x)[1]), y = expression(italic(x)[2]))
  p3 <- ggplot() +
    geom_point(data = off_pts, aes(x1, x2, colour = g), size = 1.3, alpha = 0.8) +
    geom_path(data = circles_for(off_centers), aes(x1, x2, group = g), colour = "grey25", linewidth = 0.8) +
    scale_colour_manual(values = pal3, guide = "none") +
    lims + labs(title = "Clusters anywhere", subtitle = "No single line\nexplains the cloud", x = expression(italic(x)[1]), y = expression(italic(x)[2]))
  p1 + p2 + p3
}

# The figure of Chapter 3: a kernel density estimate of bill length and depth, with three hills
peng_kde <- local({
  k <- MASS::kde2d(peng$bill_length_mm, peng$bill_depth_mm, n = 120,
                   lims = c(30, 62, 12.5, 22.5))
  list(x = k$x, y = k$y, z = k$z,
       df = expand_grid(y = k$y, x = k$x) |> mutate(z = as.vector(t(k$z))) |> dplyr::select(x, y, z))
})
fig_kde_contour <- function(colour_species = FALSE) {
  p <- ggplot() +
    geom_contour_filled(data = peng_kde$df, aes(x, y, z = z), bins = 9, alpha = 0.85) +
    scale_fill_manual(values = colorRampPalette(c("white", "#c6dbef", "#6baed6", "#2171b5", "#08306b"))(9),
                      guide = "none")
  p <- if (colour_species) {
    p + geom_point(data = peng, aes(bill_length_mm, bill_depth_mm, colour = species), size = 1.2, alpha = 0.8) +
      scale_colour_manual(values = pal_species, name = NULL) +
      theme(legend.position = c(0.86, 0.9), legend.background = element_rect(fill = "white", colour = NA))
  } else {
    p + geom_point(data = peng, aes(bill_length_mm, bill_depth_mm), colour = "grey20", size = 0.9, alpha = 0.6)
  }
  p + coord_cartesian(xlim = c(31, 61), ylim = c(13, 22), expand = FALSE) +
    labs(x = "bill length (mm)", y = "bill depth (mm)")
}


# =============================================================================
# 2  Hierarchical clustering
# =============================================================================

# Ten penguins (four Adelie, three Chinstrap and three Gentoo), numbered from the shortest bill to the longest
hc_rows <- c(72, 138, 56, 54, 273, 321, 326, 154, 159, 243)
hc_rows <- hc_rows[order(bills[hc_rows, 1])]
hc_X <- bills[hc_rows, ]; rownames(hc_X) <- seq_along(hc_rows)
hc_small <- hclust(dist(hc_X), method = "average")

# Line segments of a dendrogram. Leaf s sits at position match(s, hc$order) on the x-axis; merge i joins two
# branches (a negative number is a single observation, a positive one an earlier merge) at the height hc$height[i].
dendro_segments <- function(hc) {
  n <- length(hc$order); pos <- match(seq_len(n), hc$order); mid <- numeric(n - 1)
  x_of <- function(m) if (m < 0) pos[-m] else mid[m]
  h_of <- function(m) if (m < 0) 0 else hc$height[m]
  out <- vector("list", n - 1)
  for (i in seq_len(n - 1)) {
    a <- hc$merge[i, 1]; b <- hc$merge[i, 2]
    mid[i] <- (x_of(a) + x_of(b)) / 2
    out[[i]] <- tibble(step = i, x = c(x_of(a), x_of(b), x_of(a)), xend = c(x_of(a), x_of(b), x_of(b)),
                       y = c(h_of(a), h_of(b), hc$height[i]), yend = rep(hc$height[i], 3))
  }
  bind_rows(out)
}
# The observations in the cluster that merge i creates
merge_members <- function(hc) {
  mem <- vector("list", nrow(hc$merge))
  for (i in seq_len(nrow(hc$merge)))
    mem[[i]] <- unlist(lapply(hc$merge[i, ], function(m) if (m < 0) -m else mem[[m]]))
  mem
}
# A rounded outline around a set of points: the convex hull of small circles around each point. The radius r is a
# share of the axis ranges sx and sy, so that the circles look round on the screen.
blob <- function(P, r, sx, sy) {
  a <- seq(0, 2 * pi, length.out = 25)[-25]
  Q <- do.call(rbind, lapply(seq_len(nrow(P)), function(i) cbind(P[i, 1] + r * sx * cos(a), P[i, 2] + r * sy * sin(a))))
  h <- chull(Q); Q[c(h, h[1]), , drop = FALSE]
}
as_path <- function(d) tibble(x = as.vector(rbind(d$x, d$xend, NA)), y = as.vector(rbind(d$y, d$yend, NA)))

# Interactive: step through the merges of the ten penguins (left) and watch the dendrogram grow (right)
fig_hc_small <- function() {
  hc <- hc_small; n <- nrow(hc_X); mem <- merge_members(hc); seg <- dendro_segments(hc)
  xr <- c(35, 54.5); yr <- c(12.5, 22.5); col_new <- "#e6ab02"
  steps <- 0:(n - 1)
  cut_at <- function(s) if (s == 0) hc$height[1] / 2 else if (s < n - 1) (hc$height[s] + hc$height[s + 1]) / 2 else
    hc$height[n - 1] + 0.6
  empty <- function(s) tibble(x = -100, y = -100, s = s)     # a point outside the plot, so that every frame has data
  blobs <- purrr::map_dfr(steps, function(s) {
    cl <- cutree(hc, k = n - s)
    purrr::map_dfr(unique(cl), function(g) {
      B <- blob(hc_X[cl == g, , drop = FALSE], 0.04, diff(xr), diff(yr))
      tibble(x = c(B[, 1], NA), y = c(B[, 2], NA))
    }) |> mutate(s = s)
  })
  newest <- purrr::map_dfr(steps, function(s) if (s == 0) empty(0) else {
    B <- blob(hc_X[mem[[s]], , drop = FALSE], 0.04, diff(xr), diff(yr)); tibble(x = B[, 1], y = B[, 2], s = s) })
  done <- purrr::map_dfr(steps, function(s) if (s == 0) empty(0) else as_path(filter(seg, step <= s)) |> mutate(s = s))
  last <- purrr::map_dfr(steps, function(s) if (s == 0) empty(0) else as_path(filter(seg, step == s)) |> mutate(s = s))
  cut_line <- purrr::map_dfr(steps, function(s) tibble(x = c(0.5, n + 0.5), y = cut_at(s), s = s))
  info <- tibble(x = 0.7, y = 11.9, s = steps,
                 txt = ifelse(steps == 0, "10 clusters: each penguin alone",
                              sprintf("%d %s, merge %d at height %s mm", n - steps, ifelse(steps == n - 1, "cluster", "clusters"),
                                      steps, f2(hc$height[pmax(steps, 1)]))))
  pts <- tibble(x = hc_X[, 1], y = hc_X[, 2], lab = rownames(hc_X),
                pos = ifelse(lab == "1", "top left", ifelse(lab == "2", "bottom right", "top right")))
  scatter <- plot_ly() |>
    add_trace(data = blobs, x = ~x, y = ~y, frame = ~s, type = "scatter", mode = "lines", fill = "toself",
              fillcolor = "rgba(120,120,120,0.16)", line = list(color = "#9a9a9a", width = 1),
              hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = newest, x = ~x, y = ~y, frame = ~s, type = "scatter", mode = "lines", fill = "toself",
              fillcolor = "rgba(230,171,2,0.30)", line = list(color = col_new, width = 2),
              hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = pts, x = ~x, y = ~y, text = ~lab, type = "scatter", mode = "markers+text",
              textposition = ~pos, marker = list(color = "#333333", size = 7),
              textfont = list(size = 12, color = "#333333"), hoverinfo = "skip", showlegend = FALSE)
  dendro <- plot_ly() |>
    add_trace(data = as_path(seg), x = ~x, y = ~y, type = "scatter", mode = "lines",
              line = list(color = "#dddddd", width = 1.5), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = done, x = ~x, y = ~y, frame = ~s, type = "scatter", mode = "lines",
              line = list(color = "#333333", width = 2), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = last, x = ~x, y = ~y, frame = ~s, type = "scatter", mode = "lines",
              line = list(color = col_new, width = 4), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = cut_line, x = ~x, y = ~y, frame = ~s, type = "scatter", mode = "lines",
              line = list(color = col_a2, width = 1.5, dash = "dash"), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = info, x = ~x, y = ~y, frame = ~s, text = ~txt, type = "scatter", mode = "text",
              textposition = "middle right", textfont = list(size = 12, color = "#2c3e50"),
              hoverinfo = "skip", showlegend = FALSE)
  jump <- function(lab, s) list(label = lab, method = "animate",
                                args = list(list(as.character(s)), list(mode = "immediate",
                                            frame = list(duration = 0, redraw = TRUE), transition = list(duration = 0))))
  subplot(scatter, dendro, nrows = 1, widths = c(0.5, 0.5), margin = 0.07, titleX = TRUE, titleY = TRUE) |>
    layout(xaxis = list(title = lab_len, range = xr, zeroline = FALSE),
           yaxis = list(title = lab_dep, range = yr, zeroline = FALSE),
           xaxis2 = list(title = "penguin", tickvals = seq_len(n), ticktext = rownames(hc_X)[hc$order],
                         range = c(0.4, n + 0.6), zeroline = FALSE, showgrid = FALSE),
           yaxis2 = list(title = "distance between clusters (mm)", range = c(0, 12.4), zeroline = FALSE),
           margin = list(t = 50)) |>
    animation_opts(frame = 700, transition = 0, redraw = TRUE) |>
    animation_slider(currentvalue = list(prefix = "merges done: ")) |>
    animation_button(label = "Play") |>
    finalize_animation(list(jump("10 clusters", 0), jump("3 clusters", n - 3), jump("1 cluster", n - 1)),
                       start_frame = "0", prefix = "merges done: ")
}

# Two small clusters G and H, and four ways to measure the distance between them
link_G <- rbind(c(1.0, 1.3), c(1.7, 2.4), c(2.3, 1.5), c(1.4, 0.5))
link_H <- rbind(c(4.4, 2.1), c(5.1, 2.9), c(5.7, 1.8), c(4.9, 1.0), c(4.0, 0.9))
fig_linkage <- function(ncol = 4, base = 13) {
  pairs <- expand_grid(s = 1:nrow(link_G), t = 1:nrow(link_H)) |>
    mutate(x = link_G[s, 1], y = link_G[s, 2], xend = link_H[t, 1], yend = link_H[t, 2],
           d = sqrt((x - xend)^2 + (y - yend)^2))
  pts <- bind_rows(tibble(x = link_G[, 1], y = link_G[, 2], g = "G"), tibble(x = link_H[, 1], y = link_H[, 2], g = "H"))
  mG <- colMeans(link_G); mH <- colMeans(link_H)
  base_plot <- function(title, subtitle) ggplot() +
    geom_point(data = pts, aes(x, y, colour = g), size = 2.6) +
    annotate("text", x = c(0.55, 6.15), y = c(2.75, 2.75), label = c("G", "H"), size = 5, fontface = "bold",
             colour = pal_cl[1:2]) +
    scale_colour_manual(values = c(G = unname(pal_cl[1]), H = unname(pal_cl[2])), guide = "none") +
    coord_equal(xlim = c(0.3, 6.4), ylim = c(0.1, 3.3)) +
    labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    theme_minimal(base_size = base) +
    theme(axis.text = element_blank(), panel.grid = element_blank(), plot.title = element_text(face = "bold"))
  p1 <- base_plot("Single", "the closest pair") +
    geom_segment(data = pairs, aes(x, y, xend = xend, yend = yend), colour = "grey85", linewidth = 0.3) +
    geom_segment(data = slice_min(pairs, d), aes(x, y, xend = xend, yend = yend), colour = "grey10", linewidth = 1.2)
  p2 <- base_plot("Complete", "the farthest pair") +
    geom_segment(data = pairs, aes(x, y, xend = xend, yend = yend), colour = "grey85", linewidth = 0.3) +
    geom_segment(data = slice_max(pairs, d), aes(x, y, xend = xend, yend = yend), colour = "grey10", linewidth = 1.2)
  p3 <- base_plot("Average", sprintf("the mean of all %d pairs", nrow(pairs))) +
    geom_segment(data = pairs, aes(x, y, xend = xend, yend = yend), colour = "grey35", linewidth = 0.45)
  p4 <- base_plot("Ward", "how much W would rise") +
    geom_segment(aes(x = mG[1], y = mG[2], xend = mH[1], yend = mH[2]), colour = "grey10", linewidth = 1.2) +
    annotate("point", x = c(mG[1], mH[1]), y = c(mG[2], mH[2]), shape = 4, size = 5, stroke = 1.6,
             colour = pal_cl[1:2])
  # put the layers of points on top
  lift <- function(p) { p$layers <- c(p$layers[-1], p$layers[1]); p }
  wrap_plots(lift(p1), lift(p2), lift(p3), lift(p4), ncol = ncol)
}

# All 333 penguins: dendrograms of the four standardised measurements under four linkages, the species below
hc4 <- lapply(c(single = "single", complete = "complete", average = "average", Ward = "ward.D2"),
              function(m) hclust(dist(meas4_z), method = m))
fig_dendros <- function(methods = names(hc4), ncol = 2, base = 13) {
  # the height of a merge: the distance between the two clusters it joins, as each linkage measures it
  ylab <- list(single = "distance of the closest pair", complete = "distance of the farthest pair",
               average = "average distance of all pairs", Ward = quote(sqrt(2 %*% "rise in" ~ italic(W))))
  panels <- lapply(methods, function(m) {
    hc <- hc4[[m]]; n <- length(hc$order); top <- max(hc$height)
    sizes <- sort(table(cutree(hc, 3)), decreasing = TRUE)
    h3 <- mean(rev(hc$height)[2:3])                     # a cut between the merges that leave three and two clusters
    runs <- rle(as.character(peng$species[hc$order]))          # one bar for each run of the same species
    strip <- tibble(species = runs$values, xmax = cumsum(runs$lengths) + 0.5, xmin = xmax - runs$lengths)
    ggplot() +
      geom_segment(data = dendro_segments(hc), aes(x, y, xend = xend, yend = yend), linewidth = 0.25, colour = "grey25") +
      geom_rect(data = strip, aes(xmin = xmin, xmax = xmax, ymin = -0.115 * top, ymax = -0.025 * top, fill = species)) +
      geom_hline(yintercept = h3, linetype = 2, colour = col_a2, linewidth = 0.6) +
      scale_fill_manual(values = pal_species, name = NULL) +
      scale_x_continuous(expand = c(0.01, 0)) +
      labs(title = paste(m, "linkage"), subtitle = paste("sizes of the three clusters:", paste(sizes, collapse = ", ")),
           x = NULL, y = ylab[[m]]) +
      theme_minimal(base_size = base) +
      theme(axis.text.x = element_blank(), panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
            plot.title = element_text(face = "bold"))
  })
  wrap_plots(panels, ncol = ncol, guides = "collect") & theme(legend.position = "bottom")
}


# =============================================================================
# 3  k-means
# =============================================================================

# k-means of the bill measurements in mm, renumbered to match the species
km_bills <- local({
  set.seed(1); km <- kmeans(bills, centers = 3, nstart = 50)
  cl <- match_labels(km$cluster, peng$species)
  list(cluster = cl, centers = t(sapply(1:3, function(l) colMeans(bills[cl == l, ]))), W = km$tot.withinss)
})
km_s2 <- km_bills$W / (2 * nrow(bills))   # variance within the clusters, per measurement: the round groups of k-means
km_W <- sapply(1:6, function(k) { set.seed(1); kmeans(bills, centers = k, nstart = 50)$tot.withinss })
km_bills_z <- local({ set.seed(1); kmeans(scale(bills), centers = 3, nstart = 50)$cluster })   # standardised
# Lloyd's algorithm from 1,000 random starts (three penguins each time): the W at which each start ends
km_starts <- local({ set.seed(2); replicate(1000, kmeans(bills, centers = 3, algorithm = "Lloyd", iter.max = 50)$tot.withinss) })

nearest <- function(X, C) max.col(-sapply(seq_len(nrow(C)), function(l) colSums((t(X) - C[l, ])^2)), ties.method = "first")
W_of <- function(X, C, cl) sum((X - C[cl, ])^2)

# Lloyd's algorithm, step by step, from three Adelie penguins as starting centres
lloyd_start <- bills[c(79, 61, 98), ]
lloyd_path <- local({
  C <- lloyd_start; cl_old <- rep(0L, nrow(bills)); out <- list(list(C = C, cl = NULL, what = "start"))
  repeat {
    cl <- nearest(bills, C)
    if (all(cl == cl_old)) break
    out[[length(out) + 1]] <- list(C = C, cl = cl, what = "assign")
    C <- t(sapply(1:3, function(l) colMeans(bills[cl == l, , drop = FALSE])))
    out[[length(out) + 1]] <- list(C = C, cl = cl, what = "update")
    cl_old <- cl
  }
  out[[length(out) + 1]] <- list(C = C, cl = cl, what = "done")
  p <- match_labels(cl, peng$species)[match(1:3, cl)]   # renumber all steps as the final clusters match the species
  lapply(out, function(o) { o$C <- o$C[order(p), ]; if (!is.null(o$cl)) o$cl <- p[o$cl]; o })
})

# Interactive: Lloyd's algorithm. Shading: the region nearest to each centre.
fig_kmeans_anim <- function() {
  gx <- seq(31.125, 60.875, by = 0.25); gy <- seq(13.0625, 21.9375, by = 0.125)
  grid <- expand_grid(y = gy, x = gx)
  steps <- seq_along(lloyd_path) - 1
  shade <- purrr::map_dfr(steps, function(s) {
    C <- lloyd_path[[s + 1]]$C; tibble(x = grid$x, y = grid$y, z = nearest(as.matrix(grid[, c("x", "y")]), C), s = s)
  })
  pts <- purrr::map_dfr(steps, function(s) {
    o <- lloyd_path[[s + 1]]
    tibble(x = bills[, 1], y = bills[, 2], col = if (is.null(o$cl)) "#555555" else unname(pal_cl[o$cl]), s = s)
  })
  ctr <- purrr::map_dfr(steps, function(s) {
    C <- lloyd_path[[s + 1]]$C; tibble(x = C[, 1], y = C[, 2], col = unname(pal_cl[1:3]), s = s)
  })
  round_no <- cumsum(vapply(lloyd_path, function(o) o$what == "assign", logical(1)))
  info <- purrr::map_dfr(steps, function(s) {
    o <- lloyd_path[[s + 1]]
    txt <- switch(o$what,
      start  = "Start: three penguins serve as centres.",
      assign = sprintf("Round %d, assign: each penguin joins its nearest centre.  W = %s mm²", round_no[s + 1],
                       f0(W_of(bills, o$C, o$cl))),
      update = sprintf("Round %d, update: each centre moves to the mean of its penguins.  W = %s mm²", round_no[s + 1],
                       f0(W_of(bills, o$C, o$cl))),
      done   = sprintf("No penguin changes its centre: done.  W = %s mm²", f0(W_of(bills, o$C, o$cl))))
    tibble(x = 31.6, y = 21.5, txt = txt, s = s)
  })
  cs <- list(c(0, pal_cl[[1]]), c(1/3, pal_cl[[1]]), c(1/3, pal_cl[[2]]), c(2/3, pal_cl[[2]]), c(2/3, pal_cl[[3]]),
             c(1, pal_cl[[3]]))
  jump <- function(lab, s) list(label = lab, method = "animate",
                                args = list(list(as.character(s)), list(mode = "immediate",
                                            frame = list(duration = 0, redraw = TRUE), transition = list(duration = 0))))
  plot_ly() |>
    add_trace(data = shade, x = ~x, y = ~y, z = ~z, frame = ~s, type = "heatmap", colorscale = cs, zmin = 1, zmax = 3,
              opacity = 0.17, showscale = FALSE, hoverinfo = "skip") |>
    add_trace(data = pts, x = ~x, y = ~y, frame = ~s, type = "scatter", mode = "markers",
              marker = list(color = ~col, size = 6, opacity = 0.85), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = ctr, x = ~x, y = ~y, frame = ~s, type = "scatter", mode = "markers",
              marker = list(color = ~col, size = 17, symbol = "x", line = list(color = "black", width = 1.5)),
              hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = info, x = ~x, y = ~y, frame = ~s, text = ~txt, type = "scatter", mode = "text",
              textposition = "middle right", textfont = list(size = 13, color = "#2c3e50"),
              hoverinfo = "skip", showlegend = FALSE) |>
    layout(xaxis = list(title = lab_len, range = c(31, 61), zeroline = FALSE),
           yaxis = list(title = lab_dep, range = c(13, 22), zeroline = FALSE),
           margin = list(t = 50)) |>
    animation_opts(frame = 900, transition = 0, redraw = TRUE) |>
    animation_slider(currentvalue = list(prefix = "step ")) |>
    animation_button(label = "Play") |>
    finalize_animation(list(jump("start", 0), jump("end", max(steps))), start_frame = "0", prefix = "step ") |>
    slim_frames()
}
# Each frame of an animation repeats the coordinates of the shading grid and of the points, which never move.
# Keep in the frames only what changes (the shading and the colours), so that the page stays small.
slim_frames <- function(w) {
  w$x$frames <- lapply(w$x$frames, function(f) {
    f$data <- lapply(f$data, function(d) {
      if (identical(d$type, "heatmap")) d[c("type", "z")]
      else if (!is.null(d$marker$size) && d$marker$size == 6) { d$x <- NULL; d$y <- NULL; d } else d
    })
    f
  })
  w
}

# Within-cluster sum of squares W against the number of clusters k
fig_wss <- function(W = km_W, unit = "mm²", base = 14) {
  d <- tibble(k = seq_along(W), W = W)
  ggplot(d, aes(k, W)) + geom_line(colour = "grey50") + geom_point(size = 2.6, colour = col_a1) +
    scale_x_continuous(breaks = d$k) + scale_y_continuous(labels = function(x) f0(x), limits = c(0, NA)) +
    labs(x = "number of clusters k", y = paste0("W (", unit, ")")) + theme_minimal(base_size = base)
}


# =============================================================================
# 4  Model-based clustering
# =============================================================================

# One variable: the flipper length of the Adelie and Gentoo penguins, a mixture of two normal densities
mix1 <- Mclust(flip, G = 2, modelNames = "V", verbose = FALSE)
mix1_par <- with(mix1$parameters, list(pi = pro, mu = unname(mean), sd = sqrt(variance$sigmasq)))
mix1_agree <- agree(mix1$classification, peng_ag$species)            # most likely group = species
dens_l <- function(x, l, par = mix1_par) par$pi[l] * dnorm(x, par$mu[l], par$sd[l])   # pi_l f_l(x)
post2 <- function(x, par = mix1_par) dens_l(x, 2, par) / (dens_l(x, 1, par) + dens_l(x, 2, par))   # P(group 2 | x)
x_half <- uniroot(function(x) post2(x) - 0.5, c(195, 215))$root                    # where the two are equally likely

fig_mixture_1d <- function(base = 14) {
  xs <- seq(165, 240, by = 0.25)
  curves <- bind_rows(tibble(x = xs, y = dens_l(xs, 1), g = "1"), tibble(x = xs, y = dens_l(xs, 2), g = "2"))
  tot <- tibble(x = xs, y = dens_l(xs, 1) + dens_l(xs, 2))
  top <- ggplot() +
    geom_histogram(data = peng_ag, aes(flipper_length_mm, after_stat(density)), binwidth = 2, boundary = 0,
                   fill = "grey88", colour = "white") +
    geom_area(data = curves, aes(x, y, fill = g), position = "identity", alpha = 0.25) +
    geom_line(data = curves, aes(x, y, colour = g), linewidth = 1) +
    geom_line(data = tot, aes(x, y), linewidth = 0.9, colour = "black") +
    annotate("text", x = mix1_par$mu[1] - 6, y = 1.06 * max(dens_l(xs, 1)), label = "π[1]*f[1](x)", parse = TRUE,
             colour = pal_cl[1], size = 4.8, hjust = 1) +
    annotate("text", x = mix1_par$mu[2] + 6, y = 1.06 * max(dens_l(xs, 2)), label = "π[2]*f[2](x)", parse = TRUE,
             colour = pal_cl[3], size = 4.8, hjust = 0) +
    annotate("text", x = 232, y = 0.6 * max(tot$y), label = "f(x)", parse = TRUE, size = 4.8) +
    annotate("segment", x = 229, xend = 222, y = 0.6 * max(tot$y), yend = 0.4 * max(tot$y), linewidth = 0.4) +
    scale_fill_manual(values = c("1" = unname(pal_cl[1]), "2" = unname(pal_cl[3])), guide = "none") +
    scale_colour_manual(values = c("1" = unname(pal_cl[1]), "2" = unname(pal_cl[3])), guide = "none") +
    coord_cartesian(xlim = c(168, 237)) +
    labs(x = NULL, y = "density") + theme_minimal(base_size = base)
  bottom <- ggplot() +
    geom_hline(yintercept = c(0, 1), colour = "grey85") +
    geom_line(data = tibble(x = xs, y = post2(xs)), aes(x, y), colour = pal_cl[3], linewidth = 1.1) +
    geom_rug(data = peng_ag, aes(x = flipper_length_mm, colour = species), alpha = 0.5, length = unit(0.07, "npc")) +
    annotate("segment", x = x_half, xend = x_half, y = 0, yend = 0.5, linetype = 2, colour = "grey40") +
    annotate("text", x = x_half + 1, y = 0.42, label = paste0("equally likely at ", f0(x_half), " mm"), hjust = 0,
             size = 3.9, colour = "grey30") +
    scale_colour_manual(values = pal_species, name = NULL) +
    scale_y_continuous(breaks = c(0, 0.5, 1)) +
    coord_cartesian(xlim = c(168, 237)) +
    labs(x = "flipper length (mm)", y = "P(group 2 | x)") +
    theme_minimal(base_size = base) + theme(legend.position = "bottom")
  top / bottom + plot_layout(heights = c(2, 1.15))
}

# The EM algorithm for the two flipper groups, round by round, from a vague start: two wide, overlapping groups
em_rounds <- 40
em_path <- local({
  par <- list(pi = c(0.5, 0.5), mu = c(200, 206), sd = c(12, 12)); out <- list(par)
  for (r in seq_len(em_rounds)) {
    p2 <- post2(flip, par); w <- cbind(1 - p2, p2)                  # E-step: membership probabilities
    mu <- colSums(w * flip) / colSums(w)                              # M-step: weighted shares, means and
    dev2 <- (flip - matrix(mu, length(flip), 2, byrow = TRUE))^2      # standard deviations
    par <- list(pi = colMeans(w), mu = mu, sd = sqrt(colSums(w * dev2) / colSums(w)))
    out[[r + 1]] <- par
  }
  out
})
em_settle <- which(sapply(em_path, function(p) max(abs(p$mu - em_path[[em_rounds + 1]]$mu))) < 1)[1] - 1   # first
# round with both means within 1 mm of where EM ends
loglik_of <- function(par) sum(log(dens_l(flip, 1, par) + dens_l(flip, 2, par)))

fig_em_anim <- function(show = c(0:12, 14, 16, 18, 20, 25, 30, 40)) {
  xs <- seq(165, 240, by = 0.5)
  h <- hist(flip, breaks = seq(170, 236, by = 2), plot = FALSE)
  bars <- tibble(x = h$mids, y = h$density)
  curves <- purrr::map_dfr(show, function(r) { par <- em_path[[r + 1]]
    tibble(x = xs, g1 = dens_l(xs, 1, par), g2 = dens_l(xs, 2, par), p2 = post2(xs, par), s = r) |>
      mutate(tot = g1 + g2) })
  info <- purrr::map_dfr(show, function(r) tibble(x = 167, y = 0.0505, s = r,
    txt = paste0(if (r == 0) "Start: two wide groups in the middle." else paste0("Round ", r, "."),
                 "  ln L = ", mns(f0(loglik_of(em_path[[r + 1]]))))))
  jump <- function(lab, s) list(label = lab, method = "animate",
                                args = list(list(as.character(s)), list(mode = "immediate",
                                            frame = list(duration = 0, redraw = TRUE), transition = list(duration = 0))))
  top <- plot_ly() |>
    add_bars(data = bars, x = ~x, y = ~y, width = 2, marker = list(color = "#e3e3e3"), hoverinfo = "skip",
             showlegend = FALSE) |>
    add_trace(data = curves, x = ~x, y = ~tot, frame = ~s, type = "scatter", mode = "lines",
              line = list(color = "black", width = 2), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = curves, x = ~x, y = ~g1, frame = ~s, type = "scatter", mode = "lines",
              line = list(color = pal_cl[[1]], width = 2.5), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = curves, x = ~x, y = ~g2, frame = ~s, type = "scatter", mode = "lines",
              line = list(color = pal_cl[[3]], width = 2.5), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = info, x = ~x, y = ~y, frame = ~s, text = ~txt, type = "scatter", mode = "text",
              textposition = "middle right", textfont = list(size = 13, color = "#2c3e50"), hoverinfo = "skip",
              showlegend = FALSE)
  bottom <- plot_ly() |>
    add_trace(data = curves, x = ~x, y = ~p2, frame = ~s, type = "scatter", mode = "lines",
              line = list(color = pal_cl[[3]], width = 2.5), hoverinfo = "skip", showlegend = FALSE)
  key <- sprintf(paste0("<span style='color:%s'><b>π<sub>1</sub>f<sub>1</sub>(x)</b></span>   ",
                        "<span style='color:%s'><b>π<sub>2</sub>f<sub>2</sub>(x)</b></span>   <b>f(x)</b>"),
                 pal_cl[[1]], pal_cl[[3]])
  subplot(top, bottom, nrows = 2, heights = c(0.66, 0.34), shareX = TRUE, margin = 0.04, titleY = TRUE) |>
    layout(xaxis = list(title = "flipper length (mm)", range = c(166, 238), zeroline = FALSE),
           yaxis = list(title = "density", range = c(0, 0.053), zeroline = FALSE),
           yaxis2 = list(title = "P(group 2 | x)", range = c(-0.05, 1.05), zeroline = FALSE, tickvals = c(0, 0.5, 1)),
           annotations = list(list(text = key, x = 1, xref = "paper", xanchor = "right", y = 1.02, yref = "paper",
                                   yanchor = "bottom", showarrow = FALSE, font = list(size = 13))),
           showlegend = FALSE, margin = list(t = 60)) |>
    animation_opts(frame = 700, transition = 0, redraw = TRUE) |>
    animation_slider(currentvalue = list(prefix = "round ")) |>
    animation_button(label = "Play") |>
    finalize_animation(list(jump("start", 0), jump("round 40", 40)), start_frame = "0", prefix = "round ")
}

# Two variables: a mixture of three bivariate normals for the bill measurements (mclust chooses the model by BIC)
mix2 <- Mclust(bills, verbose = FALSE)
mix2_cl <- match_labels(mix2$classification, peng$species)
mix2_p <- match(1:3, mix2_cl[match(1:3, mix2$classification)])    # fitted group l is drawn as number mix2_p[l]
mix2_uncertain <- mix2$uncertainty > 0.2                            # the most likely group has probability < 0.8
mix2_bic <- mclustBIC(bills, G = 1:6, modelNames = c("EII", "EEE", "EVE", "VVV"), verbose = FALSE)

# Panels of the bill measurements: the species, k-means clusters, or mixture groups
bill_panel <- function(what = c("species", "kmeans", "mixture"), base = 14, subtitle = NULL) {
  what <- match.arg(what)
  sub_or <- function(default) if (is.null(subtitle)) default else subtitle
  pts <- tibble(x = bills[, 1], y = bills[, 2], species = peng$species)
  p <- ggplot() + bill_lims + labs(x = lab_len, y = lab_dep) + theme_minimal(base_size = base) +
    theme(plot.title = element_text(face = "bold"))
  if (what == "species") {
    ell <- purrr::map_dfr(levels(peng$species), function(s) {
      X <- bills[peng$species == s, ]; ell_df(colMeans(X), cov(X)) |> mutate(species = s) })
    p + geom_point(data = pts, aes(x, y, colour = species), size = 1.4, alpha = 0.75) +
      geom_path(data = ell, aes(x, y, colour = species), linewidth = 1) +
      scale_colour_manual(values = pal_species, name = NULL) +
      theme(legend.position = "top", legend.justification = "right", legend.direction = "horizontal",
            legend.margin = margin(0, 0, 0, 0), legend.box.spacing = unit(2, "pt")) +
      labs(title = "Species", subtitle = sub_or("ellipses: 80% of each species"))
  } else if (what == "kmeans") {
    grid <- expand_grid(x = seq(31.1, 60.9, by = 0.2), y = seq(13.05, 21.95, by = 0.1))
    grid$cl <- factor(nearest(as.matrix(grid), km_bills$centers))
    pts$cl <- factor(km_bills$cluster); off <- km_bills$cluster != as.integer(peng$species)
    # Round groups of one size, with the variance within the clusters: drawn, like the ellipses, to hold 80%
    circ <- purrr::map_dfr(1:3, function(l) ell_df(km_bills$centers[l, ], diag(2) * km_s2) |> mutate(cl = factor(l)))
    p + geom_raster(data = grid, aes(x, y, fill = cl), alpha = 0.16) +
      geom_point(data = pts, aes(x, y, colour = cl), size = 1.4, alpha = 0.8) +
      geom_point(data = pts[off, ], aes(x, y), shape = 21, size = 2.6, stroke = 0.7, colour = "black") +
      geom_path(data = circ, aes(x, y, colour = cl), linewidth = 1, linetype = "22") +
      annotate("point", x = km_bills$centers[, 1], y = km_bills$centers[, 2], shape = 4, size = 5, stroke = 2) +
      scale_fill_manual(values = pal_cl, guide = "none") + scale_colour_manual(values = pal_cl, guide = "none") +
      labs(title = "k-means: three clusters",
           subtitle = sub_or(sprintf("circles: round groups of one size, 80%% of each\nringed: %d penguins not in their species' cluster",
                                     sum(off))))
  } else {
    ell <- purrr::map_dfr(1:3, function(l) ell_df(mix2$parameters$mean[, l], mix2$parameters$variance$sigma[, , l]) |>
                            mutate(cl = factor(mix2_p[l])))
    pts$cl <- factor(mix2_cl)
    p + geom_point(data = pts, aes(x, y, colour = cl), size = 1.4, alpha = 0.8) +
      geom_path(data = ell, aes(x, y, colour = cl), linewidth = 1) +
      geom_point(data = pts[mix2_uncertain, ], aes(x, y), shape = 21, size = 2.6, stroke = 0.7, colour = "black") +
      scale_colour_manual(values = pal_cl, guide = "none") +
      labs(title = "Mixture: three groups",
           subtitle = sub_or(sprintf("ringed: %d penguins with P(most likely group) < 0.8", sum(mix2_uncertain))))
  }
}
fig_bills <- function(panels = c("species", "kmeans"), base = 14, ncol = 1) {   # one column: the panels share the x-axis
  ps <- lapply(panels, bill_panel, base = base)
  if (ncol == 1 && length(ps) > 1) for (i in seq_len(length(ps) - 1)) ps[[i]] <- ps[[i]] + labs(x = NULL)
  wrap_plots(ps, ncol = ncol)
}

# The fourteen shapes of mclust, drawn for three groups. Each group's ellipse has a volume (its size), a shape (how
# stretched it is: the long axis over the short one) and an orientation (its tilt). In the code, the first letter is
# the volume, the second the shape and the third the orientation: E = the same for all groups, V = varies between
# groups, I = round (shape) or along the axes (orientation).
fig_mclust_shapes <- function(base = 12, notes = c(EII = "k-means", EVE = "the bills", VEE = "four measurements")) {
  centres <- rbind(c(-1.75, -1.05), c(1.75, -1.05), c(0, 1.4))
  size  <- list(E = c(1, 1, 1), V = c(0.45, 0.95, 1.6))            # area of each group's ellipse, relative
  shape <- list(I = c(1, 1, 1), E = c(2.6, 2.6, 2.6), V = c(1.3, 2.6, 4.2))   # long axis over short axis
  tilt  <- list(I = c(0, 0, 0), E = c(35, 35, 35), V = c(35, -40, 95))        # degrees
  cols <- c(II = "round", EI = "along the axes,\nsame shape", VI = "along the axes,\nown shape",
            EE = "same shape,\nsame tilt", VE = "own shape,\nsame tilt", EV = "same shape,\nown tilt",
            VV = "own shape,\nown tilt")
  rows <- c(E = "same size", V = "own size")
  codes <- expand_grid(v = names(rows), so = names(cols)) |> mutate(code = paste0(v, so))
  ell <- purrr::map_dfr(seq_len(nrow(codes)), function(r) {
    v <- codes$v[r]; s <- substr(codes$so[r], 1, 1); o <- substr(codes$so[r], 2, 2)
    purrr::map_dfr(1:3, function(l) {
      area <- 0.38 * size[[v]][l]; ratio <- shape[[s]][l]
      M <- rotation_matrix(tilt[[o]][l] * pi / 180) %*% diag(c(sqrt(area * ratio), sqrt(area / ratio)))
      e <- ellipse_xy(M)
      tibble(x = centres[l, 1] + e[1, ], y = centres[l, 2] + e[2, ], g = factor(l))
    }) |> mutate(v = codes$v[r], so = codes$so[r])
  }) |> mutate(row = factor(rows[v], levels = rows), col = factor(cols[so], levels = cols))
  lab <- codes |> mutate(row = factor(rows[v], levels = rows), col = factor(cols[so], levels = cols),
                         note = ifelse(code %in% names(notes), notes[code], NA))
  ggplot() +
    geom_polygon(data = ell, aes(x, y, group = g, fill = g), alpha = 0.25) +
    geom_path(data = ell, aes(x, y, group = g, colour = g), linewidth = 0.8) +
    geom_text(data = lab, aes(x = -3.15, y = 3.05, label = code), hjust = 0, vjust = 1, fontface = "bold",
              size = base / 3) +
    geom_label(data = filter(lab, !is.na(note)), aes(x = 0, y = -3.05, label = note), vjust = 0,
               size = base / 3.7, label.size = 0, fill = "#fff3c4", label.padding = unit(0.12, "lines")) +
    facet_grid(row ~ col, switch = "y") +
    scale_fill_manual(values = pal_cl, guide = "none") + scale_colour_manual(values = pal_cl, guide = "none") +
    coord_equal(xlim = c(-3.2, 3.2), ylim = c(-3.1, 3.1), expand = FALSE) +
    theme_minimal(base_size = base) +
    theme(axis.text = element_blank(), axis.title = element_blank(), panel.grid = element_blank(),
          panel.background = element_rect(fill = "#fafafa", colour = "#dddddd"),
          strip.text.x = element_text(size = base * 0.85, lineheight = 0.95),
          strip.text.y.left = element_text(size = base * 0.95, face = "bold", angle = 90),
          strip.placement = "outside", panel.spacing = unit(4, "pt"))
}

# BIC against the number of groups for four shapes of the groups' ellipses
fig_bic <- function(base = 14) {
  b <- as.data.frame(unclass(mix2_bic)[, , drop = FALSE]); b$k <- seq_len(nrow(b))
  lab <- c(EII = "EII: round, all of the same size", EEE = "EEE: one ellipse for all groups",
           EVE = "EVE: same size and orientation", VVV = "VVV: each group its own ellipse")
  d <- pivot_longer(b, -k, names_to = "model", values_to = "BIC") |> mutate(model = factor(lab[model], levels = lab))
  best <- slice_max(d, BIC)
  ggplot(d, aes(k, BIC, colour = model)) + geom_line(linewidth = 0.9) + geom_point(size = 2.2) +
    geom_point(data = best, size = 5, shape = 21, stroke = 1.3, colour = "black") +
    annotate("text", x = best$k, y = best$BIC + 25, label = "largest BIC", size = 4, vjust = 0) +
    scale_colour_manual(values = c("#999999", col_a1, col_a2, "#2ca02c"), name = NULL) +
    scale_x_continuous(breaks = 1:6) + scale_y_continuous(labels = function(x) mns(f0(x))) +
    labs(x = "number of groups k", y = "BIC") +
    theme_minimal(base_size = base) + theme(legend.position = "right")
}


# =============================================================================
# 5  All four measurements
# =============================================================================

km4_raw <- local({ set.seed(1); kmeans(meas4, centers = 3, nstart = 50)$cluster })     # grams and millimetres
km4_z   <- local({ set.seed(1); kmeans(meas4_z, centers = 3, nstart = 50)$cluster })   # standardised
mix4_3  <- Mclust(meas4, G = 3, verbose = FALSE)                                      # three groups, raw units
mix4_3z <- Mclust(meas4_z, G = 3, verbose = FALSE)                                    # three groups, standardised
mix4    <- Mclust(meas4_z, verbose = FALSE)                                           # BIC chooses the number
peng_group4 <- factor(ifelse(peng$species == "Gentoo", paste("Gentoo", peng$sex), as.character(peng$species)),
                      levels = c("Adelie", "Chinstrap", "Gentoo female", "Gentoo male"))
mix4_cl <- match_labels(mix4$classification, peng_group4)

# The clusters against the species, as tables of counts
fig_units <- function(base = 13) {
  tabs <- setNames(list(km4_raw, km4_z, mix4_3$classification),       # k-means is the mixture EII
                   c("k-means (EII), grams and mm", "k-means (EII), standardised",
                     paste0("mixture ", mix4_3$modelName, ", either units")))
  panels <- lapply(names(tabs), function(nm) {
    tab <- as.data.frame(cross_tab(tabs[[nm]], peng$species))
    ggplot(tab, aes(group, factor(cluster, levels = 3:1), fill = Freq)) + geom_tile(colour = "white") +
      geom_text(aes(label = Freq, colour = Freq > 70), size = 4.6) +
      scale_fill_gradient(low = "#f7fbff", high = "#08306b", limits = c(0, 150), guide = "none") +
      scale_colour_manual(values = c("FALSE" = "black", "TRUE" = "white"), guide = "none") +
      labs(title = nm, subtitle = paste(pct(agree(tabs[[nm]], peng$species)), "in their species' cluster"),
           x = NULL, y = "cluster") +
      theme_minimal(base_size = base) + theme(panel.grid = element_blank(), plot.title = element_text(face = "bold"))
  })
  wrap_plots(panels, nrow = 1)
}

# The four mixture groups on the first two principal components of the four standardised measurements
peng_pc <- local({
  e <- eigen(cov(meas4_z), symmetric = TRUE); G <- e$vectors[, 1:2]
  if (G[3, 1] < 0) G[, 1] <- -G[, 1]                      # PC1: larger flippers to the right
  if (sum(G[1:2, 2]) < 0) G[, 2] <- -G[, 2]               # PC2: larger bills up
  list(G = G, scores = meas4_z %*% G, share = e$values[1:2] / sum(e$values))
})
fig_peng_four <- function(base = 14) {
  pal4 <- setNames(unname(pal_cl), levels(peng_group4))
  pts <- tibble(x = peng_pc$scores[, 1], y = peng_pc$scores[, 2], g = factor(levels(peng_group4)[mix4_cl],
                levels = levels(peng_group4)), sex = peng$sex)
  p_of <- match(1:4, mix4_cl[match(1:4, mix4$classification)])
  ell <- purrr::map_dfr(1:4, function(l) {
    mu <- drop(t(peng_pc$G) %*% mix4$parameters$mean[, l])
    S <- t(peng_pc$G) %*% mix4$parameters$variance$sigma[, , l] %*% peng_pc$G
    ell_df(mu, S) |> mutate(g = factor(levels(peng_group4)[p_of[l]], levels = levels(peng_group4)))
  })
  ggplot() + axes0 +
    geom_point(data = pts, aes(x, y, colour = g, shape = sex), size = 1.9, alpha = 0.8) +
    geom_path(data = ell, aes(x, y, colour = g), linewidth = 1) +
    scale_colour_manual(values = pal4, name = "group, named after\nmost of its penguins") +
    scale_shape_manual(values = c(female = 16, male = 2), name = "sex") +
    coord_equal() +
    labs(x = sprintf("first principal component (%s)", pct(peng_pc$share[1])),
         y = sprintf("second principal component (%s)", pct(peng_pc$share[2]))) +
    theme_minimal(base_size = base)
}


# =============================================================================
# 6  The US states
# =============================================================================

km_states <- local({
  set.seed(1); km <- kmeans(Z6, centers = 3, nstart = 50)
  nm <- character(3)
  nm[which.max(km$centers[, "Illiteracy"])] <- "South"
  nm[which.max(km$centers[, "Frost"])] <- "North"
  nm[nm == ""] <- "Rest"
  list(cluster = factor(nm[km$cluster], levels = c("South", "North", "Rest")),
       centers = km$centers[match(c("South", "North", "Rest"), nm), ], W = km$tot.withinss)
})
rownames(km_states$centers) <- c("South", "North", "Rest")
ward_states <- cutree(hclust(dist(Z6), method = "ward.D2"), k = 3)       # the same three clusters as k-means
avg_states  <- cutree(hclust(dist(Z6), method = "average"), k = 3)
mix_states  <- Mclust(Z6, verbose = FALSE)                               # BIC chooses the shape and the number
pal_states <- c(South = "#d95f02", North = "#1b9e77", Rest = "#7570b3", "Alaska and Nevada" = "#e7298a",
                "West" = "#e6ab02")
# Name the clusters of another method after the k-means clusters they contain
name_states <- function(cl) {
  nm <- sapply(split(rownames(Z6), cl), function(s) {
    if (all(c("AK", "NV") %in% s) && length(s) == 2) return("Alaska and Nevada")
    km <- table(km_states$cluster[rownames(Z6) %in% s])
    if (km["South"] == length(s)) "South" else if (km["North"] == length(s)) "North" else
      if ("CA" %in% s && length(s) < 10) "West" else "Rest"
  })
  factor(nm[as.character(cl)], levels = names(pal_states))
}

# PCA of the correlation matrix, with the signs of Chapter 4 (both components point towards higher income)
st_pc <- local({
  p <- prcomp(X6, scale. = TRUE)
  for (j in 1:2) if (p$rotation["Income", j] < 0) { p$rotation[, j] <- -p$rotation[, j]; p$x[, j] <- -p$x[, j] }
  p
})
fig_states_clusters <- function(base = 13) {
  d <- tibble(x = st_pc$x[, 1], y = st_pc$x[, 2], abb = rownames(X6), cl = km_states$cluster)
  share <- st_pc$sdev^2 / sum(st_pc$sdev^2)
  map <- ggplot(d, aes(x, y, colour = cl)) + axes0 +
    geom_point(size = 2) +
    geom_text_repel(aes(label = abb), size = 3.3, seed = 4, max.overlaps = 40, show.legend = FALSE) +
    scale_colour_manual(values = pal_states, name = NULL) +
    coord_equal() +
    labs(title = "The 50 states, by cluster",
         x = sprintf("PC1 score (%s of the variance)", pct(share[1])), y = sprintf("PC2 score (%s)", pct(share[2]))) +
    theme_minimal(base_size = base) + theme(legend.position = "bottom", plot.title = element_text(face = "bold"))
  ctr <- as_tibble(km_states$centers, rownames = "cl") |> pivot_longer(-cl, names_to = "var", values_to = "z") |>
    mutate(cl = factor(cl, levels = c("South", "North", "Rest")), var = factor(var, levels = rev(vars6)))
  prof <- ggplot(ctr, aes(z, var, colour = cl)) +
    geom_vline(xintercept = 0, colour = "grey70") +
    geom_point(size = 3.2) +
    scale_colour_manual(values = pal_states, guide = "none") +
    labs(title = "Cluster means", x = "standard deviations from the mean of all states", y = NULL) +
    theme_minimal(base_size = base) + theme(plot.title = element_text(face = "bold"))
  map + prof + plot_layout(widths = c(1.45, 1))
}

# Three methods on the same map: k-means, average linkage, and the mixture chosen by BIC
fig_states_methods <- function(base = 13) {
  res <- list("k-means: three clusters" = km_states$cluster,
              "Average linkage: three clusters" = name_states(avg_states),
              "Mixture (BIC): five groups" = name_states(mix_states$classification))
  d <- purrr::map_dfr(names(res), function(m) tibble(x = st_pc$x[, 1], y = st_pc$x[, 2], abb = rownames(X6),
                                                     cl = res[[m]], method = factor(m, levels = names(res))))
  ggplot(d, aes(x, y, colour = cl)) + axes0 +
    geom_point(size = 2) +
    geom_text_repel(data = filter(d, abb %in% c("AK", "NV", "CA", "FL", "HI")), aes(label = abb), size = 3.2,
                    seed = 2, show.legend = FALSE) +
    facet_wrap(~method) +
    scale_colour_manual(values = pal_states, name = NULL, drop = TRUE) +
    coord_equal() +
    labs(x = "PC1 score", y = "PC2 score") +
    theme_minimal(base_size = base) +
    theme(legend.position = "bottom", strip.text = element_text(face = "bold", size = base))
}
