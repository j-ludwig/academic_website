# =============================================================================
# introduction_figures.R  --  data and figures for the chapter "Introduction" of the ISQS 6350 course notes.
# Put this file and _common.R (shared helpers) in one folder, set it as the working
# directory, and run source("introduction_figures.R"). Every figure is a function, e.g. fig_...().
# Packages: tidyverse, plotly, patchwork, mvtnorm, MASS
# =============================================================================

source("_common.R", local = TRUE)

# ---- The airline example ------------------------------------------------------
# For each flight on one route: flight time (minutes) and fuel burned (100 kg),
# as deviations from the typical flight. Two independent shocks: headwind, load.
fuel_per_min <- 6000 / 60 / 100          # cruise fuel flow: 6 t per hour = 1 (100 kg) per minute
sd_wind      <- 1.20                     # headwind: sd of the extra flight time (minutes)
sd_load      <- 2                        # load: sd of passengers and cargo (tonnes)
min_per_t    <- 0.10                     # load: extra minutes per tonne
fuel_per_t   <- 0.65                     # load: extra fuel per tonne (100 kg)
air_A <- matrix(c(sd_wind, fuel_per_min * sd_wind,
                  min_per_t * sd_load, fuel_per_t * sd_load), 2, 2)

app_Sigma <- air_A %*% t(air_A)
app_root  <- sym_root(app_Sigma)
app_X     <- app_W %*% t(air_A)

# Structural restriction: the headwind column points in the known direction (1 min, fuel_per_min)
air_structural <- function(S) {
  u  <- c(1, fuel_per_min) / sqrt(1 + fuel_per_min^2)
  a1 <- u / sqrt(drop(t(u) %*% solve(S) %*% u))
  M  <- S - tcrossprod(a1)
  a2 <- M[, 2] / sqrt(M[2, 2])
  cbind(a1, a2, deparse.level = 0)
}
air_presets <- function(S, root) {
  P <- matrix(c(0, 1, 1, 0), 2)
  eg <- eigen(S, symmetric = TRUE)
  c("Σ½ (symmetric root)"  = phi_from_A(sym_root(S), root),
    "Cholesky"             = phi_from_A(t(chol(S)), root),
    "Cholesky (X₂ first)"  = phi_from_A(P %*% t(chol(P %*% S %*% P)), root),
    "PCA"                  = phi_from_A(eg$vectors %*% diag(sqrt(eg$values)), root),
    "Simultaneous eq. (6 t/h)" = phi_from_A(air_structural(S), root))
}
app_presets <- air_presets(app_Sigma, app_root)
app_restr   <- fa_restrictions(app_Sigma)

f2 <- function(v) sprintf("%.2f", v)
kg <- function(v) sprintf("%.0f", 100 * v)

flights <- tibble(
  flight   = seq_len(nrow(app_X)),
  time     = app_X[, 1],
  fuel     = app_X[, 2],
  headwind = app_W[, 1],
  load     = app_W[, 2]
)
targets <- rbind(c(2, 0.3), c(-0.3, 2), c(-1.6, -1.2), c(1.2, -1.6))
highlight_flights <- apply(targets, 1, function(t)
  which.min((app_W[, 1] - t[1])^2 + (app_W[, 2] - t[2])^2))

# ---- Numbers quoted in the text -------------------------------------------------
num <- list(
  corr      = cor(flights$time, flights$fuel),
  phi_struc = app_presets[["Simultaneous eq. (6 t/h)"]],
  beta      = app_Sigma[1, 2] / app_Sigma[1, 1],
  sig1      = sqrt(app_Sigma[1, 1]),
  sig2_reg  = sqrt(app_Sigma[2, 2] - app_Sigma[1, 2]^2 / app_Sigma[1, 1]),
  lam2_chol = app_Sigma[1, 2] / sqrt(app_Sigma[1, 1])
)

# ---- Population vs sample (Section: three layers) -------------------------------
col_s <- "#e67e22"   # sample
col_e <- "#16a085"   # inferred reduced form
col_h <- "#c0398b"   # estimated structural form
pop_wind <- 1.15; pop_min_per_t <- 0.12; pop_fuel_per_t <- 0.70
A_pop     <- matrix(c(pop_wind, fuel_per_min * pop_wind, pop_min_per_t * sd_load, pop_fuel_per_t * sd_load), 2, 2)
Sigma_pop <- A_pop %*% t(A_pop)
S_hat     <- cov(app_X)
A_hat     <- air_structural(S_hat)

# =============================================================================
# Static figures (ggplot)
# =============================================================================

fig_air_data <- function() {
  ggplot(flights, aes(time, fuel)) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(colour = "#6f7a8c", alpha = 0.55, size = 2) +
    geom_path(data = { e <- ellipse_xy(sym_root(cov(app_X)), 1); tibble(time = e[1, ], fuel = e[2, ]) },
              colour = "black", linewidth = 0.8) +
    coord_equal(xlim = c(-4.2, 4.2), ylim = c(-4.2, 4.2)) +
    labs(title = "200 flights on one route",
         subtitle = sprintf("Correlation %.2f (ellipse: one standard deviation)", num$corr),
         x = "flight time (min)", y = "fuel burned (100 kg)",
         caption = "Both variables are deviations from a typical flight on the route.")
}

fig_iid <- function() {
  ggplot(flights, aes(headwind, load)) +
    geom_hline(yintercept = 0, colour = "grey85") +
    geom_vline(xintercept = 0, colour = "grey85") +
    geom_point(colour = "#6f7a8c", alpha = 0.55, size = 2) +
    geom_path(data = tibble(t = seq(0, 2 * pi, length.out = 200)),
              aes(cos(t), sin(t)), colour = "black", linewidth = 0.8) +
    coord_equal(xlim = c(-3.6, 3.6), ylim = c(-3.6, 3.6)) +
    labs(title = "Two independent shocks",
         subtitle = "The cloud is round (circle: one standard deviation)",
         x = quote(epsilon[1]~~"(headwind)"), y = quote(epsilon[2]~~"(load)"))
}

fig_round_vs_tilted <- function() {
  fig_iid() + labs(title = "Independent shocks", subtitle = "round cloud") +
    fig_air_data() + labs(title = "Observed flights", subtitle = sprintf("tilted cloud, correlation %.2f", num$corr),
                          caption = NULL)
}

fig_shocks_to_data <- function() {
  hl <- flights |> dplyr::filter(flight %in% highlight_flights) |> mutate(lab = paste("flight", flight))
  p_left <- ggplot(flights, aes(headwind, load)) +
    geom_point(colour = "grey75", size = 1.6) +
    geom_point(data = hl, shape = 21, fill = "white", colour = "black", size = 3.6, stroke = 1.6) +
    geom_text(data = hl, aes(label = lab), nudge_y = 0.38, size = 3.6, fontface = "bold") +
    coord_equal(xlim = c(-3.6, 3.6), ylim = c(-3.6, 3.6)) +
    labs(title = "Shocks", subtitle = "independent", x = "headwind shock", y = "load shock")
  p_right <- ggplot(flights, aes(time, fuel)) +
    geom_point(colour = "grey75", size = 1.6) +
    geom_point(data = hl, shape = 21, fill = "white", colour = "black", size = 3.6, stroke = 1.6) +
    geom_text(data = hl, aes(label = lab), nudge_y = 0.45, size = 3.6, fontface = "bold") +
    coord_equal(xlim = c(-4.2, 4.2), ylim = c(-4.2, 4.2)) +
    labs(title = "Flights", subtitle = "X = A ε", x = "flight time (min)", y = "fuel (100 kg)")
  p_left + p_right
}

A_choices <- local({
  eg <- eigen(app_Sigma, symmetric = TRUE)
  V <- eg$vectors
  if (det(V) < 0) V[, 2] <- -V[, 2]
  if (V[1, 1] < 0) V <- -V
  list("Symmetric root" = app_root,
       "Cholesky"       = t(chol(app_Sigma)),
       "PCA"            = V %*% diag(sqrt(eg$values)),
       "Simultaneous equations" = air_A)
})

fig_useful_A <- function(which = names(A_choices), nrow = 2) {
  column_labels <- list(
    "Symmetric root" = c("closest to\ncoordinate 1", "closest to\ncoordinate 2"),
    "Cholesky"       = c("regression\nmovement", "new\ninnovation"),
    "PCA"            = c("major\nmovement", "minor\nmovement"),
    "Simultaneous equations" = c("headwind", "load"))
  choice_levels <- which
  ch <- A_choices[which]
  choice_vectors <- purrr::imap_dfr(ch, function(M, nm) {
    purrr::map_dfr(1:2, function(j) {
      v <- M[, j]; u <- v / sqrt(sum(v^2))
      tibble(method = nm, direction = factor(j), label = column_labels[[nm]][j],
             xend = v[1], yend = v[2],
             label_x = max(sqrt(sum(v^2)) + 1.3, 2.9) * u[1],
             label_y = max(sqrt(sum(v^2)) + 1.3, 2.9) * u[2])
    })
  }) |> mutate(method = factor(method, levels = choice_levels))
  choice_points <- purrr::map_dfr(choice_levels, function(nm)
    tibble(method = nm, x1 = app_X[, 1], x2 = app_X[, 2])) |>
    mutate(method = factor(method, levels = choice_levels))
  choice_ellipse <- purrr::map_dfr(choice_levels, function(nm) {
    e <- ellipse_xy(app_root, 1); tibble(method = nm, x1 = e[1, ], x2 = e[2, ])
  }) |> mutate(method = factor(method, levels = choice_levels))

  ggplot(choice_points, aes(x1, x2)) +
    geom_point(alpha = 0.18, colour = "#6f7a8c", size = 1) +
    geom_path(data = choice_ellipse, linewidth = 0.8) +
    geom_segment(data = choice_vectors,
                 aes(x = 0, y = 0, xend = xend, yend = yend, colour = direction),
                 arrow = grid::arrow(length = grid::unit(0.16, "cm")), linewidth = 1.3) +
    geom_segment(data = choice_vectors, aes(x = xend, y = yend, xend = 0.8 * label_x, yend = 0.8 * label_y,
                                            colour = direction), linewidth = 0.3, linetype = "22") +
    geom_text(data = choice_vectors, aes(label_x, label_y, label = label, colour = direction),
              size = 3.2, lineheight = 0.9, fontface = "bold") +
    scale_colour_manual(values = c("1" = col_a1, "2" = col_a2), guide = "none") +
    facet_wrap(~ method, nrow = nrow) +
    coord_equal(xlim = c(-3.6, 3.6), ylim = c(-3.6, 3.6)) +
    labs(x = "flight time (min)", y = "fuel (100 kg)")
}

# ---- Three layers --------------------------------------------------------------
lim  <- c(-4.2, 4.2)
ell  <- function(M, r = 1) { e <- ellipse_xy(M, r); tibble(x = e[1, ], y = e[2, ]) }
pts  <- tibble(x = app_X[, 1], y = app_X[, 2])
pt_alpha <- 0.6; pt_size <- 1.4
axes0 <- list(geom_hline(yintercept = 0, colour = "grey88"), geom_vline(xintercept = 0, colour = "grey88"))
arrows_of <- function(M, cols = c(col_a1, col_a2), lw = 1.3)
  geom_segment(data = tibble(xend = M[1, ], yend = M[2, ], col = cols),
               aes(x = 0, y = 0, xend = xend, yend = yend, colour = I(col)),
               linewidth = lw, arrow = arrow(length = unit(0.22, "cm"), type = "closed"))
labels_of <- function(M, labs, cols = c(col_a1, col_a2), k = 0.42) {
  P <- t(M) * (1 + k / sqrt(colSums(M^2)))
  annotate("text", x = P[, 1], y = P[, 2], label = labs, parse = TRUE, colour = cols, size = 4.6)
}
connector <- function(verb, note, col, parse_note = FALSE) {
  ggplot() +
    annotate("segment", x = 0.05, xend = 0.95, y = 0.5, yend = 0.5, colour = col, linewidth = 1.2,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
    annotate("text", x = 0.5, y = 0.62, label = verb, colour = col, size = 5, fontface = "bold") +
    annotate("text", x = 0.5, y = 0.36, label = note, colour = col, size = 3.9, parse = parse_note,
             lineheight = 0.9) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), expand = FALSE) +
    theme_void()
}
title_sf <- quote(bold("Structural form")~~epsilon)
title_rf <- quote(bold("Reduced form")~~italic(X))
title_sa <- quote(bold("Sample")~~bold(X))

fig_layers_cycle <- function() {
  bx <- c(1, 3.5, 6)
  box_txt <- tibble(
    x = rep(bx, each = 4), y = rep(c(0.33, 0.11, -0.08, -0.31), 3),
    size = rep(c(3.9, 3.4, 3.4, 4.1), 3),
    lab = c('bold("STRUCTURAL FORM")', '"model parameters and"', '"distribution of shocks"',
            'italic(A)*",  "*epsilon %~% N(0, I)',
            'bold("REDUCED FORM")', '"moments and distribution of"', '"the observed random vector"',
            'Sigma*",  "*italic(X) %~% N(0, Sigma)',
            'bold("SAMPLE")', '"200 flights (200 rows) with"', '"time and fuel (2 columns)"',
            'bold(X)'))
  ggplot() +
    annotate("rect", xmin = bx - 0.88, xmax = bx + 0.88, ymin = -0.52, ymax = 0.52,
             fill = "white", colour = "grey40", linewidth = 0.5) +
    geom_text(data = box_txt, aes(x, y, label = lab, size = I(size)), parse = TRUE) +
    annotate("segment", x = c(1.95, 4.45), xend = c(2.55, 5.05), y = 0.25, yend = 0.25,
             colour = c("black", col_s), linewidth = 1.1, arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
    annotate("text", x = c(2.25, 4.75), y = 0.55, label = c("generate", "sample"),
             colour = c("black", col_s), size = 4.6, fontface = "bold") +
    annotate("segment", x = c(2.55, 5.05), xend = c(1.95, 4.45), y = -0.25, yend = -0.25,
             colour = c("black", col_e), linewidth = 1.1, arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
    annotate("text", x = c(2.25, 4.75), y = -0.55, label = c("identify", "infer"),
             colour = c("black", col_e), size = 4.6, fontface = "bold") +
    annotate("path", x = c(6, 6, 1), y = c(-0.54, -1.05, -1.05), colour = col_h, linewidth = 1.1) +
    annotate("segment", x = 1, xend = 1, y = -1.05, yend = -0.56, colour = col_h, linewidth = 1.1,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
    annotate("label", x = 3.5, y = -1.05, label = "estimate = infer + identify", colour = col_h,
             size = 4.4, fontface = "bold", label.size = 0) +
    annotate("segment", x = 0.12, xend = 4.38, y = 0.95, yend = 0.95, colour = "grey45") +
    annotate("segment", x = c(0.12, 4.38), xend = c(0.12, 4.38), y = 0.95, yend = 0.86, colour = "grey45") +
    annotate("text", x = 2.25, y = 1.08, label = "population (unknown)", colour = "grey30", size = 4) +
    annotate("segment", x = 5.12, xend = 6.88, y = 0.95, yend = 0.95, colour = "grey45") +
    annotate("segment", x = c(5.12, 6.88), xend = c(5.12, 6.88), y = 0.95, yend = 0.86, colour = "grey45") +
    annotate("text", x = 6, y = 1.08, label = "what we observe", colour = "grey30", size = 4) +
    coord_cartesian(xlim = c(0, 7), ylim = c(-1.2, 1.2), expand = FALSE) +
    theme_void()
}

fig_layers_nature <- function() {
  s_struct <- ggplot() + axes0 +
    geom_path(data = ell(diag(2)), aes(x, y), colour = "black", linewidth = 1) +
    arrows_of(diag(2)) + labels_of(diag(2), c("e[1]", "e[2]")) +
    coord_equal(xlim = lim, ylim = lim) +
    labs(title = title_sf, subtitle = "shocks: uncorrelated,\nunit variance (a circle)",
         x = expression(epsilon[1]~"(headwind)"), y = expression(epsilon[2]~"(load)"))
  s_pop <- ggplot() + axes0 +
    geom_path(data = ell(A_pop), aes(x, y), colour = "black", linewidth = 1) +
    arrows_of(A_pop) + labels_of(A_pop, c("A*e[1]", "A*e[2]")) +
    coord_equal(xlim = lim, ylim = lim) +
    labs(title = title_rf, subtitle = "the circle moved by A:\nan ellipse, Σ = AA'",
         x = "flight time (min)", y = "fuel (100 kg)")
  s_samp <- ggplot() + axes0 +
    geom_path(data = ell(A_pop), aes(x, y), colour = "black", linewidth = 0.5) +
    geom_point(data = pts, aes(x, y), colour = col_s, alpha = pt_alpha, size = pt_size) +
    coord_equal(xlim = lim, ylim = lim) +
    labs(title = title_sa, subtitle = "data: 200 draws\nfrom the random vector",
         x = "flight time (min)", y = "fuel (100 kg)")
  s_struct + connector("generate", "italic(X)==italic(A)*epsilon", "black", parse_note = TRUE) +
    s_pop + connector("sample", "draw 200\nflights", col_s) + s_samp +
    plot_layout(widths = c(1, 0.36, 1, 0.36, 1))
}

fig_layers_analysis <- function() {
  E <- app_X %*% t(solve(A_hat))
  eps_hat <- tibble(x = E[, 1], y = E[, 2])
  a_samp <- ggplot() + axes0 +
    geom_point(data = pts, aes(x, y), colour = col_s, alpha = pt_alpha, size = pt_size) +
    coord_equal(xlim = lim, ylim = lim) +
    labs(title = title_sa, subtitle = "data: all we have",
         x = "flight time (min)", y = "fuel (100 kg)")
  a_pop <- ggplot() + axes0 +
    geom_point(data = pts, aes(x, y), colour = col_s, alpha = pt_alpha, size = pt_size) +
    geom_path(data = ell(A_pop), aes(x, y), colour = "black", linewidth = 0.5) +
    geom_path(data = ell(A_hat), aes(x, y), colour = col_e, linewidth = 1.2, linetype = "22") +
    coord_equal(xlim = lim, ylim = lim) +
    labs(title = title_rf, subtitle = "inferred: S (teal)\nthin black: true Σ (unknown)",
         x = "flight time (min)", y = "fuel (100 kg)")
  a_struct <- ggplot() + axes0 +
    geom_point(data = eps_hat, aes(x, y), colour = col_h, alpha = pt_alpha, size = pt_size) +
    geom_path(data = ell(diag(2)), aes(x, y), colour = "black", linewidth = 1) +
    coord_equal(xlim = lim, ylim = lim) +
    labs(title = title_sf, subtitle = "estimated shocks,\nwhitened via the manual",
         x = expression(hat(epsilon)[1]), y = expression(hat(epsilon)[2]))
  a_samp + connector("infer", "uses\nthe data", col_e) + a_pop +
    connector("identify", "restriction,\nno data", "black") + a_struct +
    plot_layout(widths = c(1, 0.36, 1, 0.36, 1))
}

# ---- Dimension reduction ---------------------------------------------------------
dimred <- local({
  S_tight <- matrix(c(1, 0.95, 0.95, 1.1), 2)
  X_tight <- app_W %*% t(sym_root(S_tight))
  eg_t <- eigen(cov(X_tight), symmetric = TRUE)
  v1 <- eg_t$vectors[, 1]; if (v1[1] < 0) v1 <- -v1
  s1 <- as.vector(X_tight %*% v1)
  proj <- outer(s1, v1)
  list(v1 = v1, share1 = eg_t$values[1] / sum(eg_t$values),
       df = tibble(x1 = X_tight[, 1], x2 = X_tight[, 2], p1 = proj[, 1], p2 = proj[, 2], s1 = s1))
})

fig_dimred <- function() {
  df_t <- dimred$df; v1 <- dimred$v1
  p_proj <- ggplot(df_t) +
    geom_abline(slope = v1[2] / v1[1], intercept = 0, colour = col_a1, linewidth = 0.8) +
    geom_segment(aes(x1, x2, xend = p1, yend = p2), colour = "grey70", linewidth = 0.3) +
    geom_point(aes(x1, x2), colour = "#6f7a8c", alpha = 0.6, size = 1.6) +
    coord_equal(xlim = c(-3.5, 3.5), ylim = c(-3.5, 3.5)) +
    labs(title = "Coke and Pepsi: almost one direction",
         subtitle = sprintf("Market direction: %.0f%% of the variance", 100 * dimred$share1),
         x = "Coca-Cola sales (deviation)", y = "Pepsi sales (deviation)")
  p_line <- ggplot(df_t, aes(s1, 0)) +
    geom_hline(yintercept = 0, colour = col_a1, linewidth = 0.8) +
    geom_point(colour = "#6f7a8c", alpha = 0.6, size = 1.6,
               position = position_jitter(height = 0.04, width = 0, seed = 1)) +
    scale_y_continuous(limits = c(-0.5, 0.5), breaks = NULL) +
    labs(title = "One number per week",
         subtitle = "Position along the line: the soft-drink market\n(first principal component)",
         x = "size of the soft-drink market that week", y = NULL)
  p_proj + p_line + plot_layout(widths = c(1, 1.1))
}

# ---- Factor model: four restrictions ------------------------------------------------
fig_factor_choices <- function(nrow = 2, pick = 1:4, critics = FALSE,
                               panels = c("No noise in X₁", "No noise in X₂", "Equal communalities", "Equal noise")) {
  restr_pick <- c("ψ₁ = 0 (X₁ is the factor)", "ψ₂ = 0 (X₂ is the factor)",
                  "equal communalities h₁ = h₂", "equal noise ψ₁ = ψ₂")[pick]
  panel_names <- panels[pick]
  fa_models <- lapply(restr_pick, function(nm) fa_from_split(app_Sigma, app_restr[[nm]]))
  fa_vectors <- purrr::map2_dfr(fa_models, panel_names, function(m, nm) {
    tibble(method = nm, which = c("λ", "σ1", "σ2"),
           xend = c(m$lambda[1], sqrt(m$psi[1]), 0),
           yend = c(m$lambda[2], 0, sqrt(m$psi[2])))
  }) |>
    mutate(method = factor(method, levels = panel_names),
           label_x = xend + ifelse(which == "λ", 0.3, ifelse(which == "σ1", 0.15, -0.35)),
           label_y = yend + ifelse(which == "λ", 0.25, ifelse(which == "σ1", -0.3, 0.1)))
  fa_lines <- purrr::map2_dfr(fa_models, panel_names, function(m, nm)
    tibble(method = nm, x = c(-4, 4) * m$lambda[1], y = c(-4, 4) * m$lambda[2])) |>
    mutate(method = factor(method, levels = panel_names))
  fa_pts <- purrr::map_dfr(panel_names, function(nm)
    tibble(method = nm, x1 = app_X[, 1], x2 = app_X[, 2])) |>
    mutate(method = factor(method, levels = panel_names))
  fa_ell <- purrr::map_dfr(panel_names, function(nm) {
    e <- ellipse_xy(app_root, 1); tibble(method = nm, x1 = e[1, ], x2 = e[2, ])
  }) |> mutate(method = factor(method, levels = panel_names))

  ggplot(fa_pts, aes(x1, x2)) +
    geom_line(data = fa_lines, aes(x, y), colour = col_f, linetype = 2, alpha = 0.5) +
    geom_point(alpha = 0.18, colour = "#6f7a8c", size = 1) +
    geom_path(data = fa_ell, linewidth = 0.8) +
    geom_segment(data = fa_vectors |> dplyr::filter(xend^2 + yend^2 > 1e-6),
                 aes(x = 0, y = 0, xend = xend, yend = yend, colour = which),
                 arrow = grid::arrow(length = grid::unit(0.16, "cm")), linewidth = 1.3) +
    geom_text(data = fa_vectors |> dplyr::filter(xend^2 + yend^2 > 0.01),
              aes(label_x, label_y, label = which, colour = which), size = 4.2, fontface = "bold") +
    scale_colour_manual(values = c("λ" = col_f, "σ1" = col_a1, "σ2" = col_a2), guide = "none") +
    facet_wrap(~ method, nrow = nrow) +
    coord_equal(xlim = c(-3, 3), ylim = c(-3, 3)) +
    (if (critics) labs(x = "critic 1's rating", y = "critic 2's rating")
     else labs(x = expression(x[1]~"(flight time)"), y = expression(x[2]~"(fuel)"))) +
    (if (critics) theme(axis.text = element_blank()) else NULL)
}

# ---- From factors to clusters ---------------------------------------------------------
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

# =============================================================================
# Interactive figures (plotly)
# =============================================================================

make_grid_animation <- function(A, Z, k_grid = seq(0, 1, by = 0.04), lim = 3.6, grid_max = 6,
                                highlight = integer(0), arrow_vecs = diag(2),
                                arrow_text = c("Ae<sub>1</sub>", "Ae<sub>2</sub>"),
                                xlab = "first coordinate", ylab = "second coordinate") {
  lab <- function(k) sprintf("%.2f", k)
  ticks <- seq(-grid_max, grid_max, by = 1)
  grid_paths <- function(M) {
    pieces <- c(
      lapply(ticks, function(c0) M %*% rbind(c(c0, c0), c(-grid_max, grid_max))),
      lapply(ticks, function(c0) M %*% rbind(c(-grid_max, grid_max), c(c0, c0))))
    tibble(x = unlist(lapply(pieces, function(P) c(P[1, ], NA))),
           y = unlist(lapply(pieces, function(P) c(P[2, ], NA))))
  }
  frames_grid <- purrr::map_dfr(k_grid, function(k) grid_paths(matrix_power_real(A, k)) |> mutate(k = lab(k)))
  frames_pts <- purrr::map_dfr(k_grid, function(k) {
    X <- Z %*% t(matrix_power_real(A, k))
    tibble(id = seq_len(nrow(Z)), x = X[, 1], y = X[, 2], k = lab(k))
  })
  arrows <- lapply(1:2, function(j) purrr::map_dfr(k_grid, function(k)
    arrow_geometry((matrix_power_real(A, k) %*% arrow_vecs)[, j], head_length = 0.22, head_width = 0.11,
                   shaft_width = 0.028)$arrow |> mutate(k = lab(k))))
  labels <- purrr::map_dfr(k_grid, function(k) {
    M <- matrix_power_real(A, k) %*% arrow_vecs
    tibble(x = M[1, ] + c(0.3, -0.1), y = M[2, ] + c(-0.22, 0.3),
           text = arrow_text, col = c(col_a1, col_a2), k = lab(k))
  })
  frames_hl <- purrr::map_dfr(k_grid, function(k) {
    X <- Z[highlight, , drop = FALSE] %*% t(matrix_power_real(A, k))
    tibble(id = highlight, x = X[, 1], y = X[, 2], text = paste("flight", highlight), k = lab(k))
  })
  static <- grid_paths(diag(2))

  p <- plot_ly() |>
    add_trace(data = static, x = ~x, y = ~y, type = "scatter", mode = "lines",
              line = list(color = "rgba(150,150,150,0.25)", width = 1),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(data = frames_grid, x = ~x, y = ~y, frame = ~k, type = "scatter", mode = "lines",
              line = list(color = "rgba(31,119,180,0.35)", width = 1.3),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(data = frames_pts, x = ~x, y = ~y, frame = ~k, ids = ~id,
              type = "scatter", mode = "markers",
              marker = list(size = 5, color = "#6f7a8c", opacity = 0.55),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)
  if (length(highlight)) {
    p <- p |>
      add_trace(data = frames_hl, x = ~x, y = ~y, frame = ~k, ids = ~id, text = ~text,
                type = "scatter", mode = "markers+text", textposition = "top center",
                marker = list(size = 11, color = "white", line = list(color = "black", width = 2.5)),
                textfont = list(size = 13, color = "black"),
                hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)
  }
  for (j in 1:2) p <- add_arrow_trace(p, arrows[[j]], color = c(col_a1, col_a2)[j])
  p |>
    add_trace(data = labels, x = ~x, y = ~y, frame = ~k, text = ~text,
              type = "scatter", mode = "text", textfont = list(size = 16, color = ~col),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    layout(
      xaxis = list(range = c(-lim, lim), zeroline = TRUE, zerolinecolor = "#999",
                   showgrid = FALSE, title = xlab),
      yaxis = list(range = c(-lim, lim), zeroline = TRUE, zerolinecolor = "#999",
                   showgrid = FALSE, scaleanchor = "x", title = ylab),
      margin = list(l = 60, r = 20, t = 30, b = 60),
      showlegend = FALSE) |>
    animation_opts(frame = 90, transition = 0, redraw = TRUE) |>
    animation_slider(currentvalue = list(prefix = "progress k = ")) |>
    animation_button(label = "Play")
}

fig_grid_world <- function()
  make_grid_animation(air_A, app_W, lim = 4.2, highlight = highlight_flights,
                      arrow_text = c("headwind", "load"),
                      xlab = "headwind shock → flight time (min)", ylab = "load shock → fuel (100 kg)")

fig_grid_forward <- function()
  make_grid_animation(air_A, app_W, lim = 4.2, highlight = highlight_flights,
                      xlab = "headwind → flight time", ylab = "load → fuel")

fig_grid_inverse <- function()
  make_grid_animation(solve(air_A), app_X, lim = 4.2, highlight = highlight_flights,
                      arrow_vecs = air_A, arrow_text = c("headwind", "load"),
                      xlab = "flight time → headwind", ylab = "fuel → load")

fig_rotation <- function() make_rotation_figure()
fig_factor   <- function() make_factor_figure(notation = "sigma")

# ---- MDS: the loop road ------------------------------------------------------------
t_road <- seq(0, 2 * pi, length.out = 2001)
road <- tibble(east  = 4.2 * sign(cos(t_road)) * abs(cos(t_road))^(1 / 2),
               north = 1.5 * sign(sin(t_road)) * abs(sin(t_road))^(1 / 2))
road$s <- c(0, cumsum(sqrt(diff(road$east)^2 + diff(road$north)^2)))
L_road <- max(road$s)
places <- tibble(
  place = c("Church", "Grocery", "School", "City hall", "Library", "Park", "Hospital", "Stadium"),
  frac  = c(0.03, 0.15, 0.24, 0.34, 0.50, 0.62, 0.75, 0.87)) |>
  mutate(s = frac * L_road,
         east  = approx(road$s, road$east,  s)$y,
         north = approx(road$s, road$north, s)$y)
ds <- abs(outer(places$s, places$s, "-"))
D_drive <- pmin(ds, L_road - ds); dimnames(D_drive) <- list(places$place, places$place)
D_straight <- as.matrix(dist(places[, c("east", "north")])); dimnames(D_straight) <- list(places$place, places$place)
speed   <- 30
T_loop  <- 60 * L_road / speed
s_min   <- 60 * places$s / speed
D_wait <- function(w) {
  a <- abs(outer(s_min, s_min, "-"))
  D <- pmin(a, T_loop - a + w)
  dimnames(D) <- list(places$place, places$place); D
}
w_star <- 2 * max(abs(outer(s_min, s_min, "-"))) - T_loop
w_switch <- T_loop - 2 * (T_loop - (s_min[places$place == "Stadium"] - s_min[places$place == "Grocery"]))

fig_mds <- function(height = 660) {
  pair  <- c("Stadium", "Grocery")
  c_frac <- 0.95
  site <- tibble(east  = approx(road$s, road$east,  c_frac * L_road)$y,
                 north = approx(road$s, road$north, c_frac * L_road)$y)
  w_grid <- sort(unique(c(0:36, round(w_star, 1))))
  wlab   <- function(w) sprintf("%.1f", w)
  theta <- 2 * pi * places$frac
  prev  <- cbind(cos(theta), sin(theta)) * T_loop / (2 * pi)
  mds_frames <- purrr::map_dfr(w_grid, function(w) {
    Y <- cmdscale(as.dist(D_wait(w)), k = 2)
    Y <- procrustes_align(Y, prev); prev <<- Y
    places |> dplyr::select(place) |> mutate(d1 = Y[, 1], d2 = Y[, 2], k = wlab(w))
  })
  road_mds <- mds_frames |> group_by(k) |>
    dplyr::group_modify(~ dplyr::bind_rows(.x, .x[1, ])) |> ungroup()
  pair_mds <- mds_frames |> dplyr::filter(place %in% pair)
  pr <- places |> dplyr::filter(place %in% pair)

  p_map <- plot_ly() |>
    add_trace(data = road, x = ~east, y = ~north, type = "scatter", mode = "lines", color = I("grey80"),
              line = list(color = "grey80", width = 10), hoverinfo = "skip") |>
    add_trace(x = pr$east, y = pr$north, type = "scatter", mode = "lines", color = I(col_a2),
              line = list(color = col_a2, dash = "dash", width = 1.5), hoverinfo = "skip") |>
    add_trace(data = places, x = ~east, y = ~north, type = "scatter", mode = "markers+text",
              text = ~place, color = I("#2c3e50"),
              textposition = dplyr::case_when(places$east > 3.5 ~ "middle right", places$east < -3.5 ~ "middle left",
                                              places$north > 0 ~ "top center", TRUE ~ "bottom center"),
              marker = list(color = "#2c3e50", size = 9), textfont = list(size = 12), hoverinfo = "skip") |>
    add_trace(data = site, x = ~east, y = ~north, type = "scatter", mode = "markers+text",
              text = "construction", textposition = "bottom right", color = I("#e67e22"),
              marker = list(symbol = "x", color = "#e67e22", size = 16, line = list(width = 2)),
              textfont = list(color = "#e67e22", size = 12), hoverinfo = "skip")
  p_mds <- plot_ly() |>
    add_trace(data = road_mds, x = ~d1, y = ~d2, frame = ~k, type = "scatter", mode = "lines",
              color = I("grey80"), line = list(color = "grey80", width = 6), hoverinfo = "skip") |>
    add_trace(data = mds_frames, x = ~d1, y = ~d2, frame = ~k, text = ~place, type = "scatter",
              mode = "markers+text", textposition = "middle right", color = I("#2c3e50"),
              marker = list(color = "#2c3e50", size = 9), textfont = list(size = 12),
              hovertemplate = "%{text}<extra></extra>") |>
    add_trace(data = pair_mds, x = ~d1, y = ~d2, frame = ~k, type = "scatter", mode = "lines",
              color = I(col_a2), line = list(color = col_a2, dash = "dash", width = 1.5), hoverinfo = "skip")
  p <- subplot(p_map, p_mds, widths = c(0.42, 0.58), margin = 0.05) |>
    layout(showlegend = FALSE, margin = list(t = 90),
           xaxis  = list(title = "west - east (km)", range = c(-7, 8.5), zeroline = FALSE),
           yaxis  = list(title = "south - north (km)", range = c(-5, 5), scaleanchor = "x", zeroline = FALSE),
           xaxis2 = list(title = "MDS dimension 1 (minutes)", range = c(-22, 22), zeroline = FALSE),
           yaxis2 = list(title = "MDS dimension 2", range = c(-20, 20), scaleanchor = "x2", zeroline = FALSE),
           annotations = list(
             list(x = 0, xanchor = "left", y = 1.09, xref = "paper", yref = "paper", showarrow = FALSE,
                  align = "left", font = list(size = 15),
                  text = sprintf("On the map<br><span style='font-size:12px'>Stadium to Grocery: %.1f km apart, %.1f km by road</span>",
                                 D_straight["Stadium", "Grocery"], D_drive["Stadium", "Grocery"])),
             list(x = 1, xanchor = "right", y = 1.09, xref = "paper", yref = "paper", showarrow = FALSE,
                  align = "right", font = list(size = 15),
                  text = "Rebuilt from travel times (MDS)<br><span style='font-size:12px'>grey line: the road, in driving order</span>"))) |>
    animation_opts(frame = 150, transition = 0, redraw = TRUE) |>
    animation_slider(currentvalue = list(prefix = "waiting time at the construction site: ", suffix = " minutes"))
  fig <- finalize_animation(p, list(), start_frame = wlab(0), prefix = "waiting time at the construction site: ",
                            suffix = " minutes")
  fig$x$layout$updatemenus <- head(fig$x$layout$updatemenus, -1)
  fig$height <- height
  fig
}


# =============================================================================
# Extra figures for the slides (one picture per slide)
# =============================================================================

cloud_base <- function(alpha = 0.25, lims = c(-4.2, 4.2)) {
  e <- ellipse_xy(app_root, 1)
  list(axes0,
       geom_point(data = pts, aes(x, y), colour = "#6f7a8c", alpha = alpha, size = 1.4),
       geom_path(data = tibble(x = e[1, ], y = e[2, ]), aes(x, y), colour = "black", linewidth = 0.9),
       coord_equal(xlim = lims, ylim = lims),
       labs(x = "flight time (min)", y = "fuel (100 kg)"))
}

# The two columns of A drawn as arrows on the cloud
fig_A_columns <- function() {
  lab <- tibble(x = air_A[1, ] + c(0.25, -0.15), y = air_A[2, ] + c(-0.55, 0.45),
                txt = c(sprintf("headwind\n(%s, %s)", f2(air_A[1, 1]), f2(air_A[2, 1])),
                        sprintf("load\n(%s, %s)", f2(air_A[1, 2]), f2(air_A[2, 2]))),
                col = c(col_a1, col_a2), hj = c(0, 1))
  ggplot() + cloud_base() + arrows_of(air_A, lw = 1.6) +
    geom_text(data = lab, aes(x, y, label = txt, colour = I(col), hjust = hj),
              size = 4.6, fontface = "bold", lineheight = 0.9) +
    labs(title = "Columns of A", subtitle = "what one unit of each shock does")
}

# Sigma as a picture: the ellipse's width and height are the standard deviations
fig_sigma <- function() {
  s1 <- sqrt(app_Sigma[1, 1]); s2 <- sqrt(app_Sigma[2, 2])
  ggplot() + cloud_base(alpha = 0.2) +
    annotate("segment", x = c(-s1, s1), xend = c(-s1, s1), y = -4.2, yend = 4.2,
             colour = col_a1, linetype = 2) +
    annotate("segment", y = c(-s2, s2), yend = c(-s2, s2), x = -4.2, xend = 4.2,
             colour = col_a2, linetype = 2) +
    annotate("segment", x = -s1, xend = s1, y = -3.7, yend = -3.7, colour = col_a1, linewidth = 1,
             arrow = arrow(ends = "both", length = unit(0.18, "cm"))) +
    annotate("text", x = 0, y = -3.35, label = sprintf("2σ₁  (σ₁ = %s min)", f2(s1)),
             colour = col_a1, size = 4.6, fontface = "bold") +
    annotate("segment", y = -s2, yend = s2, x = -3.7, xend = -3.7, colour = col_a2, linewidth = 1,
             arrow = arrow(ends = "both", length = unit(0.18, "cm"))) +
    annotate("text", x = -3.35, y = 0, label = sprintf("2σ₂  (σ₂ = %s)", f2(s2)), angle = 90,
             colour = col_a2, size = 4.6, fontface = "bold") +
    labs(title = "Σ as an ellipse", subtitle = sprintf("width: σ₁, height: σ₂, tilt: ρ = %s", f2(num$corr)))
}

# Several candidate A's: all arrow tips lie on the same ellipse
fig_many_A <- function(phis = c(-60, -20, 25, 70, 115)) {
  arr <- purrr::map_dfr(seq_along(phis), function(i) {
    M <- app_root %*% rotation_matrix(phis[i] * pi / 180)
    tibble(g = i, j = factor(1:2), xend = M[1, ], yend = M[2, ])
  })
  ggplot() + cloud_base(alpha = 0.15) +
    geom_segment(data = arr, aes(x = 0, y = 0, xend = xend, yend = yend, colour = j, group = interaction(g, j)),
                 linewidth = 1, alpha = 0.75, arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
    geom_point(data = arr, aes(xend, yend, colour = j), size = 2.2) +
    scale_colour_manual(values = c("1" = col_a1, "2" = col_a2), guide = "none") +
    labs(title = "Five candidates for A", subtitle = "every pair of arrows fits the same Σ")
}

# One of the four useful choices, as a single panel
fig_choice <- function(name) {
  fig_useful_A(which = name, nrow = 1) + labs(title = name) +
    theme(strip.text = element_blank(), plot.title = element_text(face = "bold"))
}

# Classical MDS of the loop road: ring (no waiting) and line (long wait), static
fig_mds_static <- function() {
  L <- 1.25 * max(abs(cmdscale(as.dist(D_wait(0)), k = 2)), abs(cmdscale(as.dist(D_wait(ceiling(w_star) + 2)), k = 2)))
  one <- function(w, ttl) {
    Y <- cmdscale(as.dist(D_wait(w)), k = 2)
    d <- places |> dplyr::select(place) |> mutate(d1 = Y[, 1], d2 = Y[, 2])
    ggplot(d, aes(d1, d2)) +
      geom_path(data = dplyr::bind_rows(d, d[1, ]), colour = "grey80", linewidth = 2.5) +
      geom_point(colour = "#2c3e50", size = 2.6) +
      (if (w == 0) geom_text(aes(label = place), size = 3.6, vjust = -0.9)
       else geom_text(aes(label = place), size = 3.6, angle = 60, hjust = 0, nudge_y = 1.2)) +
      coord_equal(xlim = c(-L, L), ylim = c(-L, L)) + labs(title = ttl, x = "MDS 1 (minutes)", y = "MDS 2")
  }
  one(0, "No waiting: a ring") + one(ceiling(w_star) + 2, sprintf("Wait > %.0f min: a line", w_star))
}
ggrepel_or_text <- function(d) geom_text(data = d, aes(label = place), size = 3.6, vjust = -0.9)

# "A world without relationships": start from the observed flights (styled like fig_air_data),
# press Play to move them back to the independent shocks. Only points and the 1-sd ellipse.
fig_world_reverse <- function(k_grid = seq(0, 1, by = 0.04), lim = 4.2) {
  lab <- function(k) sprintf("%.2f", k)
  Ainv <- solve(air_A)
  e0 <- ellipse_xy(sym_root(cov(app_X)), 1, n = 200)
  pts_f <- purrr::map_dfr(k_grid, function(k) {
    M <- matrix_power_real(Ainv, k); X <- app_X %*% t(M)
    tibble(id = seq_len(nrow(X)), x = X[, 1], y = X[, 2], k = lab(k))
  })
  ell_f <- purrr::map_dfr(k_grid, function(k) {
    E <- matrix_power_real(Ainv, k) %*% e0
    tibble(x = E[1, ], y = E[2, ], k = lab(k))
  })
  ax <- function(title) list(range = c(-lim, lim), title = title, showgrid = TRUE, gridcolor = "#ebebeb",
                             zeroline = TRUE, zerolinecolor = "#d9d9d9", tickfont = list(size = 13),
                             titlefont = list(size = 15))
  plot_ly() |>
    add_trace(data = pts_f, x = ~x, y = ~y, frame = ~k, ids = ~id, type = "scatter", mode = "markers",
              marker = list(size = 7, color = "#6f7a8c", opacity = 0.55),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(data = ell_f, x = ~x, y = ~y, frame = ~k, type = "scatter", mode = "lines",
              line = list(color = "black", width = 2.2),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    layout(xaxis = ax("flight time (min)  →  headwind shock ε<sub>1</sub>"),
           yaxis = c(ax("fuel (100 kg)  →  load shock ε<sub>2</sub>"), list(scaleanchor = "x")),
           plot_bgcolor = "white", paper_bgcolor = "white",
           margin = list(l = 60, r = 20, t = 20, b = 60), showlegend = FALSE) |>
    animation_opts(frame = 90, transition = 0, redraw = TRUE) |>
    animation_slider(currentvalue = list(prefix = "k = ")) |>
    animation_button(label = "Play")
}

# Forward version: start from the independent shocks (round), press Play to produce the observed
# flights, ending exactly at the picture of fig_air_data (same points, same 1-sd ellipse).
fig_world_forward <- function(k_grid = seq(0, 1, by = 0.04), lim = 4.2) {
  lab <- function(k) sprintf("%.2f", k)
  circ <- ellipse_xy(diag(2), 1, n = 200)
  pts_f <- purrr::map_dfr(k_grid, function(k) {
    X <- app_W %*% t(matrix_power_real(air_A, k))
    tibble(id = seq_len(nrow(X)), x = X[, 1], y = X[, 2], k = lab(k))
  })
  ell_f <- purrr::map_dfr(k_grid, function(k) {
    E <- matrix_power_real(air_A, k) %*% circ
    tibble(x = E[1, ], y = E[2, ], k = lab(k))
  })
  ax <- function(title) list(range = c(-lim, lim), title = title, showgrid = TRUE, gridcolor = "#ebebeb",
                             zeroline = TRUE, zerolinecolor = "#d9d9d9", tickfont = list(size = 13),
                             titlefont = list(size = 15))
  plot_ly() |>
    add_trace(data = pts_f, x = ~x, y = ~y, frame = ~k, ids = ~id, type = "scatter", mode = "markers",
              marker = list(size = 7, color = "#6f7a8c", opacity = 0.55),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(data = ell_f, x = ~x, y = ~y, frame = ~k, type = "scatter", mode = "lines",
              line = list(color = "black", width = 2.2),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    layout(xaxis = ax("headwind shock ε<sub>1</sub>  →  flight time (min)"),
           yaxis = c(ax("load shock ε<sub>2</sub>  →  fuel (100 kg)"), list(scaleanchor = "x")),
           plot_bgcolor = "white", paper_bgcolor = "white",
           margin = list(l = 60, r = 20, t = 20, b = 60), showlegend = FALSE) |>
    animation_opts(frame = 90, transition = 0, redraw = TRUE) |>
    animation_slider(currentvalue = list(prefix = "k = ")) |>
    animation_button(label = "Play")
}

# One choice of A drawn in the language of its own example (axes and arrow names).
# The clouds are illustrations: the flight cloud relabelled, except PCA, which uses the Coke/Pepsi cloud.
fig_choice_example <- function(name) {
  pca_A <- function(S) {
    eg <- eigen(S, symmetric = TRUE); V <- eg$vectors
    if (det(V) < 0) V[, 2] <- -V[, 2]
    if (V[1, 1] < 0) V <- -V
    V %*% diag(sqrt(eg$values))
  }
  spec <- switch(name,
    "Symmetric root" = list(X = app_X, A = app_root, labs = c("student 1's\nown shock", "student 2's\nown shock"),
                            xlab = "grade of student 1", ylab = "grade of student 2 (roommate)"),
    "Cholesky"       = list(X = app_X, A = t(chol(app_Sigma)), labs = c("education shock\n(regression line)", "earnings residual"),
                            xlab = "years of education", ylab = "earnings"),
    "PCA"            = { Xc <- as.matrix(dimred$df[, c("x1", "x2")])
                         list(X = Xc, A = pca_A(cov(Xc)), labs = c("market\n(hot weeks, holidays)", "brand\nswitching"),
                              xlab = "Coca-Cola sales", ylab = "Pepsi sales") },
    "Simultaneous equations" = list(X = app_X, A = air_A, labs = c("headwind", "load"),
                                    xlab = "flight time (min)", ylab = "fuel (100 kg)"))
  S <- cov(spec$X); e <- ellipse_xy(sym_root(S), 1)
  M <- spec$A
  lab_pos <- if (name == "Cholesky") rbind(c(2.75, 0.55), c(-1.7, 1.6)) else t(M) * (1 + 1.05 / sqrt(colSums(M^2)))
  lab <- tibble(x = lab_pos[, 1], y = lab_pos[, 2], txt = spec$labs, col = c(col_a1, col_a2))
  ggplot() + axes0 +
    geom_point(data = tibble(x = spec$X[, 1], y = spec$X[, 2]), aes(x, y), colour = "#6f7a8c", alpha = 0.25, size = 1.4) +
    geom_path(data = tibble(x = e[1, ], y = e[2, ]), aes(x, y), colour = "black", linewidth = 0.9) +
    arrows_of(M, lw = 1.6) +
    geom_label(data = lab, aes(x, y, label = txt, colour = I(col)), size = 4.3, fontface = "bold",
               lineheight = 0.9, label.size = 0, fill = alpha("white", 0.8)) +
    coord_equal(xlim = c(-4.2, 4.2), ylim = c(-4.2, 4.2)) +
    labs(title = name, x = spec$xlab, y = spec$ylab) +
    theme(axis.text = element_blank(), plot.title = element_text(face = "bold"))
}


# The food-critics version of the factor slider: same picture, axes named for the example.
fig_factor_critics <- function() {
  fig <- make_factor_figure(notation = "sigma")
  L <- fig$x$layout
  L$xaxis$title  <- list(text = "critic 1's rating x<sub>1</sub>")
  L$yaxis$title  <- list(text = "critic 2's rating x<sub>2</sub>")
  L$xaxis2$title <- list(text = "critic 1's quirks ε<sub>1</sub>")
  L$yaxis2$title <- list(text = "critic 2's quirks ε<sub>2</sub>")
  L$annotations[[1]]$text <- "<b>Ratings</b> (colour = restaurant quality f)"
  L$annotations[[2]]$text <- "<b>Quirks (ε<sub>1</sub>, ε<sub>2</sub>), quality f as colour</b>"
  fig$x$layout <- L
  fig
}

# The tempting regression: fuel on flight time, with the two "reasons for an extra minute"
fig_regression <- function() {
  b <- num$beta
  dirs <- tibble(lab = c(sprintf("headwind: +1 min,\n+%s kg", kg(fuel_per_min)),
                         sprintf("load: +%s min,\n+%s kg", f2(min_per_t * sd_load), kg(fuel_per_t * sd_load))),
                 xend = c(air_A[1, 1], air_A[1, 2]) * c(1 / air_A[1, 1], 1), yend = c(fuel_per_min, air_A[2, 2]),
                 col = c(col_a1, col_a2))
  ggplot() + cloud_base(alpha = 0.3) +
    geom_abline(slope = b, intercept = 0, colour = "grey30", linewidth = 1, linetype = 2) +
    annotate("text", x = 2.6, y = b * 2.6 + 0.9, label = sprintf("regression slope:\n%s kg per min", kg(b)),
             size = 4.2, hjust = 1, colour = "grey20", lineheight = 0.9) +
    geom_segment(data = dirs, aes(x = 0, y = 0, xend = xend, yend = yend, colour = I(col)), linewidth = 1.5,
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_label(data = dirs, aes(x = xend + c(0.3, -0.25), y = yend + c(-0.75, 0.55), label = lab, colour = I(col),
                                hjust = c(0, 1)), size = 3.9, fontface = "bold", lineheight = 0.9,
               label.size = 0, fill = alpha("white", 0.85)) +
    labs(title = "Regressing fuel on flight time", subtitle = "one slope for two very different reasons")
}

# =============================================================================
# Roadmap: from the reduced form to structural models (the chapters that follow)
# =============================================================================

fig_roadmap <- function() {
  boxes <- tibble(
    x = c(0.9, 3.2, 6.2, 6.2, 6.2, 6.2, 6.2),
    y = c(0, 0, 2.4, 1.2, 0, -1.2, -2.4),
    w = c(1.5, 1.9, 2.6, 2.6, 2.6, 2.6, 2.6),
    title = c("Sample", "Reduced form", "PCA", "MDS", "Factor analysis", "Cluster analysis", "CFA and SEM"),
    sub = c("data matrix 𝐗", "𝑋: μ, Σ", "shocks ordered by variance", "coordinates from distances",
            "few common factors + unique shocks", "discrete latent groups", "restrictions from theory"),
    fill = c("#fdebd5", "#d9f2ec", rep("#f5dfeb", 5)))
  ggplot(boxes) +
    geom_rect(aes(xmin = x - w / 2, xmax = x + w / 2, ymin = y - 0.45, ymax = y + 0.45, fill = I(fill)),
              colour = "grey40") +
    geom_text(aes(x, y + 0.17, label = title), fontface = "bold", size = 4.2) +
    geom_text(aes(x, y - 0.12, label = sub), size = 3.4) +
    annotate("segment", x = 1.7, xend = 2.2, y = 0, yend = 0, colour = col_e, linewidth = 1.1,
             arrow = arrow(length = unit(0.22, "cm"), type = "closed")) +
    annotate("text", x = 1.95, y = 0.25, label = "infer", colour = col_e, fontface = "bold", size = 4) +
    annotate("segment", x = 4.2, xend = 4.85, y = 0, yend = c(2.4, 1.2, 0, -1.2, -2.4), colour = col_h, linewidth = 0.9,
             arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    annotate("text", x = 3.2, y = 1.25, label = "identify:\nstructural models", colour = col_h,
             fontface = "bold", size = 4, lineheight = 0.9) +
    coord_cartesian(xlim = c(0, 7.7), ylim = c(-3, 3.4), expand = FALSE) +
    theme_void()
}
