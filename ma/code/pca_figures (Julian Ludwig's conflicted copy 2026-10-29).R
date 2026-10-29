# =============================================================================
# pca_figures.R  --  data and figures for the chapter "Principal Components Analysis" of the ISQS 6350
# course notes. Put this file and _common.R (shared helpers) in one folder, set it as the working
# directory, and run source("pca_figures.R"). Every figure is a function, e.g. fig_...().
# Packages: tidyverse, plotly, patchwork, mvtnorm, MASS, ggrepel
# =============================================================================

source("_common.R", local = TRUE)
suppressPackageStartupMessages(library(ggrepel))

f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
pct <- function(x, d = 0) paste0(formatC(100 * x, format = "f", digits = d), "%")

# ---- Data -----------------------------------------------------------------------
# Coke and Pepsi: weekly sales (deviations), the two-variable example of the Introduction
cola_S <- matrix(c(1, 0.95, 0.95, 1.1), 2)
cola   <- app_W %*% t(sym_root(cola_S))                 # 200 weeks, sample covariance exactly cola_S
colnames(cola) <- c("coke", "pepsi")
cola_eig <- eigen(cov(cola), symmetric = TRUE)
cola_r2  <- cola_eig$vectors[, 1]^2 * cola_eig$values[1] / diag(cov(cola))   # variance kept by the first arrow

# Flights (Introduction): flight time (minutes) and fuel (100 kg) driven by headwind and load.
# A = M D: D scales the unit-variance shocks into factors with units (headwind in km/h, load in tonnes), M holds the
# effect of one unit of each factor (the aircraft manual).
air_D <- diag(c(20, 2))                                 # sd of headwind (km/h) and of load (tonnes)
air_M <- matrix(c(0.06, 0.06,                           # per km/h of headwind: 0.06 min and 6 kg of fuel
                  0.10, 0.65), 2, 2)                    # per tonne of load: 0.1 min and 65 kg of fuel
air_A <- air_M %*% air_D                                # = the A of the Introduction
air_X <- app_W %*% t(air_A)
air_route_km <- 640; air_cruise <- 800                  # route length and cruise speed: t = 640 / (800 - w) hours
# One example per choice of A, each with factors in meaningful units. Every cloud is built from the same white
# draws app_W, turned by a different angle, so each has sample covariance exactly A A'.
white_turn <- function(deg) app_W %*% rotation_matrix(deg * pi / 180)
# Cholesky: education (years) and earnings ($1000 a year); a regression with slope 4 and residual sd 15
chol_M <- matrix(c(1, 4, 0, 1), 2, 2); chol_D <- diag(c(2.5, 15))
chol_X <- white_turn(25) %*% t(chol_M %*% chol_D)
# Symmetric root: GPAs of two randomly paired roommates; own contributions with sd 0.5, peer effect 0.2
room_A <- sym_root(matrix(c(0.26, 0.10, 0.10, 0.26), 2, 2))     # A = A': the symmetric root of S
room_D <- diag(diag(room_A)); room_M <- room_A %*% solve(room_D)  # ones on the diagonal of M
room_X <- white_turn(70) %*% t(room_A)

# US states (Chapter 2): six indicators
X6 <- state.x77[, c("Income", "Illiteracy", "Life Exp", "Murder", "HS Grad", "Frost")]
colnames(X6) <- c("Income", "Illiteracy", "LifeExp", "Murder", "HSGrad", "Frost")
st_cov <- prcomp(X6)                     # PCA of the covariance matrix
st_cor <- prcomp(X6, scale. = TRUE)      # PCA of the correlation matrix
if (st_cor$rotation["Income", 1] < 0) {  # sign convention: PC1 points towards higher income
  st_cor$rotation[, 1] <- -st_cor$rotation[, 1]; st_cor$x[, 1] <- -st_cor$x[, 1] }
if (st_cor$rotation["Income", 2] < 0) {
  st_cor$rotation[, 2] <- -st_cor$rotation[, 2]; st_cor$x[, 2] <- -st_cor$x[, 2] }
st_lam <- st_cor$sdev^2

# ---- Helpers --------------------------------------------------------------------
cloud_pts <- function(X, alpha = 0.45) geom_point(data = tibble(x = X[, 1], y = X[, 2]), aes(x, y),
                                                   colour = "#6f7a8c", alpha = alpha, size = 1.5)
axes0 <- list(geom_hline(yintercept = 0, colour = "grey88"), geom_vline(xintercept = 0, colour = "grey88"))

# =============================================================================
# 0  Factors with units: A = Gamma D
# =============================================================================

# Three layers of x_t = Gamma D eps_t for one example: the unit-variance shocks (left) are scaled by D into factors
# with units (middle), which Gamma combines into the measured variables (right). Arrows: the basis vectors (left), the
# columns of D (middle) and the columns of A = Gamma D (right), i.e. what one standard deviation of each factor does.
# One observation (orange) is followed through the layers. equal = TRUE draws the middle panel to scale (when
# the factors share the units of the data, as in PCA).
fig_md <- function(X, M, D, f_lab, x_lab, f_sub = NULL, equal = FALSE, obs = NULL) {
  A <- M %*% D; X <- scale(X, scale = FALSE); E <- X %*% t(solve(A)); Fm <- E %*% D
  if (is.null(obs)) obs <- which.min((E[, 1] - 1.3)^2 + (E[, 2] - 1.4)^2)
  col_wk <- "#e6820e"
  arrows_of2 <- function(B) geom_segment(data = tibble(xend = B[1, ], yend = B[2, ], col = c(col_a1, col_a2)),
                                         aes(0, 0, xend = xend, yend = yend, colour = I(col)), linewidth = 1.3,
                                         arrow = arrow(length = unit(0.2, "cm"), type = "closed"))
  pts <- function(P) list(geom_point(data = tibble(x = P[, 1], y = P[, 2]), aes(x, y), colour = "#6f7a8c",
                                     alpha = 0.35, size = 1.2),
                          geom_point(aes(x = P[obs, 1], y = P[obs, 2]), colour = col_wk, size = 3))
  rng <- function(P) { r <- 1.08 * max(abs(P)); c(-r, r) }
  p1 <- ggplot() + axes0 + pts(E) + arrows_of2(diag(2)) + coord_equal(xlim = c(-3.4, 3.4), ylim = c(-3.4, 3.4)) +
    labs(title = expression("Shocks "*bold(epsilon)[t]), subtitle = "variance 1, no units",
         x = expression(epsilon[1]), y = expression(epsilon[2]))
  p2 <- ggplot() + axes0 + pts(Fm) + arrows_of2(D) +
    (if (equal) coord_equal(xlim = rng(Fm), ylim = rng(Fm)) else
      coord_cartesian(xlim = rng(Fm[, 1]), ylim = rng(Fm[, 2]))) +
    labs(title = expression("Factors "*bold(f)[t] == D*bold(epsilon)[t]), subtitle = f_sub, x = f_lab[1], y = f_lab[2])
  p3 <- ggplot() + axes0 + pts(X) + arrows_of2(A) +
    (if (equal) coord_equal(xlim = rng(X), ylim = rng(X)) else
      coord_cartesian(xlim = rng(X[, 1]), ylim = rng(X[, 2]))) +
    labs(title = expression("Data "*bold(x)[t] == Gamma*bold(f)[t]*" = "*Gamma*italic(D)*bold(epsilon)[t]*" = "*italic(A)*bold(epsilon)[t]), subtitle = "deviations from the means",
         x = x_lab[1], y = x_lab[2])
  step <- function(lab) ggplot() +
    annotate("segment", x = 0.05, xend = 0.95, y = 0.5, yend = 0.5, linewidth = 1.1,
             arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
    annotate("text", x = 0.5, y = 0.62, label = lab, size = 5.5, parse = TRUE) +
    coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_void()
  p1 + step("italic(D)") + p2 + step("Gamma") + p3 + plot_layout(widths = c(1, 0.2, 1, 0.2, 1))
}

fig_md_chol <- function() fig_md(chol_X, chol_M, chol_D,
  f_lab = c("education (years)", "earnings residual ($1000)"), x_lab = c("education (years)", "earnings ($1000)"),
  f_sub = "the regressor and the residual")
fig_md_room <- function() fig_md(room_X, room_M, room_D,
  f_lab = c("student 1's own GPA", "student 2's own GPA"), x_lab = c("GPA of student 1", "GPA of student 2"),
  f_sub = "GPA before the peer effect", equal = TRUE)
fig_md_air <- function() fig_md(air_X, air_M, air_D,
  f_lab = c("headwind (km/h)", "load (tonnes)"), x_lab = c("flight time (min)", "fuel (100 kg)"),
  f_sub = "physical units")
fig_md_cola <- function() {
  G <- cola_eig$vectors; if (G[1, 1] < 0) G[, 1] <- -G[, 1]; if (det(G) < 0) G[, 2] <- -G[, 2]
  fig_md(cola, G, diag(sqrt(cola_eig$values)), f_lab = c("score 1 (market)", "score 2 (brand switching)"),
         x_lab = c("Coca-Cola sales", "Pepsi sales"), f_sub = "lengths in the units of the data", equal = TRUE)
}

# =============================================================================
# 1  The direction of greatest variance
# =============================================================================

# Slider (the Introduction's view): every A with AA' = S fits the data, A(phi) = S^(1/2) Q(phi). Left: the cloud,
# the one-sd ellipse and the two columns of A (arrows). Right: the squared length of each arrow, the variance its
# shock contributes, as phi turns. The two always add up to the total variance; PCA picks the angle at which the
# first arrow is longest.
pca_phi <- local({ g <- cola_eig$vectors[, 1]; if (g[1] < 0) g <- -g; round(atan2(g[2], g[1]) * 180 / pi, 1) })
fig_pca_arrows <- function() {
  S <- cov(cola); X <- scale(cola, scale = FALSE); root <- sym_root(S)
  phi_lab <- function(phi) sprintf("%.1f", phi)
  phi_grid <- sort(unique(c(seq(-90, 90, by = 2.5), pca_phi)))
  fk <- function(phi) factor(phi_lab(phi), levels = phi_lab(phi_grid))
  A_of <- function(phi) root %*% rotation_matrix(phi * pi / 180)
  len2 <- function(phi) colSums(A_of(phi)^2)
  cols <- c(col_a1, col_a2)
  arrows <- lapply(1:2, function(j) purrr::map_dfr(phi_grid, function(ph)
    arrow_geometry(A_of(ph)[, j], head_length = 0.2, head_width = 0.1, shaft_width = 0.022)$arrow |> mutate(k = fk(ph))))
  labs_df <- purrr::map_dfr(phi_grid, function(ph) { M <- A_of(ph); r <- sqrt(colSums(M^2))
    tibble(x = M[1, ] * (1 + 0.35 / r), y = M[2, ] * (1 + 0.35 / r), label = c("a<sub>1</sub>", "a<sub>2</sub>"),
           col = cols, k = fk(ph)) })
  ell <- ellipse_xy(root, 1)
  curve <- tibble(phi = seq(-90, 90, by = 0.5)) |>
    mutate(l1 = purrr::map_dbl(phi, ~ len2(.x)[1]), l2 = purrr::map_dbl(phi, ~ len2(.x)[2]))
  marks <- purrr::map_dfr(phi_grid, function(ph) { l <- len2(ph)
    tibble(phi = ph, l = l, col = cols, k = fk(ph)) })
  tot <- sum(diag(S))
  readout <- purrr::map_dfr(phi_grid, function(ph) { l <- len2(ph)
    tibble(x = 0, y = tot * 1.2, k = fk(ph),
           t = sprintf("‖a<sub>1</sub>‖² = %.2f (%.0f%%)    ‖a<sub>2</sub>‖² = %.2f (%.0f%%)    total = %.2f",
                       l[1], 100 * l[1] / tot, l[2], 100 * l[2] / tot, tot)) })

  p1 <- plot_ly() |>
    add_trace(x = X[, 1], y = X[, 2], type = "scatter", mode = "markers",
              marker = list(size = 5, color = "#6f7a8c", opacity = 0.4), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(x = ell[1, ], y = ell[2, ], type = "scatter", mode = "lines", line = list(color = "black", width = 2),
              hoverinfo = "skip", showlegend = FALSE)
  for (j in 1:2) p1 <- add_arrow_trace(p1, arrows[[j]], color = cols[j])
  p1 <- p1 |> add_trace(data = labs_df, x = ~x, y = ~y, frame = ~k, text = ~label, type = "scatter", mode = "text",
                        textfont = list(size = 16, color = ~col), hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)
  p2 <- plot_ly() |>
    add_trace(x = c(-90, 90), y = c(tot, tot), type = "scatter", mode = "lines",
              line = list(color = "grey", dash = "dot"), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(x = -88, y = tot * 1.07, text = "total variance tr(S)", type = "scatter", mode = "text",
              textposition = "middle right", textfont = list(size = 11, color = "grey"), hoverinfo = "skip",
              showlegend = FALSE) |>
    add_trace(data = curve, x = ~phi, y = ~l1, type = "scatter", mode = "lines", line = list(color = col_a1, width = 2.5),
              hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = curve, x = ~phi, y = ~l2, type = "scatter", mode = "lines", line = list(color = col_a2, width = 2.5),
              hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = marks, x = ~phi, y = ~l, frame = ~k, type = "scatter", mode = "markers",
              marker = list(size = 12, color = ~col), hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(data = readout, x = ~x, y = ~y, frame = ~k, text = ~t, type = "scatter", mode = "text",
              textfont = list(size = 12, color = "#2c3e50"), hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)
  jump <- list(
    list(label = "PCA: first arrow longest", method = "animate",
         args = list(list(phi_lab(pca_phi)), list(mode = "immediate", frame = list(duration = 0, redraw = TRUE),
                                                   transition = list(duration = 0)))),
    list(label = "Symmetric root", method = "animate",
         args = list(list(phi_lab(0)), list(mode = "immediate", frame = list(duration = 0, redraw = TRUE),
                                            transition = list(duration = 0)))))
  subplot(p1, p2, widths = c(0.5, 0.5), margin = 0.07, titleX = TRUE, titleY = TRUE) |>
    layout(xaxis  = list(title = "Coca-Cola sales (deviation)", range = c(-3.2, 3.2), zeroline = FALSE),
           yaxis  = list(title = "Pepsi sales (deviation)", range = c(-3.2, 3.2), scaleanchor = "x", zeroline = FALSE),
           xaxis2 = list(title = "rotation φ (degrees)", range = c(-90, 90), dtick = 45, zeroline = FALSE),
           yaxis2 = list(title = "variance carried by each arrow", range = c(0, tot * 1.3)),
           margin = list(t = 70, b = 150)) |>
    animation_opts(frame = 80, transition = 0, redraw = TRUE) |>
    animation_slider(currentvalue = list(prefix = "rotation φ = ", suffix = "°")) |>
    animation_button(x = 0, xanchor = "left", y = -0.3, yanchor = "top", label = "Play") |>
    finalize_animation(jump, start_frame = "0.0", prefix = "rotation φ = ", suffix = "°")
}

# Slider: turn a unit vector gamma = (cos phi, sin phi)'. Left: the cloud, the line, and the projections
# (scores) on it. Right: the variance of the scores, gamma' S gamma, as a function of the angle.
fig_pca_direction <- function() {
  S <- cov(cola); X <- scale(cola, scale = FALSE)
  phis <- seq(0, 175, by = 5)
  v_of <- function(phi) c(cos(phi * pi / 180), sin(phi * pi / 180))
  curve <- tibble(phi = seq(0, 180, by = 1)) |>
    mutate(v = purrr::map_dbl(phi, function(p) { g <- v_of(p); drop(t(g) %*% S %*% g) }))
  lab <- function(p) p                                               # numeric frames: slider in numeric order
  fr <- purrr::map_dfr(phis, function(p) {
    g <- v_of(p); s <- drop(X %*% g)
    tibble(frame = lab(p), x = X[, 1], y = X[, 2], px = s * g[1], py = s * g[2])
  })
  lines <- purrr::map_dfr(phis, function(p) { g <- v_of(p)
    tibble(frame = lab(p), x = c(-4 * g[1], 4 * g[1]), y = c(-4 * g[2], 4 * g[2])) })
  segs <- fr |> mutate(id = rep(seq_len(nrow(X)), length(phis)))
  seg_df <- segs |> group_by(frame) |> group_modify(~ {
    tibble(x = as.vector(rbind(.x$x, .x$px, NA)), y = as.vector(rbind(.x$y, .x$py, NA))) }) |> ungroup()
  mark <- tibble(frame = lab(phis), phi = phis,
                 v = purrr::map_dbl(phis, function(p) { g <- v_of(p); drop(t(g) %*% S %*% g) }))
  txt <- mark |> mutate(t = sprintf("angle %d°: variance of the scores = %.2f", phi, v))

  p1 <- plot_ly() |>
    add_trace(data = seg_df, x = ~x, y = ~y, frame = ~frame, type = "scatter", mode = "lines",
              line = list(color = "rgba(150,150,150,0.45)", width = 1), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = fr, x = ~x, y = ~y, frame = ~frame, type = "scatter", mode = "markers",
              marker = list(color = "rgba(111,122,140,0.55)", size = 5), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = lines, x = ~x, y = ~y, frame = ~frame, type = "scatter", mode = "lines",
              line = list(color = col_a1, width = 3), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = fr, x = ~px, y = ~py, frame = ~frame, type = "scatter", mode = "markers",
              marker = list(color = col_a1, size = 4), hoverinfo = "skip", showlegend = FALSE)
  p2 <- plot_ly() |>
    add_trace(data = curve, x = ~phi, y = ~v, type = "scatter", mode = "lines",
              line = list(color = "black", width = 2), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = mark, x = ~phi, y = ~v, frame = ~frame, type = "scatter", mode = "markers",
              marker = list(color = col_a1, size = 12), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = txt, x = I(90), y = I(max(curve$v) * 1.13), text = ~t, frame = ~frame, type = "scatter",
              mode = "text", textfont = list(size = 13, color = col_a1), hoverinfo = "skip", showlegend = FALSE)
  subplot(p1, p2, widths = c(0.5, 0.5), margin = 0.07, titleX = TRUE, titleY = TRUE) |>
    layout(xaxis  = list(title = "Coca-Cola sales (deviation)", range = c(-3.6, 3.6), zeroline = FALSE),
           yaxis  = list(title = "Pepsi sales (deviation)", range = c(-3.6, 3.6), scaleanchor = "x", zeroline = FALSE),
           xaxis2 = list(title = "angle of the line (degrees)", range = c(0, 180), dtick = 45),
           yaxis2 = list(title = "variance of the scores", range = c(0, max(curve$v) * 1.2)),
           margin = list(t = 20)) |>
    animation_opts(frame = 120, transition = 0, redraw = FALSE) |>
    animation_slider(currentvalue = list(prefix = "angle: ", suffix = "°")) |>
    animation_button(x = 0, xanchor = "left", y = -0.2, yanchor = "top", label = "Play")
}

# The first principal component line minimises the perpendicular distances; the regression line the
# vertical ones.
fig_pca_vs_ols <- function() {
  X <- scale(air_X, scale = FALSE); S <- cov(X)
  g <- eigen(S, symmetric = TRUE)$vectors[, 1]; if (g[1] < 0) g <- -g
  b <- S[1, 2] / S[1, 1]
  set.seed(3); idx <- sample(nrow(X), 45)
  d <- tibble(x = X[idx, 1], y = X[idx, 2]) |>
    mutate(s = x * g[1] + y * g[2], px = s * g[1], py = s * g[2], fy = b * x)
  base <- ggplot(d) + axes0 + coord_equal(xlim = c(-4, 4), ylim = c(-4.5, 4.5)) +
    labs(x = expression("flight time"~italic(x)[1]), y = expression("fuel"~italic(x)[2]))
  p1 <- base +
    geom_abline(slope = g[2] / g[1], intercept = 0, colour = col_a1, linewidth = 1.1) +
    geom_segment(aes(x, y, xend = px, yend = py), colour = col_a1, alpha = 0.6) +
    geom_point(aes(x, y), colour = "grey30", size = 1.6) +
    labs(title = "First principal component", subtitle = "minimises the perpendicular distances")
  p2 <- base +
    geom_abline(slope = b, intercept = 0, colour = col_a2, linewidth = 1.1) +
    geom_segment(aes(x, y, xend = x, yend = fy), colour = col_a2, alpha = 0.6) +
    geom_point(aes(x, y), colour = "grey30", size = 1.6) +
    labs(title = "Regression of fuel on time", subtitle = "minimises the vertical distances")
  p1 + p2
}

# Dropping the short arrow: x_t is predicted by eps_t1 a_1 alone. Left: each week (grey) and its prediction on the
# line of the first arrow (blue), joined by the error, which is the dropped part eps_t2 a_2. Right: predicted against
# actual sales for each brand; the closer to the 45-degree line, the better.
fig_pca_drop <- function() {
  X <- scale(cola, scale = FALSE); e <- eigen(cov(X), symmetric = TRUE)
  G <- e$vectors; if (G[1, 1] < 0) G[, 1] <- -G[, 1]
  a1 <- G[, 1] * sqrt(e$values[1]); a2 <- G[, 2] * sqrt(e$values[2])
  eps1 <- (X %*% G[, 1]) / sqrt(e$values[1])
  P <- eps1 %*% t(a1)
  d <- tibble(x = X[, 1], y = X[, 2], px = P[, 1], py = P[, 2])
  r2 <- 1 - colSums((X - P)^2) / colSums(X^2)
  p1 <- ggplot(d) + axes0 +
    geom_abline(slope = a1[2] / a1[1], intercept = 0, colour = col_a1, alpha = 0.5) +
    geom_segment(aes(x, y, xend = px, yend = py), colour = "grey55", linewidth = 0.3) +
    geom_point(aes(x, y), colour = "#6f7a8c", alpha = 0.5, size = 1.4) +
    geom_point(aes(px, py), colour = col_a1, size = 1.4) +
    coord_equal(xlim = c(-3.2, 3.2), ylim = c(-3.2, 3.2)) +
    labs(x = expression("Coca-Cola"~italic(x)[1]), y = expression("Pepsi"~italic(x)[2]),
         title = expression("Second arrow off: " * hat(bold(x))[t] == epsilon[t1] * italic(a)[1]),
         subtitle = "grey: actual, blue: predicted")
  lab <- c(coke = sprintf("Coca-Cola (R\u00b2 = %s)", f2(r2[1])), pepsi = sprintf("Pepsi (R\u00b2 = %s)", f2(r2[2])))
  d2 <- bind_rows(tibble(brand = lab[["coke"]], actual = X[, 1], pred = P[, 1]),
                  tibble(brand = lab[["pepsi"]], actual = X[, 2], pred = P[, 2]))
  p2 <- ggplot(d2, aes(actual, pred, colour = brand)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey60", linetype = "22") +
    geom_point(alpha = 0.6, size = 1.4) +
    scale_colour_manual(values = c("#b2182b", "#2166ac"), name = NULL) +
    coord_equal(xlim = c(-3.2, 3.2), ylim = c(-3.2, 3.2)) +
    labs(x = "actual sales (deviation)", y = "predicted from the first arrow",
         title = "Predicted against actual", subtitle = "dashed line: perfect") +
    theme(legend.position = c(0.3, 0.88), legend.background = element_rect(fill = "white", colour = NA))
  p1 + p2
}

# =============================================================================
# 2  Eigenvectors: PCA as a rotation (and rescaling) of the cloud
# =============================================================================

# PCA as A = Gamma D, for Coke and Pepsi, in two directions. One week (orange) is followed through all panels; the
# dashed lines show its coordinates. pca_layers() holds the shared pieces.
pca_layers <- function() {
  X <- scale(cola, scale = FALSE); e <- eigen(cov(X), symmetric = TRUE)
  G <- e$vectors; if (G[1, 1] < 0) G[, 1] <- -G[, 1]; if (det(G) < 0) G[, 2] <- -G[, 2]
  D <- diag(sqrt(e$values)); Fs <- X %*% G; E <- sweep(Fs, 2, sqrt(e$values), "/")
  list(X = X, G = G, D = D, Fs = Fs, E = E, week = which.min((E[, 1] - 1.3)^2 + (E[, 2] - 1.6)^2))
}
pca_panel <- function(P, Mell, Dcoord, week, title, sub, xl, yl, pc_lines = NULL, arrows_df = NULL, lim = c(-3.4, 3.4)) {
  col_wk <- "#e6820e"; pt <- P[week, ]
  el <- ellipse_xy(Mell, 1)
  c1 <- sum(pt * Dcoord[, 1]) * Dcoord[, 1]; c2 <- sum(pt * Dcoord[, 2]) * Dcoord[, 2]
  g <- ggplot() + axes0
  if (!is.null(pc_lines)) g <- g + geom_abline(slope = pc_lines[2, ] / pc_lines[1, ], intercept = 0,
                                               colour = c(col_a1, col_a2), alpha = 0.3)
  g <- g + cloud_pts(P, 0.35) + geom_path(data = tibble(x = el[1, ], y = el[2, ]), aes(x, y), linewidth = 0.8)
  if (!is.null(arrows_df)) g <- g +
    geom_segment(data = arrows_df, aes(0, 0, xend = x, yend = y, colour = I(col)), linewidth = 1.3,
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_label(data = filter(arrows_df, !is.na(lab)), aes(x * 1.18 + lx, y * 1.18 + ly, label = lab, colour = I(col)),
               size = 3.4, fontface = "bold", label.size = 0, fill = alpha("white", 0.8))
  g + geom_segment(aes(x = pt[1], y = pt[2], xend = c1[1], yend = c1[2]), linetype = "22", colour = col_wk) +
    geom_segment(aes(x = pt[1], y = pt[2], xend = c2[1], yend = c2[2]), linetype = "22", colour = col_wk) +
    geom_point(aes(x = pt[1], y = pt[2]), colour = col_wk, size = 3.2) +
    coord_equal(xlim = lim, ylim = lim) + labs(title = title, subtitle = sub, x = xl, y = yl)
}
pca_step <- function(lab, col = "black") ggplot() +
  annotate("segment", x = 0.05, xend = 0.95, y = 0.5, yend = 0.5, linewidth = 1.1, colour = col,
           arrow = arrow(length = unit(0.25, "cm"), type = "closed")) +
  annotate("text", x = 0.5, y = 0.63, label = lab, parse = TRUE, size = 5.5, colour = col) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1)) + theme_void()
shock_arrows <- function(B) tibble(x = B[1, ], y = B[2, ], col = c(col_a1, col_a2), lab = NA_character_, lx = 0, ly = 0)

# Building the data, x_t = Gamma D eps_t read from right to left: the round cloud of shocks (left) is stretched
# along the axes by D into the scores (middle), then rotated by Gamma into the data (right). The arrows follow the
# two shocks: the basis vectors, the columns of D, the columns of A = Gamma D.
fig_pca_forward <- function() {
  L <- pca_layers()
  p1 <- pca_panel(L$E, diag(2), diag(2), L$week, expression("Shocks "*bold(epsilon)[t]), "variance 1: round",
                  expression(epsilon[1]*" (market)"), expression(epsilon[2]*" (brand switching)"),
                  arrows_df = shock_arrows(diag(2)))
  p2 <- pca_panel(L$Fs, L$D, diag(2), L$week, expression("Scores "*bold(f)[t] == italic(D)*bold(epsilon)[t]),
                  "stretched along the axes", expression(italic(f)[1]), expression(italic(f)[2]),
                  arrows_df = shock_arrows(L$D))
  p3 <- pca_panel(L$X, L$G %*% L$D, diag(2), L$week, expression("Data "*bold(x)[t] == Gamma*bold(f)[t]*" = "*Gamma*italic(D)*bold(epsilon)[t]*" = "*italic(A)*bold(epsilon)[t]),
                  "rotated into place", "Coca-Cola (deviation)", "Pepsi (deviation)",
                  arrows_df = shock_arrows(L$G %*% L$D))
  p1 + pca_step("italic(D)") + p2 + pca_step("Gamma") + p3 + plot_layout(widths = c(1, 0.2, 1, 0.2, 1))
}

# Taking the data apart, eps_t = D^(-1) Gamma' x_t: the data (left) are rotated back by Gamma' into the scores
# (middle), and the stretch is undone by D^(-1) (right). The purple arrows follow the two variables: one unit of
# Coca-Cola and one unit of Pepsi (the axes of the data) land, after Gamma', on the rows of Gamma. The middle panel,
# scores as points and those arrows, is the biplot. The right panel is the standardised biplot: its arrows are
# the rows of A = Gamma D, so that x_ti = a_i' eps_t still holds (the stretch moves from the points to the arrows).
# Dashed lines: the week's coordinates; in the data panel they
# are its projections onto the two components.
fig_pca_inverse <- function() {
  L <- pca_layers(); col_v <- "#7b4fb5"
  var_arrows <- function(B) tibble(x = B[1, ], y = B[2, ], col = col_v, lab = c("Coca-Cola", "Pepsi"),
                                   lx = c(0.35, -0.2), ly = c(-0.25, 0.3))
  p1 <- pca_panel(L$X, L$G %*% L$D, L$G, L$week, expression("Data "*bold(x)[t]), "axes = the two variables",
                  "Coca-Cola (deviation)", "Pepsi (deviation)", pc_lines = L$G, arrows_df = var_arrows(diag(2)))
  p2 <- pca_panel(L$Fs, L$D, diag(2), L$week, expression("Scores "*bold(f)[t] == Gamma*minute*bold(x)[t]),
                  "rotated back: the biplot", expression(italic(f)[1]), expression(italic(f)[2]),
                  arrows_df = var_arrows(t(L$G)))
  p3 <- pca_panel(L$E, diag(2), diag(2), L$week, expression("Shocks "*bold(epsilon)[t] == italic(D)^{-1}*bold(f)[t]),
                  "stretch undone: standardised scores", expression(epsilon[1]), expression(epsilon[2]),
                  arrows_df = var_arrows(t(L$G %*% L$D)))                 # rows of A: the standardised biplot
  col_wk <- "#e6820e"
  p1 + pca_step("Gamma*minute", col_wk) + p2 + pca_step("italic(D)^{-1}", col_wk) + p3 +
    plot_layout(widths = c(1, 0.2, 1, 0.2, 1))
}

# What an eigenvector is: S as a movement of the plane. The slider applies the transformation gradually,
# T = (1 - k) I + k S, from the identity (k = 0) to S (k = 1). Points sit on lines through the origin (dashed). Points
# on the two eigenvector lines (blue, red) stay on their lines and are only stretched or squished, by lambda_1 and
# lambda_2 at k = 1; points on every other line (grey, the horizontal line in black) are knocked off their line.
fig_eigen_turn <- function() {
  S <- cov(cola); e <- eigen(S, symmetric = TRUE)
  ang_e <- c(pca_phi, pca_phi + 90)
  angs <- c(0, 22.5, 67.5, 90, 112.5, 157.5, ang_e)
  kind <- c("x-axis", rep("other", 5), "eig1", "eig2")
  cols <- c("x-axis" = "#222222", other = "#a6a6a6", eig1 = col_a1, eig2 = col_a2)
  r <- seq(-1.5, 1.5, by = 0.25); r <- r[r != 0]
  base <- purrr::map_dfr(seq_along(angs), function(i) {
    u <- c(cos(angs[i] * pi / 180), sin(angs[i] * pi / 180))
    tibble(line = i, kind = kind[i], x0 = r * u[1], y0 = r * u[2]) })
  ks <- seq(0, 1, by = 0.05)
  fr <- purrr::map_dfr(ks, function(k) { Tk <- (1 - k) * diag(2) + k * S
    P <- cbind(base$x0, base$y0) %*% t(Tk); mutate(base, frame = k, x = P[, 1], y = P[, 2]) })
  lines <- purrr::map_dfr(ks, function(k) purrr::map_dfr(seq_along(angs), function(i) {
    u <- c(cos(angs[i] * pi / 180), sin(angs[i] * pi / 180))
    tibble(frame = k, kind = kind[i], x = c(-3.4 * u[1], 3.4 * u[1], NA), y = c(-3.4 * u[2], 3.4 * u[2], NA)) }))
  lab <- tibble(frame = ks, t = ifelse(ks == 0, "identity: nothing moves",
                                       sprintf("blue line stretched by %.2f, red line by %.2f", 1 - ks + ks * e$values[1],
                                               1 - ks + ks * e$values[2])))
  p <- plot_ly()
  for (kd in names(cols)) {
    p <- p |>
      add_trace(data = filter(lines, kind == kd), x = ~x, y = ~y, frame = ~frame, type = "scatter", mode = "lines",
                line = list(color = cols[[kd]], width = if (kd == "other") 1 else 1.5, dash = "dash"),
                opacity = if (kd == "other") 0.5 else 0.8, hoverinfo = "skip", showlegend = FALSE) |>
      add_trace(data = filter(fr, kind == kd), x = ~x, y = ~y, frame = ~frame, type = "scatter", mode = "markers",
                marker = list(color = cols[[kd]], size = if (kd == "other") 6 else 8), hoverinfo = "skip",
                showlegend = FALSE)
  }
  p |>
    add_trace(data = lab, x = I(0), y = I(3.05), text = ~t, frame = ~frame, type = "scatter", mode = "text",
              textfont = list(size = 14), hoverinfo = "skip", showlegend = FALSE) |>
    layout(xaxis = list(title = "first coordinate", range = c(-3.4, 3.4), zeroline = FALSE),
           yaxis = list(title = "second coordinate", range = c(-3.4, 3.4), scaleanchor = "x", zeroline = FALSE),
           margin = list(t = 10)) |>
    animation_opts(frame = 150, transition = 0, redraw = FALSE) |>
    animation_slider(currentvalue = list(prefix = "how much of S is applied: ")) |>
    animation_button(x = 0, xanchor = "left", y = -0.12, yanchor = "top", label = "Play")
}


# =============================================================================
# 3  Covariance or correlation matrix?
# =============================================================================

# Slider: record Pepsi in different units (multiply by k). Left: the covariance-matrix PCA arrows swing towards
# the variable with the larger variance. Right: after standardising, the picture does not depend on k.
fig_scale_slider <- function() {
  X <- scale(cola, scale = FALSE); Z <- scale(cola)
  ks <- round(exp(seq(log(1 / 3), log(3), length.out = 41)), 3)
  lab <- function(k) sprintf("%.3f", k)
  pca_A <- function(M) { e <- eigen(cov(M), symmetric = TRUE); G <- e$vectors
    if (G[1, 1] + G[2, 1] < 0) G[, 1] <- -G[, 1]; if (G[2, 2] < 0) G[, 2] <- -G[, 2]
    list(A = G %*% diag(sqrt(e$values)), share = e$values[1] / sum(e$values), ang = atan2(G[2, 1], G[1, 1]) * 180 / pi) }
  cols <- c(col_a1, col_a2)
  pts <- purrr::map_dfr(ks, function(k) tibble(x = X[, 1], y = k * X[, 2], f = lab(k)))
  arrows <- lapply(1:2, function(j) purrr::map_dfr(ks, function(k) {
    A <- pca_A(cbind(X[, 1], k * X[, 2]))$A
    arrow_geometry(A[, j], head_length = 0.35, head_width = 0.18, shaft_width = 0.04)$arrow |> mutate(k = lab(k)) }))
  txt <- purrr::map_dfr(ks, function(k) { r <- pca_A(cbind(X[, 1], k * X[, 2]))
    tt <- sprintf("Pepsi × %.2f: PC1 at %.0f°, %.0f%% of the variance", k, r$ang, 100 * r$share)
    tibble(x = -7.6, y = 9.2, t = tt, k = lab(k)) })
  Ar <- pca_A(Z)$A; ar <- lapply(1:2, function(j) arrow_geometry(Ar[, j], head_length = 0.2, head_width = 0.1,
                                                                  shaft_width = 0.022)$arrow)
  p1 <- plot_ly() |>
    add_trace(data = pts, x = ~x, y = ~y, frame = ~f, type = "scatter", mode = "markers",
              marker = list(size = 4, color = "#6f7a8c", opacity = 0.45), hoverinfo = "skip", showlegend = FALSE)
  for (j in 1:2) p1 <- add_arrow_trace(p1, arrows[[j]], color = cols[j])
  p1 <- p1 |> add_trace(data = txt, x = ~x, y = ~y, frame = ~k, text = ~t, type = "scatter", mode = "text",
                        textposition = "middle right", textfont = list(size = 12, color = "#2c3e50"),
                        hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)
  p2 <- plot_ly() |>
    add_trace(x = Z[, 1], y = Z[, 2], type = "scatter", mode = "markers",
              marker = list(size = 4, color = "#6f7a8c", opacity = 0.45), hoverinfo = "skip", showlegend = FALSE)
  for (j in 1:2) p2 <- p2 |> add_trace(data = ar[[j]], x = ~x, y = ~y, type = "scatter", mode = "lines", fill = "toself",
                                       fillcolor = cols[j], line = list(width = 1, color = cols[j]),
                                       hoverinfo = "skip", showlegend = FALSE)
  subplot(p1, p2, widths = c(0.6, 0.4), margin = 0.07, titleX = TRUE, titleY = TRUE) |>
    layout(xaxis  = list(title = "Coca-Cola sales (deviation)", range = c(-8, 8), zeroline = FALSE),
           yaxis  = list(title = "Pepsi sales × k (deviation)", range = c(-10, 10), scaleanchor = "x", zeroline = FALSE),
           xaxis2 = list(title = "Coca-Cola (standardised)", range = c(-3.5, 3.5), zeroline = FALSE),
           yaxis2 = list(title = "Pepsi (standardised)", range = c(-3.5, 3.5), scaleanchor = "x2", zeroline = FALSE),
           annotations = list(
             list(text = "<b>Covariance matrix</b>: depends on the units", x = 0, y = 1.06, xref = "paper", yref = "paper",
                  xanchor = "left", showarrow = FALSE, font = list(size = 13)),
             list(text = "<b>Correlation matrix</b>: the same for every k", x = 0.62, y = 1.06, xref = "paper", yref = "paper",
                  xanchor = "left", showarrow = FALSE, font = list(size = 13))),
           margin = list(t = 50, b = 120)) |>
    animation_opts(frame = 100, transition = 0, redraw = TRUE) |>
    animation_slider(currentvalue = list(prefix = "Pepsi units: × ")) |>
    animation_button(x = 0, xanchor = "left", y = -0.25, yanchor = "top", label = "Play") |>
    finalize_animation(list(list(label = "original units (× 1)", method = "animate",
                                  args = list(list("1.000"), list(mode = "immediate", frame = list(duration = 0, redraw = TRUE),
                                                                  transition = list(duration = 0))))),
                       start_frame = "1.000", prefix = "Pepsi units: × ")
}

# Weights of the first component: from the covariance matrix (income dominates) and from the correlation
# matrix (all six variables count)
fig_scale <- function() {
  w_cov <- st_cov$rotation[, 1]; if (w_cov["Income"] < 0) w_cov <- -w_cov
  d <- bind_rows(
    tibble(var = names(w_cov), w = w_cov, which = sprintf("Covariance matrix: PC1 = %s of the variance",
                                                         pct(st_cov$sdev[1]^2 / sum(st_cov$sdev^2), 1))),
    tibble(var = names(w_cov), w = st_cor$rotation[, 1], which = sprintf("Correlation matrix: PC1 = %s of the variance",
                                                                        pct(st_lam[1] / sum(st_lam))))) |>
    mutate(var = factor(var, levels = rev(colnames(X6))), which = factor(which, levels = unique(which)))
  sds <- tibble(var = factor(colnames(X6), levels = rev(colnames(X6))), sd = apply(X6, 2, sd))
  p0 <- ggplot(sds, aes(sd, var)) + geom_col(fill = "grey60", width = 0.6) +
    geom_text(aes(label = formatC(sd, format = "f", digits = 1)), hjust = -0.15, size = 3.6) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.25))) +
    labs(title = "Standard deviations", x = "sd (in the variable's units)", y = NULL)
  p1 <- ggplot(d, aes(w, var, fill = w > 0)) + geom_col(width = 0.6, show.legend = FALSE) +
    geom_vline(xintercept = 0, colour = "grey40") +
    scale_fill_manual(values = c(`TRUE` = col_a1, `FALSE` = col_a2)) +
    facet_wrap(~ which, ncol = 1) + xlim(-1, 1) +
    labs(title = "Weights of the first component", x = expression("weight in"~gamma[1]), y = NULL)
  p0 + p1 + plot_layout(widths = c(0.8, 1.2))
}


# =============================================================================
# 4  How many components?
# =============================================================================
fig_scree <- function() {
  mk <- function(lam, name) tibble(data = name, j = seq_along(lam), lambda = lam, cum = cumsum(lam) / sum(lam))
  d <- mk(st_lam, "US states (6 variables)")
  ggplot(d, aes(j, lambda)) +
    geom_hline(yintercept = 1, colour = col_a2, linetype = 2) +
    geom_hline(yintercept = 0.7, colour = col_a2, linetype = 3) +
    geom_line(linewidth = 1) + geom_point(size = 2.6) +
    geom_text(aes(label = paste0(round(100 * cum), "%")), nudge_y = 0.3, size = 3.3, colour = "grey30") +
    geom_text(data = tibble(y = c(1, 0.7), lab = c("1 (Kaiser)", "0.7 (Jolliffe)")),
              aes(x = 6, y = y, label = lab), hjust = 1, vjust = c(-0.4, 1.4), colour = col_a2, size = 3.3) +
    scale_x_continuous(breaks = 1:6) +
    labs(x = "component", y = expression("eigenvalue"~lambda[j]~"(variance of the component)"),
         subtitle = "Labels: cumulative share of the total variance")
}

# =============================================================================
# 5  Reconstruction from m components
# =============================================================================
fig_reconstruct <- function(states = c("Alaska", "Hawaii", "Texas")) {
  Z <- scale(X6); G <- st_cor$rotation; Y <- st_cor$x
  rec <- function(s, m) drop(G[, 1:m, drop = FALSE] %*% Y[s, 1:m])   # x-hat = sum of the first m score * weight
  d <- purrr::map_dfr(states, function(s) bind_rows(
    tibble(state = s, var = colnames(X6), m = "actual", z = Z[s, ]),
    purrr::map_dfr(1:3, function(k) tibble(state = s, var = colnames(X6), z = rec(s, k), m = paste0("m = ", k))))) |>
    mutate(var = factor(var, levels = colnames(X6)), state = factor(state, levels = states))
  ggplot(d, aes(var, z)) +
    geom_hline(yintercept = 0, colour = "grey60") +
    geom_col(data = filter(d, m == "actual"), fill = "grey80", width = 0.7) +
    geom_point(data = filter(d, m != "actual"), aes(colour = m, shape = m), size = 2.6,
               position = position_dodge(width = 0.55)) +
    scale_colour_manual(values = c("m = 1" = col_a1, "m = 2" = "#7b4fb5", "m = 3" = col_a2), name = "reconstruction") +
    scale_shape_manual(values = c("m = 1" = 16, "m = 2" = 17, "m = 3" = 15), name = "reconstruction") +
    facet_wrap(~ state, ncol = 3) +
    labs(x = NULL, y = "standardised value (bars: actual)") +
    theme(axis.text.x = element_text(angle = 40, hjust = 1), legend.position = "top")
}


# Biplot from x_t = Gamma f_t: the observations are points at their scores f_t (the factors of A = Gamma D, in the units
# of the data), the variables are arrows from the rows of Gamma (the weights). The inner product of a point with
# an arrow gives the observation's value on that variable (exactly with all components, approximately with two),
# and distances between points approximate Euclidean distances between observations. The arrows are drawn
# arrow_scale times longer so they are visible next to the points; the top and right axes give their own scale.
# (R's biplot(prcomp(...), scale = 0) draws this version; its default, scale = 1, draws the standardised scores
# eps_t with the rows of A = Gamma D.)
biplot_gg <- function(Y, G, lam, point_labels = NULL, label_all = FALSE, arrow_scale = 1,
                      var_col = col_a2, pt_col = col_a1) {
  L <- G[, 1:2] * arrow_scale
  d <- tibble(id = rownames(Y), y1 = Y[, 1], y2 = Y[, 2])
  d$lab <- if (label_all) d$id else if (is.null(point_labels)) NA_character_ else unname(point_labels[d$id])
  a <- tibble(var = rownames(G), x = L[, 1], y = L[, 2])
  sec <- function(nm) sec_axis(~ . / arrow_scale, name = nm)
  ggplot() + axes0 +
    geom_point(data = d, aes(y1, y2), colour = pt_col, size = 1.8, alpha = 0.7) +
    geom_text_repel(data = filter(d, !is.na(lab)), aes(y1, y2, label = lab), size = 3.2, colour = pt_col,
                    seed = 4, max.overlaps = 40) +
    geom_segment(data = a, aes(0, 0, xend = x, yend = y), colour = var_col, linewidth = 1,
                 arrow = arrow(length = unit(0.2, "cm"))) +
    geom_label_repel(data = a, aes(x, y, label = var), colour = var_col, size = 3.8, fontface = "bold",
                     label.size = 0, fill = alpha("white", 0.8), seed = 7, min.segment.length = 10,
                     box.padding = 0.15, point.padding = 0) +
    scale_x_continuous(sec.axis = sec(NULL)) + scale_y_continuous(sec.axis = sec(NULL)) +
    coord_equal() +
    theme(axis.title.x.top = element_text(colour = var_col), axis.title.y.right = element_text(colour = var_col),
          axis.text.x.top = element_text(colour = var_col), axis.text.y.right = element_text(colour = var_col)) +
    labs(x = sprintf("PC1 score (%s of the variance)", pct(lam[1] / sum(lam))),
         y = sprintf("PC2 score (%s)", pct(lam[2] / sum(lam))))
}


# The states (correlation matrix): scores and the six variables
fig_states_biplot <- function() {
  Y <- st_cor$x; rownames(Y) <- state.abb
  biplot_gg(Y, st_cor$rotation, st_lam, label_all = TRUE, arrow_scale = 4)
}

# =============================================================================
# 6  PCA is one choice of A
# =============================================================================

# Flights: the PCA arrows next to the true headwind and load arrows (the manual's restriction)
fig_pca_vs_truth <- function() {
  S <- cov(air_X); e <- eigen(S, symmetric = TRUE)
  G <- e$vectors; if (G[1, 1] < 0) G[, 1] <- -G[, 1]; if (det(G) < 0) G[, 2] <- -G[, 2]
  A_pca <- G %*% diag(sqrt(e$values))
  el <- ellipse_xy(sym_root(S), 1); el <- tibble(x = el[1, ], y = el[2, ])
  panel <- function(M, labs, title) {
    lab <- tibble(x = M[1, ] * (1 + 0.9 / sqrt(colSums(M^2))), y = M[2, ] * (1 + 0.9 / sqrt(colSums(M^2))),
                  txt = labs, col = c(col_a1, col_a2))
    ggplot() + axes0 + cloud_pts(air_X, 0.25) + geom_path(data = el, aes(x, y), linewidth = 0.9) +
      geom_segment(data = tibble(xend = M[1, ], yend = M[2, ], col = c(col_a1, col_a2)),
                   aes(0, 0, xend = xend, yend = yend, colour = I(col)), linewidth = 1.5,
                   arrow = arrow(length = unit(0.22, "cm"), type = "closed")) +
      geom_label(data = lab, aes(x, y, label = txt, colour = I(col)), size = 4, fontface = "bold",
                 label.size = 0, fill = alpha("white", 0.8), lineheight = 0.9) +
      coord_equal(xlim = c(-4.2, 4.2), ylim = c(-4.2, 4.2)) +
      labs(title = title, x = "flight time (min)", y = "fuel (100 kg)")
  }
  panel(A_pca, c("PC1", "PC2"), "PCA: longest first column") +
    panel(air_A, c("headwind", "load"), "The true shocks (from the manual)")
}

# =============================================================================
# 7  Canonical correlation: education and health in the states
# =============================================================================
Zst <- scale(X6)
cc_states <- cancor(Zst[, c("Illiteracy", "HSGrad")], Zst[, c("LifeExp", "Murder")])
fig_cca <- function() {
  U <- Zst[, c("Illiteracy", "HSGrad")] %*% cc_states$xcoef
  V <- Zst[, c("LifeExp", "Murder")] %*% cc_states$ycoef
  U <- U * sqrt(nrow(U) - 1); V <- V * sqrt(nrow(V) - 1)           # canonical variates with variance 1
  s1 <- sign(cc_states$xcoef["HSGrad", 1])                           # first pair: high = more educated
  s2 <- sign(cc_states$xcoef["HSGrad", 2])                           # second pair: positive weights
  d <- tibble(abb = state.abb, u1 = s1 * U[, 1], v1 = s1 * V[, 1], u2 = s2 * U[, 2], v2 = s2 * V[, 2])
  p1 <- ggplot(d, aes(u1, v1)) + geom_point(colour = col_a1, size = 2) +
    geom_text_repel(aes(label = abb), size = 2.8, colour = "grey35", seed = 5, max.overlaps = 15) +
    labs(title = sprintf("First pair: correlation %.2f", cc_states$cor[1]),
         x = expression(italic(u)[1]*": education"), y = expression(italic(v)[1]*": health and safety"))
  p2 <- ggplot(d, aes(u2, v2)) + geom_point(colour = col_a2, size = 2) +
    geom_text_repel(aes(label = abb), size = 2.8, colour = "grey35", seed = 5, max.overlaps = 15) +
    labs(title = sprintf("Second pair: correlation %.2f", cc_states$cor[2]),
         x = expression(italic(u)[2]*": education, second index"),
         y = expression(italic(v)[2]*": health and safety, second index"))
  p1 + p2
}
