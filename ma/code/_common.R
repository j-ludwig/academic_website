# =============================================================================
# _common.R  --  shared helpers (colours, ellipses, matrix roots, animation tools)
# for the figure scripts of the ISQS 6350 course notes. Each *_figures.R script
# starts with source("_common.R"), so keep this file in the same folder.
# =============================================================================

knitr::opts_chunk$set(
  echo = TRUE,
  message = FALSE,
  warning = FALSE,
  fig.width = 8,
  fig.height = 5.5,
  out.width = "95%",
  fig.align = "center"
)

# Course colour conventions: first latent direction blue, second red,
# common factor purple, data points grey.
col_a1 <- "#1f77b4"
col_a2 <- "#d62728"
col_f  <- "#7b4fb5"
col_pt <- adjustcolor("#6f7a8c", alpha.f = 0.45)

set.seed(1234)

library(MASS)
library(tidyverse)
library(mvtnorm)
library(plotly)
library(patchwork)

theme_set(theme_minimal(base_size = 14))

wrap_plotly_title <- function(x, width = 58) {
  paste(strwrap(x, width = width), collapse = "<br>")
}

rotation_matrix <- function(theta) {
  matrix(
    c(
      cos(theta), -sin(theta),
      sin(theta),  cos(theta)
    ),
    nrow = 2,
    byrow = TRUE
  )
}

ellipse_df <- function(Sigma, mu = c(0, 0), level = 0.95, n = 400) {
  theta <- seq(0, 2 * pi, length.out = n)
  unit_circle <- rbind(cos(theta), sin(theta))

  eg <- eigen(Sigma, symmetric = TRUE)
  radius <- sqrt(qchisq(level, df = 2))

  A <- eg$vectors %*% diag(sqrt(eg$values))
  pts <- radius * A %*% unit_circle
  pts <- sweep(pts, 1, mu, FUN = "+")

  tibble(
    x1 = pts[1, ],
    x2 = pts[2, ]
  )
}

# Real fractional matrix powers for matrices with positive real eigenvalues.
matrix_power_real <- function(A, k, tol = 1e-10) {
  eg <- eigen(A)

  if (any(abs(Im(eg$values)) > tol)) {
    stop("A has complex eigenvalues; use a decomposition-specific path instead.")
  }

  values <- Re(eg$values)

  if (any(values <= 0)) {
    stop("This helper requires positive real eigenvalues.")
  }

  V <- Re(eg$vectors)

  out <- V %*%
    diag(values^k, nrow = length(values)) %*%
    solve(V)

  Re(out)
}

# Geometry for an arrow drawn as ONE closed polygon (no gaps, no NA):
#
#        h1
#   s1---b1 \
#   0         v      <- shaft rectangle flows straight into the head triangle
#   s2---b2 /
#        h2
#
# Filled with fill = "toself", it renders as a single solid arrow, so the
# shaft and head can never separate or be drawn twice.
arrow_geometry <- function(v,
                           head_length = 0.22,
                           head_width = 0.12,
                           shaft_width = 0.025) {
  v <- as.numeric(v)
  nrm <- sqrt(sum(v^2))

  if (nrm < 1e-12) {
    return(list(
      arrow = tibble(x = rep(0, 8), y = rep(0, 8)),
      label = tibble(x = 0, y = 0)
    ))
  }

  u <- v / nrm
  perp <- c(-u[2], u[1])

  hl <- min(head_length, 0.24 * nrm)
  hw <- min(head_width, 0.14 * nrm)
  sw <- min(shaft_width, 0.35 * hw)

  base <- v - hl * u

  s1 <- sw * perp            # shaft corners at the origin
  s2 <- -sw * perp
  b1 <- base + sw * perp     # shaft corners at the head
  b2 <- base - sw * perp
  h1 <- base + hw * perp     # head corners
  h2 <- base - hw * perp

  pts <- rbind(s1, b1, h1, v, h2, b2, s2, s1)

  label_pos <- v + 0.10 * u

  list(
    arrow = tibble(
      x = unname(pts[, 1]),
      y = unname(pts[, 2])
    ),
    label = tibble(
      x = label_pos[1],
      y = label_pos[2]
    )
  )
}

# Add one animated arrow as a single filled polygon trace.
# (`width` is kept for compatibility; thickness is set in arrow_geometry().)
add_arrow_trace <- function(p, data, color, width = NULL, ...) {
  add_trace(
    p,
    data = data,
    x = ~x,
    y = ~y,
    frame = ~k,
    type = "scatter",
    mode = "lines",
    fill = "toself",
    fillcolor = color,
    line = list(width = 1, color = color),
    inherit = FALSE,
    showlegend = FALSE,
    hoverinfo = "skip",
    ...
  )
}

# Build the animated arrow data for one column across all k.
arrow_frames <- function(k_grid, vec_fun) {
  purrr::map_dfr(k_grid, function(k) {
    arrow_geometry(vec_fun(k))$arrow |>
      mutate(k = sprintf("%.2f", k))
  })
}

animate_linear_map <- function(
  Z,
  map_fun,
  k_grid,
  title,
  xlab = "x<sub>1</sub>",
  ylab = "x<sub>2</sub>",
  matrix_symbol = "A"
) {
  vector_colors <- c("#1f77b4", "#d62728")

  anim_df <- purrr::map_dfr(k_grid, function(k) {
    M_k <- map_fun(k)
    X_k <- Z %*% t(M_k)

    tibble(
      id = seq_len(nrow(Z)),
      x1 = X_k[, 1],
      x2 = X_k[, 2],
      k = sprintf("%.2f", k)
    )
  })

  arrow_df <- vector("list", 2)
  label_df <- vector("list", 2)

  for (j in 1:2) {
    arrow_df[[j]] <- arrow_frames(k_grid, function(k) map_fun(k)[, j])

    label_df[[j]] <- purrr::map_dfr(k_grid, function(k) {
      M_k <- map_fun(k)
      geom <- arrow_geometry(M_k[, j])

      geom$label |>
        mutate(
          label = paste0(
            matrix_symbol,
            "<sup>", sprintf("%.2f", k), "</sup>",
            "e<sub>", j, "</sub>"
          ),
          k = sprintf("%.2f", k)
        )
    })
  }

  all_x <- c(
    anim_df$x1,
    unlist(lapply(arrow_df, `[[`, "x")),
    unlist(lapply(label_df, `[[`, "x"))
  )

  all_y <- c(
    anim_df$x2,
    unlist(lapply(arrow_df, `[[`, "y")),
    unlist(lapply(label_df, `[[`, "y"))
  )

  max_abs <- 1.08 * max(abs(c(all_x, all_y)), na.rm = TRUE)

  p <- plot_ly() |>
    add_markers(
      data = anim_df,
      x = ~x1,
      y = ~x2,
      frame = ~k,
      ids = ~id,
      marker = list(size = 6, opacity = 0.35),
      inherit = FALSE,
      showlegend = FALSE,
      hoverinfo = "skip"
    )

  for (j in 1:2) {
    p <- p |>
      add_arrow_trace(arrow_df[[j]], color = vector_colors[j], width = 4) |>
      add_trace(
        data = label_df[[j]],
        x = ~x,
        y = ~y,
        frame = ~k,
        type = "scatter",
        mode = "text",
        text = ~label,
        textposition = "middle right",
        textfont = list(
          color = vector_colors[j],
          size = 14
        ),
        inherit = FALSE,
        showlegend = FALSE,
        hoverinfo = "skip"
      )
  }

  p |>
    layout(
      title = list(
        text = wrap_plotly_title(title),
        x = 0.5,
        xanchor = "center"
      ),
      xaxis = list(
        title = list(text = xlab, standoff = 14),
        range = c(-max_abs, max_abs)
      ),
      yaxis = list(
        title = list(text = ylab, standoff = 10),
        range = c(-max_abs, max_abs),
        scaleanchor = "x",
        scaleratio = 1
      ),
      margin = list(
        l = 70,
        r = 35,
        b = 85,
        t = 105
      ),
      showlegend = FALSE
    ) |>
    animation_opts(
      frame = 70,
      transition = 0,
      redraw = TRUE
    ) |>
    animation_slider(
      currentvalue = list(prefix = "k = ")
    )
}

# ---- Helpers for latent-coordinate geometry and the two interactive figures --
# ---- Helpers for the two interactive figures ---------------------------------
col_a1 <- "#1f77b4"
col_a2 <- "#d62728"
col_pt <- adjustcolor("#6f7a8c", alpha.f = 0.45)

sigma_of <- function(sd1, sd2, rho) {
  matrix(
    c(
      sd1^2,           rho * sd1 * sd2,
      rho * sd1 * sd2, sd2^2
    ),
    2, 2,
    byrow = TRUE
  )
}

sym_root <- function(S) {
  eg <- eigen(S, symmetric = TRUE)
  eg$vectors %*% diag(sqrt(pmax(eg$values, 0))) %*% t(eg$vectors)
}

# One fixed white-noise sample whose sample covariance is exactly I.
# The observed cloud is x = Sigma^(1/2) w.
set.seed(31)
app_W_raw <- matrix(rnorm(3 * 200), ncol = 3)
app_W_raw <- scale(app_W_raw, center = TRUE, scale = FALSE)
app_W3 <- app_W_raw %*% solve(chol(cov(app_W_raw)))   # sample covariance exactly I (3 x 3)
app_W <- app_W3[, 1:2]                            # drives the observed cloud
app_V <- app_W3[, 3]                              # position along the unidentified latent direction


# Rotation phi (degrees) that turns a given factorization A (AA' = Sigma)
# into root %*% Q(phi). Column signs are flipped if needed so det(Q) = 1.
phi_from_A <- function(A, root) {
  Q <- solve(root) %*% A
  if (det(Q) < 0) Q[, 2] <- -Q[, 2]
  unname(atan2(Q[2, 1], Q[1, 1]) * 180 / pi)
}

preset_phis <- function(S, root) {
  A_chol <- t(chol(S))                                    # Cholesky, x1 first
  P <- matrix(c(0, 1, 1, 0), 2)
  A_chol2 <- P %*% t(chol(P %*% S %*% P))                 # Cholesky, x2 first
  eg <- eigen(S, symmetric = TRUE)
  A_pca <- eg$vectors %*% diag(sqrt(eg$values))           # PCA
  u <- c(1, 1) / sqrt(2)                                  # structural, eta_s = 1
  a1 <- u / sqrt(as.numeric(t(u) %*% solve(S) %*% u))
  R <- S - tcrossprod(a1)
  a2 <- if (R[1, 1] > 1e-12) R[, 1] / sqrt(R[1, 1]) else c(0, sqrt(max(R[2, 2], 0)))
  if (a2[1] > 0) a2 <- -a2

  c(
    "Σ½ (symmetric root)" = phi_from_A(sym_root(S), root),
    "Cholesky"            = phi_from_A(A_chol, root),
    "Cholesky (X₂ first)" = phi_from_A(A_chol2, root),
    "PCA"                 = phi_from_A(A_pca, root),
    "Structural (ηₛ = 1)" = phi_from_A(cbind(a1, a2), root)
  )
}

ellipse_xy <- function(M, r = 1, n = 160) {
  t <- seq(0, 2 * pi, length.out = n)
  r * M %*% rbind(cos(t), sin(t))
}

angle_between <- function(a, b) {
  acos(pmin(1, pmax(-1, sum(a * b) / sqrt(sum(a^2) * sum(b^2))))) * 180 / pi
}

wrap180 <- function(t) ((t + 180) %% 360) - 180

# Diverging colour scale for the factor value f (orange = low, purple = high).
f_pal <- colorRampPalette(c("#b35806", "#f1a340", "#d4d4d4", "#998ec3", "#542788"))(101)
f_lim <- 2.5
f_col <- function(f) f_pal[round((pmin(pmax(f, -f_lim), f_lim) + f_lim) / (2 * f_lim) * 100) + 1]
f_legend <- function(x0, y0, w, h) {
  k <- length(f_pal)
  xs <- seq(x0, x0 + w, length.out = k + 1)
  rect(xs[-(k + 1)], y0, xs[-1], y0 + h, col = f_pal, border = NA)
  rect(x0, y0, x0 + w, y0 + h, border = "grey50")
  text(c(x0, x0 + w / 2, x0 + w), y0 - 0.06 * w, c(paste0("−", f_lim), "0", f_lim), cex = 0.8, adj = c(0.5, 1))
  text(x0 + w / 2, y0 + h + 0.05 * w, "factor value f", cex = 0.85, adj = c(0.5, 0))
}

# One-factor model for a given reduced form S and communality h1 of x1.
# Normalization: lambda1 > 0 (the sign of f is not identified).
fa_from_h1 <- function(S, h1) {
  rho2 <- S[1, 2]^2 / (S[1, 1] * S[2, 2])
  h1 <- min(max(h1, rho2, 1e-9), 1)
  h2 <- rho2 / h1
  lambda <- c(sqrt(h1 * S[1, 1]), sign(S[1, 2]) * sqrt(h2 * S[2, 2]))
  psi <- pmax(c(S[1, 1] * (1 - h1), S[2, 2] * (1 - h2)), 0)
  A <- cbind(lambda, c(sqrt(psi[1]), 0), c(0, sqrt(psi[2])))
  # Null direction of A: moving z along n leaves x = Az unchanged.
  r1 <- A[1, ]; r2 <- A[2, ]
  n <- c(r1[2] * r2[3] - r1[3] * r2[2], r1[3] * r2[1] - r1[1] * r2[3], r1[1] * r2[2] - r1[2] * r2[1])
  n <- n / sqrt(sum(n^2))
  if (n[1] < 0) n <- -n
  list(h = c(h1, h2), lambda = lambda, psi = psi, A = unname(A), n = n,
       split = (1 - h1) / ((1 - h1) + (1 - h2)))
}

# Structural slider: standardized noise split s = (1 - h1) / ((1 - h1) + (1 - h2)),
# i.e. psi1/sigma1^2 relative to psi2/sigma2^2. Scale-free, so the standard deviations do not affect it.
# s is monotone (decreasing) in h1, so invert by bisection.
fa_from_split <- function(S, s) {
  rho2 <- S[1, 2]^2 / (S[1, 1] * S[2, 2])
  lo <- max(rho2, 1e-9); hi <- 1
  for (k in 1:60) {
    mid <- (lo + hi) / 2
    if (fa_from_h1(S, mid)$split > s) lo <- mid else hi <- mid
  }
  fa_from_h1(S, (lo + hi) / 2)
}

# Extra restrictions, each expressed as the noise split it selects (NA if infeasible).
fa_restrictions <- function(S) {
  s11 <- S[1, 1]; s22 <- S[2, 2]; s12 <- S[1, 2]
  rho2 <- s12^2 / (s11 * s22)
  split_at <- function(h1) {
    if (!is.finite(h1) || h1 < rho2 - 1e-9 || h1 > 1 + 1e-9) return(NA_real_)
    fa_from_h1(S, h1)$split
  }
  u_psi <- ((s11 - s22) + sqrt((s11 - s22)^2 + 4 * s12^2)) / 2   # lambda1^2 when psi1 = psi2
  c(
    "ψ₁ = 0 (X₁ is the factor)"   = 0,
    "ψ₂ = 0 (X₂ is the factor)"   = 1,
    "equal communalities h₁ = h₂"  = 0.5,
    "equal noise ψ₁ = ψ₂"          = split_at(u_psi / s11),
    "equal loadings |λ₁| = |λ₂|"   = split_at(abs(s12) / s11)
  )
}
# ---- The data, held fixed throughout the chapter -----------------------------
app_Sigma <- matrix(c(1, 0.75, 0.75, 1.25), 2, 2)
app_root <- sym_root(app_Sigma)              # symmetric square root; phi = 0 in the rotation app
app_X <- app_W %*% t(app_root)                # the observed cloud (sample covariance exactly Sigma)
app_presets <- preset_phis(app_Sigma, app_root)
app_restr <- fa_restrictions(app_Sigma)

# ---- MDS helpers ----------------------------------------------------------------
# Eigenvalues of the double-centered squared-distance matrix B = -1/2 J D^2 J.
mds_eigen <- function(D) {
  D2 <- as.matrix(D)^2
  n <- nrow(D2)
  J <- diag(n) - matrix(1 / n, n, n)
  B <- -0.5 * J %*% D2 %*% J
  eigen((B + t(B)) / 2, symmetric = TRUE)$values
}

# Rotate, reflect and translate Y so that it lines up with X as closely as possible.
procrustes_align <- function(Y, X) {
  Xc <- scale(X, scale = FALSE)
  Yc <- scale(Y, scale = FALSE)
  s <- svd(t(Yc) %*% Xc)
  R <- s$u %*% t(s$v)
  sweep(Yc %*% R, 2, colMeans(X), "+")
}

# ---- Interactive figure builders (plotly) ---------------------------------------
# plotly sorts animation frames as text ("-10" before "-100"). After building, put frames and
# slider steps back in numeric order, and add a row of buttons that jump to named frames.
finalize_animation <- function(p, jump_buttons, start_frame, prefix, suffix = "") {
  b <- plotly_build(p)
  ord <- order(as.numeric(vapply(b$x$frames, `[[`, "", "name")))
  b$x$frames <- b$x$frames[ord]
  steps <- b$x$layout$sliders[[1]]$steps
  step_ord <- order(as.numeric(vapply(steps, `[[`, "", "label")))
  b$x$layout$sliders[[1]]$steps <- steps[step_ord]
  b$x$layout$sliders[[1]]$active <- which(step_ord == 1) - 1
  b$x$layout$updatemenus <- c(b$x$layout$updatemenus, list(list(
    type = "buttons", direction = "right", showactive = FALSE,
    x = 0, xanchor = "left", y = 1.09, yanchor = "bottom",
    pad = list(t = 0, r = 4), buttons = jump_buttons
  )))
  b$x$layout$sliders[[1]]$currentvalue <- list(prefix = prefix, suffix = suffix, font = list(size = 14))
  b$x$layout$sliders[[1]]$pad <- list(t = 60)
  # Open the figure at a chosen frame (plotly otherwise starts at the first frame it built).
  b$x$layout$sliders[[1]]$active <- which(vapply(b$x$layout$sliders[[1]]$steps, `[[`, "", "label") == start_frame) - 1
  htmlwidgets::onRender(b, sprintf(
    "function(el) { Plotly.animate(el, ['%s'], {mode: 'immediate', frame: {duration: 0, redraw: true}, transition: {duration: 0}}); }",
    start_frame))
}

# ---- Interactive figure 1: rotating the columns of A --------------------------
make_rotation_figure <- function() {
  phi_lab <- function(phi) sprintf("%.1f", phi)

  phi_grid <- sort(unique(round(c(seq(-180, 175, by = 5), unname(app_presets)), 1)))
  levs <- phi_lab(phi_grid)
  fk <- function(phi) factor(phi_lab(phi), levels = levs)
  A_of <- function(phi) app_root %*% rotation_matrix(phi * pi / 180)

  # ---------- per-frame data ----------
  arrows <- lapply(1:2, function(j) purrr::map_dfr(phi_grid, function(ph)
    arrow_geometry(A_of(ph)[, j], head_length = 0.2, head_width = 0.1, shaft_width = 0.022)$arrow |>
      mutate(k = fk(ph))))

  arrow_labels <- purrr::map_dfr(phi_grid, function(ph) {
    M <- A_of(ph)
    tibble(x = M[1, ] * (1 + 0.32 / sqrt(colSums(M^2))),
           y = M[2, ] * (1 + 0.32 / sqrt(colSums(M^2))),
           label = c("a<sub>1</sub>", "a<sub>2</sub>"),
           col = c("#1f77b4", "#d62728"), k = fk(ph))
  })

  para <- purrr::map_dfr(phi_grid, function(ph) {
    M <- A_of(ph)
    cn <- M %*% cbind(c(1, 1), c(-1, 1), c(-1, -1), c(1, -1), c(1, 1))
    tibble(x = cn[1, ], y = cn[2, ], k = fk(ph))
  })

  col_lines <- lapply(1:2, function(j) purrr::map_dfr(phi_grid, function(ph) {
    v <- A_of(ph)[, j]
    tibble(x = c(-20, 20) * v[1], y = c(-20, 20) * v[2], k = fk(ph))
  }))

  latent <- purrr::map_dfr(phi_grid, function(ph) {
    z <- app_X %*% t(solve(A_of(ph)))
    tibble(id = seq_len(nrow(z)), z1 = z[, 1], z2 = z[, 2], k = fk(ph))
  })

  preset_names_short <- names(app_presets)
  readout <- purrr::map_dfr(phi_grid, function(ph) {
    M <- A_of(ph)
    ang <- acos(sum(M[, 1] * M[, 2]) / sqrt(sum(M[, 1]^2) * sum(M[, 2]^2))) * 180 / pi
    near <- which(abs(wrap180(ph - app_presets)) < 0.06 | abs(wrap180(ph - app_presets - 180)) < 0.06)
    tibble(
      x = 4.4, y = -4.35,
      text = paste0(
        "A = [", sprintf("%6.3f", M[1, 1]), "  ", sprintf("%6.3f", M[1, 2]), "]<br>",
        "    [", sprintf("%6.3f", M[2, 1]), "  ", sprintf("%6.3f", M[2, 2]), "]<br>",
        "angle(a₁, a₂) = ", sprintf("%.1f°", ang), "<br>",
        if (length(near)) paste0("<b>", preset_names_short[near[1]], "</b>") else " "
      ),
      k = fk(ph)
    )
  })

  ell1 <- ellipse_xy(app_root, 1); ell2 <- ellipse_xy(app_root, 2)
  circ <- ellipse_xy(diag(2), 1)
  e_arrows <- lapply(list(c(1, 0), c(0, 1)), function(v)
    arrow_geometry(v, head_length = 0.2, head_width = 0.1, shaft_width = 0.022)$arrow)

  # ---------- observed panel ----------
  p_obs <- plot_ly() |>
    add_trace(data = para, x = ~x, y = ~y, frame = ~k, type = "scatter", mode = "lines",
              fill = "toself", fillcolor = "rgba(44,62,80,0.06)",
              line = list(color = "rgba(44,62,80,0.35)", width = 1),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)
  for (j in 1:2) {
    p_obs <- p_obs |>
      add_trace(data = col_lines[[j]], x = ~x, y = ~y, frame = ~k, type = "scatter", mode = "lines",
                line = list(color = c("rgba(31,119,180,0.4)", "rgba(214,39,40,0.4)")[j], dash = "dash", width = 1),
                hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)
  }
  p_obs <- p_obs |>
    add_trace(data = tibble(x1 = app_X[, 1], x2 = app_X[, 2]), x = ~x1, y = ~x2,
              type = "scatter", mode = "markers",
              marker = list(size = 5, color = "#6f7a8c", opacity = 0.45),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(x = ell2[1, ], y = ell2[2, ], type = "scatter", mode = "lines",
              line = list(color = "grey", dash = "dot", width = 1),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(x = ell1[1, ], y = ell1[2, ], type = "scatter", mode = "lines",
              line = list(color = "black", width = 2),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)
  for (j in 1:2) p_obs <- add_arrow_trace(p_obs, arrows[[j]], color = c("#1f77b4", "#d62728")[j])
  p_obs <- p_obs |>
    add_trace(data = arrow_labels, x = ~x, y = ~y, frame = ~k, text = ~label,
              type = "scatter", mode = "text", textfont = list(size = 15, color = ~col),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(data = readout, x = ~x, y = ~y, frame = ~k, text = ~text,
              type = "scatter", mode = "text", textposition = "top left",
              textfont = list(size = 12, family = "Courier New, monospace", color = "#2c3e50"),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)

  # ---------- latent panel ----------
  p_lat <- plot_ly() |>
    add_trace(data = latent, x = ~z1, y = ~z2, frame = ~k, ids = ~id,
              type = "scatter", mode = "markers",
              marker = list(size = 5, color = "#6f7a8c", opacity = 0.45),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(x = circ[1, ], y = circ[2, ], type = "scatter", mode = "lines",
              line = list(color = "black", width = 2),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)
  for (j in 1:2) {
    p_lat <- p_lat |>
      add_trace(data = e_arrows[[j]], x = ~x, y = ~y, type = "scatter", mode = "lines",
                fill = "toself", fillcolor = c("#1f77b4", "#d62728")[j],
                line = list(width = 1, color = c("#1f77b4", "#d62728")[j]),
                hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)
  }
  p_lat <- p_lat |>
    add_trace(x = c(1.3, -0.25), y = c(-0.25, 1.3), text = c("e<sub>1</sub>", "e<sub>2</sub>"),
              type = "scatter", mode = "text",
              textfont = list(size = 15, color = c("#1f77b4", "#d62728")),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)

  preset_buttons <- lapply(seq_along(app_presets), function(i) list(
    label = preset_names_short[i], method = "animate",
    args = list(list(phi_lab(round(app_presets[[i]], 1))),
                list(mode = "immediate", frame = list(duration = 0, redraw = TRUE),
                     transition = list(duration = 0)))
  ))

  subplot(p_obs, p_lat, widths = c(0.55, 0.45), margin = 0.06, titleX = TRUE, titleY = TRUE) |>
    layout(
      xaxis  = list(title = "x<sub>1</sub>", range = c(-4.5, 4.5), zeroline = FALSE),
      yaxis  = list(title = "x<sub>2</sub>", range = c(-4.5, 4.5), scaleanchor = "x", zeroline = FALSE),
      xaxis2 = list(title = "ε<sub>1</sub>", range = c(-3.6, 3.6), zeroline = FALSE),
      yaxis2 = list(title = "ε<sub>2</sub>", range = c(-3.6, 3.6), scaleanchor = "x2", zeroline = FALSE),
      annotations = list(
        list(text = "<b>Observed space: X = Aε</b> (the data never move)", x = 0, y = 1.03,
             xref = "paper", yref = "paper", xanchor = "left", showarrow = FALSE, font = list(size = 13)),
        list(text = "<b>Latent space: ε = A<sup>−1</sup>X</b>", x = 0.56, y = 1.03,
             xref = "paper", yref = "paper", xanchor = "left", showarrow = FALSE, font = list(size = 13))
      ),
      margin = list(l = 50, r = 20, t = 95, b = 170),
      showlegend = FALSE
    ) |>
    animation_opts(frame = 80, transition = 0, redraw = TRUE) |>
    animation_slider(currentvalue = list(prefix = "rotation φ = ", suffix = "°")) |>
    animation_button(x = 0, xanchor = "left", y = -0.3, yanchor = "top", label = "Play") |>
    finalize_animation(preset_buttons, start_frame = "0.0", prefix = "rotation φ = ", suffix = "°")
}

# ---- Interactive figure 2: the one-factor model ---------------------------------
make_factor_figure <- function(notation = c("psi", "sigma")) {
  notation <- match.arg(notation)
  # "psi": x_i = lambda_i f + sqrt(psi_i) u_i (book notation, Chapter 5)
  # "sigma": x_i = lambda_i f + sigma_i epsilon_i (Introduction notation)
  sym <- if (notation == "psi") c(sh = "u", sd = "ψ") else c(sh = "ε", sd = "σ")
  relabel <- function(x) if (notation == "psi") x else gsub("ψ", "σ", x)
  s_lab <- function(s) sprintf("%.3f", s)
  restr_ok <- app_restr[!is.na(app_restr)]
  s_grid <- sort(unique(round(c(seq(0, 1, by = 0.025), unname(restr_ok)), 3)))
  levs <- s_lab(s_grid)
  fk <- function(s) factor(s_lab(s), levels = levs)

  models <- lapply(s_grid, function(s) fa_from_split(app_Sigma, s))
  Zs <- lapply(models, function(m) app_X %*% solve(app_Sigma) %*% m$A + outer(app_V, m$n))

  colscale <- list(c(0, "#b35806"), c(0.25, "#f1a340"), c(0.5, "#d4d4d4"),
                   c(0.75, "#998ec3"), c(1, "#542788"))

  obs_pts <- purrr::map_dfr(seq_along(s_grid), function(i)
    tibble(id = seq_len(nrow(app_X)), x1 = app_X[, 1], x2 = app_X[, 2],
           f = Zs[[i]][, 1], k = fk(s_grid[i])))
  lat_pts <- purrr::map_dfr(seq_along(s_grid), function(i)
    tibble(id = seq_len(nrow(app_X)), u1 = Zs[[i]][, 2], u2 = Zs[[i]][, 3],
           f = Zs[[i]][, 1], k = fk(s_grid[i])))

  arrow_set <- function(vec_fun) purrr::map_dfr(seq_along(s_grid), function(i)
    arrow_geometry(vec_fun(models[[i]]), head_length = 0.2, head_width = 0.1, shaft_width = 0.022)$arrow |>
      mutate(k = fk(s_grid[i])))
  arr_lam  <- arrow_set(function(m) m$lambda)
  arr_psi1 <- arrow_set(function(m) c(sqrt(m$psi[1]), 0))
  arr_psi2 <- arrow_set(function(m) c(0, sqrt(m$psi[2])))

  psi_ell <- purrr::map_dfr(seq_along(s_grid), function(i) {
    e <- ellipse_xy(diag(sqrt(models[[i]]$psi) + 1e-9), 1)
    tibble(x = e[1, ], y = e[2, ], k = fk(s_grid[i]))
  })
  fac_line <- purrr::map_dfr(seq_along(s_grid), function(i) {
    l <- models[[i]]$lambda
    tibble(x = c(-20, 20) * l[1], y = c(-20, 20) * l[2], k = fk(s_grid[i]))
  })
  lam_label <- purrr::map_dfr(seq_along(s_grid), function(i) {
    l <- models[[i]]$lambda
    tibble(x = l[1] * (1 + 0.3 / sqrt(sum(l^2))), y = l[2] * (1 + 0.3 / sqrt(sum(l^2))), k = fk(s_grid[i]))
  })

  restr_names <- names(app_restr)
  readout <- purrr::map_dfr(seq_along(s_grid), function(i) {
    m <- models[[i]]
    near <- which(!is.na(app_restr) & abs(app_restr - s_grid[i]) < 0.0006)
    tibble(
      x = 4.4, y = -4.35,
      text = paste0(
        "λ = (", sprintf("%.3f", m$lambda[1]), ", ", sprintf("%.3f", m$lambda[2]), ")<br>",
        if (notation == "psi")
          paste0("ψ = (", sprintf("%.3f", m$psi[1]), ", ", sprintf("%.3f", m$psi[2]), ")<br>")
        else
          paste0("σ = (", sprintf("%.3f", sqrt(m$psi[1])), ", ", sprintf("%.3f", sqrt(m$psi[2])), ")<br>"),
        "h = (", sprintf("%.3f", m$h[1]), ", ", sprintf("%.3f", m$h[2]), "), h₁h₂ = ", sprintf("%.3f", prod(m$h)), "<br>",
        "Var(f | x) = ", sprintf("%.3f", m$n[1]^2), "<br>",
        if (length(near)) paste0("<b>", relabel(restr_names[near[1]]), "</b>") else " "
      ),
      k = fk(s_grid[i])
    )
  })

  ell1 <- ellipse_xy(app_root, 1); ell2 <- ellipse_xy(app_root, 2)
  circ <- ellipse_xy(diag(2), 1)
  u_arrows <- lapply(list(c(1, 0), c(0, 1)), function(v)
    arrow_geometry(v, head_length = 0.2, head_width = 0.1, shaft_width = 0.022)$arrow)

  p_obs <- plot_ly() |>
    add_trace(data = fac_line, x = ~x, y = ~y, frame = ~k, type = "scatter", mode = "lines",
              line = list(color = "rgba(123,79,181,0.45)", dash = "dash", width = 1),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(data = obs_pts, x = ~x1, y = ~x2, frame = ~k, ids = ~id, type = "scatter", mode = "markers",
              marker = list(size = 6, color = ~f, colorscale = colscale, cmin = -2.5, cmax = 2.5, opacity = 0.85),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(x = ell2[1, ], y = ell2[2, ], type = "scatter", mode = "lines",
              line = list(color = "grey", dash = "dot", width = 1),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(x = ell1[1, ], y = ell1[2, ], type = "scatter", mode = "lines",
              line = list(color = "black", width = 2),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(data = psi_ell, x = ~x, y = ~y, frame = ~k, type = "scatter", mode = "lines",
              line = list(color = "#555555", dash = "dash", width = 1.5),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)
  p_obs <- add_arrow_trace(p_obs, arr_psi1, color = "#1f77b4")
  p_obs <- add_arrow_trace(p_obs, arr_psi2, color = "#d62728")
  p_obs <- add_arrow_trace(p_obs, arr_lam, color = "#7b4fb5")
  p_obs <- p_obs |>
    add_trace(data = lam_label, x = ~x, y = ~y, frame = ~k, text = "λ",
              type = "scatter", mode = "text", textfont = list(size = 16, color = "#7b4fb5"),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(data = readout, x = ~x, y = ~y, frame = ~k, text = ~text,
              type = "scatter", mode = "text", textposition = "top left",
              textfont = list(size = 12, family = "Courier New, monospace", color = "#2c3e50"),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)

  p_lat <- plot_ly() |>
    add_trace(data = lat_pts, x = ~u1, y = ~u2, frame = ~k, ids = ~id, type = "scatter", mode = "markers",
              marker = list(size = 6, color = ~f, colorscale = colscale, cmin = -2.5, cmax = 2.5, opacity = 0.85,
                            showscale = TRUE,
                            colorbar = list(title = list(text = "factor f"), len = 0.6, thickness = 12)),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(x = circ[1, ], y = circ[2, ], type = "scatter", mode = "lines",
              line = list(color = "black", width = 2),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)
  for (j in 1:2) {
    p_lat <- p_lat |>
      add_trace(data = u_arrows[[j]], x = ~x, y = ~y, type = "scatter", mode = "lines",
                fill = "toself", fillcolor = c("#1f77b4", "#d62728")[j],
                line = list(width = 1, color = c("#1f77b4", "#d62728")[j]),
                hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)
  }
  p_lat <- p_lat |>
    add_trace(x = c(1.3, -0.25), y = c(-0.25, 1.3), text = paste0(sym[["sh"]], c("<sub>1</sub>", "<sub>2</sub>")),
              type = "scatter", mode = "text",
              textfont = list(size = 15, color = c("#1f77b4", "#d62728")),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)

  short_label <- c("ψ₁ = 0 (X₁ is the factor)" = "ψ₁ = 0", "ψ₂ = 0 (X₂ is the factor)" = "ψ₂ = 0",
                   "equal communalities h₁ = h₂" = "h₁ = h₂", "equal noise ψ₁ = ψ₂" = "ψ₁ = ψ₂",
                   "equal loadings |λ₁| = |λ₂|" = "|λ₁| = |λ₂|")
  restr_buttons <- lapply(names(restr_ok), function(nm) list(
    label = relabel(short_label[[nm]]), method = "animate",
    args = list(list(s_lab(round(restr_ok[[nm]], 3))),
                list(mode = "immediate", frame = list(duration = 0, redraw = TRUE),
                     transition = list(duration = 0)))
  ))

  subplot(p_obs, p_lat, widths = c(0.55, 0.45), margin = 0.06, titleX = TRUE, titleY = TRUE) |>
    layout(
      xaxis  = list(title = "x<sub>1</sub>", range = c(-4.5, 4.5), zeroline = FALSE),
      yaxis  = list(title = "x<sub>2</sub>", range = c(-4.5, 4.5), scaleanchor = "x", zeroline = FALSE),
      xaxis2 = list(title = paste0(sym[["sh"]], "<sub>1</sub>"), range = c(-3.6, 3.6), zeroline = FALSE),
      yaxis2 = list(title = paste0(sym[["sh"]], "<sub>2</sub>"), range = c(-3.6, 3.6), scaleanchor = "x2", zeroline = FALSE),
      annotations = list(
        list(text = "<b>Observed space</b> (colour = factor value f)", x = 0, y = 1.03,
             xref = "paper", yref = "paper", xanchor = "left", showarrow = FALSE, font = list(size = 13)),
        list(text = paste0("<b>Latent space: (", sym[["sh"]], "<sub>1</sub>, ", sym[["sh"]], "<sub>2</sub>), f as colour</b>"), x = 0.56, y = 1.03,
             xref = "paper", yref = "paper", xanchor = "left", showarrow = FALSE, font = list(size = 13))
      ),
      margin = list(l = 50, r = 20, t = 95, b = 170),
      showlegend = FALSE
    ) |>
    animation_opts(frame = 120, transition = 0, redraw = TRUE) |>
    animation_slider(currentvalue = list(prefix = "relative noise in X₁ vs X₂ = ")) |>
    animation_button(x = 0, xanchor = "left", y = -0.3, yanchor = "top", label = "Play") |>
    finalize_animation(restr_buttons, start_frame = "0.500", prefix = "relative noise in X₁ vs X₂ = ")
}
