# =============================================================================
# distributions_figures.R  --  data and figures for the chapter "From Data to Distributions" of the ISQS 6350 course notes.
# Put this file and _common.R (shared helpers) in one folder, set it as the working
# directory, and run source("distributions_figures.R"). Every figure is a function, e.g. fig_...().
# Packages: tidyverse, plotly, patchwork, mvtnorm, MASS, ggrepel, palmerpenguins
# =============================================================================

source("_common.R", local = TRUE)
suppressPackageStartupMessages({
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
X6 <- as.matrix(states[, vars6]); rownames(X6) <- states$state
peng <- na.omit(penguins)
pal_species <- c(Adelie = "#1b9e77", Chinstrap = "#d95f02", Gentoo = "#7570b3")
# The flights of the first lecture (covariance of flight time and fuel)
air_A     <- matrix(c(1.2, 1.2, 0.2, 1.3), 2, 2)
air_Sigma <- air_A %*% t(air_A)

# colours for the layers (as in the Introduction lecture)
col_e <- "#16a085"   # inferred reduced form
col_h <- "#c0398b"   # estimated structural form

# ---- Numbers quoted in the text --------------------------------------------------
d2_st <- mahalanobis(X6, colMeans(X6), cov(X6))
thr6  <- qchisq(0.975, df = 6)

# =============================================================================
# Random experiments: one die, two dice
# =============================================================================

# One die: the pmf (heights) and the cdf (a staircase)
fig_die <- function() {
  pm <- tibble(x = 1:6, p = 1 / 6)
  p1 <- ggplot(pm, aes(x, p)) +
    geom_segment(aes(xend = x, yend = 0), colour = col_a1, linewidth = 1.2) +
    geom_point(colour = col_a1, size = 3.5) +
    scale_x_continuous(breaks = 1:6, limits = c(0.5, 6.5)) +
    scale_y_continuous(limits = c(0, 0.25), breaks = c(0, 1 / 6), labels = c("0", "1/6")) +
    labs(title = "Probability mass function", subtitle = expression(italic(p)(italic(x)) == italic(P)(italic(X) == italic(x))),
         x = expression("number shown"~italic(x)), y = expression(italic(p)(italic(x))))
  cd <- tibble(x0 = c(0, 1:6), x1 = c(1:6, 7), F = (0:6) / 6)
  p2 <- ggplot(cd) +
    geom_segment(aes(x = x0, xend = x1, y = F, yend = F), colour = col_a2, linewidth = 1.2) +
    geom_point(data = cd[-1, ], aes(x0, F), colour = col_a2, size = 2.8) +
    geom_point(data = cd[-7, ], aes(x1, F), colour = col_a2, size = 2.8, shape = 21, fill = "white") +
    annotate("segment", x = 3.5, xend = 3.5, y = 0, yend = 0.5, linetype = 3, colour = "grey40") +
    annotate("text", x = 0.1, y = 0.9, label = "italic(F)(3.5) == italic(P)(italic(X) <= 3.5)*' = 3/6'", parse = TRUE, hjust = 0, size = 3.8) +
    scale_x_continuous(breaks = 0:7, limits = c(0, 7)) +
    scale_y_continuous(breaks = (0:6) / 6, labels = c("0", paste0(1:5, "/6"), "1")) +
    labs(title = "Cumulative distribution function", subtitle = expression(italic(F)(italic(x)) == italic(P)(italic(X) <= italic(x))),
         x = expression(italic(x)), y = expression(italic(F)(italic(x))))
  p1 + p2
}

# Two dice: X1 = first die, X2 = sum of both dice. Joint pmf (36 equally likely outcomes) with marginals.
fig_dice_joint <- function(cond = 3) {
  jt <- expand_grid(d1 = 1:6, d2 = 1:6) |> mutate(x1 = d1, x2 = d1 + d2)
  m1 <- tibble(x1 = 1:6, p = 1 / 6)
  m2 <- jt |> count(x2) |> mutate(p = n / 36)
  p_main <- ggplot(jt, aes(x1, x2)) +
    annotate("rect", xmin = cond - 0.4, xmax = cond + 0.4, ymin = cond + 0.5, ymax = cond + 6.5,
             fill = col_a2, alpha = 0.12) +
    geom_point(size = 7, colour = col_a1, alpha = 0.85) +
    geom_text(label = "1/36", colour = "white", size = 2.3, fontface = "bold") +
    annotate("text", x = cond, y = cond + 7.1, label = sprintf("italic(X)[1] == %d", cond), parse = TRUE, colour = col_a2, size = 3.8) +
    scale_x_continuous(breaks = 1:6, limits = c(0.5, 6.5)) +
    scale_y_continuous(breaks = 2:12, limits = c(1.5, 12.5)) +
    labs(x = expression(italic(X)[1]*": first die"), y = expression(italic(X)[2]*": sum of both dice"))
  p_top <- ggplot(m1, aes(x1, p)) + geom_col(width = 0.5, fill = "grey60") +
    scale_x_continuous(limits = c(0.5, 6.5), breaks = NULL) +
    scale_y_continuous(breaks = c(0, 1 / 6), labels = c("0", "1/6")) + labs(x = NULL, y = expression(italic(P)(italic(X)[1])))
  p_right <- ggplot(m2, aes(x2, p)) + geom_col(width = 0.5, fill = "grey60") + coord_flip() +
    scale_x_continuous(limits = c(1.5, 12.5), breaks = NULL) +
    scale_y_continuous(breaks = c(0, 1 / 6), labels = c("0", "1/6")) + labs(x = NULL, y = expression(italic(P)(italic(X)[2])))
  p_top + plot_spacer() + p_main + p_right + plot_layout(ncol = 2, widths = c(4, 1.3), heights = c(1, 4))
}

# Repeating the experiment n times: relative frequencies (bubbles) approach the joint pmf (open circles)
fig_dice_sample <- function(ns = c(36, 360, 3600)) {
  set.seed(6350)
  rel <- purrr::map_dfr(ns, function(n) {
    d1 <- sample(1:6, n, TRUE); d2 <- sample(1:6, n, TRUE)
    tibble(x1 = d1, x2 = d1 + d2) |> count(x1, x2) |> mutate(f = n / !!n, n_lab = paste("n =", format(!!n, big.mark = ",", trim = TRUE)))
  }) |> mutate(n_lab = factor(n_lab, levels = paste("n =", format(ns, big.mark = ",", trim = TRUE))))
  truth <- expand_grid(x1 = 1:6, d2 = 1:6) |> mutate(x2 = x1 + d2, f = 1 / 36)
  ggplot(rel, aes(x1, x2)) +
    geom_point(data = truth, aes(size = f), shape = 21, colour = "grey40", fill = NA) +
    geom_point(aes(size = f), colour = col_a1, alpha = 0.7) +
    scale_size_area(max_size = 9, guide = "none") +
    facet_wrap(~ n_lab) +
    scale_x_continuous(breaks = 1:6) + scale_y_continuous(breaks = seq(2, 12, 2)) +
    labs(x = expression(italic(x)[1]*": first die"), y = expression(italic(x)[2]*": sum"),
         subtitle = expression("Blue: share of the"~italic(n)~"rows at each ("*italic(x)[1]*","~italic(x)[2]*"). Circles: the probability 1/36."))
}

# A continuous variable: pdf (area) and cdf for daily sales of a coffee shop
fig_cont_pdf_cdf <- function(x0 = 1500) {
  sc <- 1200 / 4
  g <- tibble(x = seq(0, 4000, length.out = 400)) |> mutate(f = dgamma(x, 4, scale = sc), F = pgamma(x, 4, scale = sc))
  F0 <- pgamma(x0, 4, scale = sc)
  p1 <- ggplot(g, aes(x, f)) +
    geom_area(data = dplyr::filter(g, x <= x0), fill = col_a1, alpha = 0.35) +
    geom_line(linewidth = 1) +
    annotate("text", x = 650, y = 0.00018, label = sprintf("area = %.2f", F0), colour = col_a1, size = 4.2) +
    scale_x_continuous(labels = scales::dollar) +
    scale_y_continuous(labels = NULL) +
    labs(title = "Probability density function", subtitle = expression("probability = area under"~italic(f)),
         x = expression("daily sales"~italic(x)), y = expression(italic(f)(italic(x))))
  p2 <- ggplot(g, aes(x, F)) +
    geom_line(colour = col_a2, linewidth = 1.1) +
    annotate("segment", x = x0, xend = x0, y = 0, yend = F0, linetype = 3) +
    annotate("segment", x = 0, xend = x0, y = F0, yend = F0, linetype = 3) +
    annotate("point", x = x0, y = F0, colour = col_a2, size = 3) +
    annotate("text", x = x0 + 120, y = F0 - 0.06, label = sprintf("italic(F)(1500) == %.2f", F0), parse = TRUE, hjust = 0, size = 4.2) +
    scale_x_continuous(labels = scales::dollar) +
    labs(title = "Cumulative distribution function",
         subtitle = expression(italic(F)(italic(x)) == italic(P)(italic(X) <= italic(x))*" = area to the left of "*italic(x)),
         x = expression("daily sales"~italic(x)), y = expression(italic(F)(italic(x))))
  p1 + p2
}

# =============================================================================
# Where the data come from: one data-generating process
# =============================================================================

set.seed(301)
coffee_n  <- 3000
coffee_k  <- 4                               # gamma shape: right-skewed sales
mu_A      <- 1200
sales_A   <- rgamma(coffee_n, shape = coffee_k, scale = mu_A / coffee_k)
m_B       <- 1200 + 0.4 * seq_len(coffee_n)
sales_B   <- rgamma(coffee_n, shape = coffee_k, scale = m_B / coffee_k)

jump_buttons <- function(labels, frames) lapply(seq_along(labels), function(i) list(
  label = labels[i], method = "animate",
  args = list(list(frames[i]), list(mode = "immediate", frame = list(duration = 0, redraw = TRUE),
                                    transition = list(duration = 0)))))
fig_coffee <- function(height = NULL) {
  n_frames <- c(10, 20, 50, 100, 200, 500, 1000, 2000, 3000)
  br_A <- seq(0, 4500, by = 150); br_B <- seq(0, 9000, by = 300)
  hist_frames <- function(x, br) purrr::map_dfr(n_frames, function(n) {
    h <- hist(pmin(x[1:n], max(br) - 1e-6), breaks = br, plot = FALSE)
    tibble(mid = h$mids, dens = h$density, n = as.character(n))
  })
  hA <- hist_frames(sales_A, br_A); hB <- hist_frames(sales_B, br_B)
  curve_A <- tibble(x = seq(0, 4500, length.out = 300)) |> mutate(y = dgamma(x, shape = coffee_k, scale = mu_A / coffee_k))
  curve_B <- purrr::map_dfr(n_frames, function(n)
    tibble(x = seq(0, 9000, length.out = 300), y = dgamma(x, shape = coffee_k, scale = m_B[n] / coffee_k),
           n = as.character(n)))
  txt <- function(v, xpos, ypos) purrr::map_dfr(n_frames, function(nn)
    tibble(x = xpos, y = ypos, n = as.character(nn),
           text = sprintf("n = %s days<br>sample mean = %s", format(nn, big.mark = ","),
                          format(round(mean(v[1:nn])), big.mark = ","))))
  pA <- plot_ly() |>
    add_bars(data = hA, x = ~mid, y = ~dens, frame = ~n, width = 150,
             marker = list(color = "rgba(111,122,140,0.6)", line = list(color = "white", width = 0.5)),
             hoverinfo = "skip", showlegend = FALSE) |>
    add_lines(data = curve_A, x = ~x, y = ~y, line = list(color = "black", width = 2.5),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_text(data = txt(sales_A, 4400, 0.00168), x = ~x, y = ~y, text = ~text, frame = ~n,
             textposition = "bottom left", textfont = list(size = 12), hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)
  pB <- plot_ly() |>
    add_bars(data = hB, x = ~mid, y = ~dens, frame = ~n, width = 300,
             marker = list(color = "rgba(111,122,140,0.6)", line = list(color = "white", width = 0.5)),
             hoverinfo = "skip", showlegend = FALSE) |>
    add_lines(data = curve_B, x = ~x, y = ~y, frame = ~n, line = list(color = col_a2, width = 2.5),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_text(data = txt(sales_B, 8800, 0.00168), x = ~x, y = ~y, text = ~text, frame = ~n,
             textposition = "bottom left", textfont = list(size = 12), hoverinfo = "skip", showlegend = FALSE, inherit = FALSE)
  fig <- subplot(pA, pB, margin = 0.06, titleX = TRUE, titleY = TRUE) |>
    layout(xaxis  = list(title = "daily sales ($), shop A", range = c(0, 4500)),
           yaxis  = list(title = "density", range = c(0, 0.0017)),
           xaxis2 = list(title = "daily sales ($), shop B", range = c(0, 9000)),
           yaxis2 = list(title = "", range = c(0, 0.0017)),
           annotations = list(
             list(text = "<b>Shop A: the same DGP every day</b>", x = 0, y = 1.02, xref = "paper", yref = "paper",
                  xanchor = "left", yanchor = "bottom", showarrow = FALSE, font = list(size = 13)),
             list(text = "<b>Shop B: a new DGP every day</b>", x = 0.53, y = 1.02, xref = "paper", yref = "paper",
                  xanchor = "left", yanchor = "bottom", showarrow = FALSE, font = list(size = 13))),
           margin = list(l = 60, r = 20, t = 90, b = 60), bargap = 0) |>
    animation_opts(frame = 700, transition = 0, redraw = TRUE) |>
    animation_button(x = 0, xanchor = "left", y = -0.2, yanchor = "top", label = "Play") |>
    animation_slider(x = 0.12, len = 0.88) |>
    finalize_animation(jump_buttons(c("10 days", "100 days", "3,000 days"), c("10", "100", "3000")),
                       start_frame = "50", prefix = "days observed: ")
  if (!is.null(height)) fig$height <- height
  fig
}
fig_coffee_mean <- function() {
  n_grid <- unique(round(10^seq(0, log10(coffee_n), length.out = 300)))
  run_df <- bind_rows(
    tibble(shop = "Shop A: same DGP", n = n_grid, running = cumsum(sales_A)[n_grid] / n_grid, target = mu_A),
    tibble(shop = "Shop B: a new DGP every day", n = n_grid, running = cumsum(sales_B)[n_grid] / n_grid,
           target = m_B[n_grid]))
  ggplot(run_df, aes(n)) +
    geom_line(aes(y = target, colour = "expected sales of the latest day"), linewidth = 1, linetype = 2) +
    geom_line(aes(y = running, colour = "sample mean so far"), linewidth = 0.9) +
    scale_colour_manual(values = c(col_a2, "black"), name = NULL) +
    scale_x_log10(labels = scales::comma) +
    scale_y_continuous(labels = scales::dollar) +
    facet_wrap(~ shop) +
    labs(x = "days observed (log scale)", y = "daily sales") +
    theme(legend.position = "bottom")
}

# =============================================================================
# One random variable: distribution and moments
# =============================================================================

sd_A <- mu_A / sqrt(coffee_k)                     # 600

fig_pmf_pdf <- function() {
  lambda <- 4
  cust <- tibble(x = 0:13, p = dpois(x, lambda)) |> mutate(shade = x <= 2)
  p_cust <- ggplot(cust, aes(x, p, fill = shade)) +
    geom_col(width = 0.6) +
    scale_fill_manual(values = c(`TRUE` = col_a2, `FALSE` = "grey65"), guide = "none") +
    scale_x_continuous(breaks = seq(0, 13, by = 2)) +
    labs(title = "Discrete: customers in ten minutes",
         subtitle = sprintf("pmf, Poisson with mean %d. Red: P(X ≤ 2) = %.2f", lambda, ppois(2, lambda)),
         x = "number of customers x", y = "p(x) = P(X = x)")
  a_lo <- 1000; b_hi <- 1500
  p_ab <- pgamma(b_hi, coffee_k, scale = mu_A / coffee_k) - pgamma(a_lo, coffee_k, scale = mu_A / coffee_k)
  dens <- tibble(x = seq(0, 4000, length.out = 400)) |> mutate(f = dgamma(x, coffee_k, scale = mu_A / coffee_k))
  p_pdf <- ggplot(dens, aes(x, f)) +
    geom_area(data = dplyr::filter(dens, x >= a_lo, x <= b_hi), fill = col_a1, alpha = 0.45) +
    geom_line(linewidth = 1) +
    annotate("text", x = 2200, y = 0.00047, label = sprintf("shaded area = %.2f\n= P(1000 ≤ X ≤ 1500)", p_ab),
             colour = col_a1, hjust = 0, size = 3.8) +
    scale_x_continuous(labels = scales::dollar) +
    labs(title = "Continuous: daily sales of shop A", subtitle = "pdf: probability = area under the curve",
         x = "daily sales x", y = "f(x)")
  p_cust + p_pdf
}

# Mean and variance, side by side: a discrete distribution (one die, pmf) and a continuous one (daily sales, pdf).
fig_moments_pair <- function() {
  mu_d <- 3.5; sd_d <- sqrt(35 / 12)
  p_d <- ggplot(tibble(x = 1:6, p = 1 / 6), aes(x, p)) +
    geom_segment(aes(xend = x, yend = 0), colour = "grey40", linewidth = 1.2) +
    geom_point(colour = "grey25", size = 3.2) +
    geom_vline(xintercept = mu_d, colour = col_a2, linewidth = 1) +
    annotate("segment", x = mu_d - sd_d, xend = mu_d + sd_d, y = 0.06, yend = 0.06,
             colour = col_a1, linewidth = 1.2, arrow = arrow(ends = "both", length = unit(0.2, "cm"))) +
    annotate("text", x = mu_d + 0.08, y = 0.235, label = "μ = 3.5", colour = col_a2, hjust = 0, size = 4.4) +
    annotate("text", x = mu_d, y = 0.085, label = "μ ± σ,  σ = 1.71", colour = col_a1, size = 4) +
    scale_x_continuous(breaks = 1:6, limits = c(0.5, 6.5)) +
    scale_y_continuous(limits = c(0, 0.25), breaks = c(0, 1 / 6), labels = c("0", "1/6")) +
    labs(title = "Discrete: one die", subtitle = "pmf: probabilities are heights",
         x = expression("number shown"~italic(x)), y = expression(italic(p)(italic(x))))
  dens <- tibble(x = seq(0, 4000, length.out = 400)) |> mutate(f = dgamma(x, coffee_k, scale = mu_A / coffee_k))
  ytop <- max(dens$f)
  p_c <- ggplot(dens, aes(x, f)) +
    geom_area(fill = "grey90") + geom_line(linewidth = 1) +
    geom_vline(xintercept = mu_A, colour = col_a2, linewidth = 1) +
    annotate("segment", x = mu_A - sd_A, xend = mu_A + sd_A, y = 0.25 * ytop, yend = 0.25 * ytop,
             colour = col_a1, linewidth = 1.2, arrow = arrow(ends = "both", length = unit(0.2, "cm"))) +
    annotate("text", x = mu_A + 60, y = 1.08 * ytop, label = "μ = $1,200", colour = col_a2, hjust = 0, size = 4.4) +
    annotate("text", x = mu_A, y = 0.36 * ytop, label = "μ ± σ,  σ = $600", colour = col_a1, size = 4) +
    scale_x_continuous(labels = scales::dollar) +
    scale_y_continuous(breaks = NULL) +
    labs(title = "Continuous: daily sales", subtitle = "pdf: probabilities are areas",
         x = expression("daily sales"~italic(x)), y = expression(italic(f)(italic(x))))
  p_d + p_c + plot_layout(widths = c(1, 1))
}
fig_sampling_var <- function() {
  set.seed(11)
  n_s <- 30; R_s <- 20
  samples <- purrr::map_dfr(1:R_s, function(r) tibble(sample = r, x = rgamma(n_s, shape = coffee_k, scale = mu_A / coffee_k)))
  summ <- samples |> group_by(sample) |> summarise(xbar = mean(x), s = sd(x))
  ggplot(samples, aes(x, factor(sample))) +
    geom_vline(xintercept = mu_A, colour = "black", linewidth = 0.8) +
    geom_point(colour = "#6f7a8c", alpha = 0.45, size = 1.2) +
    geom_segment(data = summ, aes(x = xbar - s, xend = xbar + s, y = factor(sample), yend = factor(sample)),
                 colour = col_a1, linewidth = 1, alpha = 0.7) +
    geom_point(data = summ, aes(xbar, factor(sample)), colour = col_a2, size = 3, shape = 18) +
    scale_x_continuous(labels = scales::dollar) +
    labs(x = "daily sales", y = "month (30 days each)")
}

fig_lln <- function(height = NULL) {
  set.seed(21)
  n_lln <- 2000; paths <- 6
  g <- unique(round(10^seq(0, log10(n_lln), length.out = 32)))
  runs <- purrr::map_dfr(1:paths, function(k) {
    x <- rgamma(n_lln, shape = coffee_k, scale = mu_A / coffee_k)
    tibble(path = paste("shop", k), n = g, m = (cumsum(x) / seq_len(n_lln))[g])
  })
  lln_frames <- purrr::map_dfr(g, function(nn) runs |> dplyr::filter(n <= nn) |> mutate(frame = as.character(nn)))
  band <- tibble(n = 10^seq(0, log10(n_lln), length.out = 200)) |>
    mutate(lo = mu_A - 2 * sd_A / sqrt(n), hi = mu_A + 2 * sd_A / sqrt(n))
  pal6 <- c("#1f77b4", "#d62728", "#2ca02c", "#ff7f0e", "#7b4fb5", "#8c564b")
  fig <- plot_ly() |>
    add_ribbons(data = band, x = ~n, ymin = ~lo, ymax = ~hi, fillcolor = "rgba(0,0,0,0.08)",
                line = list(color = "rgba(0,0,0,0.25)", width = 1), hoverinfo = "skip",
                name = "μ ± 2σ/√n", inherit = FALSE) |>
    add_lines(x = c(1, n_lln), y = c(mu_A, mu_A), line = list(color = "black", width = 2),
              name = "μ = 1,200", hoverinfo = "skip", inherit = FALSE) |>
    add_lines(data = lln_frames, x = ~n, y = ~m, frame = ~frame, split = ~path, color = ~path,
              colors = pal6, line = list(width = 2), hoverinfo = "skip", inherit = FALSE, showlegend = FALSE) |>
    layout(xaxis = list(title = "days observed n (log scale)", type = "log"),
           yaxis = list(title = "running sample mean ($)", range = c(0, 3000)),
           legend = list(x = 0.7, y = 0.98), margin = list(l = 70, r = 20, t = 60, b = 60)) |>
    animation_opts(frame = 120, transition = 0, redraw = TRUE) |>
    animation_button(x = 0, xanchor = "left", y = -0.2, yanchor = "top", label = "Play") |>
    animation_slider(x = 0.12, len = 0.88) |>
    finalize_animation(jump_buttons(c("n = 10", "n = 100", "n = 2,000"),
                                    as.character(g[c(which.min(abs(g - 10)), which.min(abs(g - 100)), length(g))])),
                       start_frame = as.character(g[which.min(abs(g - 10))]), prefix = "days observed: ")
  if (!is.null(height)) fig$height <- height
  fig
}

fig_clt <- function() {
  set.seed(31)
  n_clt <- c(1, 3, 10, 30)
  clt <- purrr::map_dfr(n_clt, function(n)
    tibble(n = n, xbar = colMeans(matrix(rgamma(n * 10000, shape = coffee_k, scale = mu_A / coffee_k), nrow = n))))
  clt_curve <- purrr::map_dfr(n_clt, function(n)
    tibble(n = n, x = seq(0, 4000, length.out = 300), y = dnorm(x, mu_A, sd_A / sqrt(n))))
  lab_n <- function(n) factor(paste("n =", n), levels = paste("n =", n_clt))
  ggplot(clt |> mutate(n = lab_n(n)), aes(xbar)) +
    geom_histogram(aes(y = after_stat(density)), binwidth = 40, fill = "grey65", colour = NA) +
    geom_line(data = clt_curve |> mutate(n = lab_n(n)), aes(x, y), colour = col_a1, linewidth = 1) +
    facet_wrap(~ n, nrow = 1, scales = "free_y") +
    scale_x_continuous(limits = c(0, 3500), breaks = c(0, 1200, 2400), labels = scales::dollar) +
    labs(x = "sample mean of n days", y = "density") +
    theme(axis.text.y = element_blank(), panel.grid.minor = element_blank())
}

# =============================================================================
# Random vectors
# =============================================================================

P_joint <- matrix(c(0.15, 0.10, 0.04, 0.01,
                    0.08, 0.16, 0.12, 0.04,
                    0.01, 0.05, 0.12, 0.12), nrow = 3, byrow = TRUE)
cond_mean <- drop((P_joint / rowSums(P_joint)) %*% (0:3))

fig_joint <- function() {
  v1 <- 0:2; v2 <- 0:3
  p1 <- rowSums(P_joint); p2 <- colSums(P_joint)
  jt <- expand_grid(x1 = v1, x2 = v2) |> mutate(p = as.vector(t(P_joint)))
  cm <- tibble(x1 = v1, m = cond_mean)
  p_main <- ggplot(jt, aes(x1, x2)) +
    geom_point(aes(size = p), colour = col_a1, alpha = 0.75) +
    geom_text(data = dplyr::filter(jt, p >= 0.04), aes(label = sprintf("%.2f", p)), size = 3.4,
              colour = "white", fontface = "bold") +
    geom_text(data = dplyr::filter(jt, p < 0.04), aes(label = sprintf("%.2f", p)), size = 3.4,
              colour = "black", nudge_x = 0.1, hjust = 0) +
    geom_line(data = cm, aes(x1, m), colour = col_a2, linewidth = 1.1) +
    geom_point(data = cm, aes(x1, m), colour = col_a2, size = 3, shape = 18) +
    scale_size_area(max_size = 22, guide = "none") +
    scale_x_continuous(breaks = v1, labels = c("rainy", "cloudy", "sunny"), limits = c(-0.5, 2.5)) +
    scale_y_continuous(breaks = v2, limits = c(-0.5, 3.5)) +
    labs(x = expression("weather"~italic(X)[1]), y = expression("scoops sold"~italic(X)[2]~"(hundreds)"))
  p_top <- ggplot(tibble(x1 = v1, p = p1), aes(x1, p)) + geom_col(width = 0.4, fill = "grey60") +
    scale_x_continuous(limits = c(-0.5, 2.5), breaks = NULL) + scale_y_continuous(breaks = c(0, 0.4)) +
    labs(x = NULL, y = expression(italic(P)(italic(X)[1])))
  p_right <- ggplot(tibble(x2 = v2, p = p2), aes(x2, p)) + geom_col(width = 0.4, fill = "grey60") + coord_flip() +
    scale_x_continuous(limits = c(-0.5, 3.5), breaks = NULL) + labs(x = NULL, y = expression(italic(P)(italic(X)[2])))
  p_top + plot_spacer() + p_main + p_right + plot_layout(ncol = 2, widths = c(4, 1), heights = c(1, 4))
}

# A sample of 200 flights from N(0, air_Sigma): mean vector and the one-standard-deviation ellipse
# Mean vector and covariance matrix, side by side: the joint pmf of two dice (first die, sum) and the joint pdf
# of the flights (normal). Red: the mean vector. Black: the one-standard-deviation ellipse of Sigma.
fig_mean_cov_pair <- function() {
  S_d <- 35 / 12 * matrix(c(1, 1, 1, 2), 2, 2); mu_dv <- c(3.5, 7)
  jt <- expand_grid(d1 = 1:6, d2 = 1:6) |> mutate(x1 = d1, x2 = d1 + d2)
  e_d <- ellipse_xy(sym_root(S_d), 1)
  p_d <- ggplot(jt, aes(x1, x2)) +
    geom_point(colour = "grey35", size = 2.6) +
    geom_path(data = tibble(x1 = e_d[1, ] + mu_dv[1], x2 = e_d[2, ] + mu_dv[2]), colour = "black", linewidth = 1) +
    annotate("point", x = mu_dv[1], y = mu_dv[2], colour = col_a2, size = 4.5, shape = 18) +
    annotate("text", x = mu_dv[1] + 0.35, y = mu_dv[2] - 0.6, label = "μ", colour = col_a2, size = 6) +
    scale_x_continuous(breaks = 1:6) + scale_y_continuous(breaks = 2:12) +
    labs(title = "Discrete: two dice", subtitle = "joint pmf: each point has probability 1/36",
         x = expression(italic(X)[1]*": first die"), y = expression(italic(X)[2]*": sum"))
  g <- seq(-4.5, 4.5, length.out = 120)
  grid <- expand_grid(x1 = g, x2 = g) |> mutate(f = mvtnorm::dmvnorm(cbind(x1, x2), sigma = air_Sigma))
  e_c <- ellipse_xy(sym_root(air_Sigma), 1)
  p_c <- ggplot(grid, aes(x1, x2)) +
    geom_raster(aes(fill = f), interpolate = TRUE, show.legend = FALSE) +
    scale_fill_gradient(low = "white", high = "grey45") +
    geom_path(data = tibble(x1 = e_c[1, ], x2 = e_c[2, ]), colour = "black", linewidth = 1) +
    annotate("point", x = 0, y = 0, colour = col_a2, size = 4.5, shape = 18) +
    annotate("text", x = 0.45, y = -0.6, label = "μ", colour = col_a2, size = 6) +
    coord_equal(expand = FALSE) +
    labs(title = "Continuous: flights", subtitle = "joint pdf: darker = higher density",
         x = expression("flight time"~italic(X)[1]~"(centred)"), y = expression("fuel"~italic(X)[2]~"(centred)"))
  p_d + p_c + plot_layout(widths = c(1, 1))
}
fig_conditional <- function() {
  beta <- air_Sigma[1, 2] / air_Sigma[1, 1]
  sd_c <- sqrt(air_Sigma[2, 2] - air_Sigma[1, 2]^2 / air_Sigma[1, 1])
  at <- c(-2, 0, 2)
  slices <- purrr::map_dfr(at, function(a) {
    y <- seq(beta * a - 2.8 * sd_c, beta * a + 2.8 * sd_c, length.out = 150)
    tibble(at = a, x2 = y, x1 = a + 1.3 * dnorm(y, beta * a, sd_c))
  })
  ell <- purrr::map_dfr(1:3, function(c0) { e <- ellipse_xy(sym_root(air_Sigma), c0); tibble(c = c0, x1 = e[1, ], x2 = e[2, ]) })
  ggplot() +
    geom_path(data = ell, aes(x1, x2, group = c), colour = "grey55") +
    geom_segment(data = tibble(at = at), aes(x = at, xend = at, y = beta * at - 2.8 * sd_c, yend = beta * at + 2.8 * sd_c),
                 colour = "grey65", linetype = 3) +
    geom_path(data = slices, aes(x1, x2, group = at), colour = col_a1, linewidth = 1.1) +
    geom_abline(slope = beta, intercept = 0, colour = col_a2, linewidth = 1.1) +
    annotate("label", x = 2.3, y = -4, label = sprintf("italic(E)*'['*italic(X)[2]~'|'~italic(X)[1]*' = '*italic(x)[1]*']' == %.2f*italic(x)[1]", beta),
             parse = TRUE, colour = col_a2, size = 4, label.size = 0) +
    coord_equal(xlim = c(-4.5, 4.5), ylim = c(-5, 5)) +
    labs(x = expression("flight time"~italic(X)[1]), y = expression("fuel"~italic(X)[2]))
}

# =============================================================================
# Roadmap: from the reduced form to structural models
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

# =============================================================================
# The normal distribution and density contours
# =============================================================================

fig_normal_3d <- function(height = NULL) {
  gx <- seq(-4, 4, length.out = 70)
  fz <- outer(gx, gx, function(a, b) mvtnorm::dmvnorm(cbind(a, b), sigma = air_Sigma))
  p <- plot_ly(x = gx, y = gx, z = t(fz), type = "surface", colorscale = "Viridis", showscale = FALSE,
               contours = list(z = list(show = TRUE, usecolormap = TRUE, project = list(z = TRUE)))) |>
    layout(scene = list(xaxis = list(title = "flight time X₁"), yaxis = list(title = "fuel X₂"),
                        zaxis = list(title = "density", showticklabels = FALSE),
                        camera = list(eye = list(x = 1.7, y = -1.8, z = 1.1)),
                        aspectratio = list(x = 1, y = 1, z = 0.7)),
           margin = list(l = 0, r = 0, t = 10, b = 0))
  if (!is.null(height)) p$height <- height
  p
}
# The normal hill of fig_normal_3d(), seen from above: a contour plot (elevation lines) with 200 draws
fig_normal_contour <- function(points = TRUE) {
  gx <- seq(-4, 4, length.out = 160)
  grid <- expand_grid(x = gx, y = gx) |>
    mutate(z = mvtnorm::dmvnorm(cbind(x, y), sigma = air_Sigma))
  p <- ggplot(grid, aes(x, y, z = z)) +
    geom_contour_filled(bins = 8) +
    geom_contour(bins = 8, colour = "white", linewidth = 0.3) +
    scale_fill_viridis_d(guide = "none")
  if (points) {
    set.seed(7)
    pts <- as_tibble(mvtnorm::rmvnorm(200, sigma = air_Sigma), .name_repair = ~ c("x", "y"))
    p <- p + geom_point(data = pts, aes(x, y), inherit.aes = FALSE, colour = "white", alpha = 0.7, size = 1)
  }
  p + coord_equal(xlim = c(-4, 4), ylim = c(-4, 4), expand = FALSE) +
    labs(x = expression("flight time"~italic(X)[1]), y = expression("fuel"~italic(X)[2]))
}
fig_rho_slider <- function() {
  set.seed(51)
  E <- matrix(rnorm(600), ncol = 2)
  rho_grid <- round(c(-0.95, seq(-0.9, 0.9, by = 0.1), 0.95), 2)
  rl <- function(r) sprintf("%.2f", r)
  pts_r <- purrr::map_dfr(rho_grid, function(r) {
    A <- matrix(c(1, 0, r, sqrt(1 - r^2)), 2, byrow = TRUE); X <- E %*% t(A)
    tibble(id = seq_len(nrow(X)), x1 = X[, 1], x2 = X[, 2], k = rl(r))
  })
  ell_r <- purrr::map_dfr(rho_grid, function(r) {
    S <- matrix(c(1, r, r, 1), 2)
    purrr::map_dfr(1:3, function(c0) { e <- ellipse_xy(sym_root(S), c0); tibble(x1 = c(e[1, ], NA), x2 = c(e[2, ], NA)) }) |>
      mutate(k = rl(r))
  })
  reg_r <- purrr::map_dfr(rho_grid, function(r) tibble(x1 = c(-4, 4), x2 = r * c(-4, 4), k = rl(r)))
  txt_r <- tibble(x1 = -3.9, x2 = 3.9, k = rl(rho_grid),
                  text = sprintf("E[X<sub>2</sub> | X<sub>1</sub> = x<sub>1</sub>] = %.2f x<sub>1</sub><br>Var(X<sub>2</sub> | X<sub>1</sub>) = %.2f", rho_grid, 1 - rho_grid^2))
  plot_ly() |>
    add_trace(data = pts_r, x = ~x1, y = ~x2, frame = ~k, ids = ~id, type = "scatter", mode = "markers",
              marker = list(size = 5, color = "#6f7a8c", opacity = 0.5), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = ell_r, x = ~x1, y = ~x2, frame = ~k, type = "scatter", mode = "lines",
              line = list(color = "black", width = 1.6), hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(data = reg_r, x = ~x1, y = ~x2, frame = ~k, type = "scatter", mode = "lines",
              line = list(color = col_a2, width = 2.5), hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    add_trace(data = txt_r, x = ~x1, y = ~x2, frame = ~k, text = ~text, type = "scatter", mode = "text",
              textposition = "bottom right", textfont = list(size = 13, color = "#2c3e50"),
              hoverinfo = "skip", showlegend = FALSE, inherit = FALSE) |>
    layout(xaxis = list(title = "x<sub>1</sub>", range = c(-4, 4), zeroline = FALSE),
           yaxis = list(title = "x<sub>2</sub>", range = c(-4, 4), scaleanchor = "x", zeroline = FALSE),
           margin = list(l = 60, r = 20, t = 50, b = 60)) |>
    animation_opts(frame = 150, transition = 0, redraw = TRUE) |>
    animation_button(x = 0, xanchor = "left", y = -0.12, yanchor = "top", label = "Play") |>
    animation_slider(x = 0.12, len = 0.88) |>
    finalize_animation(jump_buttons(c("ρ = −0.8", "ρ = 0", "ρ = 0.5", "ρ = 0.9"), rl(c(-0.8, 0, 0.5, 0.9))),
                       start_frame = "0.50", prefix = "correlation ρ = ")
}
peng_kde <- local({
  k <- MASS::kde2d(peng$bill_length_mm, peng$bill_depth_mm, n = 120,
                   lims = c(30, 62, 12.5, 22.5))
  list(x = k$x, y = k$y, z = k$z,
       df = expand_grid(y = k$y, x = k$x) |> mutate(z = as.vector(t(k$z))) |> dplyr::select(x, y, z))
})
fig_kde_3d <- function(height = NULL) {
  p <- plot_ly(x = peng_kde$x, y = peng_kde$y, z = t(peng_kde$z), type = "surface",
               colorscale = list(c(0, "#f7fbff"), c(0.3, "#9ecae1"), c(0.7, "#2171b5"), c(1, "#08306b")),
               showscale = FALSE,
               contours = list(z = list(show = TRUE, usecolormap = TRUE, project = list(z = TRUE)))) |>
    layout(scene = list(xaxis = list(title = "bill length X₁"), yaxis = list(title = "bill depth X₂"),
                        zaxis = list(title = "density", showticklabels = FALSE),
                        camera = list(eye = list(x = -1.2, y = -1.9, z = 1.2)),
                        aspectratio = list(x = 1.2, y = 1, z = 0.6)),
           margin = list(l = 0, r = 0, t = 10, b = 0))
  if (!is.null(height)) p$height <- height
  p
}
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
# Mahalanobis distance and the chi-square distribution
# =============================================================================

fig_chisq_share <- function() {
  set.seed(61)
  n_c <- 2000
  Xs <- mvtnorm::rmvnorm(n_c, sigma = air_Sigma)
  d2s <- mahalanobis(Xs, center = c(0, 0), cov = air_Sigma)
  band_lab <- cut(sqrt(d2s), c(0, 1, 2, 3, Inf), labels = c("d < 1", "1 to 2", "2 to 3", "d > 3"))
  shares <- round(100 * as.numeric(table(band_lab)) / n_c, 1)
  pal4 <- c("d < 1" = "#1f77b4", "1 to 2" = "#7b4fb5", "2 to 3" = "#ff7f0e", "d > 3" = "#d62728")
  ell <- purrr::map_dfr(1:3, function(c0) { e <- ellipse_xy(sym_root(air_Sigma), c0); tibble(c = c0, x1 = e[1, ], x2 = e[2, ]) })
  ggplot(tibble(x1 = Xs[, 1], x2 = Xs[, 2], band = band_lab), aes(x1, x2, colour = band)) +
    geom_point(size = 0.8, alpha = 0.6) +
    geom_path(data = ell, aes(x1, x2, group = c), colour = "black", inherit.aes = FALSE) +
    scale_colour_manual(values = pal4, name = NULL, labels = sprintf("%s (%.1f%%)", levels(band_lab), shares)) +
    coord_equal(xlim = c(-5, 5), ylim = c(-5, 5)) +
    guides(colour = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = expression("flight time"~italic(X)[1]), y = expression("fuel"~italic(X)[2])) +
    theme(legend.position = "bottom", legend.text = element_text(size = 10)) +
    guides(colour = guide_legend(nrow = 2, override.aes = list(size = 3, alpha = 1)))
}
fig_state_qq <- function() {
  d2 <- d2_st; n <- length(d2)
  qq <- tibble(id = names(d2), d2 = d2) |> arrange(d2) |>
    mutate(q = qchisq(ppoints(n), df = 6), flag = d2 > thr6)
  ggplot(qq, aes(q, d2)) +
    geom_abline(slope = 1, intercept = 0, colour = "grey40") +
    geom_hline(yintercept = thr6, linetype = 2) +
    geom_point(aes(colour = flag), size = 2.3) +
    geom_text_repel(data = qq |> dplyr::slice_max(d2, n = 4), aes(label = id),
                    size = 3.8, nudge_x = -3, direction = "y", min.segment.length = 0, seed = 1) +
    annotate("text", x = 0.3, y = thr6 + 0.8, hjust = 0, size = 3.6,
             label = sprintf("97.5%% cut-off = %.1f", thr6)) +
    scale_colour_manual(values = c(`FALSE` = "#6f7a8c", `TRUE` = col_a2), guide = "none") +
    labs(x = "χ²(6) quantile", y = "ordered squared Mahalanobis distance")
}
