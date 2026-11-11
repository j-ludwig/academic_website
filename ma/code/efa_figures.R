# =============================================================================
# efa_figures.R  --  data and figures for the chapter "Exploratory Factor Analysis" of the ISQS 6350 course
# notes. Put this file and _common.R (shared helpers) in one folder, set it as the working directory, and run
# source("efa_figures.R"). Every figure is a function, e.g. fig_...().
# Packages: tidyverse, plotly, patchwork, MASS, mvtnorm, ggrepel
# =============================================================================

source("_common.R", local = TRUE)
suppressPackageStartupMessages(library(ggrepel))

f0 <- function(x) formatC(x, format = "f", digits = 0, big.mark = ",")
f1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
f3 <- function(x) formatC(x, format = "f", digits = 3)
pct <- function(x, d = 0) paste0(formatC(100 * x, format = "f", digits = d), "%")
mns <- function(s) sub("^-", "\u2212", s)                # a proper minus sign for labels in figures
ptxt <- function(p) if (p < 1e-4) "< 0.0001" else if (p < 0.001) paste("=", formatC(p, format = "f", digits = 4)) else paste("=", f3(p))
pv <- function(p) if (p < 1e-4) "below 0.0001" else if (p < 0.001) paste("of", formatC(p, format = "f", digits = 4)) else paste("of", f3(p))

col_verb <- col_a1                                   # verbal ability (blue) and mathematical ability (red)
col_math <- col_a2
col_third <- "#e67e22"                               # the third critic
col_heat <- c(low = "#f7fbff", high = "#08306b")
axes0 <- list(geom_hline(yintercept = 0, colour = "grey88"), geom_vline(xintercept = 0, colour = "grey88"))


# =============================================================================
# 1  Critics: one common factor
# =============================================================================
# Each critic's rating of a restaurant = loading x quality + the critic's quirks (her unique shock). Critics 1 and 2
# are the two critics of the factor sliders in _common.R (variances 1 and 1.25, covariance 0.75); critics 3 and 4 join
# them here.
crit_lam <- c(0.75, 1.00, 0.60, 0.90)                # what one standard deviation of quality adds to each rating
crit_psi <- c(0.4375, 0.25, 0.64, 0.19)              # variance of each critic's quirks
crit_names <- paste("critic", 1:4)
crit_Sigma <- outer(crit_lam, crit_lam) + diag(crit_psi)
dimnames(crit_Sigma) <- list(crit_names, crit_names)

# Two critics: every pair of loadings with lambda1 * lambda2 = sigma12 fits. The four natural stories.
crit_two <- local({
  s11 <- crit_Sigma[1, 1]; s22 <- crit_Sigma[2, 2]; s12 <- crit_Sigma[1, 2]; d <- s22 - s11
  l1 <- c(sqrt(s11),                                  # critic 1 never errs: psi1 = 0
          s12 / sqrt(s22),                            # critic 2 never errs: psi2 = 0
          (s12^2 * s11 / s22)^(1 / 4),                # equally reliable: same share of variance from quality
          sqrt((-d + sqrt(d^2 + 4 * s12^2)) / 2))     # equal noise: psi1 = psi2
  tibble(story = c("critic 1 never errs", "critic 2 never errs", "equally reliable", "equal noise"),
         l1 = l1, l2 = s12 / l1,
         lx = l1 + c(-0.03, 0.05, 0.06, 0.06), ly = s12 / l1 + c(-0.07, 0.05, 0.03, 0.09), hj = c(1, 0, 0, 0))
})
# The third critic: lambda1 / lambda2 = sigma13 / sigma23, and lambda1^2 = sigma12 sigma13 / sigma23
crit_three <- with(as.list(c(s12 = crit_Sigma[1, 2], s13 = crit_Sigma[1, 3], s23 = crit_Sigma[2, 3])),
                   c(l1 = sqrt(s12 * s13 / s23), l2 = sqrt(s12 * s23 / s13), l3 = sqrt(s13 * s23 / s12)))

# Loadings of critics 1 and 2: the curve of all fits (two critics) and the third critic's line through it
fig_crit_curve <- function(panels = c("two", "three")) {
  s11 <- crit_Sigma[1, 1]; s22 <- crit_Sigma[2, 2]; s12 <- crit_Sigma[1, 2]
  ratio <- crit_Sigma[1, 3] / crit_Sigma[2, 3]
  curve_ok  <- tibble(l1 = seq(s12 / sqrt(s22), sqrt(s11), length.out = 200)) |> mutate(l2 = s12 / l1)
  curve_out <- bind_rows(tibble(l1 = seq(0.56, s12 / sqrt(s22), length.out = 60), part = "a"),
                         tibble(l1 = seq(sqrt(s11), 1.32, length.out = 60), part = "b")) |> mutate(l2 = s12 / l1)
  base <- function(title, subtitle) ggplot() +
    geom_path(data = curve_out, aes(l1, l2, group = part), colour = "#bbbbbb", linetype = "22", linewidth = 0.8) +
    geom_path(data = curve_ok, aes(l1, l2), colour = col_f, linewidth = 1.6) +
    coord_equal(xlim = c(0.4, 1.3), ylim = c(0.4, 1.3)) +
    labs(title = title, subtitle = subtitle,
         x = expression("critic 1's loading " * lambda[1]), y = expression("critic 2's loading " * lambda[2])) +
    theme(plot.title = element_text(size = 14, face = "bold"), plot.subtitle = element_text(size = 11))
  p1 <- base("Two critics", "every point on the purple curve fits the covariance matrix") +
    geom_segment(data = crit_two, aes(l1, l2, xend = lx, yend = ly), colour = "grey55", linewidth = 0.4) +
    geom_point(data = crit_two, aes(l1, l2), size = 2.8, colour = "#2c3e50") +
    geom_label(data = crit_two, aes(lx, ly, label = story, hjust = hj), size = 3.8, colour = "#2c3e50",
               label.size = 0, label.padding = unit(0.1, "lines"), fill = "white") +
    annotate("text", x = 1.27, y = 0.53, label = "dashed: a negative\nquirk variance", hjust = 1, size = 3.4,
             colour = "grey45", lineheight = 0.9)
  p2 <- base("A third critic", "a line from her covariances crosses the curve once") +
    geom_abline(slope = 1 / ratio, intercept = 0, colour = col_third, linewidth = 1.1) +
    geom_point(data = crit_two, aes(l1, l2), size = 2, colour = "#bbbbbb") +
    annotate("point", x = crit_three["l1"], y = crit_three["l2"], size = 4, colour = col_third) +
    annotate("text", x = crit_three["l1"] + 0.07, y = crit_three["l2"], hjust = 0, size = 3.8, colour = "#2c3e50",
             label = sprintf("λ₁ = %s, λ₂ = %s", f2(crit_three["l1"]), f2(crit_three["l2"]))) +
    annotate("text", x = 0.55, y = 0.52, hjust = 0, size = 3.8, colour = col_third,
             label = "the third critic's line:\nλ₁ / λ₂ = σ₁₃ / σ₂₃", lineheight = 0.9)
  if (identical(panels, "two")) p1 else if (identical(panels, "three")) p2 else p1 + p2
}

# A Heywood case with the three critics. With three variables the one-factor solution is exact,
# lambda1^2 = sigma12 sigma13 / sigma23 and so on. If critics 2 and 3 agreed much less (sigma23 = 0.25 instead of 0.60),
# critic 1's common part would exceed her whole variance. (Everitt and Hothorn, Section 5.5, show the same with
# Spearman's school marks.)
crit_S3 <- crit_Sigma[1:3, 1:3]
crit_S3h <- crit_S3; crit_S3h[2, 3] <- crit_S3h[3, 2] <- 0.25       # still a valid (positive definite) covariance matrix
solve3 <- function(S) {                                            # the exact one-factor solution for three variables
  common <- c(S[1, 2] * S[1, 3] / S[2, 3], S[1, 2] * S[2, 3] / S[1, 3], S[1, 3] * S[2, 3] / S[1, 2])   # lambda_i^2
  tibble(critic = crit_names[1:3], common = common, var = diag(S), psi = diag(S) - common, h = common / diag(S))
}
# A story with two factors that produces the changed covariances exactly: food and service, correlated 0.3. Critic 1
# judges both, critic 2 only the food (loading 1) and critic 3 only the service. With three variables two factors
# are not identified, so this is one of many such stories.
heywood_two <- local({
  S <- crit_S3h; phi <- 0.3; a2 <- 1
  b3 <- S[2, 3] / (a2 * phi)                                    # sigma23 = a2 * phi * b3
  ab <- solve(matrix(c(1, phi, phi, 1), 2), c(S[1, 2] / a2, S[1, 3] / b3))   # critic 1's two loadings
  L <- rbind(ab, c(a2, 0), c(0, b3)); Phi <- matrix(c(1, phi, phi, 1), 2)
  list(L = L, phi = phi, psi = diag(S) - diag(L %*% Phi %*% t(L)))
})
heywood_lab <- function(S, title) sprintf("%s\nσ₁₂ = %s · σ₁₃ = %s · σ₂₃ = %s", title, f2(S[1, 2]), f2(S[1, 3]), f2(S[2, 3]))
heywood_tab <- bind_rows(
  solve3(crit_S3)  |> mutate(case = heywood_lab(crit_S3, "The three critics")),
  solve3(crit_S3h) |> mutate(case = heywood_lab(crit_S3h, "Critics 2 and 3 rarely agree")))
fig_heywood <- function() {
  d <- heywood_tab |> mutate(case = factor(case, levels = unique(case)), xi = match(critic, crit_names), u = 1 - h)
  w <- 0.3
  ggplot(d) +
    geom_rect(aes(xmin = xi - w, xmax = xi + w, ymin = 0, ymax = h, fill = "communality h (quality)")) +
    geom_rect(data = filter(d, u >= 0), aes(xmin = xi - w, xmax = xi + w, ymin = h, ymax = 1,
                                            fill = "uniqueness 1 − h (quirks)")) +
    geom_hline(yintercept = 1, linetype = "22", colour = "grey30") +
    geom_text(data = filter(d, u >= 0), aes(xi, 1.07, label = sprintf("1 − h = %s", f2(u))), size = 4.2,
              colour = "#2c3e50") +
    geom_segment(data = filter(d, u < 0), aes(x = xi + w + 0.09, xend = xi + w + 0.09, y = h, yend = 1),
                 colour = "#c0392b", linewidth = 1, arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_text(data = filter(d, u < 0), aes(xi + w + 0.17, h, label = sprintf("1 − h = %s\nimpossible", mns(f2(u)))),
              hjust = 0, vjust = 1, size = 4.2, colour = "#c0392b", lineheight = 0.9) +
    scale_fill_manual(values = c(`communality h (quality)` = col_f, `uniqueness 1 − h (quirks)` = "#c9c9c9"), name = NULL) +
    scale_x_continuous(breaks = 1:3, labels = crit_names[1:3]) +
    scale_y_continuous(limits = c(0, 1.58), breaks = seq(0, 1.5, 0.5), expand = c(0, 0)) +
    facet_wrap(~ case) +
    labs(x = NULL, y = "share of the critic's rating variance",
         caption = "dashed line: the whole variance of each rating") +
    theme(legend.position = "top", strip.text = element_text(size = 12, face = "bold"),
          panel.grid.major.x = element_blank(), axis.text.x = element_text(size = 12))
}

fig_crit_split <- function() {
  heat <- function(M, title) {
    d <- as_tibble(as.table(M), .name_repair = ~ c("r", "c", "v")) |>
      mutate(r = factor(r, levels = rev(crit_names)), c = factor(c, levels = crit_names))
    ggplot(d, aes(c, r, fill = v)) + geom_tile(colour = "white", linewidth = 0.8) +
      geom_text(aes(label = ifelse(abs(v) < 1e-9, "0", f2(v)), colour = v > 0.7), size = 4) +
      scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
      scale_fill_gradient(low = col_heat["low"], high = col_heat["high"], limits = c(0, 1.3), guide = "none") +
      coord_equal() + labs(title = title, x = NULL, y = NULL) +
      theme(axis.text.x = element_text(angle = 35, hjust = 1), panel.grid = element_blank(),
            plot.title = element_text(size = 13, face = "bold", hjust = 0.5))
  }
  common <- outer(crit_lam, crit_lam); dimnames(common) <- dimnames(crit_Sigma)
  unique_ <- diag(crit_psi); dimnames(unique_) <- dimnames(crit_Sigma)
  sym <- function(s) wrap_elements(grid::textGrob(s, gp = grid::gpar(fontsize = 30, col = "#2c3e50")))
  heat(crit_Sigma, "covariance matrix Σ") + sym("=") + heat(common, "quality: λλ′") + sym("+") +
    heat(unique_, "quirks: Ψ") + plot_layout(nrow = 1, widths = c(1, 0.12, 1, 0.12, 1))
}


# =============================================================================
# 2  How many factors can the data pin down?
# =============================================================================
# Degrees of freedom of the k-factor model for p variables: p(p+1)/2 covariances minus the free parameters
# pk + p - k(k-1)/2 (loadings and unique variances, less the angles of a rotation of the factors)
efa_df <- function(p, k) ((p - k)^2 - (p + k)) / 2

fig_efa_count <- function(pmax = 10, kmax = 4) {
  d <- expand_grid(p = 2:pmax, k = 1:kmax) |> filter(k < p) |>
    mutate(df = efa_df(p, k),
           status = factor(case_when(df < 0 ~ "not identified", df == 0 ~ "just identified", TRUE ~ "testable"),
                           levels = c("not identified", "just identified", "testable")))
  ggplot(d, aes(p, k, fill = status)) + geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = df), size = 4.2, colour = "#2c3e50") +
    scale_fill_manual(values = c(`not identified` = "#f4c7c3", `just identified` = "#fde9b8", testable = "#cfe8d4"),
                      name = NULL) +
    scale_x_continuous(breaks = 2:pmax, expand = c(0, 0)) + scale_y_continuous(breaks = 1:kmax, expand = c(0, 0)) +
    coord_equal() +
    labs(x = "number of variables p", y = "number of factors k",
         subtitle = "number in each tile: covariances left over after fitting the model (degrees of freedom)") +
    theme(panel.grid = element_blank(), legend.position = "top", plot.subtitle = element_text(size = 11))
}


# =============================================================================
# 3  A sample: 200 restaurants rated by the four critics
# =============================================================================
# together = correlation of the quirks of critics 3 and 4 (they dine together and share their bad evenings)
crit_sim <- function(seed = 3, together = 0, n = 200) {
  set.seed(seed)
  f <- rnorm(n); E <- matrix(rnorm(n * 4), n); z <- rnorm(n)
  E[, 3:4] <- sqrt(1 - together) * E[, 3:4] + sqrt(together) * z
  X <- outer(f, crit_lam) + E %*% diag(sqrt(crit_psi)); colnames(X) <- crit_names
  list(X = X, quality = f)
}
crit_d   <- crit_sim()                               # the one-factor model holds
crit_X   <- crit_d$X
crit_fa  <- factanal(crit_X, factors = 1)
crit_X2  <- crit_sim(together = 0.5)$X               # critics 3 and 4 dine together
crit_fa2 <- factanal(crit_X2, factors = 1)

# One factor and four critics leave two testable restrictions: the three products s12 s34, s13 s24 and s14 s23 of
# sample covariances must be equal, up to sampling noise (in the population each equals l1 l2 l3 l4)
crit_tetrads <- function(X) {
  S <- cov(X); c(S[1, 2] * S[3, 4], S[1, 3] * S[2, 4], S[1, 4] * S[2, 3])
}
# Critic 1 switches from a 10-point to a 100-point scale (her ratings times ten). The share of each critic's rating
# variance that quality explains (communality) does not change, although her loading is ten times as large. The
# first principal component swings to her, because her variance is now a hundred times as large.
fig_units <- function() {
  alpha <- c(10, 1, 1, 1); X10 <- crit_X %*% diag(alpha)
  pc <- function(X) { p <- prcomp(X); (p$rotation[, 1] * p$sdev[1])^2 / apply(X, 2, var) }   # share kept by PC 1
  h  <- function(X) 1 - factanal(X, factors = 1)$uniquenesses                                # communalities
  lv <- c("original ratings", "critic 1 on a 100-point scale")
  pn <- c("PCA: share kept by the first component", "Factor analysis: communality")
  d <- bind_rows(
    tibble(panel = pn[1], critic = crit_names, v = pc(crit_X), when = lv[1]),
    tibble(panel = pn[1], critic = crit_names, v = pc(X10), when = lv[2]),
    tibble(panel = pn[2], critic = crit_names, v = h(crit_X), when = lv[1]),
    tibble(panel = pn[2], critic = crit_names, v = h(X10), when = lv[2])) |>
    mutate(panel = factor(panel, levels = unique(panel)), when = factor(when, levels = lv))
  dodge <- position_dodge(width = 0.78)
  ggplot(d, aes(critic, v, fill = when)) +
    geom_col(position = dodge, width = 0.72) +
    geom_text(aes(label = pct(v)), position = dodge, vjust = -0.35, size = 3.5, colour = "#2c3e50") +
    facet_wrap(~ panel) +
    scale_fill_manual(values = c("#bdbdbd", col_a1), name = NULL) +
    scale_y_continuous(limits = c(0, 1.12), breaks = seq(0, 1, 0.25), labels = scales::percent, expand = c(0, 0)) +
    labs(x = NULL, y = "share of the critic's rating variance") +
    theme(legend.position = "top", strip.text = element_text(size = 12, face = "bold"),
          panel.grid.major.x = element_blank())
}

fig_crit_tetrads <- function() {
  lab <- c("s₁₂ · s₃₄", "s₁₃ · s₂₄", "s₁₄ · s₂₃")
  d <- bind_rows(
    tibble(sample = sprintf("One quality factor\ntest of one factor: p-value %s", ptxt(crit_fa$PVAL)),
           prod = lab, v = crit_tetrads(crit_X)),
    tibble(sample = sprintf("Critics 3 and 4 dine together\ntest of one factor: p-value %s", ptxt(crit_fa2$PVAL)),
           prod = lab, v = crit_tetrads(crit_X2))) |>
    mutate(sample = factor(sample, levels = unique(sample)), prod = factor(prod, levels = lab))
  ggplot(d, aes(prod, v)) + geom_col(fill = col_f, width = 0.6) +
    geom_text(aes(label = f3(v)), vjust = -0.5, size = 4) +
    facet_wrap(~ sample) + scale_y_continuous(limits = c(0, 0.55), expand = c(0, 0)) +
    labs(x = NULL, y = "product of two sample covariances") +
    theme(strip.text = element_text(size = 11.5, face = "bold"), axis.text.x = element_text(size = 12))
}


# =============================================================================
# 4  Two factors: six tests, two abilities
# =============================================================================
# A new example with new variables: the scores of students on three verbal and three mathematics tests. Students who
# are strong in one ability tend to be strong in the other: the two factors have correlation 0.5. Each test measures
# one ability only (simple structure).
tests_items <- c("reading", "vocabulary", "grammar", "arithmetic", "algebra", "geometry")
tests_group <- setNames(rep(c("verbal", "math"), each = 3), tests_items)
tests_P <- rbind(c(0.8, 0), c(0.7, 0), c(0.6, 0), c(0, 0.8), c(0, 0.7), c(0, 0.6))   # loadings (pattern)
tests_Phi <- matrix(c(1, 0.5, 0.5, 1), 2)                                        # correlation of the factors
tests_h2 <- diag(tests_P %*% tests_Phi %*% t(tests_P))                           # communalities
tests_X <- local({
  set.seed(4); n <- 200
  F <- mvtnorm::rmvnorm(n, sigma = tests_Phi)
  X <- F %*% t(tests_P) + matrix(rnorm(n * 6), n) %*% diag(sqrt(1 - tests_h2))
  colnames(X) <- tests_items; X
})
tests_fa1 <- factanal(tests_X, factors = 1)
tests_fa  <- factanal(tests_X, factors = 2, rotation = "none")
tests_L0  <- unclass(tests_fa$loadings)             # the loadings as estimated, before any rotation

# Kaiser's varimax criterion: the spread of the squared loadings (rows scaled to communality 1) within each factor
varimax_crit <- function(L) { L2 <- (L / sqrt(rowSums(L^2)))^2; sum(colMeans(L2^2) - colMeans(L2)^2) }
tests_turn <- function(deg) tests_L0 %*% rotation_matrix(deg * pi / 180)   # loadings with the factors turned
tests_vm_deg <- optimize(function(a) -varimax_crit(tests_turn(a)), c(-90, 0), tol = 1e-8)$minimum
tests_vm <- tests_turn(tests_vm_deg)                 # varimax loadings: verbal first, math second
tests_pm <- local({                                  # promax: an oblique rotation that starts from varimax
  r <- promax(tests_L0); U <- r$rotmat
  list(L = unclass(r$loadings), Phi = solve(t(U) %*% U), axes = solve(U))   # rows of axes: the factor directions
})
tests_pm_angle <- acos(tests_pm$Phi[1, 2]) * 180 / pi
tests_fv <- factanal(tests_X, factors = 2)                         # varimax, as R computes it (stops near the top)
tests_fp <- factanal(tests_X, factors = 2, rotation = "promax")    # promax, as R computes it

# Turn the two factors: the six tests move together, and only the varimax criterion changes
fig_tests_rotate <- function(step = 2) {
  angs <- sort(unique(c(seq(-90, 90 - step, by = step), round(tests_vm_deg))))
  cols <- ifelse(tests_group == "verbal", col_verb, col_math)
  frames <- purrr::map(angs, function(a) {
    L <- tests_turn(a)
    list(pts = tibble(x = L[, 1], y = L[, 2], lab = tests_items, col = cols, k = a,
                      pos = c("middle right", "bottom left", "top left", "middle right", "top left", "middle left")),
         seg = tibble(x = as.vector(rbind(0, L[, 1], NA)), y = as.vector(rbind(0, L[, 2], NA)), k = a),
         txt = tibble(x = 0, y = 0, k = a,
                      txt = sprintf("Turning the factors changes no correlation and no communality.<br>varimax criterion: %s",
                                    f3(varimax_crit(L)))))
  })
  pts <- bind_rows(purrr::map(frames, "pts")); seg <- bind_rows(purrr::map(frames, "seg"))
  txt <- bind_rows(purrr::map(frames, "txt"))
  circ <- tibble(t = seq(0, 2 * pi, length.out = 200), x = cos(t), y = sin(t))
  jump <- function(lab, a) list(label = lab, method = "animate",
                                args = list(list(as.character(a)), list(mode = "immediate",
                                            frame = list(duration = 0, redraw = TRUE), transition = list(duration = 0))))
  txt_panel <- plot_ly() |>
    add_trace(data = txt, x = ~x, y = ~y, frame = ~k, text = ~txt, type = "scatter", mode = "text",
              textfont = list(size = 12, color = "#2c3e50"), hoverinfo = "skip", showlegend = FALSE)
  main <- plot_ly() |>
    add_trace(data = circ, x = ~x, y = ~y, type = "scatter", mode = "lines",
              line = list(color = "#dddddd", width = 1, dash = "dash"), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = seg, x = ~x, y = ~y, frame = ~k, type = "scatter", mode = "lines",
              line = list(color = "#bbbbbb", width = 1.5), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = pts, x = ~x, y = ~y, frame = ~k, text = ~lab, type = "scatter", mode = "markers+text",
              textposition = ~pos, marker = list(color = ~col, size = 9),
              textfont = list(size = 11, color = ~col), hoverinfo = "skip", showlegend = FALSE)
  subplot(txt_panel, main, nrows = 2, heights = c(0.1, 0.9), margin = 0.01, titleX = TRUE, titleY = TRUE) |>
    layout(xaxis = list(visible = FALSE, range = c(-1, 1), fixedrange = TRUE),
           yaxis = list(visible = FALSE, range = c(-1, 1), fixedrange = TRUE),
           xaxis2 = list(title = "loading on factor 1", range = c(-1.15, 1.15), zeroline = TRUE,
                         zerolinecolor = "#cccccc"),
           yaxis2 = list(title = "loading on factor 2", range = c(-1.15, 1.15), scaleanchor = "x2", zeroline = TRUE,
                         zerolinecolor = "#cccccc"),
           margin = list(t = 50)) |>
    animation_opts(frame = 120, transition = 0, redraw = TRUE) |>
    animation_slider(currentvalue = list(prefix = "factors turned by ", suffix = "°")) |>
    animation_button(label = "Play") |>
    finalize_animation(list(jump("as estimated", 0), jump("varimax", round(tests_vm_deg))),
                       start_frame = "0", prefix = "factors turned by ", suffix = "°")
}

# The same six tests, three choices of axes: as estimated, varimax (orthogonal), promax (oblique)
fig_tests_axes <- function() {
  pts <- tibble(x = tests_L0[, 1], y = tests_L0[, 2], lab = tests_items, g = tests_group)
  Qv <- rotation_matrix(tests_vm_deg * pi / 180)     # columns: the varimax factors, in the estimated coordinates
  axes_df <- function(D, labs) tibble(x = 0, y = 0, xend = D[1, ], yend = D[2, ], lab = labs, g = c("verbal", "math"),
                                      lx = 1.04 * D[1, ], ly = 1.04 * D[2, ],          # names just beyond the tips
                                      hj = (1 - D[1, ]) / 2, vj = (1 - D[2, ]) / 2)
  # The six points are the same in all three panels, so only the first panel names them.
  panel <- function(D, labs, title, subtitle, named = FALSE) {
    ax <- axes_df(D, labs)
    avoid <- bind_rows(purrr::map_dfr(seq(0.04, 1.2, by = 0.04), function(s)           # points along the arrows,
                         tibble(x = s * D[1, ], y = s * D[2, ], lab = "", g = c("verbal", "math"))),
                       tibble(x = ax$lx + 0.12 * D[1, ], y = ax$ly + 0.12 * D[2, ], lab = "", g = ax$g))  # and names
    ggplot() + axes0 +
      geom_segment(data = ax, aes(x, y, xend = xend, yend = yend, colour = g), linewidth = 1.1,
                   arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
      geom_text(data = ax, aes(lx, ly, label = lab, colour = g, hjust = hj, vjust = vj), size = 4, fontface = "bold") +
      geom_point(data = pts, aes(x, y, colour = g), size = 2.4) +
      (if (named) geom_text_repel(data = bind_rows(pts, avoid), aes(x, y, label = lab, colour = g), size = 3.6, seed = 1,
                                  box.padding = 0.3, point.padding = 0.15, force = 3, max.iter = 1e5,
                                  min.segment.length = 0.3, max.overlaps = Inf)) +
      scale_colour_manual(values = c(verbal = col_verb, math = col_math), guide = "none") +
      coord_equal(xlim = c(-0.45, 1.45), ylim = c(-0.7, 1.25)) +
      labs(title = title, subtitle = subtitle, x = "loading on factor 1 (as estimated)",
           y = "loading on factor 2 (as estimated)") +
      theme(plot.title = element_text(size = 13, face = "bold"), plot.subtitle = element_text(size = 11))
  }
  panel(diag(2), c("factor 1", "factor 2"), "As estimated", "a general factor and a contrast", named = TRUE) +
    panel(Qv, c("verbal", "math"), "Varimax", sprintf("turned by %s°, still at right angles", mns(f0(tests_vm_deg)))) +
    panel(t(tests_pm$axes), c("verbal", "math"), "Promax",
          sprintf("axes through the two groups, %s° apart", f0(tests_pm_angle)))
}


# =============================================================================
# 5  Factor scores
# =============================================================================
# Regression scores with one factor: f_hat = sum_i w_i x_i with w_i = (lambda_i / psi_i) / (1 + s),
# s = sum_i lambda_i^2 / psi_i. Here with the true loadings, for the first m critics.
crit_weights <- function(m) { l <- crit_lam[1:m]; p <- crit_psi[1:m]; (l / p) / (1 + sum(l^2 / p)) }
crit_reliability <- function(m) { tau <- sum(crit_lam[1:m]^2 / crit_psi[1:m]); tau / (1 + tau) }   # rho^2 of quality

crit_fs <- factanal(crit_X, factors = 1, scores = "regression")   # scores from the estimated model

fig_crit_scores <- function() {
  d <- purrr::map_dfr(1:4, function(m) {
    lab <- sprintf("%s\nρ² = %s", c("critic 1 alone", "critics 1 and 2", "critics 1 to 3", "all four critics")[m],
                   f2(crit_reliability(m)))
    tibble(est = drop(crit_X[, 1:m, drop = FALSE] %*% crit_weights(m)), true = crit_d$quality, panel = lab)
  }) |> mutate(panel = factor(panel, levels = unique(panel)))
  ggplot(d, aes(est, true)) + geom_abline(slope = 1, intercept = 0, colour = col_f, linetype = "22") +
    geom_point(alpha = 0.35, size = 1.2, colour = "#2c3e50") +
    facet_wrap(~ panel, nrow = 1) + coord_equal(xlim = c(-3, 3), ylim = c(-3, 3)) +
    labs(x = "estimated quality (factor score)", y = "true quality") +
    theme(strip.text = element_text(size = 11, face = "bold"))
}


# One restaurant shows why its quality cannot be computed from the ratings. Left: every quality reproduces the four
# ratings, once the quirks u_i = x_i - lambda_i f make up the difference. Right: each quality weighted by how likely it
# is together with the quirks it needs (normal shocks). This is the distribution of quality given the ratings: normal,
# centred on the regression score, with variance 1 / (1 + tau).
crit_pick <- which.min(abs(drop(crit_X %*% crit_weights(4)) - 1.2))   # a restaurant with a score of about 1.2
fig_crit_anyquality <- function(t = crit_pick) {
  x <- crit_X[t, ]; tau <- sum(crit_lam^2 / crit_psi)
  f_hat <- sum(crit_weights(4) * x); sd_f <- sqrt(1 / (1 + tau)); f_true <- crit_d$quality[t]
  fs <- seq(-1.5, 3.5, length.out = 401)
  cols <- setNames(c(col_a1, col_a2, col_third, "#16a085"), crit_names)
  quirks <- purrr::map_dfr(1:4, function(i) tibble(f = fs, u = x[i] - crit_lam[i] * fs, critic = crit_names[i]))
  ends <- quirks |> group_by(critic) |> slice_tail(n = 1) |> ungroup() |> arrange(u)
  for (j in 2:nrow(ends)) ends$u[j] <- max(ends$u[j], ends$u[j - 1] + 0.3)   # keep the labels apart
  p1 <- ggplot(quirks, aes(f, u, colour = critic)) +
    geom_hline(yintercept = 0, colour = "grey70") +
    annotate("text", x = -1.5, y = 0.12, label = "no quirk", hjust = 0, vjust = 0, size = 3.4, colour = "grey45") +
    geom_line(linewidth = 1.1) +
    geom_text(data = ends, aes(label = critic), hjust = 0, nudge_x = 0.08, size = 3.6, show.legend = FALSE) +
    scale_colour_manual(values = cols, guide = "none") +
    scale_x_continuous(limits = c(-1.5, 4.4), breaks = -1:3) +
    labs(title = "Every quality fits the four ratings",
         subtitle = "the quirk each critic would need, for each quality",
         x = "quality of the restaurant", y = "quirk needed")
  dens <- tibble(f = fs, d = dnorm(fs, f_hat, sd_f))
  top <- max(dens$d)
  p2 <- ggplot(dens, aes(f, d)) +
    geom_area(fill = col_f, alpha = 0.2) + geom_line(colour = col_f, linewidth = 1.1) +
    geom_vline(xintercept = f_hat, colour = col_f, linetype = "22") +
    geom_vline(xintercept = f_true, colour = "#2c3e50") +
    annotate("text", x = f_hat + 0.08, y = top * 1.08, hjust = 0, size = 3.7, colour = col_f,
             label = sprintf("factor score %s", f2(f_hat))) +
    annotate("text", x = f_true - 0.08, y = top * 1.08, hjust = 1, size = 3.7, colour = "#2c3e50",
             label = sprintf("true quality %s", f2(f_true))) +
    scale_x_continuous(limits = c(-1.5, 3.5), breaks = -1:3) +
    scale_y_continuous(limits = c(0, top * 1.15), expand = c(0, 0)) +
    labs(title = "But some qualities are more plausible",
         subtitle = "how likely each quality is, with the quirks it needs",
         x = "quality of the restaurant", y = "plausibility (density)")
  p1 + p2
}


# =============================================================================
# 6  The states: six indicators of Chapters 2 and 4
# =============================================================================
X6 <- state.x77[, c("Income", "Illiteracy", "Life Exp", "Murder", "HS Grad", "Frost")]
colnames(X6) <- c("Income", "Illiteracy", "LifeExp", "Murder", "HSGrad", "Frost")
st_fa1 <- factanal(X6, factors = 1)
st_fa2 <- factanal(X6, factors = 2)                  # varimax rotation (the default)
st_L <- unclass(st_fa2$loadings)

st_fa2n <- factanal(X6, factors = 2, rotation = "none")   # the two factors as estimated, before varimax

# PCA of the same six indicators (their correlation matrix), with the signs of Chapter 4: income positive on both
st_pca_A <- local({
  e <- eigen(cor(X6), symmetric = TRUE)
  A <- e$vectors[, 1:2] %*% diag(sqrt(e$values[1:2]))   # the loadings: first two columns of A = Gamma Delta
  A <- A %*% diag(sign(A[1, ]))                          # row 1 is income
  dimnames(A) <- list(colnames(X6), c("PC1", "PC2")); A
})

# One panel of loadings: an arrow for each indicator, and the circle of radius one
loading_panel <- function(L, pos, col, title, subtitle, what) {
  d <- tibble(var = rownames(L), x = L[, 1], y = L[, 2]) |> left_join(pos, by = "var")
  circ <- tibble(t = seq(0, 2 * pi, length.out = 200), x = cos(t), y = sin(t))
  ggplot() + axes0 +
    geom_path(data = circ, aes(x, y), colour = "#dddddd", linetype = "22") +
    geom_segment(data = d, aes(0, 0, xend = x, yend = y), colour = col, linewidth = 1,
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    geom_text(data = d, aes(x + dx, y + dy, label = var, hjust = hj, vjust = vj), size = 3.7, colour = "#2c3e50") +
    coord_equal(xlim = c(-1.2, 1.2), ylim = c(-1.2, 1.2)) +
    labs(title = title, subtitle = subtitle, x = paste("loading on", what, "1"), y = paste("loading on", what, "2")) +
    theme(plot.title = element_text(colour = col, face = "bold"))
}

fig_states_pca_fa <- function() {
  # where each label goes, relative to its arrow's tip: offset (dx, dy) and alignment (hj, vj)
  pos_pca <- tribble(~var,        ~dx,   ~dy,   ~hj, ~vj,
                     "Income",     0.03,  0.05,  0,   0,
                     "Illiteracy", 0,    -0.07,  0.5, 1,
                     "LifeExp",    0.05,  0,     0,   0.5,
                     "Murder",     0,     0.06,  0.5, 0,
                     "HSGrad",     0.05,  0.02,  0,   0.5,
                     "Frost",      0.04, -0.05,  0,   1)
  pos_fa  <- tribble(~var,        ~dx,   ~dy,   ~hj, ~vj,
                     "Income",    -0.04,  0.04,  1,   0.5,
                     "Illiteracy", 0,    -0.07,  0.5, 1,
                     "LifeExp",    0.06,  0.03,  0,   0.5,
                     "Murder",     0,    -0.07,  0.5, 1,
                     "HSGrad",     0,     0.07,  0.5, 0,
                     "Frost",      0.05, -0.08,  0,   1)
  loading_panel(st_pca_A, pos_pca, col_third, "PCA", "two components, correlation matrix", "component") +
    loading_panel(st_L, pos_fa, col_f, "Factor analysis", "two factors, turned by varimax", "factor")
}


# =============================================================================
# 7  Factor analysis and PCA
# =============================================================================
# For each state indicator (standardised, variance one): the share of its variance that the first two principal
# components keep, and the share that the two factors explain (its communality), each split into the first and the
# second component or factor (squared loadings; for the factors, after varimax). The totals are the squared lengths
# of its arrows in fig_states_pca_fa().
st_share <- tibble(var = colnames(X6), pca1 = st_pca_A[, 1]^2, pca2 = st_pca_A[, 2]^2,
                   fa1 = st_L[, 1]^2, fa2 = st_L[, 2]^2) |>
  mutate(pca = pca1 + pca2, fa = fa1 + fa2)
# The six tests: the share of each test's variance that two principal components claim, against the factor model
tests_pca_h2 <- local({
  e <- eigen(cor(tests_X), symmetric = TRUE)
  setNames(rowSums((e$vectors[, 1:2] %*% diag(sqrt(e$values[1:2])))^2), tests_items)
})

# Two bars for each indicator, each its whole variance (100%). The coloured part is what the two components keep
# (PCA) or what the two factors explain (factor analysis), dark for the first component or factor and light for the
# second; the grey part is the rest. The number at the end of a bar is the total of its two coloured parts.
fig_fa_pca <- function() {
  light <- function(col) grDevices::colorRampPalette(c(col, "white"))(5)[3]
  fills <- c(pca1 = col_third, pca2 = light(col_third), fa1 = col_f, fa2 = light(col_f), rest = "#e6e6e6")
  labs_m <- c(pca1 = "PCA: first component", pca2 = "PCA: second component",
              fa1 = "factor analysis: first factor (varimax)", fa2 = "factor analysis: second factor (varimax)",
              rest = "the rest of the indicator's variance")
  ord <- st_share$var[order(st_share$fa, decreasing = TRUE)]     # bottom to top: the largest unique variance on top
  pos <- tibble(var = ord, y = seq_along(ord))
  w <- 0.36                                                      # height of each bar; PCA above, factor analysis below
  bars <- st_share |> left_join(pos, by = "var") |>
    tidyr::pivot_longer(c(pca1, pca2, fa1, fa2), names_to = "part", values_to = "share") |>
    mutate(method = substr(part, 1, nchar(part) - 1), yc = ifelse(method == "pca", y + 0.2, y - 0.2)) |>
    group_by(var, method) |> mutate(x0 = cumsum(share) - share, x1 = cumsum(share)) |> ungroup()
  tots <- bars |> group_by(var, method, yc) |> summarise(total = sum(share), .groups = "drop") |>
    mutate(lab = ifelse(total > 0.99, pct(total, 1), pct(total)))
  ggplot() +
    geom_rect(data = distinct(bars, var, method, yc), aes(xmin = 0, xmax = 1, ymin = yc - w / 2, ymax = yc + w / 2,
                                                         fill = "rest")) +
    geom_rect(data = bars, aes(xmin = x0, xmax = x1, ymin = yc - w / 2, ymax = yc + w / 2, fill = part)) +
    geom_text(data = tots, aes(x = total, y = yc, label = lab, colour = method), hjust = -0.15, size = 3.9,
              fontface = "bold", show.legend = FALSE) +
    scale_fill_manual(values = fills, labels = labs_m, breaks = names(labs_m), name = NULL) +
    scale_colour_manual(values = c(pca = col_third, fa = col_f)) +
    scale_x_continuous(labels = scales::percent, breaks = seq(0, 1, 0.25), limits = c(0, 1.1), expand = c(0, 0)) +
    scale_y_continuous(breaks = pos$y, labels = pos$var, expand = c(0.02, 0.02)) +
    guides(fill = guide_legend(ncol = 2, byrow = TRUE)) +
    labs(x = "share of the indicator's variance", y = NULL) +
    theme(legend.position = "top", panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank())
}


# =============================================================================
# 8  Path diagrams
# =============================================================================
# Circles: the common factors (large) and the unique shocks (small). Boxes: the measured variables. The numbers on
# the arrows are the loadings, the numbers below the unique shocks their variances.
sub_digit <- function(i) chartr("0123456789", "\u2080\u2081\u2082\u2083\u2084\u2085\u2086\u2087\u2088\u2089", as.character(i))
fig_path <- function(model = c("critics", "tests", "heywood")) {
  model <- match.arg(model)
  phi <- NA
  if (model == "critics") {
    vars <- tibble(lab = crit_names, x = c(1, 2.45, 3.9, 5.35), psi = crit_psi, idx = 1:4)
    facs <- tibble(lab = "quality", x = mean(vars$x), col = col_f)
    arr  <- tibble(from = 1, to = 1:4, lam = crit_lam)
  } else if (model == "tests") {
    vars <- tibble(lab = tests_items, x = c(1, 2.2, 3.4, 4.9, 6.1, 7.3), psi = 1 - tests_h2, idx = 1:6)
    facs <- tibble(lab = c("verbal", "math"), x = c(2.2, 6.1), col = c(col_verb, col_math))
    arr  <- tibble(from = rep(1:2, each = 3), to = 1:6, lam = c(tests_P[1:3, 1], tests_P[4:6, 2]))
    phi <- tests_Phi[1, 2]
  } else {                                     # the two-factor story behind the Heywood case: critic 1 in the middle
    L <- heywood_two$L
    vars <- tibble(lab = crit_names[c(2, 1, 3)], x = c(1, 2.6, 4.2), psi = heywood_two$psi[c(2, 1, 3)], idx = c(2, 1, 3))
    facs <- tibble(lab = c("food", "service"), x = c(1.8, 3.4), col = c(col_verb, col_math))
    arr  <- tibble(from = c(1, 1, 2, 2), to = c(1, 2, 2, 3), lam = c(L[2, 1], L[1, 1], L[1, 2], L[3, 2]))
    phi <- heywood_two$phi
  }
  yf <- 2.45; yv <- 1; yu <- -0.1; rf <- 0.45; ru <- 0.2; bw <- 0.52; bh <- 0.22
  circle <- function(x0, y0, r, id) tibble(t = seq(0, 2 * pi, length.out = 100), x = x0 + r * cos(t),
                                           y = y0 + r * sin(t), id = id)
  fc <- purrr::map2_dfr(facs$x, seq_len(nrow(facs)), ~ circle(.x, yf, rf, .y))
  uc <- purrr::map2_dfr(vars$x, seq_len(nrow(vars)), ~ circle(.x, yu, ru, .y))
  a <- arr |> mutate(fx = facs$x[from], vx = vars$x[to], dx = vx - fx, dy = (yv + bh) - yf, len = sqrt(dx^2 + dy^2),
                     x0 = fx + rf * dx / len, y0 = yf + rf * dy / len, x1 = vx, y1 = yv + bh,
                     lx = x0 + 0.6 * (x1 - x0), ly = y0 + 0.6 * (y1 - y0), col = facs$col[from])
  p <- ggplot() +
    geom_polygon(data = fc, aes(x, y, group = id, colour = facs$col[id]), fill = "#f6f3fa", linewidth = 1) +
    geom_text(data = facs, aes(x, yf, label = lab, colour = col), size = 4.6, fontface = "bold") +
    geom_segment(data = a, aes(x0, y0, xend = x1, yend = y1, colour = col), linewidth = 0.9,
                 arrow = arrow(length = unit(0.18, "cm"), type = "closed")) +
    geom_label(data = a, aes(lx, ly, label = f2(lam), colour = col), size = 3.9, label.size = 0,
               label.padding = unit(0.12, "lines"), fill = "white") +
    geom_rect(data = vars, aes(xmin = x - bw, xmax = x + bw, ymin = yv - bh, ymax = yv + bh),
              fill = "white", colour = "#2c3e50", linewidth = 0.8) +
    geom_text(data = vars, aes(x, yv, label = lab), size = 3.9, colour = "#2c3e50") +
    geom_segment(data = vars, aes(x, yu + ru, xend = x, yend = yv - bh), colour = "grey45", linewidth = 0.7,
                 arrow = arrow(length = unit(0.15, "cm"), type = "closed")) +
    geom_polygon(data = uc, aes(x, y, group = id), fill = "#f2f2f2", colour = "grey45", linewidth = 0.7) +
    geom_text(data = vars, aes(x, yu, label = paste0("u", sub_digit(idx))), size = 3.8, colour = "grey25") +
    geom_text(data = vars, aes(x, yu - ru - 0.16, label = paste0("\u03c8 = ", f2(psi))), size = 3.6, colour = "grey30") +
    scale_colour_identity() + coord_equal(clip = "off") + theme_void()
  if (!is.na(phi)) {                            # the correlation of two factors, as a curved double arrow
    xs <- facs$x[1] + rf * cos(pi / 5); xe <- facs$x[2] - rf * cos(pi / 5); ys <- yf + rf * sin(pi / 5)
    p <- p + annotate("curve", x = xs, y = ys, xend = xe, yend = ys, curvature = -0.3, colour = "grey35",
                      linewidth = 0.8, arrow = arrow(length = unit(0.16, "cm"), type = "closed", ends = "both")) +
      annotate("label", x = mean(facs$x), y = ys + 0.3 * (xe - xs) / 2 + 0.02, label = f2(phi),
               size = 3.9, colour = "grey25", label.size = 0, fill = "white")
  }
  p
}
