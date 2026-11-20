# =============================================================================
# cfa_sem_figures.R  --  data and figures for the chapter "Confirmatory Factor Analysis and Structural Equation
# Models" of the ISQS 6350 course notes. Put this file and _common.R (shared helpers) in one folder, set it as the
# working directory, and run source("cfa_sem_figures.R"). Every figure is a function, e.g. fig_...().
# Packages: tidyverse, plotly, patchwork, MASS, mvtnorm, lavaan
# =============================================================================

source("_common.R", local = TRUE)
suppressPackageStartupMessages(library(lavaan))

f0 <- function(x) formatC(x, format = "f", digits = 0, big.mark = ",")
f1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
f3 <- function(x) formatC(x, format = "f", digits = 3)
pct <- function(x, d = 0) paste0(formatC(100 * x, format = "f", digits = d), "%")
mns <- function(s) sub("^-", "−", s)                # a proper minus sign for labels in figures
ptxt <- function(p) if (p < 1e-4) "< 0.0001" else if (p < 0.001) paste("=", formatC(p, format = "f", digits = 4)) else paste("=", f3(p))
pv <- function(p) if (p < 1e-4) "below 0.0001" else if (p < 0.001) paste("of", formatC(p, format = "f", digits = 4)) else paste("of", f3(p))
f2z <- function(x) mns(f2(ifelse(abs(x) < 0.005, 0, x)))   # two decimals, without a sign on a rounded zero
sub_digit <- function(i) chartr("0123456789", "₀₁₂₃₄₅₆₇₈₉", as.character(i))

col_verb <- col_a1                                   # verbal ability (blue) and mathematical ability (red), as in Ch 6
col_math <- col_a2
col_speed <- "#2ca02c"                               # the third ability of the 1939 tests
col_dem <- col_a1                                    # the demand curve and the demand shock
col_sup <- col_a2                                    # the supply curve and the supply shock
col_storm <- "#2c3e50"
col_grey <- "#7f8c8d"

# The fit of a lavaan model: chi-square statistic, degrees of freedom, p-value and three fit indices
fit_of <- function(fit) {
  m <- fitMeasures(fit, c("chisq", "df", "pvalue", "rmsea", "cfi", "srmr"))
  setNames(as.numeric(m), names(m))
}
# One row of the parameter table of a lavaan model, e.g. pe_of(fit, "verbal", "~~", "math")
pe_of <- function(fit, lhs, op, rhs, what = "est") {
  pe <- parameterEstimates(fit, standardized = TRUE)
  pe[pe$lhs == lhs & pe$op == op & pe$rhs == rhs, what]
}


# =============================================================================
# 1  Six tests, two abilities (the simulation of Chapter 6)
# =============================================================================
# Three verbal and three mathematics tests. Each test measures one ability only; the abilities have correlation 0.5.
# The same seed as in Chapter 6 gives the same 200 students; here the abilities themselves are kept as well.
tests_items <- c("reading", "vocabulary", "grammar", "arithmetic", "algebra", "geometry")
tests_group <- setNames(rep(c("verbal", "math"), each = 3), tests_items)
tests_P <- rbind(c(0.8, 0), c(0.7, 0), c(0.6, 0), c(0, 0.8), c(0, 0.7), c(0, 0.6))   # loadings (pattern)
tests_Phi <- matrix(c(1, 0.5, 0.5, 1), 2)                                        # correlation of the abilities
tests_h2 <- diag(tests_P %*% tests_Phi %*% t(tests_P))                           # communalities
tests_sim <- local({
  set.seed(4); n <- 200
  F <- mvtnorm::rmvnorm(n, sigma = tests_Phi)
  X <- F %*% t(tests_P) + matrix(rnorm(n * 6), n) %*% diag(sqrt(1 - tests_h2))
  colnames(X) <- tests_items; colnames(F) <- c("verbal", "math")
  list(X = X, F = F)
})
tests_X <- tests_sim$X
tests_fp <- factanal(tests_X, factors = 2, rotation = "promax")     # exploratory, with promax (Chapter 6)

# Confirmatory: each test loads on its own ability only, the abilities may correlate
tests_model <- "
  verbal =~ reading + vocabulary + grammar
  math   =~ arithmetic + algebra + geometry"
tests_cfa <- cfa(tests_model, data = as.data.frame(tests_X), std.lv = TRUE)    # scale: each ability has variance 1
tests_cfa_marker <- cfa(tests_model, data = as.data.frame(tests_X))            # scale: the first loading is 1
tests_pe <- parameterEstimates(tests_cfa, standardized = TRUE)
tests_lam <- setNames(tests_pe$est[tests_pe$op == "=~"], tests_items)          # loadings
tests_psi <- setNames(tests_pe$est[tests_pe$op == "~~" & tests_pe$lhs == tests_pe$rhs & tests_pe$lhs %in% tests_items],
                      tests_items)                                              # unique variances
tests_phi <- pe_of(tests_cfa, "verbal", "~~", "math")                          # correlation of the abilities
tests_phi_ci <- c(pe_of(tests_cfa, "verbal", "~~", "math", "ci.lower"), pe_of(tests_cfa, "verbal", "~~", "math", "ci.upper"))
tests_fit <- fit_of(tests_cfa)
tests_avg_cor <- cor(rowMeans(tests_X[, 1:3]), rowMeans(tests_X[, 4:6]))      # correlation of the two average scores
tests_resid <- residuals(tests_cfa, type = "cor")$cov                         # sample minus model correlations


# =============================================================================
# 2  The critics of Chapter 6: who dines with whom?
# =============================================================================
crit_lam <- c(0.75, 1.00, 0.60, 0.90)                # what one standard deviation of quality adds to each rating
crit_psi <- c(0.4375, 0.25, 0.64, 0.19)              # variance of each critic's quirks
crit_names <- paste("critic", 1:4)
# together = correlation of the quirks of critics 3 and 4 (they dine together and share their bad evenings)
crit_sim <- function(seed = 3, together = 0, n = 200) {
  set.seed(seed)
  f <- rnorm(n); E <- matrix(rnorm(n * 4), n); z <- rnorm(n)
  E[, 3:4] <- sqrt(1 - together) * E[, 3:4] + sqrt(together) * z
  X <- outer(f, crit_lam) + E %*% diag(sqrt(crit_psi)); colnames(X) <- crit_names
  list(X = X, quality = f)
}
crit_X2 <- crit_sim(together = 0.5)$X                # critics 3 and 4 dine together (Chapter 6)
crit_data <- setNames(as.data.frame(crit_X2), paste0("critic", 1:4))
crit_model <- "quality =~ critic1 + critic2 + critic3 + critic4"
crit_fit1  <- cfa(crit_model, crit_data, std.lv = TRUE)                                   # one factor only
crit_fit34 <- cfa(paste(crit_model, "\n  critic3 ~~ critic4"), crit_data, std.lv = TRUE)  # quirks of 3 and 4 shared
crit_fit12 <- suppressWarnings(                                                         # quirks of 1 and 2 shared:
  cfa(paste(crit_model, "\n  critic1 ~~ critic2"), crit_data, std.lv = TRUE))            # critic 4 gets a negative variance
crit_true34 <- 0.5 * sqrt(crit_psi[3] * crit_psi[4]) # the covariance of the quirks of critics 3 and 4 in the simulation


# =============================================================================
# 3  Nine tests from 1939 (Holzinger and Swineford)
# =============================================================================
# 301 seventh- and eighth-grade students of two schools. Three tests each of visual, textual and speed ability.
hs <- lavaan::HolzingerSwineford1939
hs_tests <- c(x1 = "perception", x2 = "cubes", x3 = "lozenges", x4 = "paragraphs", x5 = "sentences",
              x6 = "words", x7 = "addition", x8 = "dots", x9 = "capitals")
hs_model <- "
  visual  =~ x1 + x2 + x3
  textual =~ x4 + x5 + x6
  speed   =~ x7 + x8 + x9"
hs_cfa <- cfa(hs_model, hs, std.lv = TRUE)
hs_fit <- fit_of(hs_cfa)
hs_more <- fitMeasures(hs_cfa, c("rmsea.ci.lower", "rmsea.ci.upper", "tli", "baseline.chisq", "baseline.df"))
hs_mi <- modindices(hs_cfa, sort. = TRUE)
hs_resid <- residuals(hs_cfa, type = "cor")$cov
hs_cfa9 <- cfa(paste(hs_model, "\n  visual =~ x9"), hs, std.lv = TRUE)          # the change the data suggest
hs_fit9 <- fit_of(hs_cfa9)
hs_phi <- lavInspect(hs_cfa, "cor.lv")


# =============================================================================
# 4  A fish market: price and quantity
# =============================================================================
# Each day buyers and fishermen meet at a price (dollars per pound) and a quantity (thousands of pounds).
#   demand: price    =  9 - 0.5 quantity + demand shock           (what buyers pay for a quantity)
#   supply: quantity = -6 + 4 price - 3 storm + supply shock      (what fishermen bring at a price)
# Storms (30% of days) cut the catch by 3,000 pounds at every price; buyers do not care about the weather at sea.
mk_par <- list(d0 = 9, d1 = -0.5, s0 = -6, s1 = 4, storm = -3, sd_d = 0.5, sd_s = 1)
mk_sim <- function(n = 250, sd_d = mk_par$sd_d, sd_s = mk_par$sd_s, seed = 8, p_storm = 0.3) {
  set.seed(seed)
  storm <- rbinom(n, 1, p_storm); ed <- rnorm(n); es <- rnorm(n)
  # the price at which the two curves cross: p = 9 - 0.5 (-6 + 4p - 3 storm + sd_s es) + sd_d ed
  price <- (9 + 3 + 1.5 * storm - 0.5 * sd_s * es + sd_d * ed) / 3
  quantity <- -6 + 4 * price - 3 * storm + sd_s * es
  tibble(quantity, price, storm, ed, es)
}
mk <- mk_sim()
mk_ols <- unname(coef(lm(price ~ quantity, mk))[2])                  # regression of price on quantity
mk_centre <- mk |> group_by(storm) |> summarise(quantity = mean(quantity), price = mean(price), n = n())
mk_iv <- with(mk, cov(price, storm) / cov(quantity, storm))         # the slope through the two centres
mk_model <- "
  price    ~ quantity            # demand
  quantity ~ price + storm       # supply"
mk_sem <- sem(mk_model, mk[, c("quantity", "price", "storm")])
mk_sem_fit <- fit_of(mk_sem)
# The structural form for calm days: each shock moves price and quantity at once (columns of A, per standard deviation)
mk_A <- local({
  B <- matrix(c(0, mk_par$s1, mk_par$d1, 0), 2, byrow = TRUE)       # rows: quantity, price
  solve(diag(2) - B) %*% diag(c(mk_par$sd_s, mk_par$sd_d))        # columns: supply shock, demand shock
})


# =============================================================================
# 5  Abilities and grades: a structural equation model
# =============================================================================
# A year later the 200 students take a first-year statistics course. Their grade (points out of 100) depends on
# their mathematical ability: 6 points per standard deviation, plus luck and effort.
gr_beta_true <- 6
gr <- local({ set.seed(11); data.frame(tests_X, grade = 70 + gr_beta_true * tests_sim$F[, "math"] + rnorm(200, sd = 8)) })
gr_model <- "
  math  =~ arithmetic + algebra + geometry     # measurement
  grade ~  math                                # structure"
gr_sem <- sem(gr_model, gr, std.lv = TRUE)
gr_fit <- fit_of(gr_sem)
gr_beta <- pe_of(gr_sem, "grade", "~", "math"); gr_se <- pe_of(gr_sem, "grade", "~", "math", "se")
gr_lam_avg <- mean(tests_P[4:6, 2])                                  # the average loading of the three math tests
gr_avg <- rowMeans(tests_X[, 4:6]) / gr_lam_avg                      # the average math score, in units of the ability
gr_ols <- unname(coef(lm(gr$grade ~ gr_avg))[2])
gr_ols_true <- unname(coef(lm(gr$grade ~ tests_sim$F[, "math"]))[2]) # with the abilities themselves (never observed)
gr_h_avg <- 1 / (1 + sum(1 - tests_h2[4:6]) / 9 / gr_lam_avg^2)     # share of the average's variance due to ability
# Three models with the same fit: math -> grade, grade -> math, and a correlation only
gr_eq <- list(
  "math → grade" = gr_sem,
  "grade → math" = sem("math =~ arithmetic + algebra + geometry\n  math ~ grade", gr, std.lv = TRUE),
  "math ↔ grade" = sem("math =~ arithmetic + algebra + geometry\n  math ~~ grade", gr, std.lv = TRUE))
gr_eq_fit <- t(sapply(gr_eq, fit_of))


# =============================================================================
# 6  Industrialisation and democracy (Bollen 1989)
# =============================================================================
# 75 developing countries. Industrialisation in 1960 (three indicators), political democracy in 1960 and 1965 (four
# expert ratings each: freedom of the press, freedom of the opposition, fair elections, an effective legislature).
pd <- lavaan::PoliticalDemocracy
pd_model0 <- "
  # measurement
  ind60 =~ x1 + x2 + x3
  dem60 =~ y1 + y2 + y3 + y4
  dem65 =~ y5 + y6 + y7 + y8
  # structure
  dem60 ~ a * ind60
  dem65 ~ c * ind60 + b * dem60
  indirect := a * b
  total    := c + a * b"
pd_model <- paste(pd_model0, "
  # own shocks that may be shared: the same rating in 1960 and 1965, two ratings from the same source
  y1 ~~ y5
  y2 ~~ y4 + y6
  y3 ~~ y7
  y4 ~~ y8
  y6 ~~ y8")
pd_sem  <- sem(pd_model, pd)
pd_sem0 <- sem(pd_model0, pd)                                    # without the shared own shocks
pd_fit  <- fit_of(pd_sem); pd_fit0 <- fit_of(pd_sem0)
pd_pe <- parameterEstimates(pd_sem, standardized = TRUE)
pd_eff <- with(subset(pd_pe, label %in% c("a", "b", "c", "indirect", "total")), setNames(est, label))
pd_se  <- with(subset(pd_pe, label %in% c("a", "b", "c", "indirect", "total")), setNames(se, label))
pd_std <- with(subset(pd_pe, op == "~"), setNames(std.all, label))


# =============================================================================
# 7  Path diagrams
# =============================================================================
# Boxes: measured variables. Large circles: latent variables. Small grey circles: own shocks (the unique shock of a
# measured variable, or the shock of an equation), with their variances. One-headed arrows: effects, with their
# coefficients. Curved two-headed arrows: covariances or correlations. A missing arrow is a zero.
lat <- function(id, x, y, label = id, col = col_f)
  tibble(id = id, type = "lat", x = x, y = y, label = label, col = col, sub = "", w = NA_real_, sdx = 0, sdy = 0)
obs <- function(id, x, y, label = id, w = 0.52)
  tibble(id = id, type = "obs", x = x, y = y, label = label, col = "#2c3e50", sub = "", w = w, sdx = 0, sdy = 0)
pt <- function(id, x, y)                             # an invisible anchor, for arrows that start or end anywhere
  tibble(id = id, type = "pt", x = x, y = y, label = "", col = NA_character_, sub = "", w = NA_real_, sdx = 0, sdy = 0)
err <- function(id, x, y, label, sub = "", sdx = 0, sdy = -0.36)
  tibble(id = id, type = "err", x = x, y = y, label = label, col = "grey45", sub = sub, w = NA_real_, sdx = sdx, sdy = sdy,
         subcol = ifelse(grepl("\u2212", sub), col_a2, "grey30"))          # an impossible negative variance in red
eff <- function(from, to, label = "", col = "#2c3e50", lty = "solid", lpos = 0.6, curv = 0)
  tibble(from = from, to = to, label = label, col = col, kind = "dir", curv = curv, lty = lty, lpos = lpos)
cov2 <- function(from, to, label = "", curv = -0.3, col = "grey35", lty = "solid")
  tibble(from = from, to = to, label = label, col = col, kind = "cov", curv = curv, lty = lty, lpos = 0.5)
# Every measured variable gets its own shock below it, pointing up into its box
own_shocks <- function(vars, y = -0.1, sub = "", idx = seq_len(nrow(vars))) {
  list(nodes = err(paste0("u_", vars$id), vars$x, y, paste0("u", sub_digit(idx)), sub = sub),
       edges = eff(paste0("u_", vars$id), vars$id, col = "grey45"))
}

path_diagram <- function(nodes, edges, r_lat = 0.45, r_err = 0.2, h = 0.22, size = 3.9) {
  N <- split(nodes, nodes$id)
  bpt <- function(n, ux, uy) {                 # where a line from the centre of node n, direction (ux, uy), leaves it
    if (n$type == "lat") return(c(n$x + r_lat * ux, n$y + r_lat * uy))
    if (n$type == "err") return(c(n$x + r_err * ux, n$y + r_err * uy))
    if (n$type == "pt") return(c(n$x, n$y))
    t <- min(if (abs(ux) > 1e-9) n$w / abs(ux) else Inf, if (abs(uy) > 1e-9) h / abs(uy) else Inf)
    c(n$x + t * ux, n$y + t * uy)
  }
  rot <- function(u, a) c(cos(a) * u[1] - sin(a) * u[2], sin(a) * u[1] + cos(a) * u[2])
  geo <- purrr::pmap_dfr(edges, function(from, to, label, col, kind, curv, lty, lpos) {
    a <- N[[from]]; b <- N[[to]]; u <- c(b$x - a$x, b$y - a$y); u <- u / sqrt(sum(u^2))
    ang <- -sign(curv) * atan(2 * abs(curv))               # leave each node along the tangent of the curve
    d0 <- rot(u, ang); d1 <- rot(-u, -ang)
    p0 <- bpt(a, d0[1], d0[2]); p1 <- bpt(b, d1[1], d1[2])
    half <- sqrt(sum((p1 - p0)^2)) / 2
    m <- p0 + lpos * (p1 - p0) - curv * half * c(-u[2], u[1]) * (1 - (2 * lpos - 1)^2)
    tibble(from, to, kind, col, lty, label, curv, x0 = p0[1], y0 = p0[2], x1 = p1[1], y1 = p1[2], lx = m[1], ly = m[2])
  })
  circ <- function(n, r) purrr::pmap_dfr(n, function(id, x, y, ...)
    tibble(id = id, t = seq(0, 2 * pi, length.out = 100), x = x + r * cos(t), y = y + r * sin(t)))
  L <- filter(nodes, type == "lat"); O <- filter(nodes, type == "obs"); E <- filter(nodes, type == "err")
  straight <- filter(geo, kind == "dir", curv == 0); curved <- filter(geo, kind == "cov" | curv != 0)
  p <- ggplot() + geom_blank(data = geo, aes(lx, ly)) + scale_colour_identity() + coord_equal(clip = "off") + theme_void()
  if (nrow(L)) p <- p + geom_polygon(data = left_join(circ(L, r_lat), select(L, id, col), by = "id"),
                                     aes(x, y, group = id, colour = col), fill = "#f6f3fa", linewidth = 1) +
    geom_text(data = L, aes(x, y, label = label, colour = col), size = size + 0.7, fontface = "bold")
  if (nrow(straight)) p <- p + geom_segment(data = straight, aes(x0, y0, xend = x1, yend = y1, colour = col, linetype = lty),
                                            linewidth = 0.85, arrow = arrow(length = unit(0.17, "cm"), type = "closed"))
  for (i in seq_len(nrow(curved))) p <- p +
    annotate("curve", x = curved$x0[i], y = curved$y0[i], xend = curved$x1[i], yend = curved$y1[i],
             curvature = curved$curv[i], colour = curved$col[i], linewidth = 0.8, linetype = curved$lty[i],
             arrow = arrow(length = unit(0.16, "cm"), type = "closed", ends = if (curved$kind[i] == "cov") "both" else "last"))
  lab <- filter(geo, label != "")
  if (nrow(lab)) p <- p + geom_label(data = lab, aes(lx, ly, label = label, colour = col), size = size - 0.1,
                                     label.size = 0, label.padding = unit(0.1, "lines"), fill = "white")
  if (nrow(O)) p <- p + geom_rect(data = O, aes(xmin = x - w, xmax = x + w, ymin = y - h, ymax = y + h),
                                  fill = "white", colour = "#2c3e50", linewidth = 0.8) +
    geom_text(data = O, aes(x, y, label = label), size = size, colour = "#2c3e50")
  if (nrow(E)) p <- p + geom_polygon(data = circ(E, r_err), aes(x, y, group = id), fill = "#f2f2f2", colour = "grey45",
                                     linewidth = 0.7) +
    geom_text(data = E, aes(x, y, label = label), size = size - 0.1, colour = "grey25") +
    geom_text(data = filter(E, sub != ""), aes(x + sdx, y + sdy, label = sub, colour = subcol), size = size - 0.3)
  p + scale_linetype_identity()
}

# The six tests: the model as the analyst writes it down. Every arrow is a free parameter; the dotted zeros are not.
tests_layout <- function() {
  vars <- obs(tests_items, c(1, 2.2, 3.4, 4.9, 6.1, 7.3), 1)
  facs <- lat(c("verbal", "math"), c(2.2, 6.1), 2.45, col = c(col_verb, col_math))
  list(vars = vars, facs = facs)
}
fig_efa_cfa <- function(ncol = 1) {
  l <- tests_layout(); u <- own_shocks(l$vars, sub = paste0("\u03c8", sub_digit(1:6)))
  f_efa <- lat(c("verbal", "math"), c(2.2, 6.1), 2.45, label = c("factor 1", "factor 2"), col = col_f)
  e_efa <- bind_rows(eff(rep(c("verbal", "math"), each = 6), rep(tests_items, 2), col = col_f, lpos = 0.5), u$edges)
  lam <- paste0("\u03bb", sub_digit(1:6))
  e_cfa <- bind_rows(eff("verbal", tests_items[1:3], lam[1:3], col_verb), eff("math", tests_items[4:6], lam[4:6], col_math),
                     u$edges, cov2("verbal", "math", "\u03c6"))
  ttl <- function(t, s) labs(title = t, subtitle = s)
  p1 <- path_diagram(bind_rows(f_efa, l$vars, u$nodes), e_efa, size = 3.7) +
    ttl("Exploratory (Chapter 6)", "every factor may move every test: 12 loadings, named only afterwards")
  p2 <- path_diagram(bind_rows(l$facs, l$vars, u$nodes), e_cfa, size = 3.7) +
    ttl("Confirmatory", "each test measures one ability: 6 loadings, 6 zeros (no arrow), and a correlation")
  (if (ncol == 1) p1 / p2 else p1 | p2) &                       # stacked in the notes, side by side on slides
    theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 11))
}
# The same diagram with the estimates
fig_cfa_est <- function() {
  l <- tests_layout(); u <- own_shocks(l$vars, sub = paste0("ψ = ", f2(tests_psi)))
  e <- bind_rows(eff("verbal", tests_items[1:3], f2(tests_lam[1:3]), col_verb),
                 eff("math", tests_items[4:6], f2(tests_lam[4:6]), col_math),
                 u$edges, cov2("verbal", "math", f2(tests_phi)))
  path_diagram(bind_rows(l$facs, l$vars, u$nodes), e)
}
# The zeros pin the axes. Each test is a point (its loadings as the confirmatory model estimates them, drawn so that
# inner products are covariances, as in Chapter 6), and the two factors are two axes through the origin, with the
# angle between them set by the correlation of the abilities. The points stay; the slider turns both axes together.
# A test's loadings are its coordinates along the axes, each read off parallel to the other axis, so the dashed
# segment from a test, parallel to the other axis, to its own axis is the loading that the model fixes at zero.
cfa_turn_angle <- acos(tests_phi)                                   # angle between the axes, in radians
fig_cfa_turn <- function(angs = seq(-45, 45, by = 1), strip = TRUE, fs = 11, lims = c(-0.4, 1.25)) {
  bv <- pi / 4 - cfa_turn_angle / 2; bm <- pi / 4 + cfa_turn_angle / 2   # axes as estimated, symmetric about 45°
  verbal <- tests_group == "verbal"
  P <- cbind(tests_lam * ifelse(verbal, cos(bv), cos(bm)), tests_lam * ifelse(verbal, sin(bv), sin(bm)))
  cols <- ifelse(verbal, col_verb, col_math)
  pos <- c("bottom right", "top left", "bottom right", "middle right", "middle right", "middle right")
  frames <- purrr::map(angs, function(a) {
    th <- a * pi / 180
    v <- c(cos(bv + th), sin(bv + th)); m <- c(cos(bm + th), sin(bm + th))
    L <- t(solve(cbind(v, m), t(P)))                                 # loadings on the turned axes
    zero <- ifelse(verbal, L[, 2], L[, 1])                           # the six loadings that the model fixes at zero
    dirn <- rbind(m, v)[ifelse(verbal, 1, 2), ]                       # each segment runs parallel to the other axis
    foot <- P - zero * dirn
    ends <- c(-1.4, 1.4)
    list(ax = tibble(x = c(ends * v[1], NA, ends * m[1]), y = c(ends * v[2], NA, ends * m[2]), k = a),
         axv = tibble(x = c(0, 1.08 * v[1]), y = c(0, 1.08 * v[2]), k = a),
         axm = tibble(x = c(0, 1.08 * m[1]), y = c(0, 1.08 * m[2]), k = a),
         lab = tibble(x = 1.16 * c(v[1], m[1]), y = 1.16 * c(v[2], m[2]), txt = c("verbal", "math"),
                      col = c(col_verb, col_math), k = a),
         seg = tibble(x = as.vector(rbind(P[, 1], foot[, 1], NA)), y = as.vector(rbind(P[, 2], foot[, 2], NA)), k = a),
         txt = tibble(x = 0, y = 0, k = a,
                      txt = sprintf("Every pair of axes describes the same six points, and so the same correlations.<br>Largest of the six loadings that the model fixes at zero: %s",
                                    f2(max(abs(zero)))))) })
  g <- function(nm) bind_rows(purrr::map(frames, nm))
  pts <- tibble(x = P[, 1], y = P[, 2], lab = tests_items, col = cols, pos = pos)
  jump <- function(lab, a) list(label = lab, method = "animate",
                                args = list(list(as.character(a)), list(mode = "immediate",
                                            frame = list(duration = 0, redraw = TRUE), transition = list(duration = 0))))
  txt_panel <- plot_ly() |>
    add_trace(data = g("txt"), x = ~x, y = ~y, frame = ~k, text = ~txt, type = "scatter", mode = "text",
              textfont = list(size = 12, color = "#2c3e50"), hoverinfo = "skip", showlegend = FALSE)
  main <- plot_ly() |>
    add_trace(data = g("ax"), x = ~x, y = ~y, frame = ~k, type = "scatter", mode = "lines",
              line = list(color = "#c8c8c8", width = 1), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = g("axv"), x = ~x, y = ~y, frame = ~k, type = "scatter", mode = "lines",
              line = list(color = col_verb, width = 3), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = g("axm"), x = ~x, y = ~y, frame = ~k, type = "scatter", mode = "lines",
              line = list(color = col_math, width = 3), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = g("lab"), x = ~x, y = ~y, frame = ~k, text = ~txt, type = "scatter", mode = "text",
              textfont = list(size = fs + 2, color = ~col), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = g("seg"), x = ~x, y = ~y, frame = ~k, type = "scatter", mode = "lines",
              line = list(color = "#333333", width = 2, dash = "dash"), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = pts, x = ~x, y = ~y, text = ~lab, type = "scatter", mode = "markers+text",
              textposition = ~pos, marker = list(color = ~col, size = 9),
              textfont = list(size = fs, color = ~col), hoverinfo = "skip", showlegend = FALSE)
  plane <- function(anchor) Filter(Negate(is.null), list(title = "", range = lims, zeroline = FALSE,
                                                         showticklabels = FALSE, showgrid = FALSE, scaleanchor = anchor))
  p <- if (strip) {                                       # the readout above the plane, as in Chapter 6
    subplot(txt_panel, main, nrows = 2, heights = c(0.1, 0.9), margin = 0.01, titleX = TRUE, titleY = TRUE) |>
      layout(xaxis = list(visible = FALSE, range = c(-1, 1), fixedrange = TRUE),
             yaxis = list(visible = FALSE, range = c(-1, 1), fixedrange = TRUE),
             xaxis2 = plane(NULL), yaxis2 = plane("x2"), margin = list(t = 50))
  } else {                                                # the plane alone, for the slides
    main |> layout(xaxis = plane(NULL), yaxis = plane("x"), margin = list(t = 40))
  }
  p |>
    animation_opts(frame = 120, transition = 0, redraw = TRUE) |>
    animation_slider(currentvalue = list(prefix = "axes turned by ", suffix = "°")) |>
    animation_button(label = "Play") |>
    finalize_animation(list(jump("as estimated", 0), jump("turned by 20°", 20)),
                       start_frame = "0", prefix = "axes turned by ", suffix = "°")
}
# The four critics with shared quirks: two models with the same fit
fig_crit_shared <- function() {
  one <- function(pair, fit, title) {
    vars <- obs(paste0("critic", 1:4), c(1, 2.45, 3.9, 5.35), 1, crit_names)
    q <- lat("quality", mean(vars$x), 2.45)
    pe <- parameterEstimates(fit)
    lam <- pe$est[pe$op == "=~"]; psi <- pe$est[pe$op == "~~" & pe$lhs == pe$rhs & pe$lhs != "quality"]
    u <- own_shocks(vars, sub = paste0("ψ = ", mns(f2(psi))))
    ps <- pe$est[pe$op == "~~" & pe$lhs == paste0("critic", pair[1]) & pe$rhs == paste0("critic", pair[2])]
    e <- bind_rows(eff("quality", vars$id, f2(lam), col_f), u$edges,
                   cov2(paste0("u_critic", pair[1]), paste0("u_critic", pair[2]), f2(ps), curv = 0.55, col = col_a2))
    n <- bind_rows(q, vars, u$nodes)
    path_diagram(n, e) + labs(title = title, subtitle = sprintf("chi-square %s on %d degree of freedom, p-value %s",
                                                                f2(fitMeasures(fit, "chisq")), fitMeasures(fit, "df"),
                                                                ptxt(fitMeasures(fit, "pvalue"))))
  }
  one(c(3, 4), crit_fit34, "The quirks of critics 3 and 4 are shared") +
    one(c(1, 2), crit_fit12, "The quirks of critics 1 and 2 are shared") &
    theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 11))
}
# The nine tests of 1939: three correlated abilities (standardised estimates)
fig_hs_path <- function() {
  xs <- c(1, 2.15, 3.3, 4.75, 5.9, 7.05, 8.5, 9.65, 10.8)
  vars <- obs(paste0("x", 1:9), xs, 1, hs_tests)
  facs <- lat(c("visual", "textual", "speed"), xs[c(2, 5, 8)], 2.6, col = c(col_verb, col_math, col_speed))
  st <- standardizedSolution(hs_cfa)
  lam <- st$est.std[st$op == "=~"]
  u <- own_shocks(vars, sub = paste0("ψ = ", f2(1 - lam^2)))
  e <- bind_rows(eff(rep(facs$id, each = 3), vars$id, f2(lam), rep(facs$col, each = 3)), u$edges,
                 cov2("visual", "textual", f2(hs_phi[1, 2]), curv = -0.2), cov2("textual", "speed", f2(hs_phi[2, 3]), curv = -0.2),
                 cov2("visual", "speed", f2(hs_phi[1, 3]), curv = -0.3))
  path_diagram(bind_rows(facs, vars, u$nodes), e, size = 3.6)
}
# The fish market as a system of two equations: each curve explains one variable
fig_market_path <- function() {
  n <- bind_rows(obs("storm", 0.2, 1, "storm", w = 0.5), obs("quantity", 2.9, 1, "quantity", w = 0.62),
                 obs("price", 5.7, 1, "price", w = 0.5),
                 err("e1", 2.9, -0.15, "ε₁", sub = "supply shock"), err("e2", 5.7, -0.15, "ε₂", sub = "demand shock"))
  e <- bind_rows(eff("storm", "quantity", mns("-3"), col_sup, lpos = 0.5),
                 eff("price", "quantity", "4", col_sup, curv = 0.4, lpos = 0.5),
                 eff("quantity", "price", mns("-0.5"), col_dem, curv = 0.4, lpos = 0.5),
                 eff("e1", "quantity", col = "grey45"), eff("e2", "price", col = "grey45"))
  path_diagram(n, e) +
    annotate("text", x = 4.3, y = 2.05, label = "supply: quantity responds to price and storms", colour = col_sup, size = 3.6) +
    annotate("text", x = 4.3, y = -0.05, label = "demand: price responds to quantity", colour = col_dem, size = 3.6)
}
# Mathematical ability and the statistics grade: measurement and structure in one model
fig_grade_path <- function(truth = FALSE) {        # truth = TRUE: the values of the simulation instead of the estimates
  vars <- obs(c("arithmetic", "algebra", "geometry"), c(1.6, 2.8, 4.0), 1)
  m <- lat("math", 2.8, 2.45, col = col_math)
  g <- obs("grade", 5.7, 2.45, "grade", w = 0.48)
  if (truth) {
    lab_lam <- f1(tests_P[4:6, 2]); psi <- 1 - tests_h2[4:6]; dg <- 8^2; lab_b <- f0(gr_beta_true)
  } else {
    st <- parameterEstimates(gr_sem)
    lab_lam <- f2(st$est[st$op == "=~"]); psi <- st$est[st$op == "~~" & st$lhs %in% vars$id & st$lhs == st$rhs]
    dg <- st$est[st$op == "~~" & st$lhs == "grade" & st$rhs == "grade"]; lab_b <- f2(gr_beta)
  }
  u <- own_shocks(vars, sub = paste0("ψ = ", f2(psi)), idx = 4:6)    # tests 4 to 6 of the six tests
  eg <- err("e_g", 7.05, 2.45, "ε", sub = paste0("δ² = ", f0(dg)), sdy = -0.36)
  e <- bind_rows(eff("math", vars$id, lab_lam, col_math), u$edges, eff("math", "grade", lab_b, "#2c3e50", lpos = 0.5),
                 eff("e_g", "grade", col = "grey45"))
  path_diagram(bind_rows(m, vars, g, u$nodes, eg), e)
}
# Three models with the same fit: an arrow either way, or a correlation
fig_equivalent <- function() {
  one <- function(k) {
    vars <- obs(c("arithmetic", "algebra", "geometry"), c(1.6, 2.8, 4.0), 1)
    m <- lat("math", 2.8, 2.45, col = col_math); g <- obs("grade", 5.5, 2.45, "grade", w = 0.48)
    u <- own_shocks(vars, idx = 4:6)
    rel <- switch(k, eff("math", "grade", lpos = 0.5), eff("grade", "math", lpos = 0.5), cov2("math", "grade", curv = -0.4))
    path_diagram(bind_rows(m, vars, g, u$nodes), bind_rows(eff("math", vars$id, col = col_math), u$edges, rel), size = 3.5) +
      labs(title = names(gr_eq)[k], subtitle = sprintf("chi-square %s on %d degrees of freedom", f2(gr_eq_fit[k, "chisq"]),
                                                       gr_eq_fit[k, "df"]))
  }
  (one(1) | one(2) | one(3)) & theme(plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
                                     plot.subtitle = element_text(size = 11, hjust = 0.5))
}
# Industrialisation and democracy (standardised estimates). The own shocks of the ratings sit above them, with the
# shared ones linked by curves; the shocks of the two democracy equations sit beside the circles.
fig_pd_path <- function() {
  st <- standardizedSolution(pd_sem)
  g <- function(lhs, op, rhs) st$est.std[st$lhs == lhs & st$op == op & st$rhs == rhs]
  ytop <- obs(paste0("y", 1:8), c(0.8, 1.9, 3.0, 4.1, 5.9, 7.0, 8.1, 9.2), 3.3, paste0("y", 1:8), w = 0.4)
  xbot <- obs(paste0("x", 1:3), c(3.9, 5.0, 6.1), -0.55, paste0("x", 1:3), w = 0.4)
  facs <- lat(c("dem60", "dem65", "ind60"), c(2.45, 7.55, 5.0), c(2.0, 2.0, 0.75),
              label = c("dem60", "dem65", "ind60"), col = c(col_a1, col_f, "#e67e22"))
  ut <- err(paste0("u_", ytop$id), ytop$x, 4.3, "u")
  ub <- err(paste0("u_", xbot$id), xbot$x, -1.55, "u")
  ud <- err(c("d60", "d65"), c(1.1, 8.9), c(1.35, 1.35), c("ε₂", "ε₃"))
  e <- bind_rows(
    eff("dem60", ytop$id[1:4], f2(sapply(ytop$id[1:4], function(v) g("dem60", "=~", v))), col_a1, lpos = 0.55),
    eff("dem65", ytop$id[5:8], f2(sapply(ytop$id[5:8], function(v) g("dem65", "=~", v))), col_f, lpos = 0.55),
    eff("ind60", xbot$id, f2(sapply(xbot$id, function(v) g("ind60", "=~", v))), "#e67e22", lpos = 0.55),
    eff(ut$id, ytop$id, col = "grey45"), eff(ub$id, xbot$id, col = "grey45"), eff(ud$id, c("dem60", "dem65"), col = "grey45"),
    eff("ind60", "dem60", f2(g("dem60", "~", "ind60")), "#2c3e50", lpos = 0.5),
    eff("ind60", "dem65", f2(g("dem65", "~", "ind60")), "#2c3e50", lpos = 0.5),
    eff("dem60", "dem65", f2(g("dem65", "~", "dem60")), "#2c3e50", lpos = 0.5),
    cov2(paste0("u_y", 1:4), paste0("u_y", 5:8), curv = -c(0.16, 0.2, 0.24, 0.28), col = "grey55"),
    cov2(c("u_y2", "u_y6"), c("u_y4", "u_y8"), curv = -0.55, col = "grey55"))
  path_diagram(bind_rows(facs, ytop, xbot, ut, ub, ud), e, size = 3.6)
}
# A legend for path diagrams
fig_path_legend <- function() {
  n <- bind_rows(obs("o", 0.6, 0, "measured", w = 0.55), lat("l", 2.5, 0, "latent", col = col_f),
                 err("e", 4.1, 0, "u"), pt("a1", 5.6, 0), pt("a2", 7.0, 0), pt("c1", 8.4, -0.1), pt("c2", 9.8, -0.1))
  e <- bind_rows(eff("a1", "a2"), cov2("c1", "c2", curv = -0.5))
  path_diagram(n, e, size = 3.7) +
    annotate("text", x = c(0.6, 2.5, 4.1, 6.3, 9.1), y = -0.75, size = 3.5, colour = "grey25", lineheight = 0.95,
             label = c("a variable we\nobserve", "a variable we do\nnot observe", "a variable's\nown shock",
                       "an effect: a\ncoefficient", "a covariance or\na correlation"))
}


# =============================================================================
# 8  Figures with data
# =============================================================================
# The loadings of the six tests: the truth of the simulation, exploratory (promax, Chapter 6), confirmatory
fig_loadings <- function() {
  fp <- unclass(tests_fp$loadings)
  if (sum(fp[1:3, 1]) < sum(fp[1:3, 2])) fp <- fp[, 2:1]          # verbal first
  mk <- function(panel, L, fixed) tibble(panel = panel, test = rep(tests_items, 2), ability = rep(c("verbal", "math"), each = 6),
                                         v = c(L[, 1], L[, 2]), fixed = c(fixed[, 1], fixed[, 2]))
  zero <- tests_P == 0
  cfa_L <- cbind(c(tests_lam[1:3], 0, 0, 0), c(0, 0, 0, tests_lam[4:6]))
  lv <- c("Truth of the simulation", "Exploratory, promax (Chapter 6)", "Confirmatory")
  d <- bind_rows(mk(lv[1], tests_P, zero), mk(lv[2], fp, matrix(FALSE, 6, 2)), mk(lv[3], cfa_L, zero)) |>
    mutate(panel = factor(panel, levels = lv), test = factor(test, levels = rev(tests_items)),
           ability = factor(ability, levels = c("verbal", "math")),
           lab = ifelse(fixed & panel == lv[3], "0 fixed", ifelse(fixed, "0", f2z(v))))
  ggplot(d, aes(ability, test, fill = abs(v))) + geom_tile(colour = "white", linewidth = 1) +
    geom_text(aes(label = lab, colour = fixed), size = 4) +
    scale_fill_gradient(low = "#f7fbff", high = "#6baed6", limits = c(0, 1), guide = "none") +
    scale_colour_manual(values = c("FALSE" = "#2c3e50", "TRUE" = "grey55"), guide = "none") +
    facet_wrap(~ panel) + labs(x = NULL, y = NULL) +
    theme(panel.grid = element_blank(), strip.text = element_text(size = 12, face = "bold"),
          axis.text = element_text(size = 11.5))
}
# Residual correlations: what the model leaves unexplained, sample minus model
fig_resid <- function() {
  one <- function(R, lab, title, sub) {
    p <- nrow(R)
    d <- expand_grid(i = 1:p, j = 1:p) |> filter(j < i) |> mutate(v = R[cbind(i, j)])
    ggplot(d, aes(j, -i, fill = v)) + geom_tile(colour = "white") +
      geom_text(aes(label = f2z(v)), size = 3.3) +
      scale_fill_gradient2(low = col_a1, mid = "white", high = col_a2, limits = c(-0.2, 0.2), oob = scales::squish,
                           name = "sample\nminus model") +
      scale_x_continuous(breaks = 1:(p - 1), labels = lab[1:(p - 1)], position = "top") +
      scale_y_continuous(breaks = -(2:p), labels = lab[2:p]) +
      coord_equal() + labs(x = NULL, y = NULL, title = title, subtitle = sub) +
      theme(panel.grid = element_blank(), axis.text.x = element_text(angle = 45, hjust = 0), legend.position = "right")
  }
  one(tests_resid, tests_items, "Six simulated tests, two abilities", sprintf("chi-square %s on %d df", f1(tests_fit["chisq"]), tests_fit["df"])) +
    one(hs_resid, hs_tests, "Nine real tests from 1939, three abilities", sprintf("chi-square %s on %d df", f1(hs_fit["chisq"]), hs_fit["df"])) +
    plot_layout(guides = "collect", widths = c(1, 1.45)) &
    theme(plot.title = element_text(face = "bold", size = 13))
}
# The chi-square test: where the statistic falls in the distribution it has if the model is right
fig_chisq <- function() {
  one <- function(stat, df, title) {
    x <- seq(0, max(stat * 1.12, qchisq(0.9995, df)), length.out = 500); d <- tibble(x, y = dchisq(x, df))
    p <- pchisq(stat, df, lower.tail = FALSE)
    ggplot(d, aes(x, y)) + geom_area(data = filter(d, x >= stat), fill = col_a2, alpha = 0.3) + geom_line(linewidth = 0.8) +
      geom_vline(xintercept = stat, colour = col_a2, linewidth = 1) +
      annotate("label", x = stat, y = max(d$y) * 0.92, label = paste("statistic", f1(stat)), colour = col_a2, size = 3.8,
               hjust = ifelse(stat > max(x) * 0.6, 1.05, -0.05), label.size = 0) +
      labs(title = title, subtitle = sprintf("chi-square distribution with %d df; p-value %s", df, ptxt(p)),
           x = "chi-square statistic", y = "density")
  }
  (one(tests_fit["chisq"], tests_fit["df"], "Six simulated tests, two abilities") |
     one(hs_fit["chisq"], hs_fit["df"], "Nine real tests from 1939, three abilities")) &
    theme(plot.title = element_text(face = "bold", size = 13))
}
# The fish market: the demand curve, the supply curve on calm and on stormy days
mk_curves <- function(storm_curve = TRUE) {
  l <- list(geom_abline(intercept = mk_par$d0, slope = mk_par$d1, colour = col_dem, linewidth = 1.1),
            geom_abline(intercept = -mk_par$s0 / mk_par$s1, slope = 1 / mk_par$s1, colour = col_sup, linewidth = 1.1))
  if (storm_curve) l <- c(l, geom_abline(intercept = -(mk_par$s0 + mk_par$storm) / mk_par$s1, slope = 1 / mk_par$s1,
                                         colour = col_sup, linewidth = 0.9, linetype = "22"))
  l
}
mk_axes <- list(coord_cartesian(xlim = c(7.4, 12.6), ylim = c(3.35, 5.05), expand = FALSE),
                labs(x = "quantity (thousand pounds)", y = "price (dollars per pound)"))
# Which curve do the days trace? Only demand shifts, only supply shifts, both shift (after Working 1927)
fig_market_working <- function() {
  one <- function(d, title, sub, reg = FALSE) {
    p <- ggplot(d, aes(quantity, price)) + mk_curves(FALSE) + geom_point(colour = col_pt, size = 1.6) + mk_axes +
      annotate("text", x = 7.5, y = 5.0 - 0.06, label = "demand", colour = col_dem, hjust = 0, size = 3.8) +
      annotate("text", x = 12.5, y = 4.62 + 0.08, label = "supply", colour = col_sup, hjust = 1, size = 3.8) +
      labs(title = title, subtitle = sub)
    if (reg) p <- p + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey20", linetype = "42", linewidth = 0.9)
    p
  }
  (one(mk_sim(n = 150, sd_s = 0, p_storm = 0), "Only demand shifts", "the days lie on the supply curve") |
     one(mk_sim(n = 150, sd_d = 0, p_storm = 0), "Only supply shifts", "the days lie on the demand curve") |
     one(mk_sim(n = 150, p_storm = 0), "Both shift", sprintf("regression of price on quantity (dashed): slope %s",
                                                             mns(f2(coef(lm(price ~ quantity, mk_sim(n = 150, p_storm = 0)))[2]))), TRUE)) &
    theme(plot.title = element_text(face = "bold", size = 13), plot.subtitle = element_text(size = 10.5))
}
# Storms shift supply only, so calm and stormy days sit on the same demand curve
fig_market_storms <- function() {
  d <- mk |> mutate(day = factor(ifelse(storm == 1, "stormy day", "calm day"), levels = c("calm day", "stormy day")))
  cen <- mk_centre |> mutate(day = factor(ifelse(storm == 1, "stormy day", "calm day"), levels = levels(d$day)))
  b <- mk_iv; a <- cen$price[1] - b * cen$quantity[1]
  ggplot(d, aes(quantity, price)) + mk_curves() +
    geom_point(aes(colour = day), size = 1.6, alpha = 0.75) +
    geom_abline(intercept = a, slope = b, colour = col_storm, linewidth = 0.9, linetype = "42") +
    geom_point(data = cen, aes(fill = day), shape = 23, size = 5, colour = "black", stroke = 0.9) +
    scale_colour_manual(values = c("calm day" = "#b8c2cc", "stormy day" = col_storm), name = NULL) +
    scale_fill_manual(values = c("calm day" = "#b8c2cc", "stormy day" = col_storm), guide = "none") +
    annotate("text", x = 7.6, y = 5.0 - 0.04, label = "demand", colour = col_dem, hjust = 0, size = 3.9) +
    annotate("text", x = 12.45, y = 4.62 + 0.09, label = "supply,\ncalm days", colour = col_sup, hjust = 1, size = 3.7, lineheight = 0.9) +
    annotate("text", x = 11.3, y = 5.0, label = "supply,\nstormy days", colour = col_sup, hjust = 0, vjust = 1, size = 3.7,
             lineheight = 0.9) +
    mk_axes + labs(subtitle = sprintf("dashed: through the two centres (diamonds), slope %s; the demand curve, %s",
                                      mns(f2(mk_iv)), mns(f2(mk_par$d1)))) +
    theme(legend.position = "top", legend.justification = "left")
}
# Slide from only supply shocks to only demand shocks: the regression line turns, the line through the centres stays
fig_market_slider <- function(shares = seq(0, 1, by = 0.1)) {
  xs <- c(7.4, 12.6)
  fr <- purrr::map(shares, function(w) {
    d <- mk_sim(sd_d = mk_par$sd_d * sqrt(2 * w), sd_s = mk_par$sd_s * sqrt(2 * (1 - w)))
    k <- round(100 * w)
    reg <- coef(lm(price ~ quantity, d)); iv <- with(d, cov(price, storm) / cov(quantity, storm))
    cen <- d |> group_by(storm) |> summarise(quantity = mean(quantity), price = mean(price))
    list(pts = d |> mutate(k = k, day = ifelse(storm == 1, "stormy", "calm")),
         reg = tibble(x = xs, y = reg[1] + reg[2] * xs, k = k),
         ivl = tibble(x = xs, y = cen$price[1] + iv * (xs - cen$quantity[1]), k = k),
         cen = cen |> mutate(k = k),
         txt = tibble(x = 0, y = 0, k = k,
                      txt = sprintf("dotted: regression of price on quantity, slope %s<br>dashed: through the centres of calm (grey) and stormy (dark) days, slope %s",
                                    mns(f2(reg[2])), mns(f2(iv))))) })
  g <- function(nm) bind_rows(purrr::map(fr, nm))
  pts <- g("pts")
  line <- function(i, s) tibble(x = xs, y = i + s * xs)
  jump <- function(lab, s) list(label = lab, method = "animate",
                                args = list(list(as.character(s)), list(mode = "immediate", frame = list(duration = 0, redraw = TRUE),
                                                                        transition = list(duration = 0))))
  txt_panel <- plot_ly() |>
    add_trace(data = g("txt"), x = ~x, y = ~y, frame = ~k, text = ~txt, type = "scatter", mode = "text",
              textfont = list(size = 12, color = "#2c3e50"), hoverinfo = "skip", showlegend = FALSE)
  main <- plot_ly() |>
    add_trace(data = line(mk_par$d0, mk_par$d1), x = ~x, y = ~y, type = "scatter", mode = "lines",
              line = list(color = col_dem, width = 3), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = line(-mk_par$s0 / mk_par$s1, 1 / mk_par$s1), x = ~x, y = ~y, type = "scatter", mode = "lines",
              line = list(color = col_sup, width = 3), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = line(-(mk_par$s0 + mk_par$storm) / mk_par$s1, 1 / mk_par$s1), x = ~x, y = ~y, type = "scatter",
              mode = "lines", line = list(color = col_sup, width = 2, dash = "dash"), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = filter(pts, day == "calm"), x = ~quantity, y = ~price, frame = ~k, type = "scatter", mode = "markers",
              marker = list(color = "#b8c2cc", size = 6), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = filter(pts, day == "stormy"), x = ~quantity, y = ~price, frame = ~k, type = "scatter", mode = "markers",
              marker = list(color = col_storm, size = 6), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = g("reg"), x = ~x, y = ~y, frame = ~k, type = "scatter", mode = "lines",
              line = list(color = "#7f8c8d", width = 2.5, dash = "dot"), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = g("ivl"), x = ~x, y = ~y, frame = ~k, type = "scatter", mode = "lines",
              line = list(color = "black", width = 2, dash = "dash"), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = g("cen"), x = ~quantity, y = ~price, frame = ~k, type = "scatter", mode = "markers",
              marker = list(symbol = "diamond", size = 13, color = "white", line = list(color = "black", width = 2)),
              hoverinfo = "skip", showlegend = FALSE)
  lab <- function(x, y, t, col, anchor) list(x = x, y = y, xref = "x2", yref = "y2", text = t, showarrow = FALSE,
                                             xanchor = anchor, font = list(color = col, size = 13))
  subplot(txt_panel, main, nrows = 2, heights = c(0.1, 0.9), margin = 0.01, titleX = TRUE, titleY = TRUE) |>
    layout(xaxis = list(visible = FALSE, range = c(-1, 1), fixedrange = TRUE),
           yaxis = list(visible = FALSE, range = c(-1, 1), fixedrange = TRUE),
           xaxis2 = list(title = "quantity (thousand pounds)", range = xs, zeroline = FALSE),
           yaxis2 = list(title = "price (dollars per pound)", range = c(3.35, 5.05), zeroline = FALSE),
           margin = list(t = 50),
           annotations = list(lab(7.48, 4.9, "demand", col_dem, "left"), lab(12.55, 4.72, "supply, calm days", col_sup, "right"),
                              lab(11.3, 5.0, "supply, stormy days", col_sup, "left"))) |>
    animation_opts(frame = 300, transition = 0, redraw = TRUE) |>
    animation_slider(currentvalue = list(prefix = "share of the shocks from demand: ", suffix = "%")) |>
    animation_button(label = "Play") |>
    finalize_animation(list(jump("only supply shocks", 0), jump("half and half", 50), jump("only demand shocks", 100)),
                       start_frame = "50", prefix = "share of the shocks from demand: ", suffix = "%")
}
# The structural form of the calm market: each shock is an arrow, a column of A (here two standard deviations long)
fig_market_arrows <- function() {
  d <- mk_sim(p_storm = 0); c0 <- c(10, 4)
  arr <- tibble(lab = c("supply shock", "demand shock"), x = c0[1], y = c0[2],
                xend = c0[1] + 2 * mk_A[1, ], yend = c0[2] + 2 * mk_A[2, ], col = c(col_sup, col_dem))
  ggplot(d, aes(quantity, price)) + mk_curves(FALSE) + geom_point(colour = col_pt, size = 1.6) +
    geom_segment(data = arr, aes(x, y, xend = xend, yend = yend, colour = col), linewidth = 1.6,
                 arrow = arrow(length = unit(0.3, "cm"), type = "closed"), inherit.aes = FALSE) +
    geom_label(data = arr, aes(xend, yend, label = lab, colour = col), size = 4, label.size = 0, fill = "white",
               nudge_x = c(-0.05, 0.25), nudge_y = c(0.09, 0.07), inherit.aes = FALSE) +
    annotate("text", x = 7.5, y = 5.0 - 0.06, label = "demand", colour = col_dem, hjust = 0, size = 3.8) +
    annotate("text", x = 12.5, y = 4.62 + 0.08, label = "supply", colour = col_sup, hjust = 1, size = 3.8) +
    scale_colour_identity() + mk_axes +
    labs(subtitle = "a supply shock moves the market along the demand curve,\na demand shock along the supply curve (arrows: two standard deviations)")
}
# Regression dilution: the noise in the average test score flattens the slope
fig_atten <- function() {
  d <- tibble(math = tests_sim$F[, "math"], avg = gr_avg, grade = gr$grade)
  m <- mean(d$grade)
  one <- function(x, lab, title, sub, ref) {
    ggplot(d, aes(.data[[x]], grade)) + geom_point(colour = col_pt, size = 1.6) +
      geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = "grey20", linewidth = 1) +
      geom_abline(intercept = m - ref * mean(d[[x]]), slope = ref, colour = col_math, linewidth = 1, linetype = "42") +
      coord_cartesian(xlim = c(-3.4, 3.4), ylim = c(42, 98)) + labs(x = lab, y = "statistics grade", title = title, subtitle = sub)
  }
  (one("math", "mathematical ability (standard deviations)", "The ability itself",
       sprintf("regression slope %s points per standard deviation", f1(gr_ols_true)), gr_beta) |
      one("avg", "average math score, in units of the ability", "The average of the three math tests",
          sprintf("regression slope %s: flatter, as the score is noisy", f1(gr_ols)), gr_beta)) +
    plot_annotation(subtitle = sprintf("solid: least squares line · red dashed: the structural equation model's slope, %s", f1(gr_beta))) &
    theme(plot.title = element_text(face = "bold", size = 13))
}
