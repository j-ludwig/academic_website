# =============================================================================
# mds_figures.R  --  data and figures for the chapter "Multidimensional Scaling" of the ISQS 6350
# course notes. Put this file and _common.R (shared helpers) in one folder, set it as the working
# directory, and run source("mds_figures.R"). Every figure is a function, e.g. fig_...().
# Packages: tidyverse, plotly, patchwork, MASS, ggrepel, maps
# =============================================================================

source("_common.R", local = TRUE)
suppressPackageStartupMessages(library(ggrepel))

f0 <- function(x) formatC(x, format = "f", digits = 0, big.mark = ",")
f1 <- function(x) formatC(x, format = "f", digits = 1)
f2 <- function(x, d = 2) formatC(x, format = "f", digits = d)
pct <- function(x, d = 0) paste0(formatC(100 * x, format = "f", digits = d), "%")

# Share of the eigenvalues kept by m dimensions (with negative eigenvalues: of their absolute values or squares)
mds_share <- function(ev, m, type = c("plain", "abs", "sq")) {
  type <- match.arg(type)
  w <- switch(type, plain = ev, abs = abs(ev), sq = ev^2)
  sum(w[seq_len(m)]) / sum(w)
}
# Classical MDS with all eigenvalues: the coordinates and the eigenvalues of B = -1/2 J D^2 J
mds_fit <- function(D, k = 2) {
  m <- cmdscale(as.dist(D), k = k, eig = TRUE)
  list(points = m$points, eig = m$eig)
}
# The 2 x 2 orthogonal matrix R with Y = X R (X, Y centred): the rotation (or reflection) between two maps
orth_between <- function(X, Y) {
  s <- svd(crossprod(scale(X, scale = FALSE), scale(Y, scale = FALSE)))
  s$u %*% t(s$v)
}
axes0 <- list(geom_hline(yintercept = 0, colour = "grey88"), geom_vline(xintercept = 0, colour = "grey88"))
col_wait <- "#e67e22"        # the construction site and everything it changes
col_th <- c(drinks = "#b2182b", flights = "#2166ac", campus = "#1b7837")


# =============================================================================
# 1  The loop road of the Introduction
# =============================================================================
# Eight places on one loop road. Driving distances are measured along the road; a construction site between
# the Stadium and the Church adds a waiting time w to every trip that passes it.
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
speed   <- 30                                          # km/h
T_loop  <- 60 * L_road / speed                         # minutes for one full loop
s_min   <- 60 * places$s / speed
D_wait <- function(w) {                                # travel times (minutes) with a wait of w minutes at the site
  a <- abs(outer(s_min, s_min, "-"))
  D <- pmin(a, T_loop - a + w)
  dimnames(D) <- list(places$place, places$place); D
}
w_star <- 2 * max(abs(outer(s_min, s_min, "-"))) - T_loop   # from this wait on, nobody drives through the site
w_show <- 15                                                # the wait used in the static figures
loop_gap <- local({ D <- D_wait(0); Y <- cmdscale(as.dist(D), k = 2)   # no waiting: how far the ring is off
  e <- abs(dist(Y) - as.dist(D)); c(mean = mean(e), max = max(e), rel = mean(e / as.dist(D))) })

# The figure of the Introduction: the town (left) and the map rebuilt from the travel times (right). Drag the
# slider to add waiting time at the construction site.
fig_loop <- function(height = 660) {
  pair  <- c("Stadium", "Grocery")
  c_frac <- 0.95
  site <- tibble(east  = approx(road$s, road$east,  c_frac * L_road)$y,
                 north = approx(road$s, road$north, c_frac * L_road)$y)
  w_grid <- sort(unique(c(0:36, round(w_star, 1))))
  wlab   <- function(w) sprintf("%.1f", w)
  theta <- 2 * pi * places$frac
  prev  <- cbind(cos(theta), sin(theta)) * T_loop / (2 * pi)
  mds_frames <- purrr::map_dfr(w_grid, function(w) {
    Y <- suppressWarnings(cmdscale(as.dist(D_wait(w)), k = 2))  # classical MDS as it comes
    if (ncol(Y) < 2) Y <- cbind(Y, 0)                           # a line: the second eigenvalue is zero
    Y <- sweep(Y, 2, ifelse(colSums(Y * prev) < 0, -1, 1), "*"); prev <<- Y   # signs are arbitrary: keep them steady
    line <- max(abs(Y[, 2])) < 1e-6                             # a line: alternate the labels above and below it
    tp <- if (line) c("top center", "bottom center")[rank(Y[, 1], ties.method = "first") %% 2 + 1] else
      dplyr::case_when(abs(Y[, 1]) >= abs(Y[, 2]) & Y[, 1] >= 0 ~ "middle right", abs(Y[, 1]) >= abs(Y[, 2]) ~ "middle left",
                       Y[, 2] >= 0 ~ "top center", TRUE ~ "bottom center")   # a ring: labels point outwards
    gap <- mean(abs(dist(Y) - as.dist(D_wait(w))))               # how far the map distances are from the table
    places |> dplyr::select(place) |> mutate(d1 = Y[, 1], d2 = Y[, 2], k = wlab(w), tp = tp, gap = gap)
  })
  fit_txt <- mds_frames |> dplyr::distinct(k, gap) |>
    mutate(x = 0, y = -16, txt = ifelse(gap < 0.005, "average distance error: 0 min (exact)",
                                        sprintf("average distance error: %.1f min", gap)))
  road_mds <- mds_frames |> group_by(k) |>
    dplyr::group_modify(~ dplyr::bind_rows(.x, .x[1, ])) |> ungroup()
  pair_mds <- mds_frames |> dplyr::filter(place %in% pair)
  pr <- places |> dplyr::filter(place %in% pair)

  p_map <- plot_ly() |>
    add_trace(data = road, x = ~east, y = ~north, type = "scatter", mode = "lines", color = I("grey80"),
              line = list(color = "#cccccc", width = 10), hoverinfo = "skip") |>
    add_trace(x = pr$east, y = pr$north, type = "scatter", mode = "lines", color = I(col_a2),
              line = list(color = col_a2, dash = "dash", width = 1.5), hoverinfo = "skip") |>
    add_trace(data = places, x = ~east, y = ~north, type = "scatter", mode = "markers+text",
              text = ~place, color = I("#2c3e50"),
              textposition = dplyr::case_when(places$east > 3.5 ~ "middle right", places$east < -3.5 ~ "middle left",
                                              places$north > 0 ~ "top center", TRUE ~ "bottom center"),
              marker = list(color = "#2c3e50", size = 9), textfont = list(size = 12), hoverinfo = "skip") |>
    add_trace(data = site, x = ~east, y = ~north, type = "scatter", mode = "markers+text",
              text = "construction", textposition = "bottom right", color = I(col_wait),
              marker = list(symbol = "x", color = col_wait, size = 16, line = list(width = 2)),
              textfont = list(color = col_wait, size = 12), hoverinfo = "skip")
  p_mds <- plot_ly() |>
    add_trace(data = road_mds, x = ~d1, y = ~d2, frame = ~k, type = "scatter", mode = "lines",
              color = I("grey80"), line = list(color = "#cccccc", width = 6), hoverinfo = "skip") |>
    add_trace(data = mds_frames, x = ~d1, y = ~d2, frame = ~k, text = ~place, type = "scatter",
              mode = "markers+text", textposition = ~tp, color = I("#2c3e50"),
              marker = list(color = "#2c3e50", size = 9), textfont = list(size = 12),
              hovertemplate = "%{text}<extra></extra>") |>
    add_trace(data = fit_txt, x = ~x, y = ~y, frame = ~k, text = ~txt, type = "scatter", mode = "text",
              textfont = list(size = 12, color = col_wait), hoverinfo = "skip") |>
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

# The travel times themselves, without waiting (left) and with a wait of w_show minutes (right). Outlined in
# orange: the trips that get longer, all of them trips whose shortest way passes the construction site.
fig_loop_tables <- function(w = w_show, title = TRUE) {
  lev <- places$place
  tab <- function(D, ttl) {
    d <- as_tibble(as.table(D), .name_repair = ~ c("from", "to", "min")) |>
      mutate(from = factor(from, levels = rev(lev)), to = factor(to, levels = lev))
    d
  }
  d0 <- tab(D_wait(0)) |> mutate(panel = "No waiting")
  d1 <- tab(D_wait(w)) |> mutate(panel = sprintf("Wait of %d minutes at the site", w))
  d1$changed <- abs(d1$min - d0$min) > 1e-8
  d0$changed <- FALSE
  d <- bind_rows(d0, d1) |> mutate(panel = factor(panel, levels = unique(panel)))
  ggplot(d, aes(to, from)) +
    geom_tile(aes(fill = min), colour = "white", linewidth = 0.6) +
    geom_tile(data = filter(d, changed), fill = NA, colour = col_wait, linewidth = 1.1) +
    geom_text(aes(label = f0(min), colour = min > 18), size = 3.4) +
    scale_colour_manual(values = c(`TRUE` = "white", `FALSE` = "black"), guide = "none") +
    scale_fill_gradient(low = "#f7fbff", high = "#08306b", name = expression(italic(d)[st] ~ "(min)")) +
    facet_wrap(~ panel) + coord_equal() +
    labs(x = expression("column" ~ italic(t)), y = expression("row" ~ italic(s)),
         title = if (title) expression("Two distance matrices" ~ bold(D) * ": the travel time" ~ italic(d)[st] ~
                                         "between places" ~ italic(s) ~ "and" ~ italic(t))) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank(),
          strip.text = element_text(size = 13, face = "bold"), plot.title = element_text(size = 14))
}


# =============================================================================
# 2  A map from distances: the 48 contiguous states
# =============================================================================
# Each state's centre (longitude, latitude; datasets::state.center) is put on a flat map in km: east-west and
# south-north distances from the centre of the country. On this map the distances are Euclidean.
us_keep <- !(state.abb %in% c("AK", "HI"))
us_lon0 <- mean(state.center$x[us_keep]); us_lat0 <- mean(state.center$y[us_keep])
us_km <- function(lon, lat) cbind(east  = (lon - us_lon0) * 111.32 * cos(us_lat0 * pi / 180),
                                  north = (lat - us_lat0) * 110.57)
us48 <- tibble(abb = state.abb[us_keep], name = state.name[us_keep],
               lon = state.center$x[us_keep], lat = state.center$y[us_keep])
us_xy <- us_km(us48$lon, us48$lat)
us_center <- colMeans(us_xy)
us_X <- sweep(us_xy, 2, us_center); rownames(us_X) <- us48$abb  # centred: the mean state is at (0, 0)
us_D <- as.matrix(dist(us_X))                                 # the 48 x 48 distance matrix (km)
us_mds <- mds_fit(us_D, k = 2)
us_R <- orth_between(us_X, us_mds$points)                     # map -> MDS: rotation or reflection
us_turn <- atan2(us_R[1, 2], us_R[1, 1]) * 180 / pi           # angle between the two (degrees)
us_mirror <- det(us_R) < 0
us_err <- max(abs(as.matrix(dist(us_mds$points)) - us_D))     # largest error in any distance (km)

# Borders on the same flat map, centred like the states
us_borders <- local({
  b <- ggplot2::map_data("state")
  xy <- sweep(us_km(b$long, b$lat), 2, us_center)
  tibble(group = b$group, east = xy[, 1], north = xy[, 2])
})
us_outline <- local({                                          # the outline only, thinned (for the interactive figures)
  b <- ggplot2::map_data("usa") |> filter(region == "main")
  b <- b[seq(1, nrow(b), by = 4), ]
  sweep(us_km(b$long, b$lat), 2, us_center)
})

# One map panel: states (abbreviations) and borders, both moved by the orthogonal matrix R
us_panel <- function(R = diag(2), title = NULL, subtitle = NULL, xl = "east (km)", yl = "north (km)",
                     lim = 2900, col = "#2c3e50", ylim_frac = 0.72) {
  P <- us_X %*% R; B <- as.matrix(us_borders[, c("east", "north")]) %*% R
  ggplot() +
    geom_path(data = tibble(x = B[, 1], y = B[, 2], g = us_borders$group), aes(x, y, group = g),
              colour = "grey80", linewidth = 0.3) +
    geom_text(data = tibble(x = P[, 1], y = P[, 2], lab = us48$abb), aes(x, y, label = lab),
              size = 2.5, colour = col, fontface = "bold") +
    coord_equal(xlim = c(-lim, lim), ylim = c(-lim, lim) * ylim_frac) +
    labs(title = title, subtitle = subtitle, x = xl, y = yl) +
    theme(plot.title = element_text(size = 13, face = "bold"), plot.subtitle = element_text(size = 10))
}

# The map (left) and the configuration that classical MDS rebuilds from the 1,128 distances alone (right)
fig_us_recover <- function() {
  us_panel(title = "The map", subtitle = "state centres, km from the middle of the country") +
    us_panel(us_R, title = "Rebuilt from the distances alone", subtitle = "classical MDS of the 48 x 48 table",
             xl = "MDS dimension 1 (km)", yl = "MDS dimension 2 (km)", col = col_a1)
}

# Three answers that fit the distances equally well: as computed, turned by 90 degrees, and its mirror image
fig_us_versions <- function() {
  turn <- function(deg) rotation_matrix(deg * pi / 180)
  sub <- "every distance as in the table"
  one <- function(R, ttl) us_panel(R, title = ttl, subtitle = sub, xl = "dimension 1", yl = "dimension 2",
                                   col = col_a1, lim = 2800, ylim_frac = 1)
  one(us_R, "As computed") + one(us_R %*% turn(90), "Turned by 90°") + one(us_R %*% diag(c(1, -1)), "Mirrored") +
    plot_layout(nrow = 1)
}
# How the computed map sits relative to the real one, in words (signs of eigenvectors can differ by computer)
us_how <- local({
  P <- us_mds$points; up <- P["FL", 2] > P["ME", 2]; flip_ew <- P["WA", 1] > P["ME", 1]
  if (up && !flip_ew) paste("upside down, flipped over from top to bottom: Florida at the top, Washington and Maine",
                            "at the bottom, and the west still on the left")
  else if (up && flip_ew) "upside down, turned by half a circle: Florida at the top and Maine at the bottom left"
  else if (flip_ew) "mirrored: east and west swapped, with Maine on the left"
  else "the usual way up"
})
# What brings the computed map back to the real one
us_back <- local({
  R <- us_R %*% diag(c(1, -1))
  if (us_mirror) sprintf("Flipping the computed map over (right) happens to give back the real map, turned by only %.0f°.",
                         abs(atan2(R[2, 1], R[1, 1]) * 180 / pi))
  else sprintf("Turning the computed map by %.0f° gives back the real map.", abs(us_turn))
})
us_slider <- if (us_mirror) {
  "The computed map is a mirror image of the real one, so no angle brings back the real map: that takes a flip."
} else "At one angle the real map comes back."

# Turn the rebuilt map with the slider: every angle (and every mirror image) has the same 1,128 distances
fig_us_rotate <- function(step = 10) {
  angs <- seq(0, 360 - step, by = step)
  P0 <- us_X %*% us_R; O0 <- us_outline %*% us_R
  fr_pts <- purrr::map_dfr(angs, function(a) { P <- P0 %*% rotation_matrix(a * pi / 180)
    tibble(x = P[, 1], y = P[, 2], lab = us48$abb, k = a) })
  fr_out <- purrr::map_dfr(angs, function(a) { O <- O0 %*% rotation_matrix(a * pi / 180)
    tibble(x = O[, 1], y = O[, 2], k = a) })
  # What changes with the angle: the correlation of the two coordinates and the share of the spread on dimension 1
  fr_txt <- purrr::map_dfr(angs, function(a) { P <- P0 %*% rotation_matrix(a * pi / 180); v <- apply(P, 2, var)
    r <- sub("^-(0\\.00)$", "\\1", sprintf("%.2f", cor(P[, 1], P[, 2])))
    tibble(x = 0, y = 0, k = a,
           txt = sprintf("At every angle: the same %s distances<br>correlation of the coordinates: %s; dimension 1 carries %s of the spread",
                         f0(choose(48, 2)), r, pct(v[1] / sum(v)))) })
  jump <- function(lab, a) list(label = lab, method = "animate",
                                args = list(list(as.character(a)), list(mode = "immediate",
                                            frame = list(duration = 0, redraw = TRUE), transition = list(duration = 0))))
  txt_panel <- plot_ly() |>
    add_trace(data = fr_txt, x = ~x, y = ~y, frame = ~k, text = ~txt, type = "scatter", mode = "text",
              textfont = list(size = 12, color = "#2c3e50"), hoverinfo = "skip", showlegend = FALSE)
  map_panel <- plot_ly() |>
    add_trace(data = fr_out, x = ~x, y = ~y, frame = ~k, type = "scatter", mode = "lines",
              line = list(color = "#bbbbbb", width = 1.5), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = fr_pts, x = ~x, y = ~y, frame = ~k, text = ~lab, type = "scatter", mode = "text",
              textfont = list(size = 11, color = col_a1), hoverinfo = "skip", showlegend = FALSE)
  subplot(txt_panel, map_panel, nrows = 2, heights = c(0.1, 0.9), margin = 0.01, titleX = TRUE, titleY = TRUE) |>
    layout(xaxis = list(visible = FALSE, range = c(-1, 1), fixedrange = TRUE),
           yaxis = list(visible = FALSE, range = c(-1, 1), fixedrange = TRUE),
           xaxis2 = list(title = "dimension 1 (km)", range = c(-3000, 3000), zeroline = TRUE, zerolinecolor = "#e5e5e5"),
           yaxis2 = list(title = "dimension 2 (km)", range = c(-3000, 3000), scaleanchor = "x2", zeroline = TRUE,
                         zerolinecolor = "#e5e5e5"),
           margin = list(t = 50)) |>
    animation_opts(frame = 120, transition = 0, redraw = TRUE) |>
    animation_slider(currentvalue = list(prefix = "turned by ", suffix = "°")) |>
    animation_button(label = "Play") |>
    finalize_animation(list(jump("as computed", 0), jump("turned 90°", 90), jump("turned 180°", 180)),
                       start_frame = "0", prefix = "turned by ", suffix = "°")
}


# In three dimensions a rigid body can turn in three ways: an aeroplane can roll, pitch and yaw. None of the three
# changes a distance between two of its points. Body axes: x to the nose, y to the left wing, z up.
rot2 <- function(P, deg) { a <- deg * pi / 180; P %*% t(rbind(c(cos(a), -sin(a)), c(sin(a), cos(a)))) }
fig_attitude <- function(roll = 30, pitch = 25, yaw = 35) {
  top <- rbind(c(2.3, 0), c(2.0, -0.18), c(0.55, -0.18), c(-0.35, -2.5), c(-0.7, -2.5), c(-0.4, -0.18),
               c(-1.6, -0.15), c(-1.95, -0.95), c(-2.2, -0.95), c(-2.15, -0.1))
  top <- rbind(top, top[nrow(top):2, ] %*% diag(c(1, -1)))          # mirror the right half to the left
  top <- cbind(-top[, 2], top[, 1])                                   # seen from above: nose up, left wing left
  side_body <- rbind(c(2.3, 0.02), c(2.0, 0.17), c(-1.6, 0.17), c(-2.2, 0.08), c(-2.2, -0.06), c(-1.6, -0.14),
                     c(1.9, -0.17), c(2.3, -0.04))                    # seen from the side: nose to the right
  side_fin <- rbind(c(-1.5, 0.15), c(-2.05, 1.05), c(-2.25, 1.05), c(-2.18, 0.1))
  back_wing <- rbind(c(-2.5, 0.04), c(2.5, 0.04), c(2.5, -0.04), c(-2.5, -0.04))   # seen from behind
  back_fin <- rbind(c(-0.06, 0), c(0.06, 0), c(0.04, 1.05), c(-0.04, 1.05))
  back_stab <- rbind(c(-0.95, 0.32), c(0.95, 0.32), c(0.95, 0.27), c(-0.95, 0.27))
  t_c <- seq(0, 2 * pi, length.out = 40); body <- cbind(0.22 * cos(t_c), 0.22 * sin(t_c))
  poly <- function(P, id) tibble(x = P[, 1], y = P[, 2], id = id)
  panel <- function(shapes, deg, dots, title, sub, arc) {
    ghost  <- bind_rows(lapply(seq_along(shapes), function(i) poly(shapes[[i]], i)))
    turned <- bind_rows(lapply(seq_along(shapes), function(i) poly(rot2(shapes[[i]], deg), i)))
    D1 <- rot2(as.matrix(dots[, c("x", "y")]), deg)
    ggplot() +
      geom_polygon(data = ghost, aes(x, y, group = id), fill = NA, colour = "#bbbbbb", linetype = "22",
                   linewidth = 0.5) +
      geom_polygon(data = turned, aes(x, y, group = id), fill = "#2c3e50", colour = "#2c3e50", linewidth = 0.4) +
      geom_point(data = tibble(x = D1[, 1], y = D1[, 2], k = dots$k), aes(x, y, colour = k), size = 2.8,
                 show.legend = FALSE) +
      annotate("curve", x = arc[1], y = arc[2], xend = arc[3], yend = arc[4], curvature = arc[5],
               arrow = arrow(length = unit(0.18, "cm"), type = "closed"), colour = col_a2, linewidth = 0.7) +
      scale_colour_manual(values = c(left = "#d62728", right = "#2ca02c", nose = col_wait)) +
      coord_equal(xlim = c(-2.9, 2.9), ylim = c(-2.6, 2.6)) + theme_void() +
      labs(title = title, subtitle = sub) +
      theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
            plot.subtitle = element_text(size = 10.5, hjust = 0.5, colour = "grey30"))
  }
  p_roll  <- panel(list(back_wing, back_fin, back_stab, body), roll,
                   tibble(x = c(-2.5, 2.5), y = 0, k = c("left", "right")),
                   "Roll", "tilts the wings left or right\n(seen from behind)", c(2.6, 0.9, 2.0, 1.9, 0.3))
  p_pitch <- panel(list(side_body, side_fin), pitch, tibble(x = 2.3, y = 0, k = "nose"),
                   "Pitch", "points the nose up or down\n(seen from the side)", c(2.6, 0.5, 2.0, 1.55, 0.3))
  p_yaw   <- panel(list(top), yaw, tibble(x = 0, y = 2.3, k = "nose"),
                   "Yaw", "points the nose left or right\n(seen from above)", c(0.5, 2.55, -0.8, 2.35, 0.3))
  p_roll + p_pitch + p_yaw + plot_layout(nrow = 1)
}


# -----------------------------------------------------------------------------
# A small example for the section on classical MDS: three centred points
# -----------------------------------------------------------------------------
ex3_X <- rbind(`1` = c(-4, -1), `2` = c(2, -1), `3` = c(2, 2))   # their mean, the centre, is (0, 0)
ex3_D <- as.matrix(dist(ex3_X))                                 # distance matrix: 6, sqrt(45) and 3
ex3_B <- ex3_X %*% t(ex3_X)                                     # inner products, B = X X'

fig_inner <- function() {
  P <- tibble(id = rownames(ex3_X), x = ex3_X[, 1], y = ex3_X[, 2], b = diag(ex3_B))
  P$lab <- sprintf("italic(b)[%s%s] == %d", P$id, P$id, P$b)
  pr <- tibble(s = c(1, 2, 1), t = c(2, 3, 3)) |>
    mutate(x = ex3_X[s, 1], y = ex3_X[s, 2], xend = ex3_X[t, 1], yend = ex3_X[t, 2],
           d2 = round(ex3_D[cbind(s, t)]^2),
           lab = ifelse(sqrt(d2) == round(sqrt(d2)), sprintf("italic(d)[%d%d] == %d", s, t, round(sqrt(d2))),
                        sprintf("italic(d)[%d%d] == sqrt(%d)", s, t, d2)),
           lx = c(-1, 2.3, -1.25), ly = c(-1.38, 0.5, 0.95), hj = c(0.5, 0, 0.5))
  ggplot() + axes0 +
    geom_segment(data = pr, aes(x, y, xend = xend, yend = yend), linetype = "dashed", colour = "grey55") +
    geom_text(data = pr, aes(lx, ly, label = lab, hjust = hj), parse = TRUE, size = 4.2, colour = "grey30") +
    geom_segment(data = P, aes(0, 0, xend = x, yend = y), colour = col_a1, linewidth = 1,
                 arrow = arrow(length = unit(0.22, "cm"), type = "closed")) +
    geom_point(data = P, aes(x, y), size = 3, colour = "#2c3e50") +
    geom_text(data = P, aes(x, y, label = id), nudge_x = c(-0.3, 0.3, 0.3), nudge_y = c(0, -0.15, 0.15),
              size = 5, fontface = "bold", colour = "#2c3e50") +
    geom_text(data = P, aes(x + c(-0.3, 0.3, 0.3), y + c(-0.42, -0.55, 0.55), label = lab), parse = TRUE,
              size = 4, colour = col_a1, hjust = c(0.8, 0, 0)) +
    annotate("point", x = 0, y = 0, shape = 3, size = 4, stroke = 1.2) +
    annotate("text", x = -0.15, y = 0.32, label = "centre", hjust = 1, size = 4, colour = "grey30") +
    scale_x_continuous(breaks = -5:3) + scale_y_continuous(breaks = -2:3) +
    coord_equal(xlim = c(-5.2, 3.9), ylim = c(-1.8, 2.9)) +
    labs(x = "coordinate 1", y = "coordinate 2",
         subtitle = expression("arrows from the centre, with squared lengths" ~ italic(b)[tt] *
                                 "; dashed: distances" ~ italic(d)[st]))
}


# =============================================================================
# 3  Two steps or one: from distances to latent coordinates
# =============================================================================
# Step 1 picks an orientation for the map (the slider: any angle fits the distances). Step 2 runs PCA on those
# coordinates. The scores (right) do not move: PCA undoes whatever orientation step 1 chose.
# The roommates of Chapter 4: the covariance matrix of the two GPAs, and the standard deviations of its two principal
# components. On the map that MDS returns, the long axis of the cloud is horizontal. Random pairing gives the roommates
# no order, so swapping them must leave the cloud unchanged: their GPA axes lie at 45 degrees to the axes of the cloud,
# and the symmetric root, with its peer effect, follows from the two standard deviations alone.
room_S <- matrix(c(0.26, 0.10, 0.10, 0.26), 2, 2)
room_sd <- sqrt(eigen(room_S)$values)                          # standard deviations of the two principal components
room_peer <- (room_sd[1] - room_sd[2]) / (room_sd[1] + room_sd[2])   # the peer effect, from the map alone

# 200 pairs of roommates, with covariance matrix exactly room_S, drawn on the MDS map (the principal component scores)
room_pts <- local({ set.seed(6350); Z <- scale(matrix(rnorm(400), 200), scale = FALSE)
  Z %*% solve(chol(cov(Z))) %*% chol(room_S) })
room_map <- room_pts %*% (cbind(c(1, 1), c(1, -1)) / sqrt(2))  # long axis (both GPAs up) first, short axis second
fig_room_axes <- function() {
  E <- t(ellipse_xy(diag(room_sd), r = 2))
  ax <- tibble(xend = 1.25 / sqrt(2) * c(1, 1), yend = 1.25 / sqrt(2) * c(1, -1),
               lab = c("roommate 1's GPA", "roommate 2's GPA"), vj = c(0, 1))
  arc <- function(a0, a1, r = 0.42) { t <- seq(a0, a1, length.out = 30); tibble(x = r * cos(t), y = r * sin(t)) }
  ggplot() + axes0 +
    geom_point(data = tibble(x = room_map[, 1], y = room_map[, 2]), aes(x, y), colour = "grey60", alpha = 0.55,
               size = 1.3) +
    geom_path(data = tibble(x = E[, 1], y = E[, 2]), aes(x, y), colour = "#2c3e50", linewidth = 0.8) +
    geom_path(data = arc(0, pi / 4), aes(x, y), colour = col_a2, linewidth = 0.6) +
    geom_path(data = arc(-pi / 4, 0), aes(x, y), colour = col_a2, linewidth = 0.6) +
    annotate("text", x = 0.6 * cos(c(pi / 8, -pi / 8)), y = 0.6 * sin(c(pi / 8, -pi / 8)), label = "45°", size = 3.8,
             colour = col_a2) +
    geom_segment(data = ax, aes(0, 0, xend = xend, yend = yend), colour = col_a2, linewidth = 1.1,
                 arrow = arrow(length = unit(0.22, "cm"), type = "closed")) +
    geom_text(data = ax, aes(xend + 0.04, yend, label = lab, vjust = vj), hjust = 0, size = 4.2, colour = col_a2,
              fontface = "bold") +
    annotate("text", x = 2.15, y = 0.04, label = "long axis:\nboth GPAs rise", hjust = 1, vjust = 0, size = 3.6,
             colour = "grey35", lineheight = 0.9) +
    coord_equal(xlim = c(-1.65, 2.15), ylim = c(-1.15, 1.15)) +
    labs(x = "first dimension of the MDS map", y = "second dimension",
         title = "The MDS map of 200 pairs of roommates",
         subtitle = "No order: the two GPA axes lie at 45° to the axes of the cloud") +
    theme(plot.title = element_text(size = 13, face = "bold"), plot.subtitle = element_text(size = 12),
          axis.text = element_blank(), panel.grid = element_blank())
}

us_scores <- function(X) {                                     # PCA scores, signs fixed (Maine positive on both)
  p <- prcomp(X); Fs <- p$x
  sg <- sign(Fs["ME", ]); sweep(Fs, 2, sg, "*")
}
us_pca_R <- local({ p <- prcomp(us_X); G <- p$rotation; sg <- sign((us_X %*% G)["ME", ]); sweep(G, 2, sg, "*") })

fig_two_steps <- function(step = 15, titles = c("Step 1: coordinates, in some orientation", "Step 2: their PCA scores"),
                          slider = "orientation chosen in step 1: ") {
  angs <- seq(0, 360 - step, by = step)
  fr <- purrr::map(angs, function(a) {
    Q <- rotation_matrix(a * pi / 180)
    Xa <- us_X %*% Q; Fa <- us_scores(Xa)
    Ga <- solve(crossprod(Xa), crossprod(Xa, Fa))             # the PCA rotation for these coordinates
    Oa <- us_outline %*% Q; OF <- Oa %*% Ga
    list(pts = tibble(x = c(Xa[, 1], Fa[, 1]), y = c(Xa[, 2], Fa[, 2]), lab = rep(us48$abb, 2),
                      panel = rep(c("map", "scores"), each = nrow(Xa)), k = a),
         out = tibble(x = c(Oa[, 1], OF[, 1]), y = c(Oa[, 2], OF[, 2]),
                      panel = rep(c("map", "scores"), each = nrow(Oa)), k = a))
  })
  pts <- bind_rows(purrr::map(fr, "pts")); out <- bind_rows(purrr::map(fr, "out"))
  one <- function(pn, col) plot_ly() |>
    add_trace(data = filter(out, panel == pn), x = ~x, y = ~y, frame = ~k, type = "scatter", mode = "lines",
              line = list(color = "#bbbbbb", width = 1.2), hoverinfo = "skip", showlegend = FALSE) |>
    add_trace(data = filter(pts, panel == pn), x = ~x, y = ~y, frame = ~k, text = ~lab, type = "scatter",
              mode = "text", textfont = list(size = 9, color = col), hoverinfo = "skip", showlegend = FALSE)
  ax <- function(t) list(title = t, range = c(-3000, 3000), zeroline = TRUE, zerolinecolor = "#e5e5e5")
  subplot(one("map", "#2c3e50"), one("scores", col_a1), margin = 0.06, titleX = TRUE, titleY = TRUE) |>
    layout(xaxis = ax("coordinate 1 (km)"), yaxis = c(ax("coordinate 2 (km)"), scaleanchor = "x"),
           xaxis2 = ax("PCA score 1 (km)"), yaxis2 = c(ax("PCA score 2 (km)"), scaleanchor = "x2"),
           margin = list(t = 80),
           annotations = list(
             list(x = 0.235, y = 1.02, xref = "paper", yref = "paper", showarrow = FALSE, xanchor = "center",
                  yanchor = "bottom", font = list(size = 14), text = titles[1]),
             list(x = 0.765, y = 1.02, xref = "paper", yref = "paper", showarrow = FALSE, xanchor = "center",
                  yanchor = "bottom", font = list(size = 14), text = titles[2]))) |>
    animation_opts(frame = 150, transition = 0, redraw = TRUE) |>
    animation_slider(currentvalue = list(prefix = slider, suffix = "°")) |>
    animation_button(label = "Play") |>
    finalize_animation(list(), start_frame = "0", prefix = slider, suffix = "°")
}


# =============================================================================
# 4  Distances learned from text: a small embedding
# =============================================================================
# A made-up text written with 18 words in three themes taken from the course's examples (think of a real text with
# all other words removed). It comes in paragraphs on one theme each; within a theme some pairs of words are used
# together often, and students and pilots also drink coffee.
word_themes <- list(
  drinks  = c("coke", "pepsi", "sugar", "coffee", "tea", "cup"),
  flights = c("headwind", "fuel", "pilot", "plane", "airport", "runway"),
  campus  = c("student", "roommate", "exam", "grade", "library", "class"))
word_pairs <- list(
  drinks  = list(c("coke", "pepsi"), c("coffee", "tea"), c("cup", "coffee"), c("sugar", "coke"), c("cup", "tea")),
  flights = list(c("headwind", "fuel"), c("pilot", "plane"), c("airport", "runway"), c("fuel", "plane"),
                 c("pilot", "airport")),
  campus  = list(c("student", "roommate"), c("exam", "grade"), c("library", "class"), c("student", "exam"),
                 c("roommate", "library")))
word_vocab <- unlist(word_themes, use.names = FALSE)
word_theme <- setNames(rep(names(word_themes), lengths(word_themes)), word_vocab)
word_text <- local({
  set.seed(6350)
  sentence <- function(th) {                                  # a frequent pair, sometimes with a third word
    w <- sample(word_pairs[[th]], 1)[[1]]
    if (runif(1) < 0.5) w <- c(w, sample(setdiff(word_themes[[th]], w), 1))
    sample(w)
  }
  paragraph <- function() {                                   # 4 to 7 sentences on one theme
    th <- sample(names(word_themes), 1)
    s <- replicate(sample(4:7, 1), sentence(th), simplify = FALSE)
    if (th == "campus" && runif(1) < 0.35)
      s <- append(s, list(c("student", "coffee", sample(c("exam", "library"), 1))), sample(0:length(s), 1))
    if (th == "flights" && runif(1) < 0.25) s <- append(s, list(c("pilot", "coffee", "airport")), sample(0:length(s), 1))
    unlist(s)
  }
  unlist(replicate(300, paragraph(), simplify = FALSE))
})
# The distance between two words: wherever word i appears, count the words to the nearest word j and average over
# all those places; then take the mean of the two directions (i to j and j to i)
word_near <- function(a, b) {
  k <- findInterval(a, b); lo <- b[pmax(k, 1)]; hi <- b[pmin(k + 1, length(b))]
  mean(pmin(abs(a - lo), abs(hi - a)))
}
word_pos <- split(seq_along(word_text), factor(word_text, levels = word_vocab))
word_M <- outer(word_vocab, word_vocab,
                Vectorize(function(i, j) if (i == j) 0 else word_near(word_pos[[i]], word_pos[[j]])))
dimnames(word_M) <- list(word_vocab, word_vocab)
word_D <- (word_M + t(word_M)) / 2                            # average distance in the text, in words
word_mds <- mds_fit(word_D, k = 2)
word_Y <- word_mds$points
if (word_Y["coke", 2] > 0) word_Y[, 2] <- -word_Y[, 2]        # orientation is arbitrary: put the drinks at the bottom
if (word_Y["plane", 1] < 0) word_Y[, 1] <- -word_Y[, 1]       # and the flights on the right

# The word distances as a table (heat map), words ordered by theme
fig_word_D <- function(title = TRUE) {
  lev <- word_vocab
  d <- as_tibble(as.table(word_D), .name_repair = ~ c("w1", "w2", "d")) |>
    mutate(w1 = factor(w1, levels = rev(lev)), w2 = factor(w2, levels = lev))
  ggplot(d, aes(w2, w1, fill = d)) + geom_tile(colour = "white", linewidth = 0.4) +
    scale_fill_gradient(low = "#08306b", high = "#f7fbff", name = "words\napart") +
    coord_equal() + labs(x = NULL, y = NULL, title = if (title) expression("Distance matrix" ~ bold(D) * ": the average distance" ~
                                                               italic(d)[ij] ~ "in the text, in words")) +
    theme(axis.text.x = element_text(angle = 60, hjust = 1, colour = col_th[word_theme[lev]]),
          axis.text.y = element_text(colour = rev(col_th[word_theme[lev]])), panel.grid = element_blank())
}

# The embedding: classical MDS of the word distances. Colour shows the theme, which MDS was not told.
fig_word_map <- function() {
  d <- tibble(word = word_vocab, x = word_Y[, 1], y = word_Y[, 2], theme = word_theme[word_vocab])
  ggplot(d, aes(x, y, colour = theme)) + axes0 +
    geom_point(size = 2.4) +
    geom_text_repel(aes(label = word), size = 4.2, seed = 3, show.legend = FALSE, max.overlaps = 30,
                    box.padding = 0.35) +
    scale_colour_manual(values = col_th, name = "theme\n(not used\nby MDS)") +
    coord_equal() + labs(x = "embedding dimension 1", y = "embedding dimension 2")
}


# =============================================================================
# 5  Other distances
# =============================================================================
# The six state indicators of Chapters 2 and 4, standardised
X6 <- state.x77[, c("Income", "Illiteracy", "Life Exp", "Murder", "HS Grad", "Frost")]
colnames(X6) <- c("Income", "Illiteracy", "LifeExp", "Murder", "HSGrad", "Frost")
Z6 <- scale(X6)
st_E <- Z6 %*% solve(sym_root(cor(X6)))                       # whitened: Euclidean distance = Mahalanobis distance
st_D <- list(Euclidean = as.matrix(dist(Z6)), Manhattan = as.matrix(dist(Z6, method = "manhattan")),
             Mahalanobis = as.matrix(dist(st_E)))
st_mds <- lapply(st_D, mds_fit, k = 2)
st_pca <- prcomp(Z6)
st_prof <- st_D$Euclidean[upper.tri(st_D$Euclidean)]           # all 1,225 distances between the states' profiles
st_ex <- c("Montana", "New Hampshire", "California", "Nevada")  # far apart but alike; neighbours but different
st_pca_err <- max(abs(abs(st_mds$Euclidean$points) - abs(st_pca$x[, 1:2])))   # MDS = PCA scores

# The states mapped from three distance tables (Manhattan and Mahalanobis maps turned to match the first)
fig_states_maps <- function() {
  ref <- st_mds$Euclidean$points
  if (ref["Alaska", 2] < 0) ref[, 2] <- -ref[, 2]
  ttl <- c(Euclidean = "Euclidean", Manhattan = "Manhattan (city block)", Mahalanobis = "Mahalanobis")
  sub <- c(Euclidean = "= the PCA scores of Chapter 4", Manhattan = "sum of absolute differences",
           Mahalanobis = "relative to the shape of the cloud")
  panels <- lapply(names(st_D), function(nm) {
    Y <- if (nm == "Euclidean") ref else procrustes_align(st_mds[[nm]]$points, ref)
    d <- tibble(x = Y[, 1], y = Y[, 2], lab = state.abb)
    ggplot(d, aes(x, y, label = lab)) + axes0 + geom_text(size = 2.6, colour = "#2c3e50") +
      coord_equal() + labs(title = ttl[[nm]], subtitle = sub[[nm]], x = "dimension 1", y = "dimension 2") +
      theme(plot.title = element_text(size = 13, face = "bold"), plot.subtitle = element_text(size = 10))
  })
  wrap_plots(panels, nrow = 1)
}

# Great-circle ("as the crow flies") distances between the state centres, on the curved Earth
haversine_km <- function(lon, lat, r = 6371) {
  ph <- lat * pi / 180; la <- lon * pi / 180
  dph <- outer(ph, ph, "-"); dla <- outer(la, la, "-")
  a <- sin(dph / 2)^2 + outer(cos(ph), cos(ph)) * sin(dla / 2)^2
  2 * r * asin(pmin(sqrt(a), 1))
}
us_G <- haversine_km(us48$lon, us48$lat); dimnames(us_G) <- list(us48$abb, us48$abb)
us_mds_G <- mds_fit(us_G, k = 3)
us_bulge <- us_mds_G$points[, 3] * sign(us_mds_G$points[which.min(rowSums(us_X^2)), 3])  # dimension 3, middle up

# The map from great-circle distances (aligned with the flat map), coloured by the third MDS dimension
fig_us_sphere <- function(title = TRUE) {
  Y <- us_mds_G$points
  Y2 <- procrustes_align(Y[, 1:2], us_X)
  d <- tibble(x = Y2[, 1], y = Y2[, 2], z = us_bulge, lab = us48$abb)
  ggplot() +
    geom_path(data = us_borders, aes(east, north, group = group), colour = "grey85", linewidth = 0.3) +
    geom_label(data = d, aes(x, y, label = lab, fill = z), size = 2.5, label.size = 0, colour = "white",
               fontface = "bold", label.padding = unit(0.12, "lines")) +
    scale_fill_gradient2(low = "#2166ac", mid = "#bbbbbb", high = "#b2182b",
                         name = "dimension 3 (km)\nred: above the plane\nblue: below it") +
    coord_equal(xlim = c(-2900, 2900), ylim = c(-2100, 2100)) +
    labs(x = "dimension 1 (km)", y = "dimension 2 (km)",
         title = if (title) "Classical MDS of the great-circle distances between the state centres",
         subtitle = if (title) "position: dimensions 1 and 2 (the familiar map); colour: dimension 3, the bulge of the Earth") +
    theme(plot.title = element_text(size = 13, face = "bold"), plot.subtitle = element_text(size = 10.5))
}

# Four places evenly spaced on a loop: neighbours 1 minute apart, opposite places 2 minutes. A square with sides 1
# has diagonals of only 1.41; classical MDS gives diagonals of 2 but sides of 1.41; and no arrangement in any number
# of dimensions does better, since B and D would both have to sit halfway along the straight line from A to C.
fig_loop4 <- function(nrow = 2) {
  ttl <- theme(plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
               plot.subtitle = element_text(size = 10.5, hjust = 0.5, colour = "grey30", lineheight = 1.05))
  ang <- c(A = 180, B = 90, C = 0, D = 270) * pi / 180
  mid <- c(135, 45, 315, 225) * pi / 180                     # directions of the side midpoints A-B, B-C, C-D, D-A
  lim <- coord_equal(xlim = c(-1.45, 1.45), ylim = c(-1.45, 1.45), clip = "off")
  corners <- function(P) list(
    geom_point(data = P, aes(x, y), size = 3, colour = "#2c3e50"),
    geom_text(data = P, aes(x * (1 + 0.28 / sqrt(x^2 + y^2)), y * (1 + 0.28 / sqrt(x^2 + y^2)), label = lab),
              size = 5, fontface = "bold", colour = "#2c3e50"))
  # 1: the table, distances along the loop
  r <- 0.85
  P1 <- tibble(lab = names(ang), x = r * cos(ang), y = r * sin(ang))
  circ <- tibble(t = seq(0, 2 * pi, length.out = 300), x = r * cos(t), y = r * sin(t))
  p1 <- ggplot() + geom_path(data = circ, aes(x, y), colour = "#cccccc", linewidth = 3) + corners(P1) +
    annotate("text", x = 1.3 * r * cos(mid), y = 1.3 * r * sin(mid), label = "1 min", size = 3.8, colour = "grey35") +
    lim + theme_void() + ttl +
    labs(title = "The table", subtitle = "along the loop: neighbours 1 minute apart,\nopposite places 2 minutes")
  # 2 and 3: squares, with the distances that are wrong in red
  square <- function(a, side, diag, side_ok, title, subtitle) {
    P <- tibble(lab = names(ang), x = a * cos(ang), y = a * sin(ang))
    col <- function(ok) if (ok) "grey35" else col_a2
    ggplot() +
      annotate("segment", x = P$x, y = P$y, xend = P$x[c(2, 3, 4, 1)], yend = P$y[c(2, 3, 4, 1)],
               colour = col(side_ok), linewidth = 0.9) +
      annotate("segment", x = P$x[1:2], y = P$y[1:2], xend = P$x[3:4], yend = P$y[3:4],
               colour = col(!side_ok), linewidth = 0.9, linetype = "22") +
      annotate("text", x = (a / sqrt(2) + 0.16) * cos(mid), y = (a / sqrt(2) + 0.16) * sin(mid),
               label = side, size = 4.2, colour = col(side_ok), fontface = if (side_ok) "plain" else "bold") +
      annotate("text", x = c(0.5 * a, 0.09), y = c(0.09, -0.5 * a), label = diag, hjust = c(0.5, 0), vjust = c(0, 0.5),
               size = 4.2, colour = col(!side_ok), fontface = if (side_ok) "bold" else "plain") +
      corners(P) + lim + theme_void() + ttl + labs(title = title, subtitle = subtitle)
  }
  p2 <- square(sqrt(0.5), "1", "1.41", TRUE, "A square with sides 1", "sides 1, as in the table,\nbut diagonals of 1.41, not 2")
  p3 <- square(1, "1.41", "2", FALSE, "Classical MDS", "diagonals 2, as in the table,\nbut sides of 1.41, not 1")
  # 4: why no arrangement in any number of dimensions works
  p4 <- ggplot() +
    annotate("segment", x = -1, xend = 1, y = 0, yend = 0, colour = "#2c3e50", linewidth = 0.8) +
    annotate("point", x = c(-1, 1), y = 0, size = 3, colour = "#2c3e50") +
    annotate("point", x = 0, y = 0, size = 3.4, colour = col_a2) +
    annotate("text", x = c(-1.12, 1.12), y = 0, label = c("A", "C"), size = 5, fontface = "bold",
             colour = "#2c3e50", hjust = c(1, 0)) +
    annotate("text", x = 0, y = c(0.2, -0.2), label = c("B", "D"), size = 5, fontface = "bold", colour = col_a2) +
    annotate("text", x = c(-0.5, 0.5), y = 0.12, label = "1", size = 4.2, colour = "grey35") +
    annotate("text", x = 0, y = -0.85, size = 3.9, colour = "grey20", lineheight = 1.05,
             label = "B and D would be the same point,\nyet they are 2 minutes apart.") +
    lim + theme_void() + ttl +
    labs(title = "Why nothing works", subtitle = "A to B to C takes 1 + 1 = 2, as A to C does.\nSo B sits halfway along the line, and so does D.")
  wrap_plots(p1, p2, p3, p4, nrow = nrow)
}

# Eigenvalues of B for five distance tables: Euclidean ones have none below zero
eig_table <- function() {
  sets <- list(
    "US map, flat (Euclidean)" = us_D,
    "US, great circle" = us_G,
    "States, Manhattan" = st_D$Manhattan,
    "Loop road, no waiting" = D_wait(0),
    "Words, average distance" = word_D)
  purrr::imap_dfr(sets, function(D, nm) { ev <- mds_eigen(D)
    tibble(set = nm, j = seq_along(ev), ev = ev, rel = ev / ev[1]) }) |>
    mutate(set = factor(set, levels = names(sets)))
}
fig_eigen_types <- function() {
  d <- eig_table()
  brk <- function(l) { b <- pretty(c(1, l[2]), n = 4); c(1, b[b > 1 & b <= l[2]]) }
  ggplot(d, aes(j, rel, fill = rel < -1e-8)) + geom_col(width = 0.75) + geom_hline(yintercept = 0, colour = "grey40") +
    facet_wrap(~ set, nrow = 1, scales = "free_x") +
    scale_fill_manual(values = c(`FALSE` = col_a1, `TRUE` = col_a2), labels = c("positive", "negative"), name = NULL) +
    scale_x_continuous(breaks = brk) +
    labs(x = "eigenvalue number, largest first (all of them)", y = "eigenvalue / largest eigenvalue") +
    theme(legend.position = "top", strip.text = element_text(size = 10, face = "bold"))
}


# =============================================================================
# 6  How good is the map?
# =============================================================================
# Loop road: the eigenvalues as the waiting time grows. The second shrinks and vanishes at w_star (a line);
# the negative ones show that the travel times are not Euclidean.
fig_loop_eigen <- function() {
  ws <- seq(0, 40, by = 0.5)
  d <- purrr::map_dfr(ws, function(w) { ev <- mds_eigen(D_wait(w))
    tibble(w = w, which = c("largest", "second", "third", "smallest (most negative)"),
           ev = c(ev[1], ev[2], ev[3], min(ev[length(ev)], 0))) }) |>
    mutate(which = factor(which, levels = c("largest", "second", "third", "smallest (most negative)")))
  ggplot(d, aes(w, ev, colour = which)) + geom_hline(yintercept = 0, colour = "grey60") +
    geom_vline(xintercept = w_star, linetype = "22", colour = col_wait) +
    annotate("text", x = w_star + 0.6, y = 900, label = sprintf("%.1f min:\nnobody drives\nthrough", w_star),
             hjust = 0, size = 3.6, colour = col_wait) +
    geom_line(linewidth = 1.2) +
    scale_colour_manual(values = c(largest = col_a1, second = col_a2, third = "#7b4fb5",
                                   `smallest (most negative)` = "grey35"), name = NULL) +
    labs(x = "waiting time at the construction site (minutes)", y = expression("eigenvalue of"~bold(B))) +
    theme(legend.position = "right")
}

# Shepard diagrams: distances on the map against the distances in the table
shepard_df <- function(D, Y, nm) {
  dm <- as.matrix(dist(Y)); i <- which(upper.tri(D))
  tibble(set = nm, obs = D[i], map = dm[i])
}
fig_shepard <- function() {
  d <- bind_rows(
    shepard_df(us_D, us_mds$points, "US map (2 of 2 dimensions)"),
    shepard_df(st_D$Euclidean, st_mds$Euclidean$points, "States, 6 variables (2 dimensions)"),
    shepard_df(D_wait(0), cmdscale(as.dist(D_wait(0)), k = 2), "Loop road (2 dimensions)"),
    shepard_df(word_D, word_mds$points, "Words (2 dimensions)")) |>
    mutate(set = factor(set, levels = unique(set)))
  ggplot(d, aes(obs, map)) + geom_abline(slope = 1, intercept = 0, colour = col_a2, linetype = "22") +
    geom_point(alpha = 0.35, size = 1.1, colour = "#2c3e50") +
    facet_wrap(~ set, nrow = 1, scales = "free") +
    labs(x = "distance in the table", y = "distance on the map") +
    theme(strip.text = element_text(size = 10, face = "bold"))
}


# =============================================================================
# 7  When only the order counts: non-metric MDS
# =============================================================================
# The travel times on the loop road without waiting: classical MDS cannot reproduce them (negative eigenvalues),
# but their order is that of straight-line distances around a circle, so non-metric MDS fits them exactly
loop_D0 <- D_wait(0)
loop_cl <- cmdscale(as.dist(loop_D0), k = 2)                         # classical MDS, the start for non-metric MDS
loop_iso <- local({ set.seed(1); MASS::isoMDS(as.dist(loop_D0), y = loop_cl, k = 2, trace = FALSE) })
loop_shep <- MASS::Shepard(as.dist(loop_D0), loop_iso$points)

fig_nonmetric <- function() {
  Y <- loop_iso$points
  d <- tibble(place = places$place, x = Y[, 1], y = Y[, 2])
  p1 <- ggplot(d, aes(x, y)) + axes0 + geom_path(data = bind_rows(d, d[1, ]), colour = "#cccccc", linewidth = 2.5) +
    geom_point(size = 2.6, colour = "#2c3e50") +
    geom_text_repel(aes(label = place), size = 3.8, seed = 3, colour = "#2c3e50") + coord_equal() +
    labs(title = "Non-metric MDS of the travel times",
         subtitle = sprintf("no waiting; stress %.1f%%; grey line: the road", loop_iso$stress),
         x = "dimension 1", y = "dimension 2")
  sh <- tibble(x = loop_shep$x, y = loop_shep$y, yf = loop_shep$yf)
  p2 <- ggplot(sh, aes(x, y)) + geom_step(aes(y = yf), colour = col_a2, linewidth = 1) +
    geom_point(size = 2, colour = "#2c3e50") +
    labs(title = "Shepard diagram", subtitle = "red: the best order-preserving fit",
         x = "travel time in the table (minutes)", y = "distance on the map")
  p1 + p2 + plot_layout(widths = c(1, 1))
}
