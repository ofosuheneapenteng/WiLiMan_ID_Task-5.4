# =========================================================
# ASF Map + Animation + Ensemble CI (Full Updated Code)
# =========================================================

# -------------------------
# 0) Setup
# -------------------------
rm(list = ls())

suppressPackageStartupMessages({
  library(sf)
  library(ggplot2)
  library(rnaturalearth)
  library(rnaturalearthdata)
  library(grid)          # unit()
  library(magick)        # GIF + frame device
  library(av)            # MP4
  library(reshape2)
  library(dplyr)
  library(patchwork)
})

setwd("~/Desktop/TASK5.4/ASF_Modelling")
sf::sf_use_s2(FALSE)   # planar buffers in EPSG:25832

set.seed(42)
dir.create("Figures_Final", showWarnings = FALSE)

# -------------------------
# 1) Style
# -------------------------
style_defaults <- function(
    pal_states = c(S="#595959", E="#0072B2", I="#D55E00", D="#009E73", C="#000000"),
    zone_prot_col  = "#E41A1C", zone_prot_alpha = 0.45, zone_prot_lwd = 1.7,
    zone_surv_col  = "#4DAF4A", zone_surv_alpha = 0.30, zone_surv_lwd = 1.3,
    farm_point_size = 1.7, farm_point_alpha = 0.95,
    basemap_fill = "white", basemap_border = "grey65", basemap_lwd = 0.35,
    grid_major = "grey90", base_size = 13
){
  list(
    pal_states = pal_states,
    zone_prot_col = zone_prot_col, zone_prot_alpha = zone_prot_alpha, zone_prot_lwd = zone_prot_lwd,
    zone_surv_col = zone_surv_col, zone_surv_alpha = zone_surv_alpha, zone_surv_lwd = zone_surv_lwd,
    farm_point_size = farm_point_size, farm_point_alpha = farm_point_alpha,
    basemap_fill = basemap_fill, basemap_border = basemap_border, basemap_lwd = basemap_lwd,
    grid_major = grid_major, base_size = base_size
  )
}
STYLE <- style_defaults()

# -------------------------
# 2) Basemap + farms (EPSG:25832)
# -------------------------
den <- rnaturalearth::ne_countries(country = "Denmark", scale = 10, returnclass = "sf")
den_25832 <- st_transform(den, 25832)

# Hygiene: keep polygons only
if ("st_make_valid" %in% getNamespaceExports("sf")) {
  den_25832 <- sf::st_make_valid(den_25832)
}
den_25832 <- sf::st_collection_extract(den_25832, "POLYGON")

# Farms
N <- 1800
farms_geom <- st_sample(den_25832, size = N, type = "random", exact = TRUE)
farms_sf   <- st_as_sf(data.frame(farm_id = 1:N), geometry = farms_geom)
coords     <- st_coordinates(farms_sf)

# -------------------------
# 3) KNN network
# -------------------------
build_knn_weights <- function(coords_xy, K = 30, weight_rate = 3) {
  N <- nrow(coords_xy)
  dmat <- as.matrix(dist(coords_xy))
  W <- matrix(0, N, N)
  for (i in 1:N) {
    nn <- order(dmat[i, ])[2:(K+1)]
    W[i, nn] <- rexp(length(nn), rate = weight_rate)
  }
  diag(W) <- 0
  pmax(W, t(W))
}
W_base <- build_knn_weights(coords, K = 30, weight_rate = 3)

# -------------------------
# 4) Simulator (safe snapshot copies)
# -------------------------
run_seidc_controlled <- function(
    W_base, farms_sf,
    t_max = 140, snapshot_dt = 0.5,
    # Natural history
    sigma = 1/4, delta_0 = 1/6, kappa_cull = 1/4,
    # Transmission
    beta_m  = 0.22, rho_pre = 0.05, beta_e  = 0.015,
    omega   = 0.10, gamma   = 0.30, eta_val = 1e-5,
    # Controls & zones
    prot_radius_m = 3000, surv_radius_m = 10000,
    prot_scale = 0.30, surv_scale = 0.50,
    beta_m_post_scale = 0.25, beta_e_post_scale = 0.50, rho_post = 0.02,
    seed_id = 1L, rng_seed = 123
){
  set.seed(rng_seed)
  N <- nrow(W_base)
  S<-1L; E<-2L; I<-3L; D<-4L; C<-5L
  
  X <- rep(S, N); X[seed_id] <- E
  Eenv <- rep(0, N); t <- 0
  
  detect_time <- rep(NA_real_, N)
  first_det_time <- NA_real_
  blocked <- rep(FALSE, N)
  W_active <- W_base
  
  snap_times  <- c(0)
  snap_states <- list(as.integer(X))  # deep numeric copy
  next_snap   <- snapshot_dt
  
  recompute_active_W <- function() {
    det_ids <- which(!is.na(detect_time))
    if (length(det_ids) == 0) {
      s <- rep(1, N)
    } else {
      d_centers <- farms_sf[det_ids, , drop = FALSE]
      pts <- st_geometry(d_centers)
      ring_prot_geom <- st_buffer(st_union(pts), prot_radius_m)
      ring_surv_geom <- st_buffer(st_union(pts), surv_radius_m)
      in_prot     <- as.logical(st_intersects(farms_sf, ring_prot_geom, sparse = FALSE))
      in_surv_all <- as.logical(st_intersects(farms_sf, ring_surv_geom, sparse = FALSE))
      in_surv     <- in_surv_all & !in_prot
      s <- rep(1, N); s[in_surv] <- surv_scale; s[in_prot] <- prot_scale
    }
    W_new <- W_base * (s %o% s)
    if (any(blocked)) { W_new[blocked, ] <- 0; W_new[, blocked] <- 0 }
    pmax(W_new, t(W_new))
  }
  
  repeat {
    Svec <- as.numeric(X==S); Evec <- as.numeric(X==E)
    Ivec <- as.numeric(X==I); Dvec <- as.numeric(X==D)
    
    controls_on <- !is.na(first_det_time)
    beta_m_eff  <- if (controls_on) beta_m * beta_m_post_scale else beta_m
    beta_e_eff  <- if (controls_on) beta_e * beta_e_post_scale else beta_e
    rho_eff     <- if (controls_on) rho_post else rho_pre
    
    P_I <- as.numeric(W_active %*% Ivec)
    P_D <- as.numeric(W_active %*% Dvec)
    lambda <- beta_m_eff*P_I + rho_eff*beta_m_eff*P_D + beta_e_eff*Eenv + eta_val
    
    a_SE <- lambda * Svec; a_EI <- sigma  * Evec
    a_ID <- delta_0* Ivec; a_DC <- kappa_cull * Dvec
    
    a  <- c(a_SE, a_EI, a_ID, a_DC)
    a0 <- sum(a)
    if (a0 <= 0 || t >= t_max) break
    
    tau <- rexp(1, a0); t <- t + tau
    
    dec <- exp(-gamma * tau)
    Eenv <- Eenv*dec + (Ivec==1)*(omega/gamma)*(1 - dec)
    
    cj <- cumsum(a); j <- which(cj >= runif(1)*a0)[1]
    
    if     (j <=   N) { i <- j;        if (X[i]==S) X[i] <- E }
    else if(j <= 2*N) { i <- j -   N;  if (X[i]==E) X[i] <- I }
    else if(j <= 3*N) { i <- j - 2*N;  if (X[i]==I) {
      X[i] <- D
      if (is.na(detect_time[i])) detect_time[i] <- t
      if (is.na(first_det_time)) first_det_time <- t
      blocked[i] <- TRUE
      W_active <- recompute_active_W()
    }} else            { i <- j - 3*N; if (X[i]==D) {
      X[i] <- C; blocked[i] <- TRUE
      W_active <- recompute_active_W()
    }}
    
    if (t >= next_snap) {
      snap_times  <- c(snap_times, t)
      snap_states[[length(snap_states)+1]] <- as.integer(X + 0L)  # deep copy
      next_snap <- next_snap + snapshot_dt
    }
  }
  
  list(
    snap_times = snap_times,
    snap_states = snap_states,
    detect_time = detect_time,
    first_det_time = first_det_time,
    X_final = X
  )
}

# -------------------------
# 5) Helpers: time & zones & plotting data
# -------------------------
states_at_time <- function(sim, t) {
  idxs <- which(sim$snap_times <= t)
  if (length(idxs) == 0) return(as.integer(sim$snap_states[[1]]))
  as.integer(sim$snap_states[[max(idxs)]])
}
detected_centers_up_to <- function(sim, farms_sf, t) {
  ids <- which(!is.na(sim$detect_time) & sim$detect_time <= t)
  if (length(ids) == 0) return(farms_sf[0, , drop = FALSE])
  farms_sf[ids, , drop = FALSE]
}
snap_index_at_or_before <- function(sim, t) {
  idxs <- which(sim$snap_times <= t)
  if (length(idxs) == 0) list(idx = 1L, time_used = sim$snap_times[1])
  else { i <- max(idxs); list(idx = i, time_used = sim$snap_times[i]) }
}

# Convert POLYGON / MULTIPOLYGON sfc → data.frame for geom_polygon()
sf_to_polygon_df <- function(sfc, zone_label) {
  if (length(sfc) == 0) return(data.frame())
  sfc <- st_cast(sfc, "MULTIPOLYGON", warn = FALSE)
  sfc <- st_cast(sfc, "POLYGON", warn = FALSE)
  cc  <- st_coordinates(sfc)
  df  <- as.data.frame(cc)
  names(df)[1:2] <- c("x","y")
  grp_cols <- intersect(c("L1","L2","L3","L4"), names(df))
  df$group <- do.call(paste, c(df[grp_cols], sep = "_"))
  df$zone  <- zone_label
  df[, c("x","y","group","zone")]
}

# -------------------------
# 6) Map (zones as polygons, farms as points; stable)
# -------------------------
plot_map_with_zones <- function(sim, farms_sf, den_25832, t, style = STYLE,
                                prot_radius_m = 3000, surv_radius_m = 10000) {
  
  # States -> factor
  xt <- states_at_time(sim, t)
  state_fac <- factor(xt, levels = 1:5, labels = c("S","E","I","D","C"))
  
  # Exact snapshot time used (for title)
  info <- snap_index_at_or_before(sim, t); time_used <- info$time_used
  
  # Farms -> plain df
  xy <- st_coordinates(farms_sf)
  farms_df <- data.frame(x = xy[,1], y = xy[,2], state = state_fac)
  
  # Basemap -> polygon df
  den_df <- sf_to_polygon_df(st_geometry(den_25832), NA)
  
  # Zones
  d_centers <- detected_centers_up_to(sim, farms_sf, t)
  prot_df <- surv_df <- data.frame()
  have_zones <- nrow(d_centers) > 0
  if (have_zones) {
    pts <- st_geometry(d_centers)
    prot_sfc <- st_collection_extract(st_buffer(st_union(pts), prot_radius_m), "POLYGON")
    surv_sfc <- st_collection_extract(st_buffer(st_union(pts), surv_radius_m), "POLYGON")
    if (length(prot_sfc) > 0) prot_df <- sf_to_polygon_df(prot_sfc, "Protection (3 km)")
    if (length(surv_sfc) > 0) surv_df <- sf_to_polygon_df(surv_sfc, "Surveillance (10 km)")
    have_zones <- nrow(prot_df) + nrow(surv_df) > 0
  }
  
  p <- ggplot() +
    # basemap
    geom_polygon(data = den_df,
                 aes(x, y, group = group),
                 fill = style$basemap_fill,
                 color = style$basemap_border,
                 linewidth = style$basemap_lwd) +
    # rings
    { if (nrow(surv_df))
      geom_polygon(data = surv_df,
                   aes(x, y, group = group, fill = zone),
                   color = style$zone_surv_col,
                   alpha = style$zone_surv_alpha,
                   linewidth = style$zone_surv_lwd)
    } +
    { if (nrow(prot_df))
      geom_polygon(data = prot_df,
                   aes(x, y, group = group, fill = zone),
                   color = style$zone_prot_col,
                   alpha = style$zone_prot_alpha,
                   linewidth = style$zone_prot_lwd)
    } +
    # farms
    geom_point(data = farms_df, aes(x, y, color = state),
               size = style$farm_point_size, alpha = style$farm_point_alpha) +
    coord_equal() +
    scale_color_manual(values = style$pal_states, name = "Compartment", drop = FALSE) +
    theme_minimal(base_size = style$base_size) +
    theme(
      legend.position   = "bottom",
      legend.box        = "vertical",
      legend.title      = element_text(face = "bold"),
      legend.key.width  = unit(1.5, "lines"),
      legend.key.height = unit(1.0, "lines"),
      panel.grid.major  = element_line(color = style$grid_major),
      axis.title        = element_blank()
    ) +
    #labs(title = sprintf("ASF Zoning Snapshot — t = %.2f days", time_used))
    labs(title = sprintf("ASF simulation (Day %d)", round(time_used - first_det)))
    
  if (have_zones) {
    p <- p + scale_fill_manual(
      values = c("Protection (3 km)"    = style$zone_prot_col,
                 "Surveillance (10 km)" = style$zone_surv_col),
      name   = "Zones",
      guide  = guide_legend(override.aes = list(alpha = 1))
    )
  }
  p
}

# -------------------------
# 7) Scalebar + North arrow
# -------------------------
add_scalebar_northarrow <- function(p, den_25832,
                                    length_km = 50,
                                    pos = c("br","bl","tr","tl"),
                                    arrow = TRUE,
                                    bar_height_km = 4,
                                    text_size = 3.3) {
  pos <- match.arg(pos)
  bb <- sf::st_bbox(den_25832)
  width  <- as.numeric(bb["xmax"] - bb["xmin"])
  height <- as.numeric(bb["ymax"] - bb["ymin"])
  
  len_m_req <- length_km * 1000
  len_m     <- min(len_m_req, 0.35 * width)
  half_m <- len_m / 2
  bar_h  <- bar_height_km * 1000
  
  pad_x <- 0.02 * width
  pad_y <- 0.02 * height
  
  x0 <- if (pos %in% c("bl","tl")) bb["xmin"] + pad_x else bb["xmax"] - pad_x - len_m
  y0 <- if (pos %in% c("bl","br")) bb["ymin"] + pad_y else bb["ymax"] - pad_y - bar_h
  
  segs <- data.frame(
    xmin = c(x0, x0 + half_m),
    xmax = c(x0 + half_m, x0 + len_m),
    ymin = c(y0, y0),
    ymax = c(y0 + bar_h, y0 + bar_h),
    fill = c("black","white")
  )
  text_y  <- y0 + bar_h + (0.006 * height)
  km_lab  <- round(len_m/1000)
  lab_df  <- data.frame(
    x   = c(x0, x0 + half_m, x0 + len_m),
    y   = rep(text_y, 3),
    lab = c("0", paste0(km_lab/2, " km"), paste0(km_lab, " km"))
  )
  
  layers <- list(
    geom_rect(data = segs,
              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              inherit.aes = FALSE, color = "black",
              fill = segs$fill, linewidth = 0.35),
    geom_text(data = lab_df,
              aes(x, y, label = lab),
              inherit.aes = FALSE, size = text_size)
  )
  if (arrow) {
    xA  <- if (pos %in% c("bl","tl")) (x0 + len_m + 0.02*width) else (x0 - 0.02*width)
    yA0 <- y0
    yA1 <- y0 + 0.12 * height
    layers <- c(layers, list(
      geom_segment(aes(x = xA, xend = xA, y = yA0, yend = yA1),
                   arrow = grid::arrow(type = "closed", length = unit(6, "mm")),
                   inherit.aes = FALSE, linewidth = 0.6),
      annotate("text", x = xA, y = yA1 + 0.015*height,
               label = "N", fontface = "bold", size = text_size + 0.6)
    ))
  }
  p + layers
}

# -------------------------
# 8) Three‑panel (first detection, peak I, final)
# -------------------------
three_panel_times <- function(sim) {
  first_det <- if (all(is.na(sim$detect_time))) NA_real_ else min(sim$detect_time, na.rm = TRUE)
  ts_df <- data.frame(
    time = sim$snap_times,
    t(vapply(sim$snap_states, function(x) tabulate(x, nbins = 5L), integer(5)))
  )
  names(ts_df)[2:6] <- c("S","E","I","D","C")
  t_peakI <- ts_df$time[ which.max(ts_df$I) ]
  t_final <- max(sim$snap_times, na.rm = TRUE)
  list(first = if (is.finite(first_det)) first_det else ts_df$time[1],
       peakI = t_peakI, final = t_final)
}

save_three_panel <- function(sim, farms_sf, den_25832, style = STYLE,
                             scalebar_km = 50, add_arrow = TRUE,
                             out_file = "Figures_Final/asf_three_panel.png",
                             width = 12, height = 6.2, dpi = 300) {
  times <- three_panel_times(sim)
  p1 <- plot_map_with_zones(sim, farms_sf, den_25832, times$first, style) +
    labs(subtitle = sprintf("First detection — t = %.1f d", times$first))
  p2 <- plot_map_with_zones(sim, farms_sf, den_25832, times$peakI, style) +
    labs(subtitle = sprintf("Peak I — t = %.1f d", times$peakI))
  p3 <- plot_map_with_zones(sim, farms_sf, den_25832, times$final, style) +
    labs(subtitle = sprintf("Final — t = %.1f d", times$final))
  
  p1 <- add_scalebar_northarrow(p1, den_25832, length_km = scalebar_km, arrow = add_arrow)
  p2 <- add_scalebar_northarrow(p2, den_25832, length_km = scalebar_km, arrow = add_arrow)
  p3 <- add_scalebar_northarrow(p3, den_25832, length_km = scalebar_km, arrow = add_arrow)
  
  p_grid <- p1 + p2 + p3 + patchwork::plot_layout(ncol = 3) &
    theme(plot.title = element_text(size = style$base_size + 0.6, face = "bold"),
          plot.subtitle = element_text(size = style$base_size - 1))
  
  ggsave(out_file, plot = p_grid, width = width, height = height, dpi = dpi)
  message("Saved 3‑panel figure: ", normalizePath(out_file, winslash = "/", mustWork = FALSE))
  invisible(list(times = times, plot = p_grid))
}

# -------------------------
# 9) Export sequence (exact first detection grid) + Animation
# -------------------------
time_grid_every <- function(sim, by = 10) {
  t0 <- if (all(is.na(sim$detect_time))) 0 else min(sim$detect_time, na.rm = TRUE)
  t_max <- max(sim$snap_times, na.rm = TRUE)
  t0 + seq(0, floor((t_max - t0) / by)) * by
}

save_maps_every_10 <- function(sim, farms_sf, den_25832, style = STYLE,
                               by = 10,
                               out_dir = "Figures_Final",
                               prefix = "asf_map") {
  dir.create(out_dir, showWarnings = FALSE)
  tg <- time_grid_every(sim, by = by)
  message("Exporting maps at requested times: ", paste(sprintf("%.2f", tg), collapse = ", "))
  
  plots <- vector("list", length(tg))
  for (k in seq_along(tg)) {
    t_req <- tg[k]
    p <- plot_map_with_zones(sim, farms_sf, den_25832, t_req, style)
    info <- snap_index_at_or_before(sim, t_req)
    time_used <- info$time_used
    p <- p + labs(title = sprintf("ASF simulation (Day %d)", round(time_used - first_det)))
    plots[[k]] <- p
    fn_png <- file.path(out_dir, sprintf("%s_%03d.png", prefix, as.integer(t_req)))
    ggsave(fn_png, plot = p, width = 10, height = 7, dpi = 300)
  }
  fn_pdf <- file.path(out_dir, sprintf("%s_every_%d.pdf", prefix, by))
  grDevices::pdf(fn_pdf, width = 10, height = 7)
  for (p in plots) print(p)
  grDevices::dev.off()
  message("Saved multi‑page PDF: ", normalizePath(fn_pdf, winslash = "/", mustWork = FALSE))
  invisible(list(times = tg, plots = plots, pdf = fn_pdf))
}

animate_map_sequence <- function(sim, farms_sf, den_25832, style = STYLE,
                                 by = 10, fps = 6,
                                 out_gif = "Figures_Final/asf_anim.gif",
                                 out_mp4 = "Figures_Final/asf_anim.mp4",
                                 width = 1350, height = 900, res = 110) {
  tg <- time_grid_every(sim, by)
  if (!length(tg)) stop("No times in grid.")
  message("Animating times: ", paste(sprintf("%.2f", tg), collapse = ", "))
  
  frames <- vector("list", length(tg))
  for (i in seq_along(tg)) {
    t_req <- tg[i]
    p <- plot_map_with_zones(sim, farms_sf, den_25832, t_req, style)
    info <- snap_index_at_or_before(sim, t_req)
    time_used <- info$time_used
    p <- p + labs(title = sprintf("ASF simulation (Day %d)", round(time_used - first_det)))
    img <- magick::image_graph(width = width, height = height, res = res, bg = "white")
    print(p); dev.off()
    frames[[i]] <- img
  }
  
  delay <- max(1L, as.integer(round(100 / fps)))  # centiseconds
  eff_fps <- 100 / delay
  message(sprintf("GIF delay=%d cs → effective fps=%.2f", delay, eff_fps))
  gif <- magick::image_animate(magick::image_join(frames), delay = delay)
  magick::image_write(gif, path = out_gif)
  message("GIF saved: ", normalizePath(out_gif, winslash = "/", mustWork = FALSE))
  
  if (requireNamespace("av", quietly = TRUE)) {
    tmpdir <- tempdir()
    png_files <- file.path(tmpdir, sprintf("frame_%05d.png", seq_along(frames)))
    for (i in seq_along(frames)) magick::image_write(frames[[i]], png_files[i])
    av::av_encode_video(
      input     = png_files,
      output    = out_mp4,
      framerate = fps,
      codec     = "libx264",
      vfilter   = "format=yuv420p,scale=1280:-2",
      verbose   = FALSE
    )
    message("MP4 saved: ", normalizePath(out_mp4, winslash = "/", mustWork = FALSE))
  } else {
    warning("Install 'av' for MP4 export.")
  }
  invisible(list(times = tg, gif = out_gif, mp4 = out_mp4))
}

# -------------------------
# 10) Ensemble CI utilities (batch)
# -------------------------
scale_factor=100/N #1

compute_ts <- function(sim, scale_factor) {
  mat <- t(vapply(sim$snap_states, function(x) tabulate(as.integer(x), nbins = 5L), integer(5)))
  df  <- data.frame(time = sim$snap_times, mat)
  names(df)[2:6] <- c("S","E","I","D","C")
  df[, 2:6] <- df[, 2:6] * scale_factor
  df
}
interp_run_to_grid <- function(sim, times_grid, scale_factor) {
  ts_df <- compute_ts(sim, scale_factor)
  comps <- c("S","E","I","D","C")
  out <- data.frame(time = times_grid)
  for (cm in comps) {
    out[[cm]] <- approx(x = ts_df$time, y = ts_df[[cm]], xout = times_grid, rule = 2, ties = "ordered")$y
  }
  out_long <- reshape2::melt(out, id.vars = "time", variable.name = "Compartment", value.name = "Count")
  out_long$Compartment <- factor(out_long$Compartment, levels = comps)
  out_long
}
run_batch_extract_counts <- function(
    n_runs, times_grid, W_base, farms_sf, scale_factor,
    parallel = FALSE, seed0 = 1001, ...
){
  set.seed(seed0)
  out_list <- vector("list", n_runs)
  for (r in seq_len(n_runs)) {
    sim <- run_seidc_controlled(W_base, farms_sf, rng_seed = seed0 + r, ...)
    df  <- interp_run_to_grid(sim, times_grid, scale_factor); df$run <- r
    out_list[[r]] <- df
  }
  batch_df <- do.call(rbind, out_list)
  rownames(batch_df) <- NULL
  batch_df
}
summarise_ci_dplyr <- function(batch_df) {
  batch_df %>%
    group_by(time, Compartment) %>%
    summarise(
      lo  = quantile(Count, 0.025, na.rm = TRUE),
      med = quantile(Count, 0.50,  na.rm = TRUE),
      hi  = quantile(Count, 0.975, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(Compartment, time)
}
plot_ci_facets <- function(ci_df, style = STYLE,
                           title = "S/E/I/D/C — Median and 95% CI across runs") {
  comps <- c("S","E","I","D","C")
  ci_df <- ci_df %>% mutate(Compartment = factor(Compartment, levels = comps))
  pal   <- unname(style$pal_states[comps])
  fills <- grDevices::adjustcolor(pal, alpha.f = 0.25)
  ggplot(ci_df, aes(x = time)) +
    geom_ribbon(aes(ymin = lo, ymax = hi, fill = Compartment), color = NA, na.rm = TRUE) +
    geom_line(aes(y = med, color = Compartment), linewidth = 1.0, na.rm = TRUE) +
    facet_wrap(~ Compartment, scales = "free_y", ncol = 2) +
    scale_color_manual(values = pal, guide = "none") +
    scale_fill_manual(values = fills, guide = "none") +
    labs(title = title, x = "Time (days)", y = "Count") +
    theme_minimal(base_size = style$base_size) +
    theme(panel.grid.major = element_line(color = style$grid_major),
          strip.text = element_text(face = "bold"))
}

# =========================================================
# HOW TO RUN
# =========================================================

# A) One simulation and map
message("\n>>> Running SEIDC with EFSA-like controls ...")
sim <- run_seidc_controlled(W_base, farms_sf)

first_det <- if (all(is.na(sim$detect_time))) NA_real_ else min(sim$detect_time, na.rm = TRUE)
t_show    <- if (is.finite(first_det)) first_det + 2 else 40

p_map <- plot_map_with_zones(sim, farms_sf, den_25832, t_show, STYLE)
p_map <- add_scalebar_northarrow(p_map, den_25832, length_km = 50, pos = "br", arrow = TRUE)
print(p_map)
ggsave("Figures_Final/asf_map_final.png", p_map, width = 10, height = 7, dpi = 300)

# B) Three‑panel
save_three_panel(sim, farms_sf, den_25832, STYLE,
                 scalebar_km = 50, add_arrow = TRUE,
                 out_file = "Figures_Final/asf_three_panel.png")

# C) Sequence export (exact first detection grid)
save_maps_every_10(sim, farms_sf, den_25832, STYLE, by = 10,
                   out_dir = "Figures_Final", prefix = "asf_map_fd")

# D) Animation (GIF + MP4)
animate_map_sequence(sim, farms_sf, den_25832, STYLE,
                     by = 10, fps = 6,
                     out_gif = "Figures_Final/asf_anim.gif",
                     out_mp4 = "Figures_Final/asf_anim.mp4")

# E) Ensemble CI across runs (adjust n_runs / grid for speed)
times_grid   <- seq(0, 140, by = 1)       # coarse grid; use 0.5 for smoother
scale_factor <- 100 / N                 # or =1 for raw counts

batch_df <- run_batch_extract_counts(
  n_runs = 50, times_grid = times_grid,
  W_base = W_base, farms_sf = farms_sf,
  scale_factor = scale_factor, parallel = FALSE,
  t_max = 140, snapshot_dt = 0.5
)

ci_df_controls <- summarise_ci_dplyr(batch_df)
# p_ci  <- plot_ci_facets(ci_df, STYLE)
# print(p_ci)

p_controls <- plot_ci_facets(
  ci_df_controls,
  STYLE,
  title = "" #Projected ASF dynamics with and without EFSA-like controls"
)+
  labs(
    #     subtitle = "Median and 95% CI across stochastic simulations",
    #     x = "Time (days)",
    y = "Proportion of farms (%)"
  )

print(p_controls)

ggsave("Figures_Final/ts_all_95CI_facets.png", p_controls, width = 10, height = 6.5, dpi = 300)
ggsave("Figures_Final/ts_all_95CI_facets.pdf",  p_controls, width = 10, height = 6.5)

########################---------#######################
N=1800 # No of farms
scale_factor <- 100/N

run_seidc_no_controls <- function(
    W_base, farms_sf,
    t_max = 140,
    snapshot_dt = 0.5,
    
    sigma = 1/4,
    beta_m = 0.22,
    rho_pre = 0.05,
    beta_e = 0.015,
    
    omega = 0.10,
    gamma = 0.30,
    eta_val = 1e-5,
    
    seed_id = 1L,
    rng_seed = 123
){
  
  set.seed(rng_seed)
  
  N <- nrow(W_base)
  
  S <- 1L
  E <- 2L
  I <- 3L
  D <- 4L
  C <- 5L
  
  X <- rep(S, N)
  X[seed_id] <- E
  
  Eenv <- rep(0, N)
  t <- 0
  
  detect_time <- rep(NA_real_, N)
  first_det_time <- NA_real_
  
  snap_times <- c(0)
  snap_states <- list(as.integer(X))
  
  next_snap <- snapshot_dt
  
  repeat {
    
    Svec <- as.numeric(X == S)
    Evec <- as.numeric(X == E)
    Ivec <- as.numeric(X == I)
    
    P_I <- as.numeric(W_base %*% Ivec)
    
    lambda <-
      beta_m * P_I +
      rho_pre * beta_m * 0 +
      beta_e * Eenv +
      eta_val
    
    # ONLY S -> E -> I
    
    a_SE <- lambda * Svec
    a_EI <- sigma * Evec
    
    # NO DETECTION
    a_ID <- rep(0, N)
    
    # NO CULLING
    a_DC <- rep(0, N)
    
    a <- c(a_SE, a_EI, a_ID, a_DC)
    
    a0 <- sum(a)
    
    if(a0 <= 0 || t >= t_max)
      break
    
    tau <- rexp(1, a0)
    
    t <- t + tau
    
    dec <- exp(-gamma * tau)
    
    Eenv <- Eenv * dec +
      (Ivec == 1) * (omega/gamma) * (1 - dec)
    
    cj <- cumsum(a)
    
    j <- which(cj >= runif(1) * a0)[1]
    
    if(j <= N){
      
      i <- j
      
      if(X[i] == S)
        X[i] <- E
      
    } else if(j <= 2*N){
      
      i <- j - N
      
      if(X[i] == E)
        X[i] <- I
      
    }
    
    if(t >= next_snap){
      
      snap_times <- c(snap_times, t)
      
      snap_states[[length(snap_states)+1]] <- as.integer(X)
      
      next_snap <- next_snap + snapshot_dt
    }
  }
  
  list(
    snap_times = snap_times,
    snap_states = snap_states,
    detect_time = detect_time,
    first_det_time = first_det_time,
    X_final = X
  )
}

#N=1800 # No of farms
#scale_factor <- 100/N

run_batch_extract_counts_no_controls <- function(
    n_runs,
    times_grid,
    W_base,
    farms_sf,
    scale_factor = scale_factor,
    seed0 = 1001,
    ...
){
  
  out_list <- vector("list", n_runs)
  
  for(r in seq_len(n_runs)){
    
    sim <- run_seidc_no_controls(
      W_base,
      farms_sf,
      rng_seed = seed0 + r,
      ...
    )
    
    df <- interp_run_to_grid(
      sim,
      times_grid,
      scale_factor
    )
    
    df$run <- r
    
    out_list[[r]] <- df
  }
  
  batch_df <- do.call(rbind, out_list)
  
  rownames(batch_df) <- NULL
  
  batch_df
}

batch_df_no_controls <- run_batch_extract_counts_no_controls(
  n_runs = 50,
  times_grid = times_grid,
  W_base = W_base,
  farms_sf = farms_sf,
  scale_factor = scale_factor,
  t_max = 140,
  snapshot_dt = 0.5
)

ci_df_no_controls <- summarise_ci_dplyr(batch_df_no_controls)

p_no_controls <- plot_ci_facets(
  ci_df_no_controls,
  STYLE,
  title = "" # ASF Model Without Controls
 )+
  labs(
#     subtitle = "Median and 95% CI across stochastic simulations",
#     x = "Time (days)",
 y = "Proportion of farms (%)"
   )

print(p_no_controls)

##Both graphs

library(patchwork)
p_panel <-
  (p_no_controls + p_controls) +
  plot_layout(ncol = 2) +
  plot_annotation(
    title = "Projected ASF dynamics without and with EFSA-like controls",
    subtitle = "Median trajectories with 95% confidence intervals",
    theme = theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 14)
    )
  ) &
  theme(
    axis.title.x = element_text(size = 14),
    axis.title.y = element_text(size = 14)
  )

print(p_panel)

ggsave(
  "Figures_Final/ASF_controls_comparison_panel.png",
  p_panel,
  width = 16,
  height = 7,
  dpi = 300
)

ggsave(
  "Figures_Final/ASF_controls_comparison_panel.pdf",
  p_panel,
  width = 16,
  height = 7,
  dpi = 300
)


# ============================================================
# Spatial maps WITHOUT controls
# ============================================================

message("\n>>> Running SEI without controls ...")

sim_no_controls <- run_seidc_no_controls(
  W_base,
  farms_sf
)

# Times to plot
times_to_plot <- seq(0, 130, by = 10)

# Create output directory if needed
dir.create("Figures_Final", showWarnings = FALSE)

for (tt in times_to_plot) {
  
  p_map <- plot_map_with_zones(
    sim_no_controls,
    farms_sf,
    den_25832,
    tt,
    STYLE
  ) +
    labs(
      title = paste0(
        "ASF simulation without controls (Day ",
        tt,
        ")"
      )
    )
  
  p_map <- add_scalebar_northarrow(
    p_map,
    den_25832,
    length_km = 50,
    pos = "br",
    arrow = TRUE
  )
  
  print(p_map)
  
  ggsave(
    filename = file.path(
      "Figures_Final",
      sprintf("asf_no_controls_%03d.png", tt)
    ),
    plot = p_map,
    width = 12,
    height = 7,
    dpi = 300
  )
}

message("Finished exporting no-control maps.")

# Animation
animate_map_sequence(
  sim_no_controls,
  farms_sf,
  den_25832,
  STYLE,
  by = 10,
  fps = 6,
  out_gif = "Figures_Final/asf_no_controls.gif",
  out_mp4 = "Figures_Final/asf_no_controls.mp4"
)



