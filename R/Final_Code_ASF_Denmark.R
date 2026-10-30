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


# =========================================================
# SCENARIO COMPARISON SUITE
# =========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(tidyr)
  library(purrr)
})

# ---------- 1) Per-run metrics from a single simulation ----------
#   - peak I timing & height
#   - final C at horizon
#   - convex-hull area slope (km²/day) as spatial spread speed proxy
#
# Notes:
#   * times_grid_hull is a sparse grid (e.g., seq(0, 140, by = 5)) to speed hull computation.
#   * If there are no detections, spatial slope = 0.
#
metrics_from_sim <- function(sim, farms_sf,
                             times_grid_hull = seq(0, 140, by = 5),
                             scale_factor = 1) {
  # Time series (wide)
  ts_df <- compute_ts(sim, scale_factor = scale_factor)  # columns: time,S,E,I,D,C
  
  # Peak I
  i_row <- ts_df[which.max(ts_df$I), , drop = FALSE]
  t_peak_I <- i_row$time
  I_peak    <- i_row$I
  
  # Final C at horizon
  final_C <- ts_df$C[nrow(ts_df)]
  
  # Convex hull area (km^2) through time (detections)
  hull_area_vs_time <- function(sim, farms_sf, times) {
    sapply(times, function(tt) {
      ids <- which(!is.na(sim$detect_time) & sim$detect_time <= tt)
      if (!length(ids)) return(0)
      hull <- sf::st_convex_hull(sf::st_union(sf::st_geometry(farms_sf[ids, ])))
      as.numeric(sf::st_area(hull)) / 1e6  # m^2 -> km^2 (EPSG:25832)
    })
  }
  A_km2 <- hull_area_vs_time(sim, farms_sf, times_grid_hull)
  
  # Fit a simple slope on the portion with positive area (avoid zeros)
  pos <- which(A_km2 > 0)
  spread_slope <- if (length(pos) >= 3) {
    # fit on the last 80% of positive window to catch stable expansion
    idx <- pos[ceiling(0.2*length(pos)):length(pos)]
    if (length(idx) >= 3) {
      as.numeric(coef(lm(A_km2[idx] ~ times_grid_hull[idx]))[2])  # km^2 per day
    } else 0
  } else 0
  
  # First detection
  first_det <- if (all(is.na(sim$detect_time))) NA_real_ else min(sim$detect_time, na.rm = TRUE)
  
  data.frame(
    t_peak_I = t_peak_I,
    I_peak   = I_peak,
    C_final  = final_C,
    spread_slope_km2_day = spread_slope,
    first_det = first_det
  )
}

# ---------- 2) Run a scenario (batch) and extract metrics ----------
# Pass any model argument for run_seidc_controlled(...) via ...
run_scenario_batch <- function(name,
                               n_runs,
                               W_base, farms_sf,
                               times_grid_counts = seq(0, 140, by = 1),
                               times_grid_hull   = seq(0, 140, by = 5),
                               scale_factor = 1,
                               seed0 = 1001,
                               ...) {
  set.seed(seed0)
  
  # Reuse the existing "counts on grid" batch runner for ensemble plots if desired
  batch_df <- run_batch_extract_counts(
    n_runs      = n_runs,
    times_grid  = times_grid_counts,
    W_base      = W_base,
    farms_sf    = farms_sf,
    scale_factor = scale_factor,
    parallel     = FALSE,   # sf is not fork-safe
    seed0        = seed0,
    ...
  )
  
  # Per-run scalar metrics (simulate each run again to compute spatial slope + peaks)
  # Re-simulate with the same seeds to align with batch_df runs
  run_metrics <- vector("list", n_runs)
  for (r in seq_len(n_runs)) {
    sim_r <- run_seidc_controlled(W_base, farms_sf, rng_seed = seed0 + r, ...)
    m     <- metrics_from_sim(sim_r, farms_sf, times_grid_hull = times_grid_hull, scale_factor = scale_factor)
    m$run <- r
    run_metrics[[r]] <- m
  }
  metrics_df <- do.call(rbind, run_metrics)
  metrics_df$scenario <- name
  
  list(batch_df = batch_df, metrics_df = metrics_df)
}

# ---------- 3) Summaries across runs for each scenario ----------
# Produces tidy table with medians & IQR for key metrics.
summarise_metrics <- function(all_metrics_df) {
  all_metrics_df %>%
    group_by(scenario) %>%
    summarise(
      C_final_med = median(C_final, na.rm = TRUE),
      C_final_IQR = IQR(C_final, na.rm = TRUE),
      t_peakI_med = median(t_peak_I, na.rm = TRUE),
      I_peak_med  = median(I_peak, na.rm = TRUE),
      spread_slope_med = median(spread_slope_km2_day, na.rm = TRUE),
      spread_slope_IQR = IQR(spread_slope_km2_day, na.rm = TRUE),
      first_det_med = median(first_det, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    arrange(scenario)
}

# ---------- 4) Quick plots for comparison ----------
plot_finalC_by_scenario <- function(metrics_df, style = STYLE) {
  ggplot(metrics_df, aes(x = scenario, y = C_final, fill = scenario)) +
    geom_boxplot(outlier.alpha = 0.3) +
    scale_fill_manual(values = rep("#595959", length(unique(metrics_df$scenario))), guide = "none") +
    labs(title = "Final C distribution by scenario", x = NULL, y = "C (culled/removed)") +
    theme_minimal(base_size = style$base_size) +
    theme(panel.grid.major = element_line(color = style$grid_major),
          axis.text.x = element_text(angle = 15, hjust = 1))
}

plot_t_peakI_by_scenario <- function(metrics_df, style = STYLE) {
  ggplot(metrics_df, aes(x = t_peak_I, fill = scenario)) +
    geom_histogram(position = "identity", alpha = 0.35, bins = 30, color = "white") +
    labs(title = "Peak I timing (histograms)", x = "Time (days)", y = "Runs") +
    theme_minimal(base_size = style$base_size) +
    theme(panel.grid.major = element_line(color = style$grid_major))
}

plot_spread_slope_by_scenario <- function(metrics_df, style = STYLE) {
  ggplot(metrics_df, aes(x = scenario, y = spread_slope_km2_day, fill = scenario)) +
    geom_boxplot(outlier.alpha = 0.3) +
    scale_fill_manual(values = rep("#0072B2", length(unique(metrics_df$scenario))), guide = "none") +
    labs(title = "Spatial spread speed (km²/day) by scenario", x = NULL, y = "km² / day (convex hull slope)") +
    theme_minimal(base_size = style$base_size) +
    theme(panel.grid.major = element_line(color = style$grid_major),
          axis.text.x = element_text(angle = 15, hjust = 1))
}


# ============= SCENARIO MATRIX =============
scenarios <- tibble::tribble(
  ~name,            ~delta_0, ~prot_radius_m, ~surv_radius_m, ~beta_m_post_scale,
  "baseline",        1/6,      3000,           10000,          0.25,
  "faster_det",      1/4.5,    3000,           10000,          0.25,
  "tighter_radii",   1/6,      5000,           15000,          0.25,
  "stronger_beta",   1/6,      3000,           10000,          0.15
)

# Controls & grids for runs
n_runs_cmp         <- 20
times_grid_counts  <- seq(0, 140, by = 1)
times_grid_hull    <- seq(0, 140, by = 5)
scale_factor_batch <- 10000 / N   # keep consistent with your ensemble plots

# Run each scenario (sequentially to be sf-safe)
all_metrics <- vector("list", nrow(scenarios))
all_batches <- vector("list", nrow(scenarios))

for (i in seq_len(nrow(scenarios))) {
  scn <- scenarios[i, ]
  cat("Running scenario:", scn$name, "...\n")
  res <- run_scenario_batch(
    name  = scn$name,
    n_runs = n_runs_cmp,
    W_base = W_base, farms_sf = farms_sf,
    times_grid_counts = times_grid_counts,
    times_grid_hull   = times_grid_hull,
    scale_factor = scale_factor_batch,
    # forward the differences:
    delta_0 = scn$delta_0,
    prot_radius_m = scn$prot_radius_m,
    surv_radius_m = scn$surv_radius_m,
    beta_m_post_scale = scn$beta_m_post_scale
  )
  all_metrics[[i]] <- res$metrics_df
  all_batches[[i]] <- res$batch_df  # (optional) reuse for ci plots per scenario
}

metrics_df_all <- dplyr::bind_rows(all_metrics)

# ---- Summary table (medians & IQRs) ----
summary_tbl <- summarise_metrics(metrics_df_all)
print(summary_tbl)

# Save a CSV for your report
write.csv(summary_tbl, "Figures_Final/scenario_summary_metrics.csv", row.names = FALSE)

# ---- Quick comparison plots ----
pC  <- plot_finalC_by_scenario(metrics_df_all, STYLE)
pT  <- plot_t_peakI_by_scenario(metrics_df_all, STYLE)
pSp <- plot_spread_slope_by_scenario(metrics_df_all, STYLE)

print(pC);  ggsave("Figures_Final/scen_finalC_box.png", pC, width = 8, height = 5, dpi = 300)
print(pT);  ggsave("Figures_Final/scen_tpeak_hist.png", pT, width = 8, height = 5, dpi = 300)
print(pSp); ggsave("Figures_Final/scen_spread_slope_box.png", pSp, width = 8, height = 5, dpi = 300)

#----------------------------------------------------------------
# =========================================================
# Global Sensitivity: LHS + (manual) PRCC with bootstrap CIs
# =========================================================

suppressPackageStartupMessages({
  library(lhs)
  library(future)
  library(future.apply)
  library(dplyr)
  library(ggplot2)
  library(tibble)
  library(sf)
})

# ---------- Sanitizers ----------
sanitize_design <- function(X) {
  X_df <- as.data.frame(X, stringsAsFactors = FALSE)
  colnames(X_df) <- make.names(colnames(X_df), unique = TRUE)
  
  # Flatten lists and coerce to numeric
  X_df[] <- lapply(X_df, function(col) {
    if (is.list(col)) {
      lens <- vapply(col, length, integer(1))
      if (!all(lens == 1L)) stop("sanitize_design(): list column with varying element lengths.")
      col <- vapply(col, function(z) if (is.null(z) || length(z)==0) NA_real_ else as.numeric(z[[1]]), numeric(1))
    }
    if (inherits(col, "units")) col <- as.numeric(col)
    if (is.factor(col)) col <- as.character(col)
    suppressWarnings(as.numeric(col))
  })
  
  # Drop zero-variance columns
  keep <- vapply(X_df, function(v) sd(v, na.rm = TRUE) > 0, logical(1))
  if (!all(keep)) message("Dropping zero-variance columns: ", paste(names(X_df)[!keep], collapse = ", "))
  X_df <- X_df[, keep, drop = FALSE]
  
  # Keep only finite rows
  ok <- apply(X_df, 1L, function(r) all(is.finite(r)))
  if (!all(ok)) message("Dropping ", sum(!ok), " rows with non-finite values in X.")
  X_df[ok, , drop = FALSE]
}

sanitize_y <- function(y) {
  y_num <- if (is.list(y)) unlist(y, recursive = TRUE, use.names = FALSE) else y
  suppressWarnings(as.numeric(y_num))
}

# ---------- Manual PRCC (partial Spearman) + Bootstrap ----------
partial_spearman_prcc <- function(y, X) {
  X_df <- sanitize_design(X)
  y    <- sanitize_y(y)
  stopifnot(nrow(X_df) == length(y))
  
  # Rank-transform
  Yr <- rank(y, ties.method = "average")
  Xr <- as.data.frame(lapply(X_df, function(v) rank(v, ties.method = "average")))
  
  p <- ncol(Xr)
  prcc <- numeric(p); names(prcc) <- colnames(Xr)
  for (j in seq_len(p)) {
    others <- setdiff(seq_len(p), j)
    eY  <- resid(lm(Yr ~ as.matrix(Xr[, others, drop = FALSE])))
    eXj <- resid(lm(Xr[[j]] ~ as.matrix(Xr[, others, drop = FALSE])))
    prcc[j] <- suppressWarnings(cor(eY, eXj, method = "pearson"))
  }
  prcc[order(-abs(prcc))]
}

bootstrap_prcc <- function(y, X, B = 500, seed = 1) {
  set.seed(seed)
  X_df <- sanitize_design(X)
  y    <- sanitize_y(y)
  
  stopifnot(nrow(X_df) == length(y))
  ok <- is.finite(y) & apply(X_df, 1L, function(r) all(is.finite(r)))
  X_df <- X_df[ok, , drop = FALSE]
  y    <- y[ok]
  n    <- nrow(X_df)
  
  base <- partial_spearman_prcc(y, X_df)
  params <- names(base); P <- length(params)
  mat <- matrix(NA_real_, nrow = B, ncol = P, dimnames = list(NULL, params))
  
  for (b in seq_len(B)) {
    idx <- sample.int(n, n, replace = TRUE)
    mat[b, ] <- partial_spearman_prcc(y[idx], X_df[idx, , drop = FALSE])[params]
  }
  
  ci <- t(apply(mat, 2, function(z) stats::quantile(z, probs = c(0.025, 0.975), na.rm = TRUE)))
  out <- data.frame(
    parameter = params,
    prcc      = as.numeric(base[params]),
    low       = ci[, 1],
    high      = ci[, 2],
    row.names = NULL
  )
  out[order(-abs(out$prcc)), , drop = FALSE]
}

plot_prcc <- function(df, title = "PRCC (partial Spearman)") {
  ggplot(df, aes(x = reorder(parameter, prcc), y = prcc)) +
    geom_hline(yintercept = 0, color = "grey60") +
    geom_point(size = 2) +
    geom_errorbar(aes(ymin = low, ymax = high), width = 0.2) +
    coord_flip() +
    labs(title = title, x = NULL, y = "PRCC (rank-based)") +
    theme_minimal(base_size = STYLE$base_size) +
    theme(panel.grid.major = element_line(color = STYLE$grid_major))
}

# ---------- Metrics from a single simulation ----------
metrics_from_sim <- function(sim, farms_sf,
                             times_grid_hull = seq(0, 140, by = 5),
                             scale_factor = 1) {
  ts_df <- compute_ts(sim, scale_factor = scale_factor)  # requires your compute_ts()
  
  # Peak I
  i_row <- ts_df[which.max(ts_df$I), , drop = FALSE]
  t_peak_I <- i_row$time
  I_peak   <- i_row$I
  
  # Final C
  C_final  <- ts_df$C[nrow(ts_df)]
  
  # Convex hull spread speed (km²/day)
  hull_area_vs_time <- function(sim, farms_sf, times) {
    sapply(times, function(tt) {
      ids <- which(!is.na(sim$detect_time) & sim$detect_time <= tt)
      if (!length(ids)) return(0)
      hull <- st_convex_hull(st_union(st_geometry(farms_sf[ids, ])))
      as.numeric(st_area(hull)) / 1e6
    })
  }
  A <- hull_area_vs_time(sim, farms_sf, times_grid_hull)
  pos <- which(A > 0)
  spread_slope <- if (length(pos) >= 3) {
    idx <- pos[ceiling(0.2 * length(pos)):length(pos)]
    if (length(idx) >= 3) as.numeric(coef(lm(A[idx] ~ times_grid_hull[idx]))[2]) else 0
  } else 0
  
  data.frame(
    t_peak_I = t_peak_I,
    I_peak   = I_peak,
    C_final  = C_final,
    spread_slope_km2_day = spread_slope,
    row.names = NULL
  )
}

# ---------- Main: run_sensitivity ----------
# Arguments:
# - W_base, farms_sf: your existing objects
# - param_ranges: tibble/data.frame with columns: name, min, max
# - n_samples: LHS samples; B: bootstrap for PRCC CIs
# - t_max, snapshot_dt: simulator setup; scale_factor for output scaling
# - workers: #processes for multisession; outputs saved to Figures/ if save=TRUE
run_sensitivity <- function(W_base, farms_sf,
                            param_ranges = tibble::tribble(
                              ~name,              ~min,    ~max,
                              "beta_m",            0.15,    0.30,
                              "beta_e",            0.005,   0.030,
                              "sigma",             1/6,     1/3,
                              "delta_0",           1/8,     1/4,
                              "kappa_cull",        1/7,     1/3,
                              "beta_m_post_scale", 0.15,    0.50,
                              "beta_e_post_scale", 0.30,    0.80,
                              "prot_radius_m",     2000,    6000,
                              "surv_radius_m",     8000,   16000
                            ),
                            n_samples   = 200,
                            B           = 500,        # bootstrap reps
                            t_max       = 140,
                            snapshot_dt = 0.5,
                            times_grid_hull = seq(0, 140, by = 5),
                            scale_factor = 10000 / nrow(farms_sf),
                            workers     = max(1L, parallel::detectCores() - 1L),
                            seed        = 123,
                            save        = TRUE,
                            out_dir     = "Figures_Final") {
  
  dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
  
  # 1) LHS design
  set.seed(seed)
  U    <- lhs::randomLHS(n_samples, nrow(param_ranges))
  mins <- param_ranges$min
  maxs <- param_ranges$max
  X    <- sweep(U, 2, (maxs - mins), `*`)
  X    <- sweep(X, 2, mins, `+`)
  colnames(X) <- param_ranges$name
  X_df <- as.data.frame(X)
  
  # 2) Parallel plan (sf-safe)
  plan(multisession, workers = workers)
  
  # 3) Run simulations (one per sample)
  run_one <- function(i) {
    pars <- as.list(X_df[i, , drop = FALSE])
    args <- c(list(
      W_base = W_base, farms_sf = farms_sf,
      t_max = t_max, snapshot_dt = snapshot_dt,
      rng_seed = 10000 + i
    ), pars)
    sim_i <- do.call(run_seidc_controlled, args)  # your simulator
    m     <- metrics_from_sim(sim_i, farms_sf, times_grid_hull = times_grid_hull, scale_factor = scale_factor)
    m$sample <- i
    m
  }
  
  metrics_list <- future.apply::future_lapply(seq_len(n_samples), run_one, future.seed = TRUE)
  plan(sequential)
  
  metrics_df <- bind_rows(metrics_list)
  # Re-align design to metrics row order (sample ascending)
  metrics_df <- metrics_df[order(metrics_df$sample), ]
  X_df       <- X_df[metrics_df$sample, , drop = FALSE]
  
  # 4) Outcomes for PRCC
  y_Cfinal <- metrics_df$C_final
  y_tpeak  <- metrics_df$t_peak_I
  y_Ipeak  <- metrics_df$I_peak
  y_spread <- metrics_df$spread_slope_km2_day
  
  # 5) PRCC (manual partial Spearman) with bootstrap CIs
  prcc_Cfinal <- bootstrap_prcc(y_Cfinal, X_df, B = B, seed = seed + 1)
  prcc_tpeak  <- bootstrap_prcc(y_tpeak,  X_df, B = B, seed = seed + 2)
  prcc_Ipeak  <- bootstrap_prcc(y_Ipeak,  X_df, B = B, seed = seed + 3)
  prcc_spread <- bootstrap_prcc(y_spread, X_df, B = B, seed = seed + 4)
  
  # 6) Plots
  pC <- plot_prcc(prcc_Cfinal, "PRCC — Final C")
  pT <- plot_prcc(prcc_tpeak,  "PRCC — Time to peak I")
  pI <- plot_prcc(prcc_Ipeak,  "PRCC — Peak I size")
  pS <- plot_prcc(prcc_spread, "PRCC — Spread speed (km²/day)")
  
  if (save) {
    readr::write_csv(bind_cols(sample = metrics_df$sample, X_df, metrics_df),
                     file.path(out_dir, "sensitivity_samples_and_metrics.csv"))
    readr::write_csv(prcc_Cfinal, file.path(out_dir, "prcc_finalC.csv"))
    readr::write_csv(prcc_tpeak,  file.path(out_dir, "prcc_tpeakI.csv"))
    readr::write_csv(prcc_Ipeak,  file.path(out_dir, "prcc_Ipeak.csv"))
    readr::write_csv(prcc_spread, file.path(out_dir, "prcc_spread.csv"))
    
    ggsave(file.path(out_dir, "prcc_finalC.png"), pC, width = 7, height = 5, dpi = 300)
    ggsave(file.path(out_dir, "prcc_tpeakI.png"), pT, width = 7, height = 5, dpi = 300)
    ggsave(file.path(out_dir, "prcc_Ipeak.png"),  pI, width = 7, height = 5, dpi = 300)
    ggsave(file.path(out_dir, "prcc_spread.png"), pS, width = 7, height = 5, dpi = 300)
  }
  
  list(
    design  = X_df,
    metrics = metrics_df,
    prcc = list(C_final = prcc_Cfinal,
                t_peakI = prcc_tpeak,
                I_peak  = prcc_Ipeak,
                spread  = prcc_spread),
    plots = list(C_final = pC, t_peakI = pT, I_peak = pI, spread = pS)
  )
}

sens <- run_sensitivity(
  W_base, farms_sf,
  n_samples   = 200,     # start with 100–200 (minutes)
  B           = 500,     # bootstrap reps for CIs; bump to 1000 later
  t_max       = 140,
  snapshot_dt = 0.5,
  workers     = max(1, parallel::detectCores() - 1),
  save        = TRUE,
  out_dir     = "Figures_Final"
)

# Inspect top drivers for final C:
sens$prcc$C_final

# Show the PRCC plot:
print(sens$plots$C_final)


#---------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr); library(tibble); library(glue); library(readr); library(tidyr)
})

# Map parameter codes to descriptive labels (edit as you wish)
param_dictionary <- tribble(
  ~name,               ~Parameter, ~Description,
  "beta_m",            "β_m",      "Contact transmission",
  "beta_e",            "β_e",      "Environmental transmission",
  "sigma",             "σ",        "Incubation rate",
  "delta_0",           "δ_0",      "Detection rate",
  "kappa_cull",        "κ",        "Removal/clearance rate",
  "omega",             "ω",        "Environmental shedding",
  "gamma",             "γ",        "Environmental decay",
  "rho_pre",           "ρ",        "Detected-farm transmission fraction",   # pre-detection fraction
  "beta_m_post_scale", "β_m(post)", "Post-detection contact multiplier",
  "beta_e_post_scale", "β_e(post)", "Post-detection environmental multiplier",
  "prot_radius_m",     "R_prot",   "Protection radius (m)",
  "surv_radius_m",     "R_surv",   "Surveillance radius (m)"
)

# Build the table
# - baseline: named numeric vector of YOUR central values (the ones used in your sim)
# - prior_ranges: tibble with columns name, min, max (lower/upper limits you want to show)  OR
# - empirical: pass sens$design to use empirical 2.5%–97.5% limits instead of prior min/max
# - digits: number of decimals for formatting
make_parameter_table <- function(baseline, prior_ranges = NULL, empirical = NULL, digits = 3) {
  
  # Turn baseline into a tibble
  base_tbl <- enframe(baseline, name = "name", value = "value")
  
  # Decide where the limits come from
  if (!is.null(empirical)) {
    # Empirical limits from samples: 2.5%–97.5%
    lim_tbl <- empirical |>
      summarize(across(everything(),
                       list(lo = ~quantile(.x, 0.025, na.rm = TRUE),
                            hi = ~quantile(.x, 0.975, na.rm = TRUE)),
                       .names = "{.col}_{.fn}")) |>
      pivot_longer(everything(),
                   names_to = c("name","bound"),
                   names_sep = "_",
                   values_to = "val") |>
      pivot_wider(names_from = bound, values_from = val)
  } else if (!is.null(prior_ranges)) {
    # Prior limits (min/max)
    lim_tbl <- prior_ranges |>
      transmute(name, lo = .data$min, hi = .data$max)
  } else {
    stop("Provide either empirical = sens$design (for 2.5%–97.5%) or prior_ranges (for min–max).")
  }
  
  # Merge and format
  out <- base_tbl |>
    left_join(lim_tbl, by = "name") |>
    left_join(param_dictionary, by = "name") |>
    mutate(
      value_fmt = formatC(value, format = "f", digits = digits),
      lo_fmt    = formatC(lo,    format = "f", digits = digits),
      hi_fmt    = formatC(hi,    format = "f", digits = digits),
      `Values (lower and upper limits)` = glue("{value_fmt} [{lo_fmt} – {hi_fmt}]")
    ) |>
    select(Parameter, Description, `Values (lower and upper limits)`) |>
    arrange(factor(Parameter, levels = unique(param_dictionary$Parameter)))
  
  out
}


# Your simulator defaults (baseline central values) — edit if yours differ
baseline <- c(
  beta_m = 0.22,
  beta_e = 0.015,
  sigma = 1/4,        # 0.25
  delta_0 = 1/6,      # 0.167
  kappa_cull = 1/4,   # 0.25
  omega = 0.10,
  gamma = 0.30,
  rho_pre = 0.05,
  beta_m_post_scale = 0.25,
  beta_e_post_scale = 0.50,
  prot_radius_m = 3000,
  surv_radius_m = 10000
)

# The prior bounds you used for LHS (example — edit if you used different)
param_ranges <- tibble::tribble(
  ~name,              ~min,    ~max,
  "beta_m",            0.052,   0.149,   # <-- REPLACE with your desired lower/upper
  "beta_e",            0.004,   0.044,
  "sigma",             0.088,   0.517,
  "delta_0",           0.091,   0.449,
  "kappa_cull",        0.086,   0.399,
  "omega",             0.037,   0.234,
  "gamma",             0.134,   0.805,
  "rho_pre",           0.001,   0.104,
  "beta_m_post_scale", 0.15,    0.50,
  "beta_e_post_scale", 0.30,    0.80,
  "prot_radius_m",     2000,    6000,
  "surv_radius_m",     8000,   16000
)

tab_prior <- make_parameter_table(
  baseline = baseline,
  prior_ranges = param_ranges,   # <= show min–max
  empirical = NULL,
  digits = 3
)
tab_prior
readr::write_csv(tab_prior, "Figures_Final/parameter_table_prior_limits.csv")

#----------------------

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(tibble); library(glue); library(readr)
})

# -- empirical quantiles per column, no name-splitting (safe with underscores) --
compute_empirical_limits <- function(design_df, candidates, probs = c(0.025, 0.975)) {
  cols <- intersect(candidates, colnames(design_df))
  if (!length(cols)) return(tibble(name = character(), lo = numeric(), hi = numeric()))
  
  rows <- lapply(cols, function(nm) {
    v <- suppressWarnings(as.numeric(design_df[[nm]]))
    v <- v[is.finite(v)]
    if (length(v) == 0) {
      tibble(name = nm, lo = NA_real_, hi = NA_real_)
    } else {
      qs <- stats::quantile(v, probs = probs, names = FALSE, na.rm = TRUE)
      if (length(qs) == 1L) qs <- rep(qs, 2L)  # constant vector safety
      tibble(name = nm, lo = qs[1], hi = qs[2])
    }
  })
  bind_rows(rows)
}

# -- blend empirical with prior fallback; pretty-format; units-aware digits --
make_parameter_table_empirical_or_prior <- function(baseline,
                                                    design_df,
                                                    dictionary,
                                                    prior_ranges,         # tibble: name, min, max
                                                    digits_rates = 3,
                                                    digits_radius = 0) {
  
  base_tbl <- tibble(name = names(baseline), value = as.numeric(baseline))
  
  # 1) empirical limits for columns present
  emp_tbl <- compute_empirical_limits(design_df, candidates = base_tbl$name)
  
  # 2) prior bounds (min/max)
  prior_tbl <- prior_ranges %>%
    transmute(name, lo_prior = .data$min, hi_prior = .data$max)
  
  # 3) merge and fallback: empirical if present, else prior
  merged <- base_tbl %>%
    left_join(emp_tbl,   by = "name") %>%
    left_join(prior_tbl, by = "name") %>%
    mutate(
      lo_use = dplyr::coalesce(lo, lo_prior),
      hi_use = dplyr::coalesce(hi, hi_prior)
    ) %>%
    left_join(dictionary, by = "name")
  
  # 4) format with units-aware precision
  radius_names <- c("prot_radius_m","surv_radius_m")
  
  fmt_val <- function(nm, x) {
    if (nm %in% radius_names) formatC(x, format = "f", digits = digits_radius)
    else                      formatC(x, format = "f", digits = digits_rates)
  }
  
  out <- merged %>%
    mutate(
      value_fmt = mapply(fmt_val, name, value),
      lo_fmt    = mapply(fmt_val, name, lo_use),
      hi_fmt    = mapply(fmt_val, name, hi_use),
      `Values (lower and upper limits)` = glue("{value_fmt} [{lo_fmt} – {hi_fmt}]")
    ) %>%
    select(Parameter, Description, `Values (lower and upper limits)`) %>%
    arrange(match(Parameter, dictionary$Parameter))
  
  out
}

# Your baseline values (edit κ baseline to 0.125 or 0.25 as you prefer)
baseline <- c(
  beta_m = 0.22, beta_e = 0.015,
  sigma = 0.30, delta_0 = 0.167, kappa_cull = 0.25,
  omega = 0.10, gamma = 0.30, rho_pre = 0.05,
  beta_m_post_scale = 0.25, beta_e_post_scale = 0.50,
  prot_radius_m = 3000, surv_radius_m = 10000
)

# The prior limits you want to show as fallback for missing empirical columns
param_ranges <- tibble::tribble(
  ~name,              ~min,    ~max,
  "beta_m",            0.052,   0.149,
  "beta_e",            0.004,   0.044,
  "sigma",             0.088,   0.517,
  "delta_0",           0.091,   0.449,
  "kappa_cull",        0.086,   0.399,
  "omega",             0.037,   0.234,
  "gamma",             0.134,   0.805,
  "rho_pre",           0.001,   0.104,
  "beta_m_post_scale", 0.15,    0.50,
  "beta_e_post_scale", 0.30,    0.80,
  "prot_radius_m",     2000,    6000,
  "surv_radius_m",     8000,   16000
)

# Your dictionary from earlier (param_dictionary)
# ... (assumed already defined in your session)

# Build the blended table
tab_emp_blended <- make_parameter_table_empirical_or_prior(
  baseline   = baseline,
  design_df  = sens$design,        # empirical source
  dictionary = param_dictionary,   # pretty names/descriptions
  prior_ranges = param_ranges,     # fallback min/max for missing params
  digits_rates  = 3,               # decimals for rates/fractions
  digits_radius = 0                # integers for meter radii
)

tab_emp_blended
readr::write_csv(tab_emp_blended, "Figures_Final/parameter_table_empirical_limits_blended.csv")

#----#

df <- readr::read_csv("Figures_Final/sensitivity_samples_and_metrics.csv", show_col_types = FALSE)
names(df)

# Known outcome columns (we'll pick the ones that exist)
known_metrics <- c("t_peak_I","I_peak","C_final","spread_slope_km2_day",
                   "peak_I","final_size","final_leakage_rate")

metric_cols <- intersect(names(df), known_metrics)
param_cols  <- setdiff(names(df), metric_cols)

# Reconstruct sens object
sens <- list(
  design  = as.data.frame(df[, param_cols, drop = FALSE]),
  metrics = as.data.frame(df[, metric_cols, drop = FALSE])
)

str(sens$design)
str(sens$metrics)

library(dplyr)
library(tidyr)
library(ggplot2)

# assuming you have a 'sens' list with sens$design (parameters) and sens$metrics (outcomes)
design  <- sens$design
metrics <- sens$metrics

# Compute Spearman correlations: outcomes × parameters
rho_mat <- sapply(metrics, function(y) {
  sapply(design, function(x) suppressWarnings(cor(x, y, method = "spearman", use = "pairwise.complete.obs")))
})

rho_long <- as.data.frame(rho_mat) |>
  tibble::rownames_to_column(var = "Outcome") |>
  pivot_longer(-Outcome, names_to = "Parameter", values_to = "rho")

# Optional: order axes
# rho_long$Outcome   <- factor(rho_long$Outcome,   levels = c("peak_I","final_size","final_leakage_rate"))
# rho_long$Parameter <- factor(rho_long$Parameter, levels = c("beta_m","beta_e","rho","sigma","delta_0","kappa","omega","gamma","eta_val","K"))

# Build the plot
p <- ggplot(rho_long, aes(x = Parameter, y = Outcome, fill = rho)) +
  geom_tile(color = "white") +
  scale_fill_gradient2(
    low = "#D55E00", mid = "white", high = "#0072B2",
    midpoint = 0, name = "Spearman rho"
  ) +
  labs(title = "Global Sensitivity (Spearman correlations)", x = NULL, y = NULL) +
  theme_minimal(base_size = 14) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid = element_blank())

# PNG
ggsave("Figures_Final/sensitivity_spearman_heatmap.png",
       plot = p, width = 10, height = 4, dpi = 300)

# PDF (vector, great for publications)
ggsave("Figures_Final/sensitivity_spearman_heatmap.pdf",
       plot = p, width = 10, height = 4, device = cairo_pdf)

#------------------------------------------------------
# Load required packages
library(tidyverse)
library(reshape2)
library(pheatmap)

#---------------------------------------------------------
# 0. Create output folder if it does not already exist
#---------------------------------------------------------
outdir <- "Figures_Final"
if (!dir.exists(outdir)) dir.create(outdir)

# 1. Load the data
df <- read_csv("Figures_Final/sensitivity_samples_and_metrics.csv")
#View(sensitivity_samples_and_metrics)

# 2. Remove the unwanted sample columns
df_clean <- df %>% 
  select(-starts_with("sample"))

# 3. Define model outputs
outputs <- c("C_final", "I_peak", "spread_slope_km2_day", "t_peak_I")

# 4. Select parameter columns
param_cols <- setdiff(colnames(df_clean), outputs)

# 5. Compute Spearman correlation matrix (parameters × outputs)
corr_matrix <- df_clean %>%
  select(all_of(c(param_cols, outputs))) %>%
  cor(method = "spearman") %>%
  as.data.frame() %>%
  rownames_to_column("parameter") %>%
  filter(parameter %in% param_cols) %>%
  column_to_rownames("parameter") %>%
  select(all_of(outputs))

# 6. Save correlation matrix to CSV (in Figures_Final)
write.csv(corr_matrix,
          file.path(outdir, "spearman_param_to_outcomes.csv"),
          row.names = TRUE)

# 7. Plot and save heatmap (PNG + PDF in Figures_Final)
pheatmap(corr_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize = 10,
         fontsize_row = 10,
         fontsize_col = 10,
         main = "Global Sensitivity of Model Outputs to Parameters (Spearman ρ)",
         filename = file.path(outdir, "sensitivity_spearman_heatmap_clean.png"),
         width = 8, height = 6)

pheatmap(corr_matrix,
         color = colorRampPalette(c("blue", "white", "red"))(50),
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize = 10,
         fontsize_row = 10,
         fontsize_col = 10,
         main = "Global Sensitivity of Model Outputs to Parameters (Spearman ρ)",
         filename = file.path(outdir, "sensitivity_spearman_heatmap_clean.pdf"),
         width = 8, height = 6)
#-----------------------------------------------


# Packages
library(tidyverse)
library(pheatmap)

# ------------------------------------------------------------
# Paths
# ------------------------------------------------------------
infile <- "sensitivity_samples_and_metrics.csv"
outdir <- "Figures_Final"
if (!dir.exists(outdir)) dir.create(outdir, recursive = TRUE)

# ------------------------------------------------------------
# Load and clean
# ------------------------------------------------------------
df <- read_csv(infile, show_col_types = FALSE)

# Drop artefact index columns (e.g., sample...1, sample...15, etc.)
df_clean <- df %>% select(-starts_with("sample"))

# Define outputs and derive parameter set
outputs <- c("C_final", "I_peak", "spread_slope_km2_day", "t_peak_I")
outputs <- outputs[outputs %in% names(df_clean)]
param_cols <- setdiff(names(df_clean), outputs)

# ------------------------------------------------------------
# Spearman correlations: parameters (rows) × outputs (cols)
# ------------------------------------------------------------
corr_mat <- df_clean %>%
  select(all_of(c(param_cols, outputs))) %>%
  cor(method = "spearman") %>%
  as.data.frame() %>%
  rownames_to_column("name") %>%
  filter(name %in% param_cols) %>%
  column_to_rownames("name") %>%
  select(all_of(outputs))

# Save coefficients
write.csv(corr_mat,
          file.path(outdir, "spearman_param_to_outcomes.csv"),
          row.names = TRUE)

# ------------------------------------------------------------
# Clustered heatmap
# ------------------------------------------------------------
# Palette: blue-white-red with centered 0
pal <- colorRampPalette(c("#2166AC", "white", "#B2182B"))(101)

# Main plot (PNG)
pheatmap(corr_mat,
         color = pal,
         breaks = seq(-1, 1, length.out = 102),
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize = 9,
         fontsize_row = 9,
         fontsize_col = 9,
         legend_labels = "Spearman ρ",
         main = "Global Sensitivity of Model Outputs to Parameters (Spearman ρ)",
         filename = file.path(outdir, "sensitivity_spearman_heatmap_clean.png"),
         width = 7.5, height = 5.5)

# Also save as vector PDF
pheatmap(corr_mat,
         color = pal,
         breaks = seq(-1, 1, length.out = 102),
         clustering_distance_rows = "euclidean",
         clustering_method = "complete",
         clustering_distance_cols = "euclidean",
         display_numbers = TRUE,
         number_format = "%.2f",
         fontsize = 9,
         fontsize_row = 9,
         fontsize_col = 9,
         legend_labels = "Spearman ρ",
         main = "Global Sensitivity of Model Outputs to Parameters (Spearman ρ)",
         filename = file.path(outdir, "sensitivity_spearman_heatmap_clean.pdf"),
         width = 7.5, height = 5.5)
