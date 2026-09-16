###################################################################################
# Functions for PlanetScope Flowering Phenology 
#
# Workflow:
# 1. Detect and remove outliers in VI and FI time series
# 2. Detect annual green-up phenophases
# 3. Define phenology-based fitting and detection windows
# 4. Smooth VI and FI within the fitting window
# 5. Identify and characterize candidate flowering peaks
# 6. Select the final flowering peak
###################################################################################


#----------------------------------------------------------------------------------
# Outlier detection
#----------------------------------------------------------------------------------

## Detect temporal outliers using a MAD-based second-difference criterion
detect_outliers_mad <- function(x, dates, z, maxd, n_repeat) {
  out <- rep(FALSE, length(x))
  
  for (iter in seq_len(n_repeat)) {
    
    keep <- !out & !is.na(x)
    
    x_sub <- x[keep]
    dates_sub <- dates[keep]
    
    if (length(x_sub) < 3) break
    
    t <- as.numeric(dates_sub)
    
    i1 <- 1:(length(x_sub) - 2)
    i2 <- 2:(length(x_sub) - 1)
    i3 <- 3:length(x_sub)
    
    di <- (x_sub[i2] - x_sub[i1]) - (x_sub[i3] - x_sub[i2])
    
    gap_ok <- (t[i2] - t[i1] <= maxd) &
      (t[i3] - t[i2] <= maxd)
    
    m <- median(di[gap_ok], na.rm = TRUE)
    mad <- median(abs(di[gap_ok] - m), na.rm = TRUE)
    th <- z * mad / 0.6745
    
    out_sub <- rep(FALSE, length(x_sub))
    out_sub[which(gap_ok & (di < m - th | di > m + th)) + 1] <- TRUE
    
    original_idx <- which(keep)
    new_out_idx <- original_idx[out_sub]
    
    if (length(new_out_idx) == 0) break
    
    out[new_out_idx] <- TRUE
  }
  
  out
}

## Detect residual outliers from deviations around a smoothing spline
detect_outliers_spline <- function(date, value, spar,
                                   z_thresh, n_repeat) {
  df_iter <- data.frame(date = date, value = value)
  all_outliers <- as.Date(character())
  
  for (r in seq_len(n_repeat)) {
    valid_idx <- which(!is.na(df_iter$value))
    if (length(valid_idx) < 4) break
    
    fit <- smooth.spline(
      x = as.numeric(df_iter$date[valid_idx]),
      y = df_iter$value[valid_idx],
      spar = spar
    )
    
    pred <- predict(fit, x = as.numeric(df_iter$date))$y
    resid <- df_iter$value - pred
    
    med <- median(resid, na.rm = TRUE)
    s <- mad(resid, center = med, constant = 1.4826, na.rm = TRUE)
    if (!is.finite(s) || s == 0) break
    
    z <- (resid - med) / s
    freq_outliers <- df_iter$date[abs(z) > z_thresh]
    
    if (length(freq_outliers) == 0) break
    
    all_outliers <- unique(c(all_outliers, freq_outliers))
    df_iter <- df_iter[!df_iter$date %in% freq_outliers, ]
    
    if (nrow(df_iter) < 20) break
  }
  
  sort(as.Date(all_outliers))
}


#----------------------------------------------------------------------------------
# Green-up phase detection
#----------------------------------------------------------------------------------

## Detect green-up phenophases from the annual VI time series
detect_greenup <- function(dfvi, dffi, phenopar) {

  spar <- phenopar$AsplineSpar
  thresholds <- phenopar$gup_threshes
  bufferdays <- phenopar$AsplineBuffer
  
  df <- dfvi
  
  ## Detect coincident VI and FI outliers using MAD-based filtering
  vi_ois <- detect_outliers_mad(x=dfvi$value,
                                dates=dfvi$date,
                                z=phenopar$VspikeThresh,
                                maxd=phenopar$maxDistance,
                                n_repeat=phenopar$spike_n_repeat)
  
  fi_ois <- detect_outliers_mad(x=dffi$value,
                                dates=dffi$date,
                                z=phenopar$FspikeThresh,
                                maxd=phenopar$maxDistance,
                                n_repeat=phenopar$spike_n_repeat)

  ois <- vi_ois & fi_ois
  
  df$value_cleaned <- df$value
  if (any(ois)) {
    df$value_cleaned[ois] <- NA
  }
  
  yyyy <- unique(lubridate::year(df$date))
  full_dates <- tibble(
    date = seq(
      as.Date(paste0(yyyy, "-01-01")) - bufferdays,
      as.Date(paste0(yyyy + 1, "-01-01")) + bufferdays - 1,
      by = "1 day"
    )
  )
  df_full <- full_dates %>%
    left_join(df, by = "date") %>%
    arrange(date)

  df_full <- df_full %>%
    mutate(
      across(
        -c(date, value, value_cleaned),
        ~ first(na.omit(.x))
      ),
      doy = yday(date),
      dos= 1:nrow(df_full)
    )
  
  
  if (sum(!is.na(df_full$value_cleaned)) < 5) {
    df_full$value_smooth <- NA_real_
    result <- tibble(stage = c("greenup-onset", "mid-greenup", "maturity"),
                     date = as.Date(NA),
                     doy = NA_integer_,
                     value_smooth = NA_real_)
    return(list(phase = result, VI = df_full, amp = NA_real_, max = NA_real_, min = NA_real_))
  }
  
  ## Estimate and apply the dormant-season baseline
  gs_idxs <- which(month(df_full$date) >= 3 & month(df_full$date) <= 11)
  gs_vals <- df_full$value_cleaned[gs_idxs]
  gs_vals <- gs_vals[gs_vals>0]
  dormant_val <- as.numeric(
    quantile(gs_vals, probs = phenopar$dormantQuantile, na.rm = TRUE)
  )
  df_full$value_cleaned[!is.na(df_full$value_cleaned) & df_full$value_cleaned < dormant_val] <- dormant_val
  df_full$value_cleaned[month(df_full$date) < 3 | month(df_full$date) > 11] <- dormant_val

  ## Smooth the annual VI time series using a smoothing spline
  fit <- smooth.spline(x = df_full$dos[!is.na(df_full$value_cleaned)],
                       y = df_full$value_cleaned[!is.na(df_full$value_cleaned)],
                       spar = spar)
  y_smooth <- predict(fit, x = df_full$dos)$y

  ## Calculate VI amplitude and green-up threshold values
  y_min <- min(y_smooth[gs_idxs], na.rm = TRUE)
  y_max <- max(y_smooth[gs_idxs], na.rm = TRUE)
  amp   <- y_max - y_min
  thresh_vals <- y_min + amp * thresholds
  
  y_smooth <- y_smooth[df_full$date>=min(df$date) & df_full$date<=max(df$date)]
  df_full <- df_full[df_full$date>=min(df$date) & df_full$date<=max(df$date),]
  
  ## Identify green-up dates 
  min_idx <- which.min(y_smooth[1:180])
  rise_idxs <- sapply(thresh_vals, function(th) {
    idx <- which(y_smooth[min_idx:length(y_smooth)] >= th)[1]
    if (!is.na(idx)) idx <- idx + min_idx - 1
    if (is.na(idx)) NA_integer_ else idx
  })
  
  ## Store green-up phenophase dates and corresponding VI values
  result <- tibble(
    stage     = c("greenup-onset", "mid-greenup", "maturity"),
    date      = df_full$date[rise_idxs],
    doy       = df_full$doy[rise_idxs],
    value_smooth  = y_smooth[rise_idxs]
  )
  
  df_full <- df_full %>%
    mutate(value_smooth = y_smooth)
  
  return(list(phase=result, VI=df_full, amp=amp, max=y_max, min=y_min))
  
}


#----------------------------------------------------------------------------------
# Peak, valley, and increase detection
#----------------------------------------------------------------------------------

## Identify local maxima in a smoothed time series
find_peaks <- function(y) {
  n <- length(y)
  peaks <- c()
  for (i in 2:(n - 1)) {
    if (!is.na(y[i - 1]) && !is.na(y[i + 1]) && !is.na(y[i])) {
      if (y[i] > y[i - 1] && y[i] > y[i + 1]) {
        peaks <- c(peaks, i)
      }
    }
  }
  return(peaks)
}

## Identify local minima in a smoothed time series
find_valleys <- function(y) {
  n <- length(y)
  valleys <- c()
  for (i in 2:(n - 1)) {
    if (!is.na(y[i - 1]) && !is.na(y[i + 1]) && !is.na(y[i])) {
      if (y[i] < y[i - 1] && y[i] < y[i + 1]) {
        valleys <- c(valleys, i)
      }
    }
  }
  return(valleys)
}

## Identify the start of a sustained increase in the smoothed time series
find_increase <- function(y, thr, pitw) {
  
  if (length(y) < 5) return(NA_integer_)
  
  diffs <- diff(y)

  posdiff_idx <- which(diffs>0)
  
  if(all(diff(posdiff_idx)<pitw+2)) {
    inc_idx <- posdiff_idx
  } else {
    inc_st_idx <- max(which(diff(posdiff_idx)>=pitw+2))+ 1
    inc_idx <- posdiff_idx[inc_st_idx:length(posdiff_idx)]
  }
  inc_idx <- sort(inc_idx)
  inc_diffs <- diffs[inc_idx]
  thr_diff <- as.numeric(quantile(inc_diffs, probs = thr))

  inc_st_idx <- inc_idx[which(inc_diffs >= thr_diff)[1]]

  if (is.na(inc_st_idx)) {
    return(NA_integer_)
  } else {
    return(inc_st_idx)
    
  }
}


#----------------------------------------------------------------------------------
# Within-window smoothing
#----------------------------------------------------------------------------------

## Smooth VI/FI time series within the fitting window
smooth_in_window <- function(df,
                             k_interp = 15,
                             min_obs_interp = 15,
                             k_sg = 15,
                             p = 2,
                             spar = NULL,
                             max_gap = 10) {

  ## Construct a complete daily time series within the fitting window
  full_dates <- tibble(date = seq(min(df$date), max(df$date), by = "1 day"))

  df <- full_dates %>%
    left_join(df, by = "date") %>%
    arrange(date)

  df <- df %>%
    mutate(
      across(
        -c(date, value, value_cleaned),
        ~ first(na.omit(.x))
      ),
      doy = yday(date)
    )

  y2 <- df$value_cleaned
  x2 <- df$doy

  valid_idx0 <- which(!is.na(y2))

  if (length(valid_idx0) < 10) {
    df$value_interp <- NA_real_
    df$value_filled <- NA_real_
    df$value_filtered <- NA_real_
    df$value_smooth <- NA_real_
    df$long_gap <- is.na(y2)
    return(df)
  }

  ## Identify long data gaps that should not be bridged by smoothing
  is_na <- is.na(y2)
  r <- rle(is_na)

  gap_len <- inverse.rle(list(
    lengths = r$lengths,
    values = ifelse(r$values, r$lengths, 0)
  ))

  long_gap <- is_na & gap_len >= max_gap
  df$long_gap <- long_gap

  ## Fill short gaps using the local median when sufficient observations are available
  y2_interp <- y2
  na_idx <- which(is.na(y2_interp))

  for (i in na_idx) {
    left  <- max(1, i - floor(k_interp / 2))
    right <- min(length(y2_interp), i + floor(k_interp / 2))

    win <- y2[left:right]
    win <- win[!is.na(win)]

    if (length(win) >= min_obs_interp) {
      y2_interp[i] <- median(win)
    }
  }

  df$value_interp <- y2_interp
  
  ## Fill remaining gaps by linear interpolation and nearest-value extension
  y2_fill <- zoo::na.approx(
    y2_interp,
    x = x2,
    na.rm = FALSE
  )

  y2_fill <- zoo::na.locf(y2_fill, na.rm = FALSE)
  y2_fill <- zoo::na.locf(y2_fill, fromLast = TRUE, na.rm = FALSE)

  df$value_filled <- y2_fill

  if (any(is.na(y2_fill))) {
    df$value_filtered <- NA_real_
    df$value_smooth <- NA_real_
    return(df)
  }
  
  ## Apply Savitzky–Golay filtering to the gap-filled time series
  if (k_sg %% 2 == 0) k_sg <- k_sg + 1
  if (p >= k_sg) stop("p must be smaller than k_sg")

  y2_hat <- signal::sgolayfilt(y2_fill, p = p, n = k_sg)

  ## Restore long gaps before spline smoothing
  y2_hat_for_spline <- y2_hat
  y2_hat_for_spline[long_gap] <- NA_real_

  df$value_filtered <- y2_hat_for_spline

  ## Apply spline smoothing to the SG-filtered time series
  if (!is.null(spar)) {

    valid_idx <- which(!is.na(y2_hat_for_spline) & !is.na(x2))

    if (length(valid_idx) < 10) {
      df$value_smooth <- NA_real_
    } else {
      fit2 <- smooth.spline(
        x = x2[valid_idx],
        y = y2_hat_for_spline[valid_idx],
        spar = spar
      )

      df$value_smooth <- predict(fit2, x = x2)$y
    }

  } else {
    df$value_smooth <- df$value_filtered
  }
  
  ## Return the smoothed within-window time series
  return(df)
}


#----------------------------------------------------------------------------------
# Flowering phenology detection
#----------------------------------------------------------------------------------

## Detect flowering phenology using green-up stages and FI peak dynamics
detect_flower <- function(dfvi, dffi, sg, phenopar) {
  
  spar      <- phenopar$WsplineSpar           # for within-window FI spline smoothing
  sparv     <- phenopar$AsplineSpar           # for annual VI smoothing
  
  hsfws     <- phenopar$half_fitWsize         # half-width of fitting window
  hdfws     <- phenopar$half_detectWsize      # half-width of detection window
  
  vmaxthr  <- phenopar$VImaxThresh            # for detecting vegetation
  vampthr   <- phenopar$VIampThresh           # for detecting deciduous trees
  
  gpthr     <- phenopar$gup_threshes          # for detecting green-phase
  fthr      <- phenopar$FIampThresh           # for detecting flowering peak
  fincthr   <- phenopar$Finc_Fdiff_pc         # for finding start point
  pitw      <- phenopar$Finc_pitWmax          # for finding start point
  
  ## Extract green-up phenophases and annual VI metrics
  gdetecs <- detect_greenup(dfvi, dffi, phenopar)
  green_date <- gdetecs$phase
  vamp <- gdetecs$amp
  vmax <- gdetecs$max
  vmin <- gdetecs$min
  
  ## Remove coincident VI and FI outliers using MAD-based detection
  vi_ois <- detect_outliers_mad(x=dfvi$value,
                                dates=dfvi$date,
                                z=phenopar$VspikeThresh,
                                maxd=phenopar$maxDistance,
                                n_repeat=phenopar$spike_n_repeat)
  fi_ois <- detect_outliers_mad(x=dffi$value,
                                dates=dffi$date,
                                z=phenopar$FspikeThresh,
                                maxd=phenopar$maxDistance,
                                n_repeat=phenopar$spike_n_repeat)
  ois <- vi_ois & fi_ois
  
  dfvi$value_cleaned <- dfvi$value
  dffi$value_cleaned <- dffi$value
  if (any(ois)) {
    dfvi$value_cleaned[ois] <- NA
    dffi$value_cleaned[ois] <- NA
  }
  
  ## Set species-specific fitting windows based on green-up phenophases
  if (str_detect(sg, "Py")) {
    ref_date <- green_date$date[1]    # onset_greenup
    fitw <- c(ref_date-hsfws, ref_date+hsfws)
  } else if (str_detect(sg, "Rp")) {
    ref_date <- green_date$date[2]    # mid_greenup
    fitw <- c(ref_date-hsfws, ref_date+hsfws)
  } else if (str_detect(sg, "Cc")) {
    ref_date <- green_date$date[3]    # end_greenup
    fitw <- c(green_date$date[1], green_date$date[1]+(hsfws*2))

  }
  if (is.na(ref_date)) {
    return(list(
      vi = NA,
      fi = NA,
      peak = tibble(
        peak_num=NA,
        ref_date=NA,
        peak_date=NA,
        peak_value=NA,
        start_date=NA,
        start_value=NA,
        end_date=NA,
        end_value=NA,
        peak_amp=NA,
        noise_level=NA,
        period=NA
      )
    ))
  }
  subfi <- dffi %>%
    filter(date>fitw[1] & date<fitw[2])
  subvi <- dfvi %>%
    filter(date>fitw[1] & date<fitw[2])
  
  ## Remove residual FI outliers within the fitting window
  rm_date <- detect_outliers_spline(date=subfi$date,
                                    value=subfi$value_cleaned,
                                    spar=spar,
                                    z_thresh=phenopar$z_thresh,
                                    n_repeat=phenopar$n_repeat)

  if (length(rm_date)>0) {
    subfi[subfi$date %in% rm_date,]$value_cleaned <- NA
    subvi[subvi$date %in% rm_date,]$value_cleaned <- NA
  }

  ## Smooth FI and VI within the fitting window
  subfi <- smooth_in_window(df=subfi, spar=spar)
  subvi <- smooth_in_window(df=subvi, spar=spar)

  noise <- subfi$value_cleaned - subfi$value_smooth
  noise_lv <- median(abs(noise - median(noise, na.rm = TRUE)), na.rm = TRUE)     # noise level
  
  ## Exclude pixels with insufficient signals or unreliable green-up phenology
  if (all(is.na(subfi$value_smooth)) || 
      all(is.na(subvi$value_smooth)) || 
      vmax < vmaxthr  ||                                   # for detecting vegetation
      vamp < vampthr  ||                                   # for detecting deciduous trees
      is.na(green_date$date[3]) ||                         # incomplete annual time series
      as.numeric(substr(green_date$date[3],6,7)) >= 7      # unstable/late green-up
  ) {
    return(list(
      vi = subvi,
      fi = subfi,
      peak = tibble(
        peak_num=NA,
        ref_date = ref_date,
        peak_date = NA,
        peak_value = NA,
        start_date = NA,
        start_value = NA,
        end_date = NA,
        end_value = NA,
        peak_amp=NA,
        noise_level=NA,
        period=NA
      )
    ))
  }
  
  ## Identify candidate flowering peaks from the smoothed FI time series
  peaks_idx <- find_peaks(subfi$value_smooth)
  peak_dates <- subfi$date[peaks_idx]
  
  ## Apply species-specific phenological constraints to candidate peaks
  if (str_detect(sg, "Cc")) {
    # Cc: retain peaks after canopy maturity
    keep <- which(peak_dates > green_date$date[3])
    peaks_idx <- peaks_idx[keep]
    peak_dates <- peak_dates[keep]
    if (length(peak_dates)!=0) {
      # valley idx after mid-greenup & closer to maturity
      valleys_idx <- find_valleys(subfi$value_smooth)
      valleys_idx <- valleys_idx[
        abs(as.numeric(subfi$date[valleys_idx] - green_date$date[3])) <
          abs(as.numeric(subfi$date[valleys_idx] - green_date$date[2]))
      ]
      valid_peaks <- sapply(peaks_idx, function(p) {
        any(valleys_idx < p)
      })
      peaks_idx <- peaks_idx[valid_peaks]
      peak_dates <- peak_dates[valid_peaks]
    }
  } else if (str_detect(sg, "Rp")) {
    # Rp: retain peaks between green-up onset and maturity
    keep <- which(peak_dates < green_date$date[3] & peak_dates > green_date$date[1])
    peaks_idx <- peaks_idx[keep]
    peak_dates <- peak_dates[keep]
    if (length(peak_dates)!=0) {
      # valley idx after greenup onset & closer to mid-greenup
      valleys_idx <- find_valleys(subfi$value_smooth)
      valleys_idx <- valleys_idx[
        abs(as.numeric(subfi$date[valleys_idx] - green_date$date[2])) <
          abs(as.numeric(subfi$date[valleys_idx] - green_date$date[1]))
      ]
      valid_peaks <- sapply(peaks_idx, function(p) {
        any(valleys_idx < p)
      })
      peaks_idx <- peaks_idx[valid_peaks]
      peak_dates <- peak_dates[valid_peaks]
    }
  } else  { 
    # Py: retain peaks before mid-green-up
    keep <- which(peak_dates < green_date$date[2])
    peaks_idx <- peaks_idx[keep]
    peak_dates <- peak_dates[keep]
  }
  
  ## Check whether any valid candidate peaks remain
  if (is.null(peak_dates) ||
      length(peak_dates) == 0 ||
      all(is.na(peak_dates))) {
    peak <- NULL
  } else {
    ## Store valid candidate peak dates and values
    peak <- tibble(
      peak_num = seq(1, length(peak_dates)),
      ref_date = rep(ref_date,length(peak_dates)),
      peak_date = peak_dates,
      peak_value = subfi$value_smooth[peaks_idx]
    )
  }
  
  if (is.null(peak)) {
    best_peak <- tibble(
      peak_num=NA,
      ref_date = ref_date,
      peak_date = NA,
      peak_value = NA,
      start_date = NA,
      start_value = NA,
      end_date = NA,
      end_value = NA,
      peak_amp=NA,
      noise_level=noise_lv,
      period=NA
    )
  } else {
    
    ## Characterize each candidate peak
    peak_info <- data.frame()
    
    for (pn in 1:nrow(peak)) {
      apeak <- peak[pn,] 
      
      # Identify the flowering increase start point
      peak_i <- which(subfi$date == apeak$peak_date)
      prior_idx <- 1:(peak_i)
      
      start_i <- find_increase(y=subfi$value_smooth[prior_idx], thr=fincthr, pitw=pitw)
      if (!is.na(start_i) && start_i == 1) {
        start_i <- find_increase(y=subfi$value_smooth[prior_idx], thr=fincthr, pitw=0)
      }
      start <- if (!is.na(start_i)) subfi$date[start_i] else NA
      start_val <- if (!is.na(start_i)) subfi$value_smooth[start_i] else NA
      
      # Identify the flowering end point
      if (!is.na(start_val)) {
        after_idx <- (peak_i+1):nrow(subfi)
        diffs <- subfi$value_smooth[after_idx] - start_val
        rel_end_i <- which(diffs < 0)[1]
        end_i <- after_idx[rel_end_i]
        end <- subfi$date[end_i]
        end_val <- subfi$value_smooth[end_i]
      } else {
        end <- NA
        end_val <- NA
      }
      
      apeak_info <- apeak %>%
        mutate(
          start_date = as.Date(start),
          start_value = start_val,
          end_date = as.Date(end),
          end_value = end_val
        )
      datinperiod <- subfi$value_cleaned[subfi$date>=apeak_info$peak_date-14 & subfi$date<=apeak_info$peak_date+15]
      datnuminperiod <- sum(!is.na(datinperiod))
      apeak_info$data_num <- datnuminperiod
      peak_info <- rbind(peak_info, apeak_info)
      
    }
    
    ## Calculate flowering metrics for each candidate peak
    peak_info <- peak_info %>% 
      mutate(peak_amp=peak_value-start_value, 
             diff_ref=abs(peak_date-ref_date),
             period=end_date-start_date,
             noise_level=noise_lv)
    
    ## Restrict candidate peaks to the flowering detection window
    if (str_detect(sg, "Cc")) {
      peak_info <- peak_info %>%
        filter(peak_date>green_date$date[3] & peak_date<=green_date$date[3]+(hdfws*2)) 
    } else {
      peak_info <- peak_info %>%
        filter(peak_date>ref_date-hdfws & peak_date<=ref_date+hdfws) 
    }

    ## Select the final flowering peak
    ## Retain near-maximum peaks and select the one closest to the reference date
    if (nrow(peak_info)>=1) {
      
      max_amp <- max(peak_info$peak_amp , na.rm = TRUE)
      if (str_detect(sg, "Py")) {
        best_peak <- peak_info %>%
          filter(peak_amp/max_amp >= 0.9 | peak_value==max(peak_value)) %>%
          arrange(diff_ref) %>%
          slice(1)
      } else {
        best_peak <- peak_info %>%
          filter(peak_amp/max_amp >= 0.9) %>%
          arrange(diff_ref) %>%
          slice(1)
      }

      best_peak <- best_peak %>% dplyr::select(-data_num, -diff_ref) %>% filter(peak_amp>=fthr)

    }
    
    ## Return NA when no valid flowering peak is selected
    if (nrow(peak_info)==0 || !exists("best_peak") || nrow(best_peak)==0) {
      best_peak <- tibble(
        peak_num=NA,
        # ref_date = greenphase,
        ref_date = ref_date,
        peak_date = NA,
        peak_value = NA,
        start_date = NA,
        start_value = NA,
        end_date = NA,
        end_value = NA,
        peak_amp=NA, 
        noise_level=noise_lv,
        period=NA
      )
    }

  }
  
  return(list(
    vi = subvi,
    fi = subfi,
    peak = best_peak
  ))
  
}

