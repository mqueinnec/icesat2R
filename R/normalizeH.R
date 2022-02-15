#' Normalize photons height
#'
#' Use local polynomial regression fitting with degree 0 (moving average) to interpolate ground surface and normalize photons using the interpolated surface
#'
#'Segments with continuous ground points are first identified and interpolation is performaed at these segment level. THe continuity between ground point is determined by the gap_dist threshold: if two groudn points are more than gap_dist along-distance apart, a gap is detected and ground surface will be interpolated for each gap. Additionally, segments with less than min_ph_gap photons in them are discarded. Therefore, not all photons in the original df are guaranteed to be normalized.
#'
#' @param df Data frame
#' @param x Name of df column indicating photon x coordinate
#' @param y Name of df column indicating photon y coordinate
#' @param h x Name of df column indicating photon height
#' @param dh min and max height thresholds used to flag unsual photons height
#' @param dx min and max horizontal distance
#'
#' @export

normalizeH <- function(df,
                       h = "h",
                       beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                       gap_dist = 10,
                       min_ph_gap = 20,
                       returnGaps = TRUE,
                       keeps_cols,
                       dist_along = "cum_dist_along",
                       segment_id = "seg_id",
                       segment_idx = "seg_idx",
                       ph_id = "ph_id",
                       beam_n = "beam",
                       ATL08_class = "ph_class") {

  # Check if columns to keep are in df

  if(!missing(keep_cols)) {
    if(!is.character(keep_cols)){
      stop("keep_cols must be a character vector")
    }else{
      kc <- keep_cols[keeps_cols %in% colnames(df)]
      if(length(kc) > 0) {
        if(length(kc) != length(keep_columns)) {
          warning(sprintf("%s if not in df",keep_cols[!keep_cols %in% kc]))
        }
        dat <- df[,c(h, dist_along, segment_id, segment_idx, beam_n, ph_id, ATL08_class, kc)]
        colnames(dat) <- c("h", "dist_along", "seg_id", "seg_idx", "beam", "ph_id", "ph_class", kc)
      }else{
        warning("keep columns are not in df")
        dat <- df[,c(h, dist_along, segment_id, segment_idx, beam_n, ph_id, ATL08_class)]
        colnames(dat) <- c("h", "dist_along", "seg_id", "seg_idx", "beam", "ph_id", "ph_class")
      }
    }
  }else{
    dat <- df[,c(h, dist_along, segment_id, segment_idx, beam_n, ph_id, ATL08_class)]
    colnames(dat) <- c("h", "dist_along", "seg_id", "seg_idx", "beam", "ph_id", "ph_class")
  }

  unique_beam <- unique(dat$beam)

  if (any(! beam %in% unique_beam)) {
    missing_beams <- beam[! beam %in% unique_beam]
    message(sprintf("The following requested beams are not in df and will be discarded: %s",paste(missing_beams, collapse = ", ")))
    beam = beam[beam %in% unique_beam]
  }

  if (length(beam) == 0 ) {
    stop("No valid beam")
  }else{
    unique_beam <- unique_beam[unique_beam %in% beam]
  }


  out <- list()
  for (n in unique_beam) {
    message(sprintf("Working on beam %s", n))

    dat_beam <- dplyr::filter(dat, beam == n)
    dat_beam_ground <- dat_beam %>%
      dplyr::filter(ph_class == "ground") %>%
      dplyr::arrange(dist_along)

    if (NROW(dat_beam_ground) > 0) {
      # Find ground gaps in data > threshold

      pair_dist <- dat_beam_ground[2:NROW(dat_beam_ground), "dist_along"] - dat_beam_ground[1:(NROW(dat_beam_ground) - 1),"dist_along"]
      gaps_end <- c(which(pair_dist > gap_dist) ,NROW(dat_beam_ground))
      gaps_start <- c(1, gaps_end[1:(length(gaps_end) - 1 )] + 1)


      gaps_length <- gaps_end - gaps_start + 1
      keep <- gaps_length > min_ph_gap #Only keep groups with more than min_ph_gap photons

      gaps_length <- gaps_length[keep]
      gaps_start <- gaps_start[keep]
      gaps_end <- gaps_end[keep]

      if (returnGaps) {
        out_gap <- list()
      }else{
        out_gap <- data.frame()
      }

      if (length(gaps_length) == 0 ) {
        warning("Not enough ground points for beam %s to interpolate any ground surface")
      }else{
        for (i in 1:length(gaps_length)) {

          # Ground surface interpolation
          dat_gap_ground <- dat_beam_ground[gaps_start[i]:gaps_end[i],]

          calcSSE <- function(x){
            loessMod <- try(loess(h ~ dist_along, data=dat_gap_ground, span=x, degree = 0), silent=T)
            res <- try(loessMod$residuals, silent=T)
            if(class(res)!="try-error"){
              if((sum(!is.na(res)) > 0)){
                sse <- sum(res^2)
              }else{
                sse <- 99999
              }
            }else{
              sse <- 99999
            }
            return(sse)
          }

          opt <- optimize(calcSSE, lower = 0.0001, upper = 0.2)

          loessMod <- loess(h ~ dist_along, data=dat_gap_ground, span=opt$minimum, degree = 0)

          # Normalize height
          dat_gap_all <- dplyr::filter(dat_beam, dist_along >= dat_beam_ground[gaps_start[i],"dist_along"] & dist_along <= dat_beam_ground[gaps_end[i],"dist_along"])

          dat_gap_all$h_norm <- dat_gap_all$h - predict(loessMod, dat_gap_all$dist_along)

          if (returnGaps) {
            out_gap[[i]] <- dat_gap_all
          }else{
            out_gap <- rbind(out_gap, dat_gap_all)
          }
        }
    }
      out[[n]] <- out_gap
    }else{
      warning("Not ground points for beam %s to interpolate any ground surface")
    }
  }
  return(out)
}
