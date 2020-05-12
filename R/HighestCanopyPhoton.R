#' Detect the highest canopy photon within moving window
#'
#' @param df
#' @param window
#'
#' @export

HighestPhotonWindow <- function(df,
                                window = 0.7,
                                h_norm = "h_norm",
                                beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                                gap_dist = 10,
                                min_ph_gap = 20,
                                keep_cols,
                                dist_along = "dist_along",
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
        dat <- df[,c(h_norm, dist_along, beam_n, ph_id, ATL08_class, kc)]
        colnames(dat) <- c("h_norm", "dist_along", "beam", "ph_id", "ph_class", kc)
      }else{
        warning("keep columns are not in df")
        dat <- df[,c(h_norm, dist_along, beam_n, ph_id, ATL08_class)]
        colnames(dat) <- c("h_norm", "dist_along", "beam", "ph_id", "ph_class")
      }
    }
  }else{
    dat <- df[,c(h_norm, dist_along, beam_n, ph_id, ATL08_class)]
    colnames(dat) <- c("h_norm", "dist_along", "beam", "ph_id", "ph_class")
  }

  if (is.factor(dat$ph_id)) dat$ph_id <- as.character(dat$ph_id)

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
    dat_beam <- dplyr::filter(dat, beam == n)

    inter <- findInterval(dat_beam$dist_along,seq(from = min(dat_beam$dist_along),to = max(dat_beam$dist_along),by = window))
    dat_inter <- split(dat_beam, inter)

    ph_id_top <- unlist(lapply(dat_inter, function(x) {
      #dat <- dplyr::filter(x, ph_class == "canopy" | ph_class == "top of canopy")
      dat <- dplyr::filter(x, ! ph_class %in% "noise")
      if (NROW(dat) > 0) {
        max_h <- dat[which.max(dat$h_norm),"ph_id"]
      }else{
        max_h <- NA
      }
      max_h
      }))
    ph_id_top <- ph_id_top[!is.na(ph_id_top)]
    dat_beam$top_ph <- FALSE
    dat_beam[dat_beam$ph_id %in% ph_id_top, "top_ph"] <- TRUE

    out[[n]] <- dat_beam

  }
  return(out)
}
