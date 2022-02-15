#' Calculates vegetation fill index
#'
#' @param df
#' @param th
#' @param top_ph
#' @param ph_class
#'
#' @export
#'

vegetationFill <- function(df,
                           th,
                           dist_along = "dist_along",
                           h_norm = "h_norm",
                           top_ph = "top_ph",
                           ph_class = "ph_class") {

  dat <- df[,c(h_norm, dist_along, top_ph, ph_class)]
  colnames(dat) <- c("h_norm", "dist_along", "top_ph", "ph_class")

  dat <- dplyr::filter(dat,  ! ph_class %in% "noise" & top_ph == TRUE)

  if (missing(th)) {
    th <- max(dat$h_norm)
  }else{
    if (max(dat$h_norm) < th) {
      warning(sprintf("th of %0.2f m higher than max height (%0.2f m)", th, max(dat$h_norm)))
      #th <- max(dat$h_norm)
    }
  }

  if (th > 0) {
    dat[dat$h_norm >= th, "h_norm"] <- th
    dat <- dplyr::arrange(dat, dist_along)

    sum_h <- dat[2:NROW(dat), "h_norm"] + dat[1:(NROW(dat) - 1), "h_norm"]
    dist <- dat[2:NROW(dat), "dist_along"] - dat[1:(NROW(dat) - 1), "dist_along"]
    l <- max(dat$dist_along) - min(dat$dist_along)

    VF <- 0.5*sum(dist*sum_h) / (th * l)
  } else {
    VF <- 0
  }

  return(VF)
}
