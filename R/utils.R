#' Check if file follows conventional ATLAS file naming
#'
#'
#'@export

check_ATL_naming <- function(file, product) {

  if(!is.character(file)) stop("File must be a character")


  if(missing(product)) {
    grepl("ATL[0-9]{2}_[0-9]{14}_[0-9]{8}_[0-9]{3}_[0-9]{2}", tools::file_path_sans_ext(basename(file)))
  }else{
    product <- toupper(product)
    grepl(paste0(product,"_[0-9]{14}_[0-9]{8}_[0-9]{3}_[0-9]{2}"), tools::file_path_sans_ext(basename(file)))
  }
}

#' Extract conventional ATLAS file naming
#'
#'
#'@export

extract_ATL_ID <- function(file, product) {

  ID <- character()
  date <- character()
  time <- character()
  rgt <- character()

  for (f in 1:length(file)) {
    if(check_ATL_naming(file[f], product)) {
      ID[f] <- stringr::str_extract(file[f],"[0-9]{14}_[0-9]{8}_[0-9]{3}_[0-9]{2}")
      date[f] <- substr(ID[f],1,8)
      time[f] <- substr(ID[f],9,14)
      rgt[f] <- substr(ID[f],16,19)
    }else{
      warning(sprintf("Input file name has been modified from original: %s", file[f]))
      ID[f] <- NA
      date[f] <- NA
      time[f] <- NA
      rgt[f] <- NA
    }
  }
  data.frame(ID = ID,
             date = date,
             time = time,
             Date = as.Date(paste0(date, time),"%Y%m%d%H%M%S"),
             rgt = rgt)
}

#' Find UTM zone
#' @export

find_UTM_zone <- function(longitude, latitude) {

  # Special zones for Svalbard and Norway
  if (latitude >= 72.0 && latitude < 84.0 )
    if (longitude >= 0.0  && longitude <  9.0)
      return(31);
  if (longitude >= 9.0  && longitude < 21.0)
    return(33)
  if (longitude >= 21.0 && longitude < 33.0)
    return(35)
  if (longitude >= 33.0 && longitude < 42.0)
    return(37)

  (floor((longitude + 180) / 6) %% 60) + 1
}

#' Find UTM zone
#' @export

find_UTM_hemisphere <- function(latitude) {

  ifelse(latitude > 0, "N", "S")
}

#' UTM projstring
#' @export

UTM_proj <- function(longitude, latitude) {

  zone <- find_UTM_zone(longitude, latitude)
  hemisphere <- find_UTM_hemisphere(latitude)

  proj <- paste("+proj=utm",
                paste0("+zone=",zone),
                paste0("+",hemisphere),
                "+ellps=WGS84",
                "+units=m")

  return(proj)
}

#' Rotation matrix ATL03
#' @export

get_rot_matrix <- function(xy_start, xy_end) {

  x1 = xy_start[1]
  x2 = xy_end[1]
  y1 = xy_start[2]
  y2 = xy_end[2]

  deltaX = x2 - x1
  deltaY = y2 - y1

  theta = atan2(deltaX, deltaY)

  phi = pi/2 - theta

  R_mat <- matrix(c(cos(phi), -sin(phi), sin(phi), cos(phi)), ncol = 2, byrow = FALSE)

  xRotPt = x1
  yRotPt = y1

  return(list(R_mat = R_mat,
              xRotPt = xRotPt,
              yRotPt = yRotPt))
}

#' Convert x/y to along and across distance
#' @export


get_along_distance <- function(x_in, y_in, R_mat, xRotPt, yRotPt) {

  xTranslated = x_in - xRotPt
  yTranslated = y_in - yRotPt

  xyTranslated_mat <- matrix(c(xTranslated, yTranslated), byrow = TRUE, nrow = 2)

  measRot_mat = R_mat %*% xyTranslated_mat

  return(list(along_distance = measRot_mat[1,],
         across_distance = measRot_mat[2,]))

}

#' Find pairs of ATL03 and ATL08 files
#' @param atl03_flist Vector of character listing the file path to ATL03 HDF files
#' @param atl08_flist Vector of character listing the file path to ATL08 HDF files
#' @export

find_atl_pairs <- function(atl03_flist,
                           atl08_flist) {

  atl_pattern = "(ATL[0-9]{2})_([0-9]{14}_[0-9]{8}_[0-9]{3}_[0-9]{2})"

  check_pattern <- lapply(list(atl03_flist, atl08_flist), function(x) {
    out <- stringr::str_match_all(x, atl_pattern)
    out <- data.frame(do.call(rbind, out))
    colnames(out) <- c("basename", "product", "id")
    out$fname <- x
    return(out)
  })

  if(!all(check_pattern[[1]]$product == "ATL03")) {
    stop("Not all files in atl03_flist are ATL03 files")
  }

  if(!all(check_pattern[[2]]$product == "ATL08")) {
    stop("Not all files in atl08_flist are ATL08 files")
  }

  atl03_df <- check_pattern[[1]]
  colnames(atl03_df) <- stringr::str_c(colnames(atl03_df), "atl03", sep = "_")

  atl08_df <- check_pattern[[2]]
  colnames(atl08_df) <- stringr::str_c(colnames(atl08_df), "atl08", sep = "_")


  joined_list <- as.list(dplyr::inner_join(atl03_df, atl08_df,
                                    by = c("id_atl03" = "id_atl08")))

  out <- mapply(function(x, y) list(ATL03 = x, ATL08 = y),
                x = as.list(joined_list$fname_atl03),
                y = as.list(joined_list$fname_atl08), SIMPLIFY = FALSE)

  return(out)
}
