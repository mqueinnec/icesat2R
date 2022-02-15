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
