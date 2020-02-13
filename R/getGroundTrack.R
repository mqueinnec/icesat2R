#' Get the ground track of ATLAS beams based on ATL03 data product
#'
#'@param ATL03 Path to ATL03 HDF5 file (.h5 extension)
#'@param beam Character. List of beams for which groudn track needs to be extracted
#'@param rgt Logical. Should RGT number be added as a field to output feature?
#'@param cycle Logical. Should cycle number be added as a field to output feature?
#'@param segment Logical. Should segment number be added as a field to output feature?
#'@param date Logical. Should date be added as a field to output feature?
#'@param time Logical. Should time be added as a field to output feature?
#'@param points Logical. Should each reference point be returned in a SpatialPointsDataFrame?
#'@param lines Logical. Should lines features between all reference points be returned in a SpatialLinesDataFrame?
#'@param ... Parameters passed to writeOGR (e.g. overwrite_layer = TRUE)

getGroundTrack <- function(ATL03,
                           beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                           rgt = TRUE,
                           cycle = TRUE,
                           segment = TRUE,
                           date = TRUE,
                           time = TRUE,
                           points = FALSE,
                           lines = TRUE,
                           write = TRUE,
                           outDir,
                           filename_points,
                           filename_lines,
                           ...) {

  if(missing(outDir)) outDir <- ""
  if (write & (gsub(" ", "", outDir) == "" | is.null(outDir))) outDir <- getwd()

  #Check ATL03 file extension

  if(!is.character(ATL03)) {
    stop("ATL03 must be a HDF5 file name (extension h5")
  }else{
    ATL03 <- normalizePath(ATL03)
    if(!file.exists(ATL03)) stop("ATL03 filename does not exist")
    if(tools::file_ext(ATL03) != "h5") stop("ATL03 must be a HDF5 file (h5 extension")
  }

  #Check validity of beams

  if (any(beam %in% c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r") == FALSE)) {
    stop("beam must be gt1l, gt1r, gt2l, gt2r, gt3l or gt3r")
  }


  test_format <- grepl("ATL03_[0-9]{14}_[0-9]{8}_[0-9]{3}_[0-9]{2}", tools::file_path_sans_ext(basename(ATL03)))

  if(!test_format) {
    warning("ATL03 does not follow the conventional ICESat-2 file naming. rgt, cycle, segment, dat and time (if set to TRUE) will not be written as attributes")
    rgt <- FALSE
    cycle <- FALSE
    segment = FALSE
    date <- FALSE
    time <- FALSE

    h5file_filter <- tools::file_path_sans_ext(basename(ATL03))
  }else{
    h5file_filter <- stringr::str_extract(basename(ATL03), pattern = "A.*5") #Keep only file name format
  }

  if(missing(filename_points)) filename_points = ""
  if (write & (gsub(" ", "", filename_points) == "" | is.null(filename_points)) ) {
    filename_points <- paste0(tools::file_path_sans_ext(h5file_filter),"_points")
  }

  if(missing(filename_lines)) filename_lines = ""
  if (write & (gsub(" ", "", filename_lines) == "" | is.null(filename_lines)) ) {
    filename_lines <- paste0(tools::file_path_sans_ext(h5file_filter),"_lines")
  }



  h5list <- rhdf5::h5ls(ATL03)

  required <- c(paste0(beam,"/geolocation/reference_photon_lat"),
                paste0(beam,"/geolocation/reference_photon_lon"),
                paste0(beam,"/geolocation/segment_ph_cnt"),
                paste0(beam,"/geolocation/segment_id"))

  check <- paste0("/",required) %in% paste(h5list$group, h5list$name, sep = "/")

  if (any(check == FALSE)) {
    stop(paste0("Missing dataset in ATL03: ", paste(required[!check], collapse = ",")))
  }

  dat_ref_ph <- data.frame() #Stores all reference photons, grouped per beam
  att_table <- data.frame() #Attribute table for lines

  for ( n in beam) {
    ref_lat <- as.numeric(rhdf5::h5read(file = ATL03, name = paste0(n,"/geolocation/reference_photon_lat")))
    ref_lon <- as.numeric(rhdf5::h5read(file = ATL03, name = paste0(n,"/geolocation/reference_photon_lon")))
    segment_ph_cnt <- as.numeric(rhdf5::h5read(file = ATL03, name = paste0(n,"/geolocation/segment_ph_cnt")))
    segment_id <- as.character(rhdf5::h5read(file = ATL03, name = paste0(n,"/geolocation/segment_id")))

    temp_df <- data.frame(lon = ref_lon, lat = ref_lat, beam = n,segID = segment_id, segCount = segment_ph_cnt)
    temp_df <- temp_df[temp_df$segCount > 0, ]

    dat_ref_ph <- rbind(dat_ref_ph, temp_df)


    if (!all(segment_ph_cnt == 0)){
      toAdd <- data.frame(beam = n)
      rownames(toAdd) <- n
      att_table<- rbind(att_table, toAdd)
    }
  }

  if (rgt) dat_ref_ph$rgt <- substr(h5file_filter, start = 22, stop = 25)
  if(cycle) dat_ref_ph$cycle <- substr(h5file_filter, start = 26, stop = 27)
  if(segment) dat_ref_ph$segment <- substr(h5file_filter, start = 28, stop = 29)
  if(date) dat_ref_ph$date <- substr(h5file_filter, start = 7, stop = 14)
  if(time) dat_ref_ph$time <- substr(h5file_filter, start = 15, stop = 20)

  if (rgt) att_table$rgt <- substr(h5file_filter, start = 22, stop = 25)
  if(cycle) att_table$cycle <- substr(h5file_filter, start = 26, stop = 27)
  if(segment) att_table$segment <- substr(h5file_filter, start = 28, stop = 29)
  if(date) att_table$date <- substr(h5file_filter, start = 7, stop = 14)
  if(time) att_table$time <- substr(h5file_filter, start = 15, stop = 20)

  ref_spdf <- sp::SpatialPointsDataFrame(coords = dat_ref_ph[,c("lon", "lat")], data = dat_ref_ph, proj4string = CRS("+init=epsg:4326"))

  paths <- sp::split(ref_spdf, ref_spdf[["beam"]])

  test <- lapply(paths, sp::Line)
  test2 <- sapply(names(test), function(x){Lines(list(test[[x]]), x)})
  test3 <- sp::SpatialLines(test2, proj4string = CRS("+init=epsg:4326"))

  lines_df <- sp::SpatialLinesDataFrame(test3, data = att_table)

  if (write) {
    if(points) {
      rgdal::writeOGR(ref_spdf, dsn = outDir, layer = filename_points, driver = "ESRI Shapefile", ...)
    }
    if(lines) {
      rgdal::writeOGR(lines_df, dsn = outDir, layer = filename_lines, driver = "ESRI Shapefile", ...)
    }
  }


  return(list(if(lines)lines = lines_df,
              if(points) points = ref_spdf))

}
