#' Get the ground track of ATLAS beams based on ATL03 ot ATL08 data product
#'
#'@param file Path to ATL03 or ATL08 HDF5 file (.h5 extension)
#'@param product Either "ATL03" or "ATL08"
#'@param beam Character. List of beams for which groudn track needs to be extracted
#'@param rgt Logical. Should RGT number be added as a field to output feature?
#'@param cycle Logical. Should cycle number be added as a field to output feature?
#'@param segment Logical. Should segment number be added as a field to output feature?
#'@param date Logical. Should date be added as a field to output feature?
#'@param time Logical. Should time be added as a field to output feature?
#'@param points Logical. Should each reference point be returned in a SpatialPointsDataFrame?
#'@param lines Logical. Should lines features between all reference points be returned in a SpatialLinesDataFrame?
#'@param ... Parameters passed to writeOGR (e.g. overwrite_layer = TRUE)
#'
#'@export

getGroundTrack <- function(file,
                           product,
                           beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                           beam_strength = c("weak", "strong"),
                           points = FALSE,
                           lines = TRUE,
                           df = TRUE,
                           write = TRUE,
                           outDir,
                           filename_points,
                           filename_lines,
                           ...) {

  if(missing(outDir)) outDir <- ""
  if (write & (gsub(" ", "", outDir) == "" | is.null(outDir))) outDir <- getwd()

  #Check file file extension

  if(!is.character(file)) {
    stop("file must be a HDF5 file name (extension h5")
  }else{
    file <- normalizePath(file)
    if(!file.exists(file)) stop("file  does not exist")
    if(tools::file_ext(file) != "h5") stop("file must be a HDF5 file (h5 extension")
  }

  #Check validity of beams

  if (any(beam %in% c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r") == FALSE)) {
    stop("beam must be gt1l, gt1r, gt2l, gt2r, gt3l or gt3r")
  }

  if (any(!beam_strength %in% c("weak", "strong"))) {
    stop("beam_strength must contain weak, beam or both ")
  }

  # Check required datasets

  h5list <- rhdf5::h5ls(file)

  if (product == "ATL03") {
    required <- list(ref_ph_lat = paste0(beam,"/geolocation/reference_photon_lat"),
                     ref_ph_lon = paste0(beam,"/geolocation/reference_photon_lon"),
                     seg_ph_cnt = paste0(beam,"/geolocation/segment_ph_cnt"),
                     seg_id = paste0(beam,"/geolocation/segment_id"),
                     start_region = "ancillary_data/start_region",
                     sc_orient = "orbit_info/sc_orient",
                     rgt = "orbit_info/rgt",
                     cycle = "orbit_info/cycle_number",
                     date_start = "ancillary_data/data_start_utc",
                     date_end = "ancillary_data/data_end_utc")
  }else if (product == "ATL08") {
    required <- list(lat = paste0(beam,"/land_segments/latitude"),
                  lon = paste0(beam,"/land_segments/longitude"),
                  start_region = "ancillary_data/start_region",
                  sc_orient = "orbit_info/sc_orient",
                  rgt = "orbit_info/rgt",
                  cycle = "orbit_info/cycle_number",
                  date_start = "ancillary_data/data_start_utc",
                  date_end = "ancillary_data/data_end_utc")
  }

  # Retrieve ancillary and metadata

  names_metadata <- c("start_region", "sc_orient", "rgt", "cycle", "date_start", "date_end")
  required_metadata <- as.character(unlist(required[names(required) %in% names_metadata]))
  check_metadata <-  paste0("/", required_metadata) %in% paste(h5list$group, h5list$name, sep = "/")


  if (any(check_metadata == FALSE)) {
    warning(paste0("Missing metadata in file: ", paste(required_metadata[!check_metadata], collapse = ","),". Attributes will be set to NA"))
    eval(parse(text = paste0(names_metadata[!check_metadata], "<- NA"))) #Initiate to NA
  }

  ## Strong/weak beams
  if ("sc_orient" %in% names_metadata[check_metadata]) {
    sc_orient <- as.numeric(rhdf5::h5read(file = file, name = required[["sc_orient"]]))

    if(sc_orient == 0) { # Backward
      list_strength <- list(weak = c("gt1r", "gt2r", "gt3r"),
                            strong = c("gt1l", "gt2l", "gt3l"))
    }else{ #Forward
      list_strength <- list(weak = c("gt1l", "gt2l", "gt3l"),
                            strong = c("gt1r", "gt2r", "gt3r"))
    }
  }else{
    list_strength <- NULL
  }

  ## Orbit direction
  if ("start_region" %in% names_metadata[check_metadata]){
    start_region <- as.numeric(rhdf5::h5read(file = file, name = required[["start_region"]]))

    if (start_region %in% c(1, 2, 3, 12, 13, 14)) {
      orbit_dir = "ascending"
    }else if (start_region %in% c(5,6,7,8,9,10)) {
        orbit_dir = "descending"
    }else{
      orbit_dir  = "poles"
    }
  }else{
    orbit_dir <- NA
  }

  if("cycle" %in% names_metadata[check_metadata]) cycle <- as.numeric(rhdf5::h5read(file = file, name = required[["cycle"]]))
  if("rgt" %in% names_metadata[check_metadata]) rgt <- as.numeric(rhdf5::h5read(file = file, name = required[["rgt"]]))
  if("date_start" %in% names_metadata[check_metadata]) date_start <- as.character(rhdf5::h5read(file = file, name = required[["date_start"]]))
  if("date_end" %in% names_metadata[check_metadata]) date_end <- as.character(rhdf5::h5read(file = file, name = required[["date_end"]]))

  # Check for photon data
  required_c <- as.character(unlist(required))
  required_c <- required_c[!required_c %in% required_metadata]

  check <- paste0("/", required_c) %in% paste(h5list$group, h5list$name, sep = "/")

  if (any(check == FALSE)) {
    beam = beam[!(beam %in% unique(stringr::str_extract(required_c[!check],"gt..")))]
    warning(paste0("Missing dataset in file: ", paste(required_c[!check], collapse = ","),". Ignoring these beams"))
  }

  # -- Deprecated-- Retrieve information from file, not name
  # if (product == "ATL03") {
  #   test_format <- icesat2R::check_ATL_naming(file, product = "ATL03")
  # }else if (product == "ATL08"){
  #   test_format <- icesat2R::check_ATL_naming(file, product = "ATL08")
  #   }
  # if(!test_format) {
  #   warning("file does not follow the conventional ICESat-2 file naming. rgt, cycle, segment, dat and time (if set to TRUE) will not be written as attributes")
  #   rgt <- FALSE
  #   cycle <- FALSE
  #   segment = FALSE
  #   date <- FALSE
  #   time <- FALSE
  #
  #   h5file_filter <- tools::file_path_sans_ext(basename(file))
  # }else{
  #   h5file_filter <- stringr::str_extract(basename(file), pattern = "A.*5") #Keep only file name format
  # }

  # Create file to save outputs if not provided
  if(missing(filename_points)) filename_points = ""
  if (write & (gsub(" ", "", filename_points) == "" | is.null(filename_points)) ) {
    filename_points <- paste0(tools::file_path_sans_ext(basename(file)),"_points")
  }

  if(missing(filename_lines)) filename_lines = ""
  if (write & (gsub(" ", "", filename_lines) == "" | is.null(filename_lines)) ) {
    filename_lines <- paste0(tools::file_path_sans_ext(basename(file)),"_lines")
  }


  # Filter by beam strength
  if(length(beam_strength) == 1 & !is.null(list_strength)) {

    keep_beams <- beam %in% list_strength[[beam_strength]]

    if (sum(keep_beams) == 0) stop(sprintf("Selected beams are not %s",beam_strength))
    if (sum(keep_beams) != length(beam)) message(sprintf("Keeping only the following %s beams: %s", beam_strength, paste(beam[keep_beams], collapse = ",")))

    beam <- beam[keep_beams]
  }

  dat_ref_ph <- data.frame() #Stores all reference photons, grouped per beam
  att_table <- data.frame() #Attribute table for lines

  for ( n in beam) {

    #Beam strength
    if(!is.null(list_strength)) {
      if (n %in% list_strength[["weak"]]) {
        beam_strength = "weak"
      }else if (n %in% list_strength[["strong"]]) {
        beam_strength <- "strong"
      }
    }else{
      beam_strength <- NA
    }

    if (product == "ATL03") {
      ref_lat <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/geolocation/reference_photon_lat")))
      ref_lon <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/geolocation/reference_photon_lon")))
      segment_ph_cnt <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/geolocation/segment_ph_cnt")))
      segment_id <- as.character(rhdf5::h5read(file = file, name = paste0(n,"/geolocation/segment_id")))

      temp_df <- data.frame(lon = ref_lon, lat = ref_lat, beam = n,  beam_strength = beam_strength, segID = segment_id, segCount = segment_ph_cnt)
      temp_df <- temp_df[temp_df$segCount > 0, ]

      dat_ref_ph <- rbind(dat_ref_ph, temp_df)


      if (!all(segment_ph_cnt == 0)){
        toAdd <- data.frame(beam = n, beam_strength = beam_strength)
        rownames(toAdd) <- n
        att_table<- rbind(att_table, toAdd)
      }

    }else if (product == "ATL08") {
      ref_lat <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/land_segments/latitude")))
      ref_lon <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/land_segments/longitude")))

      temp_df <- data.frame(lon = ref_lon, lat = ref_lat, beam = n, beam_strength = beam_strength)
      dat_ref_ph <- rbind(dat_ref_ph, temp_df)

      toAdd <- data.frame(beam = n, beam_strength = beam_strength)
      rownames(toAdd) <- n
      att_table<- rbind(att_table, toAdd)
    }
  }

  dat_ref_ph$rgt <- rgt
  dat_ref_ph$odir <- orbit_dir
  dat_ref_ph$cycle <- cycle
  dat_ref_ph$segment <- start_region
  dat_ref_ph$date_start <- date_start
  dat_ref_ph$date_end <- date_end


  att_table$rgt <- rgt
  att_table$odir <- orbit_dir
  att_table$cycle <- cycle
  att_table$segment <- start_region
  att_table$date_start <- date_start
  att_table$date_end <- date_end

  ref_spdf <- sp::SpatialPointsDataFrame(coords = dat_ref_ph[,c("lon", "lat")], data = dat_ref_ph, proj4string = sp::CRS("+init=epsg:4326"))

  paths <- sp::split(ref_spdf, ref_spdf[["beam"]])

  test <- lapply(paths, sp::Line)
  test2 <- sapply(names(test), function(x){sp::Lines(list(test[[x]]), x)})
  test3 <- sp::SpatialLines(test2, proj4string = sp::CRS("+init=epsg:4326"))

  lines_df <- sp::SpatialLinesDataFrame(test3, data = att_table)

  if (write) {
    if(points) {
      rgdal::writeOGR(ref_spdf, dsn = outDir, layer = filename_points, driver = "ESRI Shapefile", ...)
    }
    if(lines) {
      rgdal::writeOGR(lines_df, dsn = outDir, layer = filename_lines, driver = "ESRI Shapefile", ...)
    }
  }

  return(list(if(lines)lines = lines_df else lines = NULL,
              if(points) points = ref_spdf else points = NULL,
              if(df) df = dat_ref_ph else df = NULL))

}
