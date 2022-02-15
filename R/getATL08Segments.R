#' Get the segment of each ATL08 point
#'
#'@param file Path to ATL03 or ATL08 HDF5 file (.h5 extension)
#'@param beam Character. List of beams for which segments to be extracted
#'@param rgt Logical. Should RGT number be added as a field to output feature?
#'@param cycle Logical. Should cycle number be added as a field to output feature?
#'@param segment Logical. Should segment number be added as a field to output feature?
#'@param date Logical. Should date be added as a field to output feature?
#'@param time Logical. Should time be added as a field to output feature?
#'@param ... Parameters passed to writeOGR (e.g. overwrite_layer = TRUE)
#'
#'@export
#'@importFrom dplyr %>%


getATL08Segments <- function(file,
                             beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                             beam_strength = c("weak", "strong"),
                             crs_proj = 26917,
                             odir) {

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


  required <- list(lat = paste0(beam,"/land_segments/latitude"),
                     lon = paste0(beam,"/land_segments/longitude"),
                     start_region = "ancillary_data/start_region",
                     sc_orient = "orbit_info/sc_orient",
                     rgt = "orbit_info/rgt",
                     cycle = "orbit_info/cycle_number",
                     date_start = "ancillary_data/data_start_utc",
                     date_end = "ancillary_data/data_end_utc")


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

   cycle <- as.numeric(rhdf5::h5read(file = file, name = required[["cycle"]]))
   rgt <- as.integer(rhdf5::h5read(file = file, name = required[["rgt"]]))
   rgt <- sprintf("%04d", rgt)
  date_start <- as.character(rhdf5::h5read(file = file, name = required[["date_start"]]))
  date_end <- as.character(rhdf5::h5read(file = file, name = required[["date_end"]]))

  # Check for photon data
  required_c <- as.character(unlist(required))
  required_c <- required_c[!required_c %in% required_metadata]

  check <- paste0("/", required_c) %in% paste(h5list$group, h5list$name, sep = "/")

  if (any(check == FALSE)) {
    beam = beam[!(beam %in% unique(stringr::str_extract(required_c[!check],"gt..")))]
    warning(paste0("Missing dataset in file: ", paste(required_c[!check], collapse = ","),". Ignoring these beams"))
  }

  # Filter by beam strength
  if(length(beam_strength) == 1 & !is.null(list_strength)) {

    keep_beams <- beam %in% list_strength[[beam_strength]]

    if (sum(keep_beams) == 0) stop(sprintf("Selected beams are not %s",beam_strength))
    if (sum(keep_beams) != length(beam)) message(sprintf("Keeping only the following %s beams: %s", beam_strength, paste(beam[keep_beams], collapse = ",")))

    beam <- beam[keep_beams]
  }

  #Initiate ref_pts and and seg_line
  ref_pts <- list()
  seg_lines <- list()

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

    ref_lat <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/land_segments/latitude")))

    #Make sure there are at least 3 points
    if (length(ref_lat) > 2) {
      ref_lon <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/land_segments/longitude")))

      segid_beg <- as.integer(rhdf5::h5read(file = file, name = paste0(n,"/land_segments/segment_id_beg")))
      segid_end <- as.integer(rhdf5::h5read(file = file, name = paste0(n,"/land_segments/segment_id_end")))

      id <- paste(rgt,segid_beg, segid_end, sep = "_")

      attributes_df <- data.frame(lon = ref_lon,
                                  lat = ref_lat,
                                  seg_beg = segid_beg,
                                  seg_end = segid_end,
                                  seg_id = id,
                                  beam = n,
                                  beam_strength = beam_strength,
                                  rgt = rgt,
                                  cycle = cycle,
                                  dir = orbit_dir,
                                  date_st = date_start,
                                  date_end = date_end)

      ref_pts[[n]] <- sf::st_as_sf(attributes_df,
                                   coords = c("lon", "lat"),
                                   crs = 4326)

      ref_pts[[n]] <- sf::st_transform(ref_pts[[n]], crs = crs_proj)

      ref_buffer <- sf::st_buffer(ref_pts[[n]], dist = units::set_units(50, m))

      ref_line <- ref_pts[[n]] %>%
        dplyr::group_by() %>%
        dplyr::summarise(do_union = FALSE) %>%
        sf::st_cast("LINESTRING")

      seg_lines[[n]] <- sf::st_intersection(ref_line, ref_buffer)

      # Remove first and last segment (half-segments)

      seg_lines[[n]] <- seg_lines[[n]][c(-1, -NROW(seg_lines[[n]])),]
    }
  }

  # Bind beams
  ref_pts <- do.call(rbind, ref_pts)
  seg_lines <- do.call(rbind, seg_lines)

    #Save to SHP
    sf::st_write(seg_lines,
             dsn = odir,
             layer = paste0(tools::file_path_sans_ext(basename(file)),
             "_segments"),
             driver = "ESRI Shapefile",
             delete_layer = TRUE)

    sf::st_write(ref_pts,
             dsn = odir,
             layer = paste0(tools::file_path_sans_ext(basename(file)),
                            "_pts"),
             driver = "ESRI Shapefile",
             delete_layer = TRUE)

}
