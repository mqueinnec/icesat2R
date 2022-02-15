#' Extract individual photons from an ATL03 file
#'
#' Returns a data.frame with lat, lon and height of each photon and their associated segement and unique ID
#'
#'@param file
#'@param beam
#'@param beam_strength
#'@param lat_range
#'@param seg_filter Vector of integers. Filter by 20 m segments of the ATL03 product
#'@param id_sep
#'
#'@export
#'@importFrom dplyr %>%

getATL03ph <- function(file,
                       beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                       beam_strength = c("weak", "strong"),
                       lat_range,
                       seg_filter,
                       id_sep = "_") {

  if(!is.character(file)) {
    stop("file must be a HDF5 file name (extension h5)")
  }else{
    file <- normalizePath(file)
    if(!file.exists(file)) stop("file  does not exist")
    if(tools::file_ext(file) != "h5") stop("file must be a HDF5 file (h5 extension)")
  }

  if(!is.character(id_sep)) {
      stop("id_sep must be a character")
  }

  # Filter beams to analyze
  if (any(!beam_strength %in% c("weak", "strong"))) {
    stop("beam_strength must contain weak, beam or both ")
  }


  if(length(beam_strength) == 1) {
      # Check beams to select
      sc_orient <- as.numeric(rhdf5::h5read(file = file, name = "orbit_info/sc_orient"))

      if(sc_orient == 0) { # Backward
        list_strength <- list(weak = c("gt1r", "gt2r", "gt3r"),
                              strong = c("gt1l", "gt2l", "gt3l"))
      }else{ #Forward
        list_strength <- list(weak = c("gt1l", "gt2l", "gt3l"),
                              strong = c("gt1r", "gt2r", "gt3r"))
      }

      keep_beams <- beam %in% list_strength[[beam_strength]]

      if (sum(keep_beams) == 0) stop(sprintf("Selected beams are not %s",beam_strength))
      if (sum(keep_beams) != length(beam)) message(sprintf("Keeping only the following %s beams: %s", beam_strength, paste(beam[keep_beams], collapse = ",")))

      beam <- beam[keep_beams]

  }

  # Setting up min and max lat
  if(!missing(lat_range)) {
    if (length(lat_range) != 2) stop("lat_range must have a min and and max value")
    if(lat_range[1] >= lat_range[2]) stop("First element of lat_range must be lower than second element")
  }else{
    lat_range <- c(NA, NA)
  }

  # Filter by segments
  if(!missing(seg_filter)) {
    if (!is(seg_filter, "integer")) stop("seg_filter must be a vector of integer")
  }else{
    seg_filter <- NA
  }

  # Check that the necessary data is available
  h5list <- rhdf5::h5ls(file)

  required <- c(paste0(beam,"/geolocation/segment_ph_cnt"),
                paste0(beam,"/geolocation/segment_id"),
                paste0(beam,"/heights/lat_ph"),
                paste0(beam,"/heights/lat_ph"),
                paste0(beam,"/heights/h_ph"),
                paste0(beam,"/heights/dist_ph_along"),
                paste0(beam,"/heights/dist_ph_across"),
                paste0(beam,"/heights/delta_time"),
                paste0(beam, "/geolocation/reference_photon_index"),
                paste0(beam, "/geolocation/segment_length"))

  check <- paste0("/",required) %in% paste(h5list$group, h5list$name, sep = "/")

  if (any(check == FALSE)) {
    missing_beams <- beam[(beam %in% unique(stringr::str_extract(required[!check],"gt..")))]
    warning(sprintf("Missing datasets in file: %s. Ignoring following beams: %s ", paste(required[!check], collapse = ","), paste(missing_beams, collapse = ", ") ))
    beam <- beam[!beam %in% missing_beams]
  }

  # Retrieve RGT, region and cycle
  rgt <- as.numeric(rhdf5::h5read(file = file, name = "orbit_info/rgt"))
  cycle <- as.numeric(rhdf5::h5read(file = file, name = "orbit_info/cycle_number"))
  region <- as.numeric(rhdf5::h5read(file = file, name = "ancillary_data/start_region"))

  out <- list() #Initiate list to return
  for (n in beam) {

    lat_ph <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/heights/lat_ph")))
    if (any(is.na(lat_range))) {
      lat_range <- c(min(lat_ph), max(lat_ph))
    }
    # Get index of data to keep
    lat_idx <- which(lat_ph >= lat_range[1] & lat_ph<=lat_range[2])

    if (length(lat_idx) > 0) {

      lat_ph <- lat_ph[lat_idx]

      lon_ph <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/heights/lon_ph"),
                                         index = list(lat_idx)))
      h_ph <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/heights/h_ph"),
                                       index = list(lat_idx)))

      segcnt <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/geolocation/segment_ph_cnt")))
      segID <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/geolocation/segment_id")))
      ref_ph_idx <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/geolocation/reference_photon_index")))
      seg_length <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/geolocation/segment_length")))
      # We initiate the first photon along distance to 0
      seg_length <- c(0, seg_length[1:(length(segID) - 1)])
      # seg_length <- seg_length[1:length(segID)]
      cum_seg_length <- cumsum(seg_length)
      # And remove segments with no photons
      #cum_seg_length <- cum_seg_length[segcnt > 0]


      dist_along <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/heights/dist_ph_along")))
      dist_along <- dist_along[lat_idx]
      dist_along <- dist_along - dist_along[1]
      dist_across <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/heights/dist_ph_across")))
      dist_across <- dist_across[lat_idx]

      # Time tag (ATLAS shot) for each photon
      delta_time <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/heights/delta_time")))
      delta_time <- delta_time[lat_idx]

      #PCE frame
      pce_mframe <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/heights/pce_mframe_cnt")))
      pce_mframe <- pce_mframe[lat_idx]

      # Pulse ID
      ph_id_pulse <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/heights/ph_id_pulse")))
      ph_id_pulse <- ph_id_pulse[lat_idx]

      # Photon signal confidence
      signal_conf <- matrix(as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/heights/signal_conf_ph"))), ncol = 5, byrow = TRUE)

      signal_conf_land <- signal_conf[,1]
      signal_conf_land <- signal_conf_land[lat_idx]

      #SegID of each photon
      seg_id <- unlist(mapply(rep, segID, segcnt))
      seg_id <- seg_id[lat_idx]
      #Remove segments with no photons and not within lat bounds from cum_seg_length
      cum_seg_length <- cum_seg_length[which(segID %in% unique(seg_id))]
      #Photon index in segment
      seg_idx <- unlist(lapply(segcnt, function(x) {
        if (x > 0) {
          seq(from = 1, to = x, by = 1)
        }}))
      seg_idx <- seg_idx[lat_idx]
      #Unique photon ID
      ph_id <- paste(rgt,cycle,region,n,seg_id, seg_idx, sep = id_sep)

      #Identify reference photon in segment
      ref_ph <- unlist(mapply(function(x, y) {
        idx <- rep(FALSE, x)
        idx[y] <- TRUE
        idx }, segcnt, ref_ph_idx))
      ref_ph <- ref_ph[lat_idx]

      # Cumulative along distance
      cum_dist_along <- as.vector(unlist(mapply(function(x, y) x + y, split(dist_along, seg_id), as.list(cum_seg_length))))

      #Quick check to see if everything consistent

      if (any(length(ph_id) != c(length(lat_ph), length(lon_ph), length(h_ph)))) {
        stop("The length of photon ID doesn't match the number of photons")
      }

      df_beam <- data.frame(lon = lon_ph,
                              lat = lat_ph,
                              h = h_ph,
                              dist_along_seg = dist_along,
                              cum_dist_along = cum_dist_along,
                              dist_across = dist_across,
                              delta_time = delta_time,
                            pce_mframe = pce_mframe,
                            ph_id_pulse = ph_id_pulse,
                             signal_conf_land = signal_conf_land,
                              beam = n,
                              rgt = rgt,
                              cycle = cycle,
                              region = region,
                              seg_id = seg_id,
                              seg_idx = seg_idx,
                              ph_id = ph_id,
                              ref_ph = ref_ph)

      # Make spatial object and transform coordinates
      # df_beam <- sf::st_as_sf(df_beam,
      #                         coords = c("lon", "lat"),
      #                         crs = 4326, remove = FALSE)
      #
      #
      # df_beam <- sf::st_transform(df_beam, crs = UTM_proj(longitude = df_beam$lon[1], latitude = df_beam$lat[1]))
      #
      # UTM_coords <- sf::st_coordinates(df_beam)
      # colnames(UTM_coords) <- c("Easting", "Northing")
      #
      # df_beam <- cbind(df_beam, UTM_coords)

      # Filter by segment
      if(!any(is.na(seg_filter))) {
        df_beam <- df_beam %>%
          filter(seg_id %in% seg_filter)
      }
      out[[n]] <- df_beam
    }else{
      warning(sprintf("No photons to return within specified lat_range for beam %s", n))
    }

  }
  return(out)
}
