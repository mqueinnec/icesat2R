#' Extract individual classified photons from an ATL08 file
#'
#' Returns a data.frame with lat, lon and height of each photon and their associated segement and unique ID
#'@param file
#'@param beam
#'@param beam_strength
#'@param lat_range
#'@param seg_filter Vector of integers. Filter by 20 m segments of the ATL03 product
#'@param id_sep
#'@param flags Character. ATL08 segement-level flag to return. All photons within the segments will be assigned the same flag. Impemented flags are: "snow", "night", "msw".
#'
#'@export

getATL08ph <- function(file,
                       beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                       beam_strength = c("weak", "strong"),
                       lat_range,
                       seg_filter,
                       id_sep = "_",
                       flags = c("snow", "night", "msw")) {

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

  # Check flags
  validFlags <- c("snow", "night", "msw")
  if (any(!flags %in% validFlags)) {
    skipFlags <- which(!flags %in% validFlags)
    warning(sprintf("Following flags not implemented and skipped: %s", paste(flags[skipFlags], collapse = ",")))

    flags <- flags[-skipFlags]
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

  required <- c(paste0(beam,"/signal_photons/ph_segment_id"),
                paste0(beam,"/signal_photons/classed_pc_indx"),
                paste0(beam,"/signal_photons/classed_pc_flag"))

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

    # Find all 100-m segement center points within lat_range
    lat_center <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/land_segments/latitude")))
    if (any(is.na(lat_range))) {
      lat_range <- c(min(lat_center), max(lat_center))
    }
    lat_idx <- which(lat_center >= lat_range[1] & lat_center <= lat_range[2])

    if (length(lat_idx) != 0) {
      beg_seg <- as.integer(rhdf5::h5read(file = file, name = paste0(n,"/land_segments/segment_id_beg")))
      end_seg <- as.integer(rhdf5::h5read(file = file, name = paste0(n,"/land_segments/segment_id_end")))

      all_seg <- beg_seg[lat_idx[1]]:end_seg[lat_idx[length(lat_idx)]]

      # Filter by segment ID
      if(!any(is.na(seg_filter))) {
        all_seg <- all_seg[all_seg %in% seg_filter]
      }

      if(length(all_seg) == 0) {
        warning("No segments match lat range or seg_filter")
      }

      segID <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/signal_photons/ph_segment_id")))
      keep_ph <- which(segID %in% all_seg)
      segID <- segID[keep_ph]

      ph_idx <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/signal_photons/classed_pc_indx"), index = list(keep_ph)))
      ph_flag <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/signal_photons/classed_pc_flag"), index = list(keep_ph)))
      ph_flag <- factor(ph_flag, levels = c(0,1,2,3), labels = c("noise", "ground", "canopy", "top of canopy"))

      ph_id <- paste(rgt,cycle,region,n,segID, ph_idx, sep = id_sep)

      #Assign flags
      flags_dat <- list()
      for( f in flags) {
        if ( f == "snow") {
          required <- paste0(n,"/land_segments/segment_snowcover")
          check <- paste0("/",required) %in% paste(h5list$group, h5list$name, sep = "/")
          flags_dat$snow_flag <- as.integer(rhdf5::h5read(file = file, name = required))
        } else if (f == "night") {
          required <- paste0(n,"/land_segments/night_flag")
          check <- paste0("/",required) %in% paste(h5list$group, h5list$name, sep = "/")
          flags_dat$night_flag <- as.integer(rhdf5::h5read(file = file, name = required))
        } else if (f == "msw") {
          required <- paste0(n,"/land_segments/msw_flag")
          check <- paste0("/",required) %in% paste(h5list$group, h5list$name, sep = "/")
          flags_dat$msw_flag <- as.integer(rhdf5::h5read(file = file, name = required))
        }
      }

      count <- 1
      flags_df <- data.frame()
      for (s in seq_len(length(beg_seg))) {
        seg_temp <- list(seg_id = beg_seg[s]:end_seg[s])
        flags_temp <- lapply(flags_dat, function(x) x[count])
        flags_df <- rbind(flags_df, as.data.frame(do.call(cbind, c(seg_temp, flags_temp))))
        count <- count + 1
      }

      out[[n]] <- data.frame(rgt = rgt,
                             cycle = cycle,
                             region = region,
                             beam = n,
                             ph_id = ph_id,
                             ph_class = ph_flag,
                             seg_id = segID,
                             seg_idx = ph_idx)

      out[[n]] <- merge(out[[n]], flags_df, by = "seg_id")


    }else{
      warning(sprintf("No photons to return within specified lat_range for beam %s", n))
    }
  }

  return(out)
}
