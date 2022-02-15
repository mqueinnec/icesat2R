#' Read ATL03 h5 file as a data frame
#'
#' @param file Path to h5 file
#' @param beam Character vector indicating beams to process
#' @param beam_strength Character vector indicating the strength of beams to process
#' @param lat_range Numeric vector. Lower and upper latitude to return
#' @param odir_atl03 Character. Output directory of ATL03 product
#'
#' @export
#'

read_ATL03 <- function(file,
                       beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                       beam_strength = c("weak", "strong"),
                       lat_range,
                       odir_atl03) {

  # Check file input
  if(!is.character(file) | !tools::file_ext(file) == "h5") {
    stop("file must be a HDF5 file name (extension h5)")
  }else{
    file <- normalizePath(file)
    if(!file.exists(file)) stop("file  does not exist")
  }

  # Filter beams to analyze
  if (any(!beam_strength %in% c("weak", "strong"))) {
    stop("beam_strength must contain weak, beam or both ")
  }

  # Check beams to select
  sc_orient <- as.numeric(rhdf5::h5read(file = file, name = "orbit_info/sc_orient"))

  if(sc_orient == 0) { # Backward
    list_strength <- list(weak = c("gt1r", "gt2r", "gt3r"),
                          strong = c("gt1l", "gt2l", "gt3l"))
  }else{ #Forward
    list_strength <- list(weak = c("gt1l", "gt2l", "gt3l"),
                          strong = c("gt1r", "gt2r", "gt3r"))
  }

  if(length(beam_strength) == 1) {

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
    warning(sprintf("Missing datasets in file: %s. Ignoring following beams: %s. ", paste(required[!check], collapse = ","), paste(missing_beams, collapse = ", ")))
    beam <- beam[!beam %in% missing_beams]
  }

  # Retrieve RGT, region and cycle
  rgt <- as.numeric(rhdf5::h5read(file = file, name = "orbit_info/rgt"))
  cycle <- as.numeric(rhdf5::h5read(file = file, name = "orbit_info/cycle_number"))
  region <- as.numeric(rhdf5::h5read(file = file, name = "ancillary_data/start_region"))

  out <- list() #Initiate list to return

  for (n in beam) {

    beam_strength_n <- names(list_strength)[unlist(lapply(list_strength, function(x) n %in% x))]

    lat_ph <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/heights/lat_ph")))
    if (any(is.na(lat_range))) {
      lat_range <- c(min(lat_ph), max(lat_ph))
    }
    # Get index of data to keep
    lat_idx <- lat_ph >= lat_range[1] & lat_ph<=lat_range[2]

    if (sum(lat_idx) > 0) {

      # Geolocation

      trg_fields <-  h5list %>%
        dplyr::filter(stringr::str_detect(group, stringr::str_c(n, "/geolocation$")),
                      otype == "H5I_DATASET") %>%
        dplyr::mutate(full_name = stringr::str_c(group, name, sep = "/"))

      trg_fields_list <- lapply(1:length(trg_fields$name), function(x) {
        out <- rhdf5::h5read(file = file,
                             name = trg_fields$full_name[x])



        if(trg_fields$dclass[x] == "FLOAT") {
          if(length(dim(out)) > 1) {
            out <- as.data.frame(apply(out, 1, as.numeric))
            colnames(out) <- stringr::str_c("X", 1:ncol(out))

          }else{
            out <- as.numeric(out)
          }
        }else if (trg_fields$dclass[x] == "INTEGER") {
          if(length(dim(out)) > 1) {
            out <- as.data.frame(apply(out, 1, as.integer))
            colnames(out) <- stringr::str_c("X", 1:ncol(out))
          }else{
            out <- as.integer(out)
          }
        }

        return(out)

      })

      names(trg_fields_list) <- trg_fields$name

      geolocation_df <- as.data.frame(trg_fields_list)

      geolocation_df <- geolocation_df %>%
        dplyr::mutate(temp_length = c(0, segment_length[1:(nrow(.) - 1)]),
               cum_segment_length = cumsum(temp_length)) %>%
        dplyr::select(!temp_length) %>%
        dplyr::relocate(cum_segment_length, .after = segment_length)

      # Heights

      trg_fields <-  h5list %>%
        dplyr::filter(stringr::str_detect(group, stringr::str_c(n, "/heights$")),
                      otype == "H5I_DATASET") %>%
        dplyr::mutate(full_name = stringr::str_c(group, name, sep = "/"))

      trg_fields_list <- lapply(1:length(trg_fields$name), function(x) {
        out <- rhdf5::h5read(file = file,
                             name = trg_fields$full_name[x])



        if(trg_fields$dclass[x] == "FLOAT") {
          if(length(dim(out)) > 1) {
            out <- as.data.frame(apply(out, 1, as.numeric))
            colnames(out) <- stringr::str_c("X", 1:ncol(out))

          }else{
            out <- as.numeric(out)
          }
        }else if (trg_fields$dclass[x] == "INTEGER") {
          if(length(dim(out)) > 1) {
            out <- as.data.frame(apply(out, 1, as.integer))
            colnames(out) <- stringr::str_c("X", 1:ncol(out))
          }else{
            out <- as.integer(out)
          }
        }

        return(out)

      })

      names(trg_fields_list) <- trg_fields$name

      heights_df <- as.data.frame(trg_fields_list)

      # Join geolocation and height

      seg_idx <- unlist(lapply(geolocation_df$segment_ph_cnt, function(x) {
        if (x > 0) {
          seq(from = 1, to = x, by = 1)
        }}))

      seg_id_cum <- unlist(mapply(rep, geolocation_df$segment_id, geolocation_df$segment_ph_cnt))

      heights_df$ph_segment_id <- seg_id_cum
      heights_df$classed_pc_indx <- seg_idx

      heights_df <- inner_join(heights_df, select(geolocation_df, !delta_time), by = c("ph_segment_id" = "segment_id"))


      # ----

      atl03_out <- cbind(data.frame(rgt = rgt,
                                    cycle = cycle,
                                    region = region,
                                    beam = n,
                                    beam_strength = beam_strength_n),
                         heights_df)

      colnames(atl03_out) <- stringr::str_replace(colnames(atl03_out),
                                                     pattern = "\\.",
                                                     replacement = "_")

      # keep only within lat_range
      atl03_out <- atl03_out[lat_idx, ]

      # Along distance

      atl03_out <- atl03_out %>%
        mutate(cum_dist_along = dist_ph_along + cum_segment_length) %>%
        relocate(cum_dist_along, .after = dist_ph_along)


      if(!missing(odir_atl03)) {
        if(!dir.exists(odir_atl03)) {
          dir.create(odir_atl03)
        }

        readr::write_csv(atl03_out,
                         file = file.path(odir_atl03, stringr::str_c(tools::file_path_sans_ext(basename(file)), "_", toupper(n), "_", toupper(beam_strength_n),".csv")))
      }

      out[[n]] <- atl03_out
    }else{
      warning(sprintf("No photons to return within specified lat_range for beam %s", n))
    }

  }

  return(out)
}
