#' Classify ATL03 photons to ATL08 classes
#'
#' @param atl03_h5 Path to ATL03 HDF file
#' @param atl08_h5 Path to ATL08 file
#' @param lat_range Optional. Vector indicating the range of latitude to process
#' @param beam Character vector indicating beams to process
#' @param beam_strength Character vector indicating the strength of beams to process
#' @param odir Directory where ouput files will be written
#'
#' @export

classify_ATL03 <- function(atl03_h5,
                           atl08_h5,
                           lat_range,
                           beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                           beam_strength = c("weak", "strong"),
                           odir) {

  if (!is.character(atl03_h5) | !tools::file_ext(atl03_h5) == "h5") {
    stop("atl03_h5 must be a path to a h5 file")
  }

  if (!is.character(atl08_h5) | !tools::file_ext(atl08_h5) == "h5") {
    stop("atl08_h5 must be a path to a h5 file")
  }

  # Check that both files are consistent
  atl_pattern = "(ATL[0-9]{2})_([0-9]{14}_[0-9]{8}_[0-9]{3}_[0-9]{2})"

  check_pattern <- stringr::str_match_all(c(basename(atl03_h5),
                                            basename(atl08_h5)),
                                          atl_pattern)

  if(!check_pattern[[1]][, 2] == "ATL03") {
    warning("Looks like atl03_h5 is not an ATL03 file (based on file name)")
  }

  if(!check_pattern[[2]][, 2] == "ATL08") {
      warning("Looks like atl08_h5 is not an ATL08 file (based on file name)")
    }

  if(!check_pattern[[1]][, 3] == check_pattern[[2]][, 3]) {
    warning("Looks like atl03_h5 and atl08_h5 are not from the same acquisition (based on file name)")
  }

  # Filter beams to analyze
  if (any(!beam_strength %in% c("weak", "strong"))) {
    stop("beam_strength must contain weak, beam or both ")
  }

  # Check beams to select
  suppressWarnings(sc_orient <- as.numeric(rhdf5::h5read(file = atl08_h5, name = "orbit_info/sc_orient")))

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

  atl03_h5list <- suppressWarnings(rhdf5::h5ls(atl03_h5))
  atl08_h5list <- suppressWarnings(rhdf5::h5ls(atl08_h5))

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

  check <- paste0("/",required) %in% paste(atl03_h5list$group, atl03_h5list$name, sep = "/")


  if (any(check == FALSE)) {
    missing_beams <- beam[(beam %in% unique(stringr::str_extract(required[!check],"gt..")))]
    warning(sprintf("Missing datasets in file: %s. Ignoring following beams: %s ", paste(required[!check], collapse = ","), paste(missing_beams, collapse = ", ") ))
    beam <- beam[!beam %in% missing_beams]
  }

  required <- c(paste0(beam,"/signal_photons/ph_segment_id"),
                paste0(beam,"/signal_photons/classed_pc_indx"),
                paste0(beam,"/signal_photons/classed_pc_flag"))

  check <- paste0("/",required) %in% paste(atl08_h5list$group, atl08_h5list$name, sep = "/")


  if (any(check == FALSE)) {
    missing_beams <- beam[(beam %in% unique(stringr::str_extract(required[!check],"gt..")))]
    warning(sprintf("Missing datasets in file: %s. Ignoring following beams: %s ", paste(required[!check], collapse = ","), paste(missing_beams, collapse = ", ") ))
    beam <- beam[!beam %in% missing_beams]
  }

  # Retrieve RGT, region and cycle
  suppressWarnings(rgt <- as.integer(rhdf5::h5read(file = atl08_h5, name = "orbit_info/rgt")))
  suppressWarnings(cycle <- as.integer(rhdf5::h5read(file = atl08_h5, name = "orbit_info/cycle_number")))
  suppressWarnings(region <- as.integer(rhdf5::h5read(file = atl08_h5, name = "ancillary_data/start_region")))

  out <- list() #Initiate list to return

  # Iterate through all beams
  for (n in beam) {

    print(sprintf("Working on beam %s", n))

    beam_strength_n <- names(list_strength)[unlist(lapply(list_strength, function(x) n %in% x))]

    # ATL03 data

    ## Index to keep within lat range

    suppressWarnings(lat_ph <- rhdf5::h5read(file = atl03_h5,
                            name = stringr::str_c(n, "/heights/lat_ph")))

    if (is.na(lat_range[1])) {
      lat_range[1] <- min(lat_ph)
    }

    if (is.na(lat_range[2])) {
      lat_range[2] <- max(lat_ph)
    }

    lat_idx <- which(lat_ph >= min(lat_range) & lat_ph <= max(lat_range))


    if (length(lat_idx) != 0) {

      # ATL03

      ## Heights group

      suppressWarnings(trg_fields <-  atl03_h5list %>%
                         dplyr::filter(stringr::str_detect(group, stringr::str_c(n, "/heights$")),
                                       otype == "H5I_DATASET") %>%
                         dplyr::mutate(full_name = stringr::str_c(group, name, sep = "/")) %>%
                         tidyr::separate(dim, into = c("dim1", "dim2"), sep = " x "))

      trg_fields_list <- lapply(1:length(trg_fields$name), function(x) {

        if (any(is.na(c(trg_fields$dim1[x], trg_fields$dim2[x])))) {
          trg_index <- list(lat_idx)
        }else{
          trg_index <- list(1:trg_fields$dim1[x], lat_idx)
        }

        suppressWarnings(out <- rhdf5::h5read(file = atl03_h5,
                             name = trg_fields$full_name[x],
                             index = trg_index))



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

      ## Geolocation group

      # Observations to keep

      # unique_segments <- unique(heights_df$ph_segment_id)
      #
      # segment_id <- rhdf5::h5read(atl03_h5,
      #                             name = stringr::str_c(n, "/geolocation/segment_id"))
      #
      # seg_idx <- which(segment_id %in% unique_segments)

      suppressWarnings(trg_fields <-  atl03_h5list %>%
                         dplyr::filter(stringr::str_detect(group, stringr::str_c(n, "/geolocation$")),
                                       otype == "H5I_DATASET") %>%
                         dplyr::mutate(full_name = stringr::str_c(group, name, sep = "/")) %>%
                         tidyr::separate(dim, into = c("dim1", "dim2"), sep = " x "))

      trg_fields_list <- lapply(1:length(trg_fields$name), function(x) {

        suppressWarnings(out <- rhdf5::h5read(file = atl03_h5,
                             name = trg_fields$full_name[x]))


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

      # ### Cumulated segment length
      # geolocation_df <- geolocation_df %>%
      #   dplyr::mutate(temp_length = c(0, segment_length[1:(nrow(.) - 1)]),
      #                 cum_segment_length = cumsum(temp_length)) %>%
      #   dplyr::select(!temp_length) %>%
      #   dplyr::relocate(cum_segment_length, .after = segment_length)

      # Add segment_id and photon count to height_df

      seg_idx <- unlist(lapply(geolocation_df$segment_ph_cnt, function(x) {
        if (x > 0) {
          seq(from = 1, to = x, by = 1)
        }}))

      seg_id_cum <- unlist(mapply(rep, geolocation_df$segment_id, geolocation_df$segment_ph_cnt))

      heights_df$ph_segment_id <- seg_id_cum[lat_idx]
      heights_df$classed_pc_indx <- seg_idx[lat_idx]

      # heights_df <- inner_join(heights_df, dplyr::select(geolocation_df, !delta_time), by = c("ph_segment_id" = "segment_id"))

      # Find UTM zone from first photon record

      zone <- find_UTM_zone(heights_df$lon_ph[1], heights_df$lat_ph[1])

      hemi <- find_UTM_hemisphere(heights_df$lat_ph[1])

      if (hemi == "N") {
        if(stringr::str_length(zone) == 1) zone = stringr::str_c("0", zone, sep = "")
        epsg_trg <- as.integer(stringr::str_c("326", zone, sep = ""))
      }else{
        if(stringr::str_length(zone) == 1) zone = stringr::str_c("0", zone, sep = "")
        epsg_trg <- as.integer(stringr::str_c("327", zone, sep = ""))
      }

      # Project to UTM Zone

      heights_sf <- sf::st_as_sf(heights_df,
                                 coords = c("lon_ph", "lat_ph"),
                                 remove = FALSE)
      sf::st_crs(heights_sf) <- 4326

      heights_sf_proj <- sf::st_transform(heights_sf, crs = epsg_trg)

      heights_sf_proj <- heights_sf_proj %>%
        dplyr::mutate(data.frame(sf::st_coordinates(heights_sf_proj))) %>%
        dplyr::rename(Easting = X,
               Northing = Y)


      # Find rotation matrix and rotation points

      rot_mat_list <- get_rot_matrix(xy_start = c(heights_sf_proj$Easting[1], heights_sf_proj$Northing[1]),
                                     xy_end = c(heights_sf_proj$Easting[nrow(heights_sf_proj)],
                                                heights_sf_proj$Northing[nrow(heights_sf_proj)]))

      # Calculate along and across track

      xy_rot <- get_along_distance(heights_sf_proj$Easting,
                                   heights_sf_proj$Northing,
                                   rot_mat_list$R_mat,
                                   rot_mat_list$xRotPt,
                                   rot_mat_list$yRotPt)

      heights_sf_proj <- cbind(heights_sf_proj, data.frame(xy_rot))

      # ATL08 data

      # Filter based on segment id in ATL03

      unique_segments <- unique(heights_sf_proj$ph_segment_id)

      suppressWarnings(atl08_seg_beg <- rhdf5::h5read(atl08_h5,
                                     name = stringr::str_c(n, "/land_segments/segment_id_beg")))

      suppressWarnings(atl08_seg_end <- rhdf5::h5read(atl08_h5,
                                     name = stringr::str_c(n, "/land_segments/segment_id_end")))

      all_segments <- atl08_seg_beg[1]:atl08_seg_end[length(atl08_seg_end)]

      seg_idx <- which(atl08_seg_beg >= unique_segments[1] & atl08_seg_end <= unique_segments[length(unique_segments)])

      suppressWarnings(trg_fields <-  atl08_h5list %>%
        dplyr::filter(stringr::str_detect(group, stringr::str_c(n, "/land_segments$")),
                      otype == "H5I_DATASET") %>%
        dplyr::mutate(full_name = stringr::str_c(group, name, sep = "/")) %>%
        tidyr::separate(dim, into = c("dim1", "dim2"), sep = " x "))

      trg_fields_list <- lapply(1:length(trg_fields$name), function(x) {

        if (any(is.na(c(trg_fields$dim1[x], trg_fields$dim2[x])))) {
          trg_index <- list(seg_idx)
        }else{
          trg_index <- list(1:trg_fields$dim1[x], seg_idx)
        }

        suppressWarnings(out <- rhdf5::h5read(file = atl08_h5,
                             name = trg_fields$full_name[x],
                             index = trg_index))

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

      land_df <- as.data.frame(trg_fields_list)

      # Canopy products

      suppressWarnings(trg_fields <-  atl08_h5list %>%
        dplyr::filter(stringr::str_detect(group, stringr::str_c(n, "/land_segments/canopy$")),
                      otype == "H5I_DATASET") %>%
        dplyr::mutate(full_name = stringr::str_c(group, name, sep = "/"))  %>%
        tidyr::separate(dim, into = c("dim1", "dim2"), sep = " x "))

      trg_fields_list <- lapply(1:length(trg_fields$name), function(x) {

        if (any(is.na(c(trg_fields$dim1[x], trg_fields$dim2[x])))) {
          trg_index <- list(seg_idx)
        }else{
          trg_index <- list(1:trg_fields$dim1[x], seg_idx)
        }

        suppressWarnings(out <- rhdf5::h5read(file = atl08_h5,
                             name = trg_fields$full_name[x],
                             index = trg_index))

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

      canopy_df <- as.data.frame(trg_fields_list)

      # Terrain products

      suppressWarnings(trg_fields <-  atl08_h5list %>%
                         dplyr::filter(stringr::str_detect(group, stringr::str_c(n, "/land_segments/terrain$")),
                                       otype == "H5I_DATASET") %>%
                         dplyr::mutate(full_name = stringr::str_c(group, name, sep = "/")) %>%
                         tidyr::separate(dim, into = c("dim1", "dim2"), sep = " x "))

      trg_fields_list <- lapply(1:length(trg_fields$name), function(x) {

        if (any(is.na(c(trg_fields$dim1[x], trg_fields$dim2[x])))) {
          trg_index <- list(seg_idx)
        }else{
          trg_index <- list(1:trg_fields$dim1[x], seg_idx)
        }

        suppressWarnings(out <- rhdf5::h5read(file = atl08_h5,
                             name = trg_fields$full_name[x],
                             index = trg_index))

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

      terrain_df <- as.data.frame(trg_fields_list)

      # Combine info from all groups

      atl08_df <- cbind(data.frame(cycle = cycle,
                                   region = region,
                                   beam = n,
                                   beam_strength = beam_strength_n),
                        land_df,
                        canopy_df,
                        terrain_df)

      colnames(atl08_df) <- stringr::str_replace(colnames(atl08_df),
                                                 pattern = "\\.",
                                                 replacement = "_")

      atl08_sf <- sf::st_as_sf(atl08_df,
                               coords = c("longitude", "latitude"),
                               remove = FALSE,
                               crs = sf::st_crs(4326))

      atl08_sf_proj <- sf::st_transform(atl08_sf, epsg_trg)

      atl08_sf_proj <- dplyr::mutate(atl08_sf_proj, data.frame(sf::st_coordinates(geometry))) %>%
        dplyr::rename(Easting = X,
               Northing = Y)

      # Along distance of ATL08 center points

      xy_rot <- get_along_distance(atl08_sf_proj$Easting,
                                   atl08_sf_proj$Northing,
                                   rot_mat_list$R_mat,
                                   rot_mat_list$xRotPt,
                                   rot_mat_list$yRotPt)

      atl08_sf_proj <- cbind(atl08_sf_proj, data.frame(xy_rot))


      # ATL03 classified photons

      suppressWarnings(ph_segment_id <- rhdf5::h5read(atl08_h5,
                                     name = stringr::str_c(n, "/signal_photons/ph_segment_id")))

      seg_idx <- which(ph_segment_id %in% unique_segments)


      suppressWarnings(trg_fields <-  atl08_h5list %>%
                         dplyr::filter(stringr::str_detect(group, stringr::str_c(n, "/signal_photons$")),
                                       otype == "H5I_DATASET") %>%
                         dplyr::mutate(full_name = stringr::str_c(group, name, sep = "/")) %>%
                         tidyr::separate(dim, into = c("dim1", "dim2"), sep = " x "))

      trg_fields_list <- lapply(1:length(trg_fields$name), function(x) {

        if (any(is.na(c(trg_fields$dim1[x], trg_fields$dim2[x])))) {
          trg_index <- list(seg_idx)
        }else{
          trg_index <- list(1:trg_fields$dim1[x], seg_idx)
        }

        suppressWarnings( out <- rhdf5::h5read(file = atl08_h5,
                             name = trg_fields$full_name[x],
                             index = trg_index))

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

      signal_ph_df <- as.data.frame(trg_fields_list)

      signal_ph_df <- signal_ph_df %>%
        dplyr::mutate(classed_pc_flag = factor(classed_pc_flag,
                                        levels = c(0, 1, 2, 3),
                                        labels = c("noise", "ground", "canopy", "top of canopy")))

      # Classify ATL03 product

      classified_atl03 <- dplyr::left_join(heights_sf_proj,
                                            dplyr::select(signal_ph_df, !delta_time),
                                            by = c("ph_segment_id", "classed_pc_indx"))

      # There might be some photons referenced in ATL03 but not in ATL08 (starting segment)

      if (!missing(odir)) {

        if (!dir.exists(odir)) {
          dir.create(odir)
        }

        # Write ATL08

        readr::write_csv(sf::st_drop_geometry(atl08_sf_proj),
                         file = file.path(odir, stringr::str_c(tools::file_path_sans_ext(basename(atl08_h5)), "_", toupper(n), "_", toupper(beam_strength_n),".csv")))

        sf::st_write(atl08_sf_proj,
                     dsn = file.path(odir, stringr::str_c(tools::file_path_sans_ext(basename(atl08_h5)), "_", toupper(n), "_", toupper(beam_strength_n),".gpkg")),
                     quiet = TRUE,
                     delete_layer = TRUE)

        # Write Classified ATL03

        readr::write_csv(sf::st_drop_geometry(classified_atl03),
                         file = file.path(odir, stringr::str_c(tools::file_path_sans_ext(basename(atl03_h5)), "_", toupper(n), "_", toupper(beam_strength_n),".csv")))


      }

      out[[n]] <- list(ATL08 = sf::st_drop_geometry(atl08_sf_proj),
                       ATL08_sf = atl08_sf_proj,
                       ATL03 = classified_atl03)

    }else{
      warning(sprintf("No photons to return within specified lat_range for beam %s", n))
    }
  }

  rhdf5::h5closeAll()
  return(out)
}

