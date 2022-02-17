#' Read ATL08 h5 file as a data frame
#'
#' @param atl08_h5 Path to h5 file
#' @param beam Character vector indicating beams to process
#' @param beam_strength Character vector indicating the strength of beams to process
#' @param lat_range Numeric vector. Lower and upper latitude to return
#' @param odir Character. Output directory
#' @param atl03_pc Logical indicating if the classification of ATL03 photons should be returned
#'
#' @export

read_ATL08 <- function(atl08_h5,
                       beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                       beam_strength = c("weak", "strong"),
                       lat_range,
                       odir,
                       atl03_pc = FALSE) {

  # Check file input
  if (!is.character(atl08_h5) | !tools::file_ext(atl08_h5) == "h5") {
    stop("atl08_h5 must be a path to a h5 file")
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

  # Check that the necessary data is available
  atl08_h5list <- suppressWarnings(rhdf5::h5ls(atl08_h5))

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

  for ( n in beam) {

    beam_strength_n <- names(list_strength)[unlist(lapply(list_strength, function(x) n %in% x))]

    print(sprintf("Working on beam %s (%s)", n, beam_strength_n))

    # Find all 100-m segment center points within lat_range
    lat_center <- as.numeric(rhdf5::h5read(file = atl08_h5, name = paste0(n,"/land_segments/latitude")))

    if (is.na(lat_range[1])) {
      lat_range[1] <- min(lat_center)
    }

    if (is.na(lat_range[2])) {
      lat_range[2] <- max(lat_center)
    }

    # Index of points to keep
    lat_idx <- which(lat_center >= lat_range[1] & lat_center <= lat_range[2])

    if (sum(lat_idx) != 0) {


      # Land segment flags and info
      suppressWarnings(trg_fields <-  atl08_h5list %>%
                         dplyr::filter(stringr::str_detect(group, stringr::str_c(n, "/land_segments$")),
                                       otype == "H5I_DATASET") %>%
                         dplyr::mutate(full_name = stringr::str_c(group, name, sep = "/")) %>%
                         tidyr::separate(dim, into = c("dim1", "dim2"), sep = " x "))

      trg_fields_list <- lapply(1:length(trg_fields$name), function(x) {

        if (any(is.na(c(trg_fields$dim1[x], trg_fields$dim2[x])))) {
          trg_index <- list(lat_idx)
        }else{
          trg_index <- list(1:trg_fields$dim1[x], lat_idx)
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
          trg_index <- list(lat_idx)
        }else{
          trg_index <- list(1:trg_fields$dim1[x], lat_idx)
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
          trg_index <- list(lat_idx)
        }else{
          trg_index <- list(1:trg_fields$dim1[x], lat_idx)
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

      # ATL08 output

      atl08_out <- cbind(data.frame(cycle = cycle,
                              region = region,
                              beam = n,
                              beam_strength = beam_strength_n),
                         land_df,
                         canopy_df,
                         terrain_df)

      colnames(atl08_out) <- stringr::str_replace(colnames(atl08_out),
                                                  pattern = "\\.",
                                                  replacement = "_")

      atl08_sf <- sf::st_as_sf(atl08_out,
                               coords = c("longitude", "latitude"),
                               remove = FALSE,
                               crs = sf::st_crs(4326))
      # Project to UTM Zone

      zone <- find_UTM_zone(atl08_sf$longitude[1], atl08_sf$latitude[1])

      hemi <- find_UTM_hemisphere(atl08_sf$latitude[1])

      if (hemi == "N") {
        if(stringr::str_length(zone) == 1) zone = stringr::str_c("0", zone, sep = "")
        epsg_trg <- as.integer(stringr::str_c("326", zone, sep = ""))
      }else{
        if(stringr::str_length(zone) == 1) zone = stringr::str_c("0", zone, sep = "")
        epsg_trg <- as.integer(stringr::str_c("327", zone, sep = ""))
      }


      atl08_sf_proj <- sf::st_transform(atl08_sf, epsg_trg)

      atl08_sf_proj <- atl08_sf_proj %>%
        dplyr::mutate(data.frame(sf::st_coordinates(atl08_sf_proj))) %>%
        dplyr::rename(Easting = X,
                      Northing = Y)

      # Find rotation matrix and rotation points

      rot_mat_list <- get_rot_matrix(xy_start = c(atl08_sf_proj$Easting[1], atl08_sf_proj$Northing[1]),
                                     xy_end = c(atl08_sf_proj$Easting[nrow(atl08_sf_proj)],
                                                atl08_sf_proj$Northing[nrow(atl08_sf_proj)]))

      # Calculate along and across track

      xy_rot <- get_along_distance(atl08_sf_proj$Easting,
                                   atl08_sf_proj$Northing,
                                   rot_mat_list$R_mat,
                                   rot_mat_list$xRotPt,
                                   rot_mat_list$yRotPt)

      atl08_sf_proj <- cbind(atl08_sf_proj, data.frame(xy_rot))


      # ATL03 classified photons
      if(atl03_pc) {

        beg_seg <- atl08_sf_proj$segment_id_beg
        end_seg <- atl08_sf_proj$segment_id_end

        unique_segments <- beg_seg[1]:end_seg[length(end_seg)]

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

        atl03_cl_out <- cbind(data.frame(cycle = cycle,
                                         region = region,
                                         beam = n,
                                         beam_strength = beam_strength_n),
                              signal_ph_df)

        colnames(atl03_cl_out) <- stringr::str_replace(colnames(atl03_cl_out),
                                                       pattern = "\\.",
                                                       replacement = "_")

      }

      # Clean column order

      atl08_sf_proj <- atl08_sf_proj %>%
        dplyr::mutate(UTM_zone = zone,
               UTM_hemi = hemi) %>%
        dplyr::relocate(rgt, .before = cycle) %>%
        dplyr::relocate(c(latitude, longitude,
                          UTM_zone, UTM_hemi,
                          Easting, Northing,
                          along_distance, across_distance), .after = beam_strength)

      if(!missing(odir)) {
        if(!dir.exists(odir)) {
          dir.create(odir)
        }

        readr::write_csv(sf::st_drop_geometry(atl08_sf_proj),
                         file = file.path(odir, stringr::str_c(tools::file_path_sans_ext(basename(atl08_h5)), "_", toupper(n), "_", toupper(beam_strength_n),".csv")))

        sf::st_write(atl08_sf_proj,
                     dsn = file.path(odir, stringr::str_c(tools::file_path_sans_ext(basename(atl08_h5)), "_", toupper(n), "_", toupper(beam_strength_n),".gpkg")),
                     quiet = TRUE,
                     delete_layer = TRUE)

        if(atl03_pc) {
          readr::write_csv(atl03_cl_out,
                           file = file.path(odir, stringr::str_c(tools::file_path_sans_ext(basename(atl08_h5)), "_", toupper(n), "_", toupper(beam_strength_n),"_ATL03_class.csv")))
        }
      }

      if(atl03_pc) {
        out[[n]] <- list(ATL08 = sf::st_drop_geometry(atl08_sf_proj),
                         ATL08_sf = atl08_sf_proj,
                         ATL03_cl = atl03_cl_out)
      }else{
        out[[n]] <- list(ATL08 = sf::st_drop_geometry(atl08_sf_proj),
                         ATL08_sf = atl08_sf_proj)
      }

    }else{
      warning(sprintf("No photons to return within specified lat_range for beam %s", n))
    }
  }
  return(out)
}
