#' Read ATL08 h5 file as a data frame
#'
#' @param file Path to h5 file
#' @param beam Character vector indicating beams to process
#' @param beam_strength Character vector indicating the strength of beams to process
#' @param lat_range Numeric vector. Lower and upper latitude to return
#' @param odir Character. Output directory
#'
#' @export

read_ATL08 <- function(file,
                       beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
                       beam_strength = c("weak", "strong"),
                       lat_range,
                       odir) {

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
  rgt <- as.integer(rhdf5::h5read(file = file, name = "orbit_info/rgt"))
  cycle <- as.integer(rhdf5::h5read(file = file, name = "orbit_info/cycle_number"))
  region <- as.integer(rhdf5::h5read(file = file, name = "ancillary_data/start_region"))

  out <- list() #Initiate list to return

  for ( n in beam) {

    beam_strength_n <- names(list_strength)[unlist(lapply(list_strength, function(x) n %in% x))]

    # Find all 100-m segment center points within lat_range
    lat_center <- as.numeric(rhdf5::h5read(file = file, name = paste0(n,"/land_segments/latitude")))

    if (any(is.na(lat_range))) {
      lat_range <- c(min(lat_center), max(lat_center))
    }

    # Index of points to keep
    #lat_idx <- which(lat_center >= lat_range[1] & lat_center <= lat_range[2])
    lat_idx <- lat_center >= lat_range[1] & lat_center <= lat_range[2]

    if (sum(lat_idx) != 0) {
      beg_seg <- as.integer(rhdf5::h5read(file = file,
                                          name = paste0(n,"/land_segments/segment_id_beg")))
      end_seg <- as.integer(rhdf5::h5read(file = file,
                                          name = paste0(n,"/land_segments/segment_id_end")))

      all_seg <- beg_seg[which(lat_idx)[1]]:end_seg[which(lat_idx)[length(which(lat_idx))]]

      # Land segment flags and info
      trg_fields <-  h5list %>%
        dplyr::filter(stringr::str_detect(group, stringr::str_c(n, "/land_segments$")),
                      otype == "H5I_DATASET") %>%
                        dplyr::mutate(full_name = stringr::str_c(group, name, sep = "/"))

      trg_fields_list <- lapply(1:length(trg_fields$name), function(x) {
        out <- rhdf5::h5read(file = file,
                      name = trg_fields$full_name[x])

        if(trg_fields$dclass[x] == "FLOAT") {
          if(length(dim(out)) > 1) {
            out <- as.data.frame(apply(out, 1, as.numeric))
            colnames(out) <- str_c("X", 1:ncol(out))

          }else{
            out <- as.numeric(out)
          }
        }else if (trg_fields$dclass[x] == "INTEGER") {
          if(length(dim(out)) > 1) {
            out <- as.data.frame(apply(out, 1, as.integer))
            colnames(out) <- str_c("X", 1:ncol(out))
          }else{
            out <- as.integer(out)
          }
        }

        return(out)

      })

      names(trg_fields_list) <- trg_fields$name

      land_df <- as.data.frame(trg_fields_list)

      # Canopy products

      trg_fields <-  h5list %>%
        dplyr::filter(stringr::str_detect(group, stringr::str_c(n, "/land_segments/canopy$")),
                      otype == "H5I_DATASET") %>%
        dplyr::mutate(full_name = stringr::str_c(group, name, sep = "/"))

      trg_fields_list <- lapply(1:length(trg_fields$name), function(x) {
        out <- rhdf5::h5read(file = file,
                             name = trg_fields$full_name[x])



        if(trg_fields$dclass[x] == "FLOAT") {
          if(length(dim(out)) > 1) {
            out <- as.data.frame(apply(out, 1, as.numeric))
            colnames(out) <- str_c("X", 1:ncol(out))

          }else{
            out <- as.numeric(out)
          }
        }else if (trg_fields$dclass[x] == "INTEGER") {
          if(length(dim(out)) > 1) {
            out <- as.data.frame(apply(out, 1, as.integer))
            colnames(out) <- str_c("X", 1:ncol(out))
          }else{
            out <- as.integer(out)
          }
        }

        return(out)

      })

      names(trg_fields_list) <- trg_fields$name

      canopy_df <- as.data.frame(trg_fields_list)

      # Terrain products

      trg_fields <-  h5list %>%
        dplyr::filter(stringr::str_detect(group, stringr::str_c(n, "/land_segments/terrain$")),
                      otype == "H5I_DATASET") %>%
        dplyr::mutate(full_name = stringr::str_c(group, name, sep = "/"))

      trg_fields_list <- lapply(1:length(trg_fields$name), function(x) {
        out <- rhdf5::h5read(file = file,
                             name = trg_fields$full_name[x])



        if(trg_fields$dclass[x] == "FLOAT") {
          if(length(dim(out)) > 1) {
            out <- as.data.frame(apply(out, 1, as.numeric))
            colnames(out) <- str_c("X", 1:ncol(out))

          }else{
            out <- as.numeric(out)
          }
        }else if (trg_fields$dclass[x] == "INTEGER") {
          if(length(dim(out)) > 1) {
            out <- as.data.frame(apply(out, 1, as.integer))
            colnames(out) <- str_c("X", 1:ncol(out))
          }else{
            out <- as.integer(out)
          }
        }

        return(out)

      })

      names(trg_fields_list) <- trg_fields$name

      terrain_df <- as.data.frame(trg_fields_list)

      # ATL08 output

      atl08_out <- cbind(data.frame(rgt = rgt,
                              cycle = cycle,
                              region = region,
                              beam = n,
                              beam_strength = beam_strength_n),
                         land_df,
                         canopy_df,
                         terrain_df)

      colnames(atl08_out) <- stringr::str_replace(colnames(atl08_out),
                                                  pattern = "\\.",
                                                  replacement = "_")

      atl08_out <- atl08_out[lat_idx, ]

      atl08_out_sf <- sf::st_as_sf(atl08_out,
                               coords = c("longitude", "latitude"),
                               remove = FALSE,
                               crs = sf::st_crs(4326))

      # ATL03 classified photons

      trg_fields <-  h5list %>%
        dplyr::filter(stringr::str_detect(group, stringr::str_c(n, "/signal_photons")),
                      otype == "H5I_DATASET") %>%
        dplyr::mutate(full_name = stringr::str_c(group, name, sep = "/"))

      trg_fields_list <- lapply(1:length(trg_fields$name), function(x) {
        out <- rhdf5::h5read(file = file,
                             name = trg_fields$full_name[x])



        if(trg_fields$dclass[x] == "FLOAT") {
          if(length(dim(out)) > 1) {
            out <- as.data.frame(apply(out, 1, as.numeric))
            colnames(out) <- str_c("X", 1:ncol(out))

          }else{
            out <- as.numeric(out)
          }
        }else if (trg_fields$dclass[x] == "INTEGER") {
          if(length(dim(out)) > 1) {
            out <- as.data.frame(apply(out, 1, as.integer))
            colnames(out) <- str_c("X", 1:ncol(out))
          }else{
            out <- as.integer(out)
          }
        }

        return(out)

      })

      names(trg_fields_list) <- trg_fields$name

      signal_ph_df <- as.data.frame(trg_fields_list)

      signal_ph_df <- signal_ph_df %>%
        mutate(classed_pc_flag = factor(classed_pc_flag,
                                        levels = c(0, 1, 2, 3),
                                        labels = c("noise", "ground", "canopy", "top of canopy")))

      # Filter from lat range
      signal_ph_df <- signal_ph_df %>%
        dplyr::filter(ph_segment_id %in% dplyr::all_of(all_seg))


      atl03_cl_out <- cbind(data.frame(cycle = cycle,
                                    region = region,
                                    beam = n,
                                    beam_strength = beam_strength_n),
                         signal_ph_df)

      colnames(atl03_cl_out) <- stringr::str_replace(colnames(atl03_cl_out),
                                                  pattern = "\\.",
                                                  replacement = "_")

      if(!missing(odir)) {
        if(!dir.exists(odir)) {
          dir.create(odir)
        }

        readr::write_csv(atl08_out,
                         file = file.path(odir_atl08, str_c(tools::file_path_sans_ext(basename(file)), "_", toupper(n), "_", toupper(beam_strength_n),".csv")))

        sf::st_write(atl08_out_sf,
                     dsn = file.path(odir_atl08, str_c(tools::file_path_sans_ext(basename(file)), "_", toupper(n), "_", toupper(beam_strength_n),".gpkg")),
                     quiet = TRUE,
                     delete_layer = TRUE)

        readr::write_csv(atl03_cl_out,
                         file = file.path(odir, str_c(tools::file_path_sans_ext(basename(file)), "_", toupper(n), "_", toupper(beam_strength_n),"_ATL03_class.csv")))

      }

      out[[n]] <- list(ATL08 = atl08_out,
                       ATL08_sf = atl08_out_sf,
                       ATL03_cl = atl03_cl_out)
    }else{
      warning(sprintf("No photons to return within specified lat_range for beam %s", n))
    }
  }
  return(out)
}
