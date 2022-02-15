#' Combine ATL08 shp files
#'
#' @param directory
#' @param out_dir
#' @param out_type
#'
#' @export


combine_atl08_shp <- function(directory,
                              out_dir = directory,
                              out_type = "gpkg") {


  all_files <- list.files(directory,
                          pattern = ".shp$",
                          full.names = TRUE)

  un_beams <- unique(stringr::str_match(basename(all_files), "gt[0-9]{1}[l|r]{1}")[, 1])

  for (beam in un_beams) {
    flist <- all_files[stringr::str_detect(all_files, beam)]

    # Open and rbind all layers
    shp_list <- lapply(1:length(flist), function(x) {


      subdir <- ifelse(stringr::str_detect(basename(flist[x]), "land_segments_canopy"), "canopy",
                       ifelse(stringr::str_detect(basename(flist[x]), "land_segments_terrain"), "terrain", ""))

      if (subdir == "canopy") {
        param <- stringr::str_match(basename(flist[x]),
                                    pattern = "^gt[0-9]{1}[r|l]{1}_land_segments_canopy_(.*).shp")[1, 2]
      }else if (subdir == "terrain"){
        param <- stringr::str_match(basename(flist[x]),
                                    pattern = "^gt[0-9]{1}[r|l]{1}_land_segments_terrain_(.*).shp")[1, 2]
      }else if (subdir == ""){
        param <- stringr::str_match(basename(flist[x]),
                                    pattern = "^gt[0-9]{1}[r|l]{1}_land_segments_(.*).shp")[1, 2]
      }

      df <- sf::st_drop_geometry(sf::st_read(flist[x], quiet = TRUE))

      if(ncol(df) == 3) {

        if(subdir != "") {
          colnames(df)[1] <- paste(subdir, param, sep = "_")
        }else{
          colnames(df)[1] <- param
        }


      } else if(any(stringr::str_detect(colnames(df), "X[0-9]+"))) {

        if(param %in% c("canopy_h_metrics", "canopy_h_metrics_abs")) {

          colnames(df)[which(stringr::str_detect(colnames(df), "X[0-9]+"))] <- paste0(subdir, "_", param, "_zq", seq(10, 95, by = 5))

        }else{
          colnames(df)[which(stringr::str_detect(colnames(df), "X[0-9]+"))] <- paste0(stringr::str_match(basename(flist[x]),
                                                                                                "segments_(.*).shp")[1, 2], "_", colnames(df)[which(stringr::str_detect(colnames(df), "X[0-9]+"))])
        }
      }

      return(df)

    })

    multi_join <- Reduce(
      function(x, y, ...) merge(x, y, by = c("latitude", "longitude")),
      shp_list
    )

    multi_join$beam <- beam

    multi_join_sf <- sf::st_as_sf(multi_join,
                              coords = c("longitude", "latitude"),
                              crs = 4326,
                              remove = FALSE)


    if (out_type == "rds") {
      saveRDS(multi_join_sf,
              file = file.path(out_dir,
                               paste0(str_match(all_files[1], "processed_ATL08_[0-9]{14}_[0-9]{8}_[0-9]{3}_[0-9]{2}"), "_", beam, ".rds")))
    }else if (out_type == "gpkg"){
      sf::st_write(multi_join_sf,
               dsn = file.path(out_dir,
                               paste0(str_match(all_files[1], "processed_ATL08_[0-9]{14}_[0-9]{8}_[0-9]{3}_[0-9]{2}"), "_", beam, ".gpkg")))
    }else{
      sf::st_write(multi_join_sf,
               dsn = file.path(out_dir,
                               paste0(str_match(all_files[1], "processed_ATL08_[0-9]{14}_[0-9]{8}_[0-9]{3}_[0-9]{2}"), "_", beam, ".shp")))
    }

  }
}





