#' Generate LINESTRING from ATL08 center points
#'
#' @param atl08_df Data.frame or sf object containing the center points coordinates of ATL08 100-m segments. Can be generated from \link[icesat2R]{read_ATL08}
#' @param out_dir Directory where output will be written (as gpkg file)
#' @param tr_length Length of transect (DO NOT CHANGE)
#'
#' @export

make_transect <- function(atl08_sf,
                          out_dir,
                          tr_length = 100) {

  if (is.character(atl08_sf)) {
    if (tools::file_ext(atl08_sf) == "h5") {

      atl08_dat <- read_ATL08(atl08_sf)
      atl08_sf <- atl08_dat$ATL08_sf

    }else if (tools::file_ext(atl08_sf) == "csv") {

      atl08_sf <- readr::read_csv(atl08_sf)

      atl08_sf <- sf::st_as_sf(atl08_sf,
                               coords = c("longitude", "latitude"),
                               crs = 4326)


    }else if (tools::file_ext(atl08_sf) == "gpkg"){

      atl08_sf <- sf::st_read(atl08_sf, quiet = TRUE)
    }

  }else{
    if (!inherits(atl08_sf, "sf")) {
      atl08_sf <- sf::st_as_sf(atl08_sf,
                               coords = c("longitude", "latitude"),
                               crs = 4326, )
    }
  }


  lat_lon <- sf::st_coordinates(atl08_sf)

  utm_zone <- find_UTM_zone(longitude = lat_lon[1, 1],
                latitude = lat_lon[1, 2])
  utm_hem <- find_UTM_hemisphere(latitude = lat_lon[1, 2])

  trg_crs <- rgdal::make_EPSG() %>%
    dplyr::filter(note == paste0("NAD83(CSRS) / UTM zone ", utm_zone, utm_hem))

  atl08_sf_proj <- sf::st_transform(atl08_sf, trg_crs$code)

  xy <- sf::st_coordinates(atl08_sf_proj)

  # Calculate sin and cos of theta.
  dx <- diff(xy[, "X"])
  dy <- diff(xy[, "Y"])

  cos_theta = dx/sqrt(dx^2+dy^2)
  sin_theta = dy/sqrt(dx^2+dy^2)

  # Calculate start and end point coordinates.
  #print("Calculate the start- and end- point coordinates...")
  x_start = xy[, "X"][-nrow(xy)]-tr_length/2*cos_theta
  y_start = xy[, "Y"][-nrow(xy)]-tr_length/2*sin_theta
  x_end = xy[, "X"][-nrow(xy)]+tr_length/2*cos_theta
  y_end = xy[, "Y"][-nrow(xy)]+tr_length/2*sin_theta

  tr_sfg_list <- list()

  for(i in 1:length(x_start)) {
    tr_sfg_list <- c(tr_sfg_list,
                     list(sf::st_linestring(x = matrix(c(x_start[i], y_start[i], x_end[i], y_end[i]),
                                                   byrow = TRUE,
                                                   ncol = 2),
                                   dim = "XY")))


  }


  tr_sfc <- sf::st_sfc(tr_sfg_list)

  tr_sf <- sf::st_set_geometry(atl08_sf_proj[1:(nrow(atl08_sf_proj)-1), ], tr_sfc)
  sf::st_crs(tr_sf) <- trg_crs$code

  if (!missing(out_dir)) {

    if(!dir.exists(out_dir)) dir.create(out_dir)

    sf::st_write(tr_sf,
                 dsn = file.path(out_dir,
                                 paste0(tools::file_path_sans_ext(basename(in_shp)), ".gpkg")),
                 driver = "GPKG",
                 delete_layer = TRUE)
  }


  return(tr_sf)
}


test <- unlist(mapply(function(x, y) x:y, x = tr_sf$segment_id_beg, y = tr_sf$segment_id_end, SIMPLIFY = FALSE))

test_split <- sf::st_line_sample(tr_sf,
                                 sample = c(0, 0.25, 0.5, 0.75, 1))

test_split <- st_cast(test_split, "POINT")

df <- data.frame(segment_id = test)
df$geometry <- test_split
st_geometry(df) <- test_split


seg_sfg_list <- list()

for(i in 1:length(x_start)) {
  seg_sfg_list <- c(seg_sfg_list,
                   list(sf::st_linestring(x = matrix(c(x_start[i], y_start[i], x_end[i], y_end[i]),
                                                     byrow = TRUE,
                                                     ncol = 2),
                                          dim = "XY")))


}

