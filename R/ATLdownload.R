#' Downloads ATLAS data products with a curl request
#'
#' curl needs to be installed on your machine and a token available to download data
#'
#'@param short_name Short name ID of the dataset. For example ATL03 or ATL08
#'@param bounding_box Bounding box used as a search filter. lower_left_long, lower_left_lat, upper_right_long, upper_righ_lat
#'@param bbox Bounding box used for spatial subsetting. Must match the coordinates entered for bounding_box. If left NULL, automatically copies bounding_box values.
#'@param shapefile An ESRI Shapefile used for spatial subsetting. Output files will be cropped to the shapefile extent.
#'@param time Temporal range filter
#'@param token Token required to access data using Eathdata login credentials
#'@param version Dataset version number. Default is 002
#'@param coverage Subsets dataset parameters (based on HDF5 dataset structure)
#'@param projection Projection of the output files.
#'@param tempDir Full path to temporary folder where download request and eventually clipping shapefile will be saved
#'@param dwlDir Full path to folder where data will be downloaded
#'
#'@export

ATLdownload <- function(short_name,
                        bounding_box,
                        bbox,
                        shapefile,
                        time,
                        token,
                        version,
                        coverage,
                        projection,
                        tempDir,
                        dwlDir,
                        email = "yes") {


  #Checking inputs
  if(missing(shapefile)) {
    doSHP <- FALSE
  }else{
    doSHP <- TRUE
  }

  if(missing(tempDir) | is.null(tempDir) | tempDir == ""){
    tempDir <- getwd()
  }else if (!is.character(tempDir)){
    stop("tempDir must be character")
  }else{
    if(!file.exists(tempDir)) dir.create(tempDir)
    tempDir <- normalizePath(tempDir)
  }

  if(missing(dwlDir) | is.null(dwlDir) | dwlDir == ""){
    dwlDir <- getwd()
  }else if (!is.character(dwlDir)){
    stop("dwlDir must be character")
  }else{
    if(!file.exists(dwlDir)) dir.create(dwlDir)
    dwlDir <- normalizePath(dwlDir)
  }


  if(doSHP) {
    if(is.character(shapefile)) {
      if(!file.exists(shapefile)) stop("File name of shapefile does not exist")
      shapefile <- rgdal::readOGR(shapefile, verbose = FALSE, stringsAsFactors = F)
    }
  }


  if(!is.character(token))
    stop("token must be a character")

  # Create a unique ID for temp files
  tempID <- basename(tempfile("", tmpdir = "", fileext = ""))

  #Key-Value-Pairs

  kvp_def <- list(short_name = short_name,
                  bounding_box = if(missing(bounding_box)) NULL else bounding_box,
                  bbox = if(missing(bbox)) NULL else bbox,
                  time = if(missing(time)) NULL else time,
                  version = if(missing(version)) NULL else version,
                  coverage = if(missing(coverage)) NULL else coverage,
                  projection = if(missing(projection)) NULL else projection,
                  email = email,
                  token = token
                  )

  if(doSHP) {
    #If sf object transformn to sp
    if (is(shapefile, "sf")){
      if(sf::st_geometry_type(shapefile) == "POLYGON") {
        shapefile <- sf::as_Spatial(shapefile)
      }else{
        stop("The sf shapefile object must have a POLYGON geometry type")
      }
    }
    #Check if SpatialPolygonsDataFrame
    if(!is(shapefile, "SpatialPolygonsDataFrame")) stop("shapefile must be a polygon")
    if(length(shapefile) == 0 | length(shapefile) == 0) stop("shapefile must contain only one feature")
    if(is.null(raster::crs(shapefile))){
      stop("shapefile must have a defined CRS")
    } else {
      shapefile_crs <- sp::proj4string(shapefile)
    }
    if (!raster::compareCRS(shapefile_crs, sp::CRS("+init=epsg:4326"))) {
      message("Projecting shapefile to WGS84 (EPSG:4326)")
      shapefile <- sp::spTransform(shapefile, sp::CRS("+init=epsg:4326"))
    }

    # #Get vertices coordinates
    # poly <- shapefile@polygons[[1]]@Polygons[[1]]@coords
    # #Concatenate in one string
    # poly <- paste0(apply(format(poly), 1, paste, collapse = ","), collapse = ",")
    # poly <- gsub(" ", "", poly)
    # Get bbox in counter clockwise order
    box <- sp::bbox(shapefile)
    poly <- c(box[1,1], box[2,1], box[1,2], box[2,1], box[1,2], box[2,2],
              box[1,1], box[2,2],box[1,1], box[2,1])
    poly <- format(round(poly,3), nsmall = 3)
    poly <- paste0(as.character(poly), collapse = ",")
    poly <- gsub(" ", "", poly)

    rgdal::writeOGR(shapefile, layer = paste0("clipping_",tempID), driver = "ESRI Shapefile", dsn = tempDir)
    utils::zip(zipfile = file.path(tempDir, paste0("clipping_",tempID)), files = list.files(file.path(tempDir), pattern = tempID, full.names = T), flags = "a -tzip" ,zip = "C:/Program Files/7-Zip/7Z")
    shp_name <- file.path(tempDir, paste0("clipping_",tempID,".zip"))

    kvp_def$bbox <- NULL
    kvp_def$bounding_box <- NULL
    kvp_def$polygon <- poly

    kvp_call_shp <- paste0("shapefile=@",shp_name)
  }

  kvp_call <- kvp_def[!unlist(lapply(kvp_def, is.null))]

  kvp_call <- paste(names(kvp_call), kvp_call, sep = "=", collapse = "&")

  if (doSHP) {
    dwl_call <- paste0("curl -O -J -F ","\"", kvp_call_shp,"\"", " \"","https://n5eil02u.ecs.nsidc.org/egi/request?", kvp_call,"\"")
  }else{
    dwl_call <- paste0("curl -O -J ","\"","https://n5eil02u.ecs.nsidc.org/egi/request?", kvp_call, "\"")
  }

  #Save the command in the bat file
  dwlDrive <- strsplit(dwlDir,"\\\\")[[1]][1]

  tmpfile <- file.path(tempDir,paste0("request_", tempID,".bat"))

  fileConn <- file(tmpfile)
  writeLines(c(dwlDrive, paste0("cd ", dwlDir), dwl_call, "pause"), fileConn, sep = "\n")
  close(fileConn)

  #Execute .bat file
  shell.exec(tmpfile)

}
