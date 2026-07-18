library(terra)
library(sf)
library(stats)
library(utils)
library(tools)
library(httr)
library(parallel)



parallelism <- function(params, iterover, iterwith) {
    cl <- parallel::makeCluster(detectCores() - 2)
    parallel::clusterExport(cl, params, envir = environment())
    parallel::parLapply(cl, iterover, iterwith)
    parallel::stopCluster(cl)
}



download_tile <- function(wfs_url, tile, typename, output) {
    tile_id <- names(tile)[1]
    tile_url <- paste0(
        wfs_url,
        "?SERVICE=WFS&VERSION=2.0.0&REQUEST=GetFeature&",
        "TYPENAMES=", typename, "&",
        "FEATUREID=", tile_id, "&",
        "OUTPUTFORMAT=application/gml+json"
    )

    temp_file <- tempfile(fileext = ".json")
    on.exit(file.remove(temp_file))
    httr::GET(tile_url, httr::write_disk(temp_file, overwrite = TRUE))

    tile_data <- sf::st_read(temp_file, quiet = TRUE)
    tile_raster <- terra::rast(terra::vect(tile_data, geom = sf::st_geometry(tile_data)))
    output_file <- file.path(output, paste0("MNT_", tile_id, ".tif"))
    terra::writeRaster(tile_raster, output_file, filetype = "GTiff", overwrite = TRUE)

    file.remove(temp_file)
}



process_nc <- function(url, roi, output) {
    filename <- tools::file_path_sans_ext(url)
    output_path <- file.path(output, paste0(filename, ".tif"))

    temp_file <- tempfile(fileext = ".nc")
    download.file(url, destfile = temp_file, mode = "wb", quiet = TRUE)

    nc_data <- ncdf4::nc_open(temp_file)
    lon <- ncdf4::ncvar_get(nc_data, "lon")
    lat <- ncdf4::ncvar_get(nc_data, "lat")

    grid <- terra::rast(nrows = length(lat), ncols = length(lon), ext = c(min(lon), max(lon), min(lat), max(lat)))
    terra::crs(raster) <- "EPSG:2154"
    raster <- terra::crop(grid, roi, snap = "out")

    interest_vars_raw <- c("prtotAdjust", "tasAdjust", "tasmaxAdjust", "tasminAdjust", "hussAdjust", "sfcWindAdjust", "rldsAdjust", "rsdsAdjust", "evspsblpotAdjust")
    interest_vars_mod <- c("prec", "t_mean", "t_max", "t_min", "spec_hum", "wind", "ir_rad", "vis_rad", "etp")
    var_mapping <- stats::setNames(interest_vars_mod, interest_vars_raw)
    var_names <- names(nc_data$var)
    for (var in var_names) {
        if (var %in% interest_vars_raw) {
        var_data <- ncdf4::ncvar_get(nc_data, var)
        var_data_2d <- aperm(var_data[1, , ], c(2, 1))
        var_layer <- terra::rast(nrows = nrow(var_data_2d), ncols = ncol(var_data_2d), ext = ext(raster), crs = crs(raster))
        var_layer <- setValues(var_layer, as.matrix(var_data_2d))
        raster <- terra::addLayer(raster, var_layer)
        names(raster)[nlyr(raster)] <- var_mapping[[var]]
        }
    }

    terra::writeRaster(raster, output_path, overwrite = TRUE)

    ncdf4::nc_close(nc_data)
    file.remove(temp_file)
}