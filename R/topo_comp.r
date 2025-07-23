#' topo_comp
#'
#' This function loads the spruce forest mask and computes topographic data for the study area using the provided DEM.
#'
#' @param dem A raster (preferably a SpatRaster in EPSG:2154) representing elevation, at a 25m-pixel spatial resolution, covering the entire study area plus a 10-pixel buffer (ie 250m). It has to be located within the French alpine arc (see README.md for further information).
#' @return A raster stack containing the spruce forest mask and each topographic raster in a separate band.
#' @examples
#' \dontrun{
#'  dem <- terra::rast("dem_roi.tif")
#'  topography <- topo_comp(dem)
#' }
#' @export
topo_comp <- function(dem) {
    cat("Calculating the topography raster...\n")
    if (!inherits(dem, c("SpatRaster", "RasterLayer"))) {
        stop("!! ERROR - The object passed as `dem` must be an elevation raster (a georeferenced matrix of numerical values).")
    }
    if (class(dem)[1] != "SpatRaster") {
        try(dem <- terra::rast(dem))
        cat("<!> Warning - The input DEM has been converted to a SpatRaster (terra format).\n")
    }
    if (terra::crs(dem) != "EPSG:2154") {
        try(dem <- terra::project(dem, "EPSG:2154"))
        cat("<!> Warning - The input DEM has been reprojected to Lambert 93 (EPSG:2154).\n")
    }
    if (terra::res(dem)[1] != 25) {
        try(grid <- terra::rast(ext = terra::ext(dem), resolution = 25, crs = terra::crs(dem)))
        try(dem <- terra::resample(dem, grid, method = "near", threads = TRUE))
        cat("<!> Warning - The input DEM has been resampled to a 25m-pixel spatial resolution.\n")
    }
    if (terra::nlyr(dem) > 1) {
        try(dem <- dem[[1]])
        cat("<!> Warning - The input DEM contains more than one elevation layer. Only the first one is considered.\n")
    }
    if (any(is.na(terra::values(dem))) == TRUE) {
        try(dem[is.na(dem)] <- 0)
        cat("<!> Warning - The input DEM contains pixels with NA values. They have been assigned a value of 0.\n")
    }

    bbox <- terra::vect(terra::ext(dem))

    spruce_mask <- terra::rast(system.file("extdata/spruce_mask.tif", package = "b2spm"))
    terra::crs(spruce_mask) <- "EPSG:2154"
    if (terra::xmin(terra::ext(bbox)) >= terra::xmin(terra::ext(spruce_mask)) &&
        terra::xmax(terra::ext(bbox)) <= terra::xmax(terra::ext(spruce_mask)) &&
        terra::ymin(terra::ext(bbox)) >= terra::ymin(terra::ext(spruce_mask)) &&
        terra::ymax(terra::ext(bbox)) <= terra::ymax(terra::ext(spruce_mask))) {
            spruce_mask_bbox <- terra::crop(spruce_mask, bbox)
            spruce_mask_bbox <- terra::resample(spruce_mask_bbox, dem, method = "near", threads = TRUE)

            dem <- terra::clamp(dem, lower = 0, upper = 4810)
            dem <- round(dem, 0)

            slope <- terra::terrain(dem, v = "slope", unit = "radians")
            slope[is.na(slope)] <- 0
            slope <- terra::clamp(slope, lower = 0, upper = pi / 2)
            if (any(abs(terra::values(slope)) > 1.1)) {
                cat("<!> Warning - The computed topography contains pixels with a slope higher than 1.1rad (ie more than 200m of elevation gain in 25m of distance).\n")
                cat("Make sure that no DRIAS points fall on these pixels, or expect to get aberrant phenological results.\n")
            }

            aspect <- terra::terrain(dem, v = "aspect", unit = "radians")
            aspect[is.na(aspect)] <- 0.5 * pi # as if the slope was null
            aspect <- terra::clamp(aspect, lower = 0.5 * pi, upper = 2 * pi)

            tpi <- terra::terrain(dem, v = "TPI")
            tpi[is.na(tpi)] <- 0
            if (any(abs(terra::values(tpi)) > 200)) {
              cat("<!> Warning - The computed topography contains pixels with shoulders higher than 200m (ie ridges above 200m of the surrounding topography in a 25m-radius).\n")
              cat("Make sure that no DRIAS points fall on these pixels, or expect to get aberrant phenological results.\n")
            }

            topography <- c(spruce_mask_bbox, dem, slope, aspect, tpi)
            names(topography) <- c("spruce_forests", "alt", "slope", "aspect", "tpi")
            terra::varnames(topography) <- c("spruce_forests", "alt", "slope", "aspect", "tpi")
        } else {
            stop("!! ERROR - The input DEM is not within the validity range of the phenological model.")
        }

    cat("Topography raster -- OK\n")
    gc()

    return(topography)
}