library(terra)



dem_checker <- function(dem) {
    cat("Checking DEM...\n")
    if (!inherits(dem, c("SpatRaster", "RasterLayer"))) {
        stop("<!> ERROR - The object passed as `dem` must be an elevation raster (a georeferenced matrix of numerical values).\n")
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
    if (terra::xmin(terra::ext(bbox)) < terra::xmin(terra::ext(spruce_mask)) ||
        terra::xmax(terra::ext(bbox)) > terra::xmax(terra::ext(spruce_mask)) ||
        terra::ymin(terra::ext(bbox)) < terra::ymin(terra::ext(spruce_mask)) ||
        terra::ymax(terra::ext(bbox)) > terra::ymax(terra::ext(spruce_mask))) {
            stop("<!> ERROR - The input DEM is not within the validity range of the phenological model.")
    }

    gc()

    return(dem)
}



topo_spruce <- function(dem) {
    cat("Calculating spruce forests...\n")

    bbox <- terra::vect(terra::ext(dem))
    spruce_mask <- terra::rast(system.file("extdata/spruce_mask.tif", package = "b2spm"))
    terra::crs(spruce_mask) <- "EPSG:2154"
    
    spruce_mask_bbox <- terra::crop(spruce_mask, bbox)
    spruce_mask_bbox <- terra::resample(spruce_mask_bbox, dem, method = "near", threads = TRUE)
    terra::names(aspect) <- "spruce_forests"
    terra::varnames(aspect) <- "spruce_forests"

    cat("Spruce forests -- OK\n")
    gc()

    return(spruce_mask_bbox)
}



topo_alt <- function(dem) {
    cat("Calculating altitude...\n")

    alt <- terra::clamp(dem, lower = 0, upper = 4810)
    alt <- round(alt, 0)

    terra::names(alt) <- "alt"
    terra::varnames(alt) <- "alt"

    cat("Altitude -- OK\n")
    gc()

    return(alt)
}



topo_slope <- function(dem) {
    cat("Calculating slope...\n")

    slope <- terra::terrain(dem, v = "slope", unit = "radians")
    slope[is.na(slope)] <- 0
    slope <- terra::clamp(slope, lower = 0, upper = pi / 2)
    if (any(abs(terra::values(slope)) > 1.1)) {
        cat("<!> Warning - The computed topography contains pixels with a slope higher than 1.1rad (ie more than 200m of elevation gain in 25m of distance).\n")
        cat("Make sure that no DRIAS points fall on these pixels, or expect to get aberrant phenological results.\n")
    }

    terra::names(slope) <- "slope"
    terra::varnames(slope) <- "slope"

    cat("Slope -- OK\n")
    gc()

    return(slope)
}



topo_aspect <- function(dem) {
    cat("Calculating aspect...\n")

    aspect <- terra::terrain(dem, v = "aspect", unit = "radians")
    aspect[is.na(aspect)] <- 0.5 * pi # as if the slope was null
    aspect <- terra::clamp(aspect, lower = 0.5 * pi, upper = 2 * pi)

    terra::names(aspect) <- "aspect"
    terra::varnames(aspect) <- "aspect"

    cat("Aspect -- OK\n")
    gc()

    return(aspect)
}



topo_tpi <- function(dem) {
    cat("Calculating TPI...\n")

    tpi <- terra::terrain(dem, v = "TPI")
    tpi[is.na(tpi)] <- 0
    if (any(abs(terra::values(tpi)) > 200)) {
        cat("<!> Warning - The computed topography contains pixels with shoulders higher than 200m (ie ridges above 200m of the surrounding topography in a 25m-radius).\n")
        cat("Make sure that no DRIAS points fall on these pixels, or expect to get aberrant phenological results.\n")
    }

    terra::names(tpi) <- "tpi"
    terra::varnames(tpi) <- "tpi"

    cat("TPI -- OK\n")
    gc()

    return(tpi)
}