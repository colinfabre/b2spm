library(terra)
library(sf)
library(httr)



geopf_download <- function(type, roi, output) {
    if (!dir.exists(output)) {
        dir.create(output, recursive = TRUE)
    }

    if (!inherits(roi, "SpatVector")) {
        roi <- terra::vect(roi)
    }

    wfs_url <- "https://data.geopf.fr/wfs/ows"
    if(type == "mnt") {
        typename <- "IGNF_MNT-LIDAR-HD:dalle"
    }
    if(type == "mnh") {
        typename <- "IGNF_MNH-LIDAR-HD:dalle"
    }
    if(type == "bdforetv2") {
        typename <- "LANDCOVER.FORESTINVENTORY.V2:formation_vegetale"
    }

    bbox <- terra::ext(roi)
    query <- paste0(
        wfs_url,
        "?SERVICE=WFS&VERSION=2.0.0&REQUEST=GetFeature&",
        "TYPENAMES=", typename, "&",
        "BBOX=", paste(bbox, collapse = ","), "&",
        "SRSNAME=EPSG:2154&",
        "OUTPUTFORMAT=application/json"
    )

    response <- httr::GET(query)
    httr::stop_for_status(response)
    tiles <- terra::vect(sf::st_read(httr::content(response, "text"), quiet = TRUE))
    tiles_intersect <- tiles[intersect(tiles, roi)]

    params <- c("wfs_url", "feature_type", "output", "download_tile")
    iterover <- seq_len(nrow(tiles_intersect))
    iterwith <- function(i) {
        download_tile(wfs_url, tiles_intersect[i, ], typename, output)
    }
    parallelism(params, iterover, iterwith)
}