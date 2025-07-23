#' drias_fetcher
#'
#' This function fetches the DRIAS database corresponding to the provided DEM and year, and formats it into a dataframe compatible with the B2SPM pipeline.
#'
#' @param topography The raster stack returned by the topo_comp() function.
#' @param year Year of analysis. Currently stored databases are 2050, 2075 and 2100, for the RCP8.5 scenario (cf. README).
#' @return A data.frame containing the columns: `id`, `X93`, `Y93`, `date`, `doy`, `tmin`, `tmax`, `tmean`, `tot_pr`, `spec_hum`, `vis_solrad`, `ir_solrad`, `wind`, `pet`.
#' @examples
#' \dontrun{
#'  drias_table <- drias_fetcher(topography, 2050)
#' }
#' @export
drias_fetcher <- function(topography, year) {
    cat("Fetching the corresponding online DRIAS database...\n")
    if (!year %in% c(2050, 2075, 2100)) {
        stop("!! ERROR - 'year' must be equal either to 2050, 2075 or 2100.")
    }

    cat("Fetching the online database corresponding to the provided year of analysis can take up to 60s.\n")
    cat("Please wait...\n")
    drias_points <- tryCatch({
        terra::vect(paste0("https://github.com/colinfabre/b2spm_database/raw/refs/heads/main/alpine_arc_", year, ".gpkg"))},
        error = function(e) {stop("!! ERROR -  Impossible to fetch the online database corresponding to the provided year. Please check your internet connection.\n")})
    cat("The online database has been correctly fetched.\n")
    terra::crs(drias_points) <- "EPSG:27572"
    drias_points <- terra::project(drias_points, "EPSG:2154")

    roi <- terra::vect(terra::ext(topography))
    terra::crs(roi) <- "EPSG:2154"

    drias_points_roi <- terra::crop(drias_points, roi)
    drias_table <- as.data.frame(drias_points_roi)
    drias_table <- cbind(drias_table[, 1], terra::geom(drias_points_roi)[, c("x", "y")], drias_table[, 2:11])
    names(drias_table) <- c("id", "X93", "Y93", "date", "tmin", "tmax", "tmean", "tot_pr", "spec_hum", "vis_solrad", "ir_solrad", "wind", "pet")

    stations <- unique(drias_table$id)
    years <- unique(as.numeric(substr(drias_table$date, 7, 10)))
    results <- list()

    for (station in stations) {
        station_data <- drias_table[drias_table$id == station, ]

        for (year in years) {
            
            year_mask <- tryCatch({
                as.numeric(substr(station_data$date, 7, 10)) == year},
                error = function(e) {stop("!! ERROR -  Impossible to extract and convert records dates for station ", station, ".\n")})
            doy_val <- seq_len(sum(year_mask))
            station_data$doy[year_mask] <- doy_val
        }

        results[[as.character(station)]] <- station_data
    }

    drias_table <- do.call(rbind, results)
    drias_table <- cbind(drias_table[, 1:4], drias_table[, "doy"], drias_table[, 5:13])
    names(drias_table) <- c("id", "X93", "Y93", "date", "doy", "tmin", "tmax", "tmean", "tot_pr", "spec_hum", "vis_solrad", "ir_solrad", "wind", "pet")

    drias_table$tmin <- drias_table$tmin - 273.15
    drias_table$tmax <- drias_table$tmax - 273.15
    drias_table$tmean <- drias_table$tmean - 273.15
    drias_table$tot_pr <- drias_table$tot_pr * 86400
    drias_table$pet <- pmax(drias_table$pet * 86400, 0)

    cat("Fetching DRIAS table -- OK\n")
    gc()
    
    return(drias_table)
}