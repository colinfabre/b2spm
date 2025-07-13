#' hsc
#'
#' This function calculates the hydraulic stress index of spruce forests according to the water balance P-ETP and the AWC_max (Available Water Content) map (PIEDALLU, C. et al. Cartographie de la Réserve Utile Maximale en Eau des sols forestiers de France, UMR SILVAE, 2012).
#'
#' @param drias_table The DRIAS table processed by the ppc() function.
#' @return A data.frame with the columns: `id`, `X93`, `Y93`, `hsi` (maximal hydraulic stress index throughout the year).
#' @examples
#' \dontrun{
#'  drias_table <- hsc(drias_table)
#' }
#' @export
hsc <- function(drias_table) {
    cat("Computing hydraulic stress of spruce stands regarding soil available water content...\n")
    awc_max <- terra::rast(system.file("extdata/awc_max.tif", package = "b2spm"))
    terra::crs(awc_max) <- "EPSG:2154"

    stations <- unique(drias_table$id)

    results <- lapply(stations, function(station) {
        station_data <- drias_table[drias_table$id == station, ]
        station_point <- terra::vect(station_data, geom = c("X93", "Y93"), crs = "EPSG:2154")[1]
        station_point <- station_point[1]

        station_awc_max <- terra::extract(awc_max, station_point)[, 2]
        station_awc <- 0
        station_stress <- 0
        station_drought <- 0
        station_data$awc <- NA
        station_data$hsi <- NA

        for (doy in station_data$doy) {
            station_wb <- station_data$tot_pr[station_data$doy == doy] - station_data$pet[station_data$doy == doy]
            if (station_wb >= 0) {
                station_awc <- min(station_awc + station_wb, station_awc_max)
            } else {
                station_awc <- max(station_awc + station_wb, 0)
            }

            station_data[station_data$doy == doy, "awc"] <- station_awc

            if (station_wb < 0 && station_awc == 0) {
                station_stress <- station_stress + 1
                station_drought <- station_drought + 1
            } else {
                station_drought <- 0
            }

            stress_curve <- function(x, ths, max) {
                if (x <= ths) {
                    return(1 + 0.5 * (x / ths))
                } else {
                    k <- -log(0.01) / (max - ths)
                    return(1.5 + 0.5 * (1 - exp(-k * (x - ths))))
                }
            }
            hsi <- max(stress_curve(station_stress, ths = 10, max = 21), stress_curve(station_drought, ths = 5, max = 10))

            station_data[station_data$doy == doy, "hsi"] <- hsi
        }

        return(station_data)
    })

    temp <- do.call(rbind, results)
    hsi_table <- aggregate(hsi ~ id + X93 + Y93, data = temp, FUN = max, na.rm = TRUE)
    hsi_table$hsi <- round(hsi_table$hsi, 2)
    names(hsi_table) <- c("id", "X93", "Y93", "hsi")

    cat("Computing hydraulic stress -- OK\n")
    gc()

    return(hsi_table)
}