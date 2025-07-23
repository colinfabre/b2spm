#' ppc
#'
#' This function calculates the photoperiod (time of incident light between sunrise and sunset) according to the latitude, doy and surrounding topography for the cast shadows.
#'
#' @param drias_table The DRIAS table processed by the phloem_rm() function.
#' @param topography The raster stack returned by the topo_comp() function.
#' @return The updated DRIAS table with additional columns `cshd` containing the cast-shawoding factor, and `pp_cshd` containing the duration of the photoperiod for each day moderated by that factor.
#' @examples
#' \dontrun{
#'  drias_table <- ppc(drias_table, topography)
#' }
#' @export
ppc <- function(drias_table, topography) {
    cat("Computing the cast-shadowing-moderated photoperiod using FOURIER's Series...\n")
    cat("<!> Warning - This step can be quite long, depending on your computing resources and the size of your study area.\n")

    hour <- c(12)
    months <- 1:12
    astro_table <- expand.grid(month = months, hour = hour)
    astro_table$doy <- floor(seq(from = 1, to = 365, length.out = 12 * 24))[astro_table$month]
    astro_table$gamma <- 2 * pi / 365 * (astro_table$doy - 1 + astro_table$hour / 24)
    astro_table$delta <- with(astro_table, 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) - 0.006758 * cos(2 * gamma) + 0.000907 * sin(2 * gamma) - 0.002697 * cos(3 * gamma) + 0.00148 * sin(3 * gamma))
    astro_table$hour_angle <- pi * (astro_table$hour - 12) / 12

    stations <- unique(drias_table$id)
    results <- lapply(stations, function(station) {
        station_data <- drias_table[drias_table$id == station, ]
        station_point <- terra::vect(station_data, geom = c("X93", "Y93"), crs = "EPSG:2154")[1]
        station_point_4326 <- terra::project(station_point, "EPSG:4326")
        lat <- round(terra::crds(station_point_4326)[2], 5)
        lat_rad <- lat * pi / 180

    for (month in months) {
        cum_cshd <- terra::rast(topography$alt)
        terra::values(cum_cshd) <- 0
        i <- 0
        precomputed_rows <- astro_table[astro_table$month == month & astro_table$hour %in% 12, ]

        for (hour in 12) {
            row <- precomputed_rows[precomputed_rows$hour == hour, ]
            delta <- row$delta
            hour_angle <- row$hour_angle
            solar_angle <- pmax(asin(sin(lat_rad) * sin(delta) + cos(lat_rad) * cos(delta) * cos(hour_angle)) * 180 / pi, 0)

            if (solar_angle > 0) {
            solar_azimuth <- atan2(-cos(delta) * sin(hour_angle), sin(delta) * cos(lat_rad) - cos(delta) * sin(lat_rad) * cos(hour_angle))
            solar_azimuth <- (solar_azimuth * 180 / pi) %% 360
            terra::values(cum_cshd) <- terra::values(cum_cshd) + terra::values(terra::shade(topography$slope, topography$aspect, angle = solar_angle, direction = solar_azimuth))
            cum_cshd <- terra::clamp(cum_cshd, lower = 0, values = TRUE)
            i <- i + 1
            }
        }

        mean_cshd <- round(cum_cshd / i, 2)
        station_cshd <- terra::extract(mean_cshd, station_point)[, 2]
        row12 <- astro_table[astro_table$month == month & astro_table$hour == 12, ]
        pp <- (24 / pi) * acos(-tan(lat_rad) * tan(row12$delta))

        station_data[station_data$month == month, "cshd"] <- station_cshd
        station_data[station_data$month == month, "pp_cshd"] <- round(pp * station_cshd, 1)
    }
    gc()

    return(station_data)
  })

    cat("Computing photoperiod -- OK\n")
    gc()
    return(do.call(rbind, results))
}