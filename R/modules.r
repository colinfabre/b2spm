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
    cat("===== TOPOGRAPHY COMPUTER =====\n")
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

    cat("== TOPOGRAPHY COMPUTER -- OK ==\n")
    cat("===============================\n")
    gc()

    return(topography)
}

#' drias_reader
#'
#' This function reads a text file containing climatic data from DRIAS and formats it into a dataframe compatible with the B2SPM pipeline.
#'
#' @param drias_txt_path The path to the text file containing DRIAS climatic data. The data must be comma-separated, the date must be in DD/MM/YYYY format, and it must include tmin (K), tmax (K), tmean (K), tot_pr (kg/m2/s), spec_hum (kg/kg), vis_solrad (W/m2), ir_solrad (W/m2), wind (m/s) and pet (kg/m2/s).
#' @param smoothing Whether the data should be averaged around the central year (must provide an odd number of years) or not. Simulated climate data are usually averaged with a 10-year range on either side of the central year to analyze.
#' @return A data.frame containing the columns: `id`, `X93`, `Y93`, `date`, `doy`, `tmin`, `tmax`, `tmean`, `tot_pr`, `spec_hum`, `vis_solrad`, `ir_solrad`, `wind`, `pet`.
#' @examples
#' \dontrun{
#'  drias_table <- drias_reader("chablais_2050.txt")
#' }
#' @export
drias_reader <- function(drias_txt_path, smoothing = FALSE) {
    cat("===== DRIAS READER =====\n")
    if (!inherits(drias_txt_path, c("character"))) {
        stop("!! ERROR - 'drias_txt_path' must be a valid string reaching to the DRIAS .txt file.")
    }
    if (!is.logical(smoothing)) {
        stop("!! ERROR - 'smoothing' must be TRUE or FALSE.")
    }

    drias_table <- utils::read.table(drias_txt_path, sep = ",", row.names = NULL)
    names(drias_table) <- c("id", "X_LambII", "Y_LambII", "date", "tmin", "tmax", "tmean", "tot_pr", "spec_hum", "vis_solrad", "ir_solrad", "wind", "pet")

    drias_points <- terra::vect(drias_table, geom = c("X_LambII", "Y_LambII"), crs = "EPSG:27572")
    drias_points_2154 <- terra::project(drias_points, "EPSG:2154")
    coords_2154 <- terra::crds(drias_points_2154)
    drias_table <- cbind(drias_table[, 1], coords_2154, drias_table[, 4:13])
    names(drias_table) <- c("id", "X93", "Y93", "date", "tmin", "tmax", "tmean", "tot_pr", "spec_hum", "vis_solrad", "ir_solrad", "wind", "pet")

    stations <- unique(drias_table$id)
    years <- unique(as.numeric(substr(drias_table$date, 7, 10)))
    results <- list()

    for (station in stations) {
        station_data <- drias_table[drias_table$id == station, ]

        for (year in years) {
            
            year_mask <- tryCatch({
                as.numeric(substr(station_data$date, 7, 10)) == year},
                error = function(e) {stop(cat(paste0("!! ERROR -  Impossible to extract and convert records dates for station ", station, ".\n")))})
            doy_val <- seq_len(sum(year_mask))
            station_data$doy[year_mask] <- doy_val
        }

        results[[as.character(station)]] <- station_data
    }

    temp <- do.call(rbind, results)
    
    if (length(years) > 1) {
        if (smoothing && (length(years) %% 2 == 0)) {
            stop("!! ERROR - 'smoothing = TRUE' needs an odd number of years in the provided record.")
        } 
        if (!smoothing && (length(years) %% 2 == 1)) {
            stop("!! ERROR - 'smoothing = FALSE' while there is an odd number of years in the provided record. Fix 'smoothing = TRUE' or provide a single-year record.")
        } 
        if (!smoothing && (length(years) %% 2 == 0)) {
            stop("!! ERROR - 'smoothing = FALSE' and there is an even number of years in the provided record. Fix 'smoothing = TRUE' and provide a record with an odd number of years, or fix 'smoothinng = FALSE' and provide a single-year record.")
        }
        if (smoothing && (length(years) %% 2 == 1)) {
            pack <- (length(years) - 1) / 2
            yta <- years[pack + 1]

            mean_values <- stats::aggregate(temp[, 6:14], by = list(doy = temp$doy), FUN = mean)
            yta_mask <- as.numeric(substr(drias_table$date, 7, 10)) == yta
            doy_indices <- match(drias_table$doy[yta_mask], mean_values$doy)
            drias_table[yta_mask, 6:14] <- mean_values[doy_indices, -1]
        }
    } else {drias_table <- temp}

    drias_table <- cbind(drias_table[, 1:4], drias_table[, "doy"], drias_table[, 5:13])
    names(drias_table) <- c("id", "X93", "Y93", "date", "doy", "tmin", "tmax", "tmean", "tot_pr", "spec_hum", "vis_solrad", "ir_solrad", "wind", "pet")

    drias_table$tmin <- drias_table$tmin - 273.15
    drias_table$tmax <- drias_table$tmax - 273.15
    drias_table$tmean <- drias_table$tmean - 273.15
    drias_table$tot_pr <- drias_table$tot_pr * 86400
    drias_table$pet <- pmax(drias_table$pet * 86400, 0)

    cat("== DRIAS READER -- OK ==\n")
    cat("========================\n")
    gc()

    return(drias_table)
}

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
    cat("===== DRIAS FETCHER =====\n")
    if (!year %in% c(2050, 2075, 2100)) {
        stop("!! ERROR - 'year' must be equal either to 2050, 2075 or 2100.")
    }

    cat("Fetching the online database corresponding to the provided year of analysis can take up to 60s.\n")
    cat("Please wait...\n")
    drias_points <- tryCatch({
        terra::vect(paste0("https://github.com/colinfabre/b2spm_database/raw/refs/heads/main/alpine_arc_", year, ".gpkg"))},
        error = function(e) {stop("!! ERROR -  Impossible to fetch the online database corresponding to the provided year. Please check your internet connection.\n")})
    cat("The online database was correctly fetched.\n")
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

    cat("== DRIAS FETCHER -- OK ==\n")
    cat("=========================\n")
    gc()
    
    return(drias_table)
}

#' phloem_rm
#'
#' This function uses a constrained non-linear radiative model calibrated for spruce to calculate the temperatures beneath the phloem regulating bark beetle development.
#'
#' @param drias_table The DRIAS table processed either by the drias_reader() or the drias_fetcher() function.
#' @return The updated DRIAS table with additional columns `tmin_phloem`, `tmax_phloem`, and `tmean_phloem` containing the temperatures beneath the phloem.
#' @examples
#' \dontrun{
#'  drias_table <- phloem_rm(drias_table)
#' }
#' @export
phloem_rm <- function(drias_table) {
    cat("===== RADIATIVE MODEL FOR UNDER-PHLOEM TEMPERATURE CALCULATION =====\n")

    beta_vis_solrad <- 0.02
    beta_ir_solrad <- 0.03
    lambda <- 0.2
    x <- 5
    h0 <- 10
    gamma_wind <- 0.1
    gamma_spec_hum <- 0.3
    delta <- 0.2
    b <- 1.3

    drias_table$tmin_phloem <- round(drias_table$tmin + ((beta_vis_solrad * exp(-lambda * x) * drias_table$vis_solrad + beta_ir_solrad * exp(-lambda * x) * drias_table$ir_solrad)/(h0 * (1 + gamma_wind * exp(delta * drias_table$wind) * (1 + gamma_spec_hum * (drias_table$spec_hum)^b)))), 2)
    drias_table$tmax_phloem <- round(drias_table$tmax + ((beta_vis_solrad * exp(-lambda * x) * drias_table$vis_solrad + beta_ir_solrad * exp(-lambda * x) * drias_table$ir_solrad)/(h0 * (1 + gamma_wind * exp(delta * drias_table$wind) * (1 + gamma_spec_hum * (drias_table$spec_hum)^b)))), 2)
    drias_table$tmean_phloem <- round(drias_table$tmean + ((beta_vis_solrad * exp(-lambda * x) * drias_table$vis_solrad + beta_ir_solrad * exp(-lambda * x) * drias_table$ir_solrad)/(h0 * (1 + gamma_wind * exp(delta * drias_table$wind) * (1 + gamma_spec_hum * (drias_table$spec_hum)^b)))), 2)

    drias_table$tmin_phloem <- pmin(drias_table$tmin + 2.5, drias_table$tmin_phloem)
    drias_table$tmax_phloem <- pmin(drias_table$tmax + 2.5, drias_table$tmax_phloem)
    drias_table$tmean_phloem <- pmin(drias_table$tmean + 2.5, drias_table$tmean_phloem)

    drias_table <- cbind(drias_table[, 1:5], drias_table[, "tmin"], drias_table[, "tmin_phloem"], drias_table[, "tmax"], drias_table[, "tmax_phloem"], drias_table[, "tmean"], drias_table[, "tmean_phloem"], drias_table[, 9:14])
    names(drias_table) <- c("id", "X93", "Y93", "date", "doy", "tmin", "tmin_phloem", "tmax", "tmax_phloem", "tmean", "tmean_phloem", "tot_pr", "spec_hum", "vis_solrad", "ir_solrad", "wind", "pet")

    cat("== RADIATIVE MODEL FOR UNDER-PHLOEM TEMPERATURE CALCULATION -- OK ==\n")
    cat("====================================================================\n")
    gc()

    return(drias_table)
}

#' ppc
#'
#' This function calculates the photoperiod (time of incident light between sunrise and sunset) according to the latitude, doy and surrounding topography for the cast shadows.
#'
#' @param drias_table The DRIAS table processed by the phloem_rm() function.
#' @param topography The raster stack returned by the topo_comp() function.
#' @param precision The level of computation precision for the cast shadowing. `precision = 1` means the cast shadows are computed as a constant through a 15-days window. `precision = 2` means the cast shadows are daily computed at 12.00 (minimum shadowing). `precision = 3` means the cast shadows are daily computed throughout the course of the sun (at 06.00, 09.00, 12.00, 15.00 and 18.00); in that case, expect the computation time to be multiplied by 2.5 from `precision = 2`.
#' @return The updated DRIAS table with additional column `photoperiod` containing the duration of the photoperiod for each day.
#' @examples
#' \dontrun{
#'  drias_table <- ppc(drias_table, topography)
#' }
#' @export
ppc <- function(drias_table, topography, precision = 1) {
    cat("===== PHOTOPERIOD COMPUTER =====\n")
    if (precision != 1 && precision != 2 && precision != 3) {
        stop("!! ERROR - 'precision' must be set to 1, 2 or 3.")
    }

    stations <- unique(drias_table$id)

    results <- lapply(stations, function(station) {
        station_data <- drias_table[drias_table$id == station, ]
        station_point <- terra::vect(station_data, geom = c("X93", "Y93"), crs = "EPSG:2154")[1]
        station_point_4326 <- terra::project(station_point, "EPSG:4326")

        lat <- round(terra::crds(station_point_4326)[2], 5)
        lat_rad <- round(lat * pi / 180, 5)

        if (precision == 1) {
            hours <- c(12)
            windows_centers <- 8:(max(station_data$doy) - 7)

            for (center in windows_centers) {
                window_days <- (center - 7):(center + 7) 

                cum_cshd <- terra::rast(topography$alt)
                terra::values(cum_cshd) <- 0
                i <- 0

                for (hour in hours) {
                    gamma <- 2 * pi / 365 * (center - 1 + hour / 24)
                    delta <- 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) -
                            0.006758 * cos(2 * gamma) + 0.000907 * sin(2 * gamma) -
                            0.002697 * cos(3 * gamma) + 0.00148 * sin(3 * gamma)

                    hour_angle <- pi * (hour - 12) / 12

                    solar_angle <- pmax(asin(sin(lat_rad) * sin(delta) + cos(lat_rad) * cos(delta) * cos(hour_angle)) * 180 / pi, 0)

                    if (solar_angle > 0) {
                        solar_azimuth <- atan2(
                            -cos(delta) * sin(hour_angle),
                            sin(delta) * cos(lat_rad) - cos(delta) * sin(lat_rad) * cos(hour_angle)
                        )
                        solar_azimuth <- (solar_azimuth * 180 / pi) %% 360

                        terra::values(cum_cshd) <- terra::values(cum_cshd) + terra::values(terra::shade(slope = topography$slope, aspect = topography$aspect, angle = solar_angle, direction = solar_azimuth))

                        cum_cshd <- terra::clamp(cum_cshd, lower = 0, values = TRUE)
                        i <- i + 1
                    }
                }

                mean_cshd <- round(cum_cshd / i, 2)
                station_cshd <- terra::extract(mean_cshd, station_point)[, 2]

                gamma <- 2 * pi / 365 * (center - 1 + 12 / 24)
                delta <- 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) -
                        0.006758 * cos(2 * gamma) + 0.000907 * sin(2 * gamma) -
                        0.002697 * cos(3 * gamma) + 0.00148 * sin(3 * gamma)

                pp <- (24 / pi) * acos(-tan(lat_rad) * tan(delta))

                station_data[station_data$doy %in% window_days, "pp_cshd"] <- round(pp * station_cshd, 1)
            }

        } else {
            if (precision == 2) {
                hours <- c(12)
            } else if (precision == 3) {
            hours <- c(6, 9, 12, 15, 18)
            }

            for (doy in station_data$doy) {
                cum_cshd <- terra::rast(topography$alt)
                terra::values(cum_cshd) <- 0
                i <- 0

                for (hour in hours) {
                    gamma <- 2 * pi / 365 * (doy - 1 + hour / 24)
                    delta <- 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) - 0.006758 * cos(2 * gamma) + 0.000907 * sin(2 * gamma) - 0.002697 * cos(3 * gamma) + 0.00148 * sin(3 * gamma)

                    hour_angle <- pi * (hour - 12) / 12

                    solar_angle <- pmax(asin(sin(lat_rad) * sin(delta) + cos(lat_rad) * cos(delta) * cos(hour_angle)) * 180 / pi, 0)
                    
                    if (solar_angle > 0) {
                        solar_azimuth <- atan2(-cos(delta) * sin(hour_angle), sin(delta) * cos(lat_rad) - cos(delta) * sin(lat_rad) * cos(hour_angle))
                        solar_azimuth <- (solar_azimuth * 180 / pi) %% 360
                        
                        terra::values(cum_cshd) <- terra::values(cum_cshd) + terra::values(terra::shade(slope = topography$slope, aspect = topography$aspect, angle = solar_angle, direction = solar_azimuth))
                        cum_cshd <- terra::clamp(cum_cshd, lower = 0, values = TRUE)

                        i <- i + 1
                    } else {
                        next
                    }
                }

                mean_cshd <- round(cum_cshd / i, 2)
                station_cshd <- terra::extract(mean_cshd, station_point)[, 2]

                gamma <- 2 * pi / 365 * (doy - 1 + 12 / 24)
                delta <- 0.006918 - 0.399912 * cos(gamma) + 0.070257 * sin(gamma) - 0.006758 * cos(2 * gamma) + 0.000907 * sin(2 * gamma) - 0.002697 * cos(3 * gamma) + 0.00148 * sin(3 * gamma)
                pp <- (24 / pi) * acos(-tan(lat_rad) * tan(delta))

                station_data[station_data$doy == doy, "pp_cshd"] <- round(pp * station_cshd, 1)
            }
        }

        return(station_data)
    })

    cat("== PHOTOPERIOD COMPUTER -- OK ==\n")
    cat("================================\n")
    gc()

    return(do.call(rbind, results))
}

#' awakening
#'
#' This function calculates the awakening day (`awakening`) of adult bark beetles, weighted by the topographic modulator CSI_adj.
#'
#' @param drias_table The DRIAS table processed by the phloem_rm() function.
#' @param topography The raster stack returned by the topo_comp() function.
#' @return A data.frame with the columns: `id`, `X93`, `Y93`, `awakening_doy`, and `CSI_adj` (adjusted CSI index).
#' @examples
#' \dontrun{
#'  awakening_table <- awakening(drias_data, topography)
#' }
#' @export
awakening <- function(drias_table, topography) {
    cat("===== AWAKENING CALCULATION =====\n")
    
    stations <- unique(drias_table$id)

    results <- lapply(stations, function(station) {
        station_data <- drias_table[drias_table$id == station, ]
        station_data$dd_cumsum <- cumsum(pmax(0, station_data$tmean_phloem))
        awakening_day <- pmax(min(station_data$doy[station_data$dd_cumsum >= 100]), 0)

        window <- awakening_day - 30
        start_day <- max(1, window)
        station_data <- station_data[station_data$doy >= start_day & station_data$doy <= awakening_day, ]
        station_point <- terra::vect(station_data, geom = c("X93", "Y93"), crs = "EPSG:2154")[1]
        station_data$alt <- terra::extract(topography$alt, station_point)[, 2]
        station_data$aspect <- terra::extract(topography$aspect, station_point)[, 2]

        csi_table <- stats::aggregate(cbind(vis_solrad, ir_solrad) ~ id, data = station_data, sum, na.rm = TRUE)
        csi_table$csi <- csi_table$vis_solrad + csi_table$ir_solrad
        topo_table <- stats::aggregate(cbind(alt, aspect) ~ id, data = station_data, stats::median, na.rm = TRUE)
        csi_table <- merge(csi_table, topo_table, by = "id")
        csi_table$csi_adj <- pmax(round(csi_table$csi * cos(csi_table$aspect * pi / 180) * (1 - (csi_table$alt / 3500)), 2), 0)
        csi_table$csi_adj_log <- log1p(csi_table$csi_adj)

        awakening_day <- ceiling(awakening_day * exp(0.05 * csi_table$csi_adj_log))
        if (awakening_day > 365) {
            awakening_day <- 0
            cat(paste0("<!> Warning - awakening() returned a 0-doy for the station ", station, ", meaning adult bark-beetles won't awake this year.\n"))
            cat("<!> Warning - It's caused by inadapted climate conditions.\n")
        }

        coords <- station_data[1, c("X93", "Y93")]
        awakening_table <- data.frame(id = station, X93 = coords$X93, Y93 = coords$Y93, awakening_doy = awakening_day, csi_adj = csi_table$csi_adj)
        names(awakening_table) <- c("id", "X93", "Y93", "awakening_doy", "csi_adj")
        return(awakening_table)
    })

    cat("== AWAKENING CALCULATION -- OK ==\n")
    cat("=================================\n")
    gc()
    
    return(do.call(rbind, results))
}

#' swarming
#'
#' This function calculates the swarming day (`swarming`) of bark beetles after their awakening (`awakening`).
#'
#' @param drias_table The DRIAS table processed by the phloem_rm() function.
#' @param awakening_table The data.frame returned by the awakening() function.
#' @param topography The raster stack returned by the topo_comp() function.
#' @return A data.frame with the columns: `id`, `X93`, `Y93`, `awakening_doy`, and `swarming_doy`.
#' @examples
#' \dontrun{
#'  swarming_table <- swarming(drias_data, awakening_table, topography)
#' }
#' @export
swarming <- function(drias_table, awakening_table, topography) {
    cat("===== SWARMING CALCULATION =====\n")

    stations <- unique(drias_table$id)

    results <- lapply(stations, function(station) {
        station_data <- drias_table[drias_table$id == station, ]
        awakening_day <- awakening_table$awakening_doy[awakening_table$id == station]

        swarming_day <- max(min(station_data$doy[station_data$doy > awakening_day & station_data$tmean >= 15 & station_data$tmax <= 30 & station_data$tot_pr == 0 & station_data$wind < 3]), awakening_day + 1)
        if (awakening_day == 0) {
            swarming_day <- 0
            cat(paste0("<!> Warning - swarming() returned a 0-doy for station ", station, " meaning adult bark-beetles won't swarm this year.\n"))
            cat("<!> Warning - It's caused by an impossible awakening day (check awakening's results)\n")
        }
        if (swarming_day > 365) {
            swarming_day <- 0
            cat(paste0("<!> Warning - swarming() returned a 0-doy for station ", station, " meaning adult bark-beetles won't swarm this year.\n"))
            cat("<!> Warning - It's caused by inadapted climate conditions.\n")
        }

        coords <- station_data[1, c("X93", "Y93")]
        swarming_table <- data.frame(id = station, X93 = coords$X93, Y93 = coords$Y93, swarming_doy = swarming_day)
        names(swarming_table) <- c("id", "X93", "Y93", "swarming_doy")
        return(swarming_table)
    })

    cat("== SWARMING CALCULATION -- OK ==\n")
    cat("================================\n")
    gc()
    
    return(do.call(rbind, results))
}

#' maturing
#'
#' This function calculates the maturing date (`maturing`) of bark beetle larvae after swarming (`swarming`), weighted by the topographic modulators CSIadj and MDI.
#'
#' @param drias_table The DRIAS table processed by the phloem_rm() function.
#' @param swarming_table The data.frame returned by the swarming() function.
#' @param topography The raster stack returned by the topo_comp() function.
#' @return A data.frame with the columns: `id`, `X93`, `Y93`, `maturing_doy`, `CSI_adj` and `MDI`.
#' @examples
#' \dontrun{
#'  maturing_table <- maturing(drias_data, swarming_table, topography)
#' }
#' @export
maturing <- function(drias_table, swarming_table, topography) {
    cat("===== MATURING CALCULATION =====\n")

    stations <- unique(drias_table$id)

    results <- lapply(stations, function(station) {
        station_data <- drias_table[drias_table$id == station, ]
        swarming_day <- swarming_table$swarming_doy[swarming_table$id == station]
        station_data <- station_data[station_data$doy >= swarming_day, ]
        station_data$dev_cumsum <- cumsum(station_data$tmean_phloem)

        maturing_day <- ifelse(length(station_data$doy[station_data$dev_cumsum >= 143]) > 0, min(station_data$doy[station_data$dev_cumsum >= 143]), 0)

        if (swarming_day != 0) {
            window <- swarming_day - 30
            start_day <- max(1, window)
            topomod_data <- station_data[station_data$doy >= start_day & station_data$doy <= swarming_day, ]
            station_point <- terra::vect(topomod_data, geom = c("X93", "Y93"), crs = "EPSG:2154")[1]
            topomod_data$alt <- terra::extract(topography$alt, station_point)[, 2]
            topomod_data$aspect <- terra::extract(topography$aspect, station_point)[, 2]

            csi_table <- stats::aggregate(cbind(vis_solrad, ir_solrad) ~ id, data = topomod_data, sum, na.rm = TRUE)
            csi_table$csi <- csi_table$vis_solrad + csi_table$ir_solrad
            topo_table <- stats::aggregate(cbind(alt, aspect) ~ id, data = topomod_data, stats::median, na.rm = TRUE)
            csi_table <- merge(csi_table, topo_table, by = "id")
            csi_table$csi_adj <- pmax(round(csi_table$csi * cos(csi_table$aspect * pi / 180) * (1 - (csi_table$alt / 3500)), 2), 0)
            csi_table$csi_adj_log <- log1p(csi_table$csi_adj)

            mdi_table <- stats::aggregate(cbind(tot_pr, tmean) ~ id, data = topomod_data, sum, na.rm = TRUE)
            mdi_table$pr_t <- mdi_table$tot_pr / mdi_table$tmean
            topo_table <- stats::aggregate(cbind(alt, aspect) ~ id, data = topomod_data, stats::median, na.rm = TRUE)
            mdi_table <- merge(mdi_table, topo_table, by = "id")
            mdi_table$mdi <- pmax(round(mdi_table$pr_t * cos(mdi_table$aspect * pi / 180) * (1 - (mdi_table$alt / 3500)), 2), 0)
            mdi_table$mdi_log <- log1p(mdi_table$mdi)

            maturing_day <- ceiling(maturing_day * exp(0.005 * csi_table$csi_adj_log * mdi_table$mdi_log))
        } else {
            station_point <- terra::vect(station_data, geom = c("X93", "Y93"), crs = "EPSG:2154")[1]
            station_point <- terra::project(station_point, terra::crs(topography))
            station_data$alt <- terra::extract(topography$alt, station_point)[, 2]
            station_data$aspect <- terra::extract(topography$aspect, station_point)[, 2]

            csi_table <- stats::aggregate(cbind(vis_solrad, ir_solrad) ~ id, data = station_data, sum, na.rm = TRUE)
            csi_table$csi <- csi_table$vis_solrad + csi_table$ir_solrad
            topo_table <- stats::aggregate(cbind(alt, aspect) ~ id, data = station_data, stats::median, na.rm = TRUE)
            csi_table <- merge(csi_table, topo_table, by = "id")
            csi_table$csi_adj <- pmax(round(csi_table$csi * cos(csi_table$aspect * pi / 180) * (1 - (csi_table$alt / 3500)), 2), 0)
            csi_table$csi_adj_log <- log1p(csi_table$csi_adj)

            mdi_table <- stats::aggregate(cbind(tot_pr, tmean) ~ id, data = station_data, sum, na.rm = TRUE)
            mdi_table$pr_t <- mdi_table$tot_pr / mdi_table$tmean
            topo_table <- stats::aggregate(cbind(alt, aspect) ~ id, data = station_data, stats::median, na.rm = TRUE)
            mdi_table <- merge(mdi_table, topo_table, by = "id")
            mdi_table$mdi <- pmax(round(mdi_table$pr_t * cos(mdi_table$aspect * pi / 180) * (1 - (mdi_table$alt / 3500)), 2), 0)
            mdi_table$mdi_log <- log1p(mdi_table$mdi)

            maturing_day <- 0
            cat(paste0("<!> Warning - maturing() returned a 0-doy for station ", station, " meaning bark-beetle pupae won't mature this year.\n"))
            cat("<!> Warning - It's caused by an impossible swarming day.\n")
        }

        if (maturing_day > 365) {
            maturing_day <- 0
            cat(paste0("<!> Warning - maturing() returned a 0-doy for station ", station, " meaning bark-beetle pupae won't mature this year.\n"))
            cat("<!> Warning - It's caused by inadapted climate conditions.\n")
        }

        coords <- station_data[1, c("X93", "Y93")]
        maturing_table <- data.frame(id = station, X93 = coords$X93, Y93 = coords$Y93, maturing_doy = maturing_day, csi_adj = csi_table$csi_adj, mdi = mdi_table$mdi)
        names(maturing_table) <- c("id", "X93", "Y93", "maturing_doy", "csi_adj", "mdi")
        return(maturing_table)
    })
    
    cat("== MATURING CALCULATION -- OK ==\n")
    cat("================================\n")
    gc()
    
    return(do.call(rbind, results))
}

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
    cat("===== HYDRAULIC STRESS COMPUTER =====\n")
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

    cat("== HYDRAULIC STRESS COMPUTER -- OK ==\n")
    cat("=====================================\n")
    gc()

    return(hsi_table)
}

#' kpi
#'
#' This function spatializes phenological indicators (`awakening_doy`, `swarming_doy`, `maturing_doy`) and hydraulic stress index (`hsi`) using IDW (nmax = 8 and idp = 2).
#'
#' @param awakening_table The data.frame returned by the awakening() function.
#' @param swarming_table The data.frame returned by the swarming() function.
#' @param maturing_table The data.frame returned by the maturing() function.
#' @param hsi_table The data.frame returned by the hsc() function.
#' @param topography The raster stack returned by the topo_comp() function.
#' @return A raster stack containing the interpolated phenological indicators and hydraulic stress index, restricted to pure spruce forests.
#' @examples
#' \dontrun{
#'  spat_ind <- kpi(awakening_table, swarming_table, maturing_table, hsi_table, topography)
#' }
#' @export
kpi <- function(awakening_table, swarming_table, maturing_table, hsi_table, topography) {
    cat("===== INDICATORS SPATIALISATION =====\n")

    data <- Reduce(function(x, y) merge(x, y, by = c("id", "X93", "Y93"), all = TRUE), list(awakening_table, swarming_table, maturing_table, hsi_table))
    data <- data[, c(1:4, 6:7, 10)]
    data <- data["awakening_doy" != 0, ]
    data_sf <- sf::st_as_sf(data, coords = c("X93", "Y93"), crs = 2154)
    vars <- c("awakening_doy", "swarming_doy", "maturing_doy", "hsi")

    grid <- terra::rast(terra::ext(topography), crs = terra::crs(topography), resolution = 250)
    grid_df <- as.data.frame(grid, xy = TRUE, na.rm = FALSE)
    grid_sf <- sf::st_as_sf(grid_df, coords = c("x", "y"), crs = 2154)
    spat_ind <- terra::rast()

    idw_spationer <- function(var) {
        formula <- stats::as.formula(paste(var, "~1"))

        idw_result <- gstat::idw(formula = formula, locations = data_sf, newdata = grid_sf, idp = 2, nmax = 8, debug.level = -1)
        idw_df <- cbind(sf::st_coordinates(idw_result), idw_result$var1.pred)
        colnames(idw_df) <- c("x", "y", var)

        idwed_ind <- terra::rast(idw_df, type = "xyz", crs = terra::crs(grid))
        idwed_ind <- terra::resample(idwed_ind, grid, method = "near")
        idwed_ind[is.na(terra::resample(topography$spruce_forests, grid, method = "near"))] <- NA
        idwed_ind[idwed_ind < 0] <- 0
        if (var != "hsi") {
            terra::values(idwed_ind) <- ceiling(terra::values(idwed_ind)) 
        }
        names(idwed_ind) <- var

        return(idwed_ind)
    }

    idwed_vars <- lapply(vars, idw_spationer)
    spat_ind <- do.call(c, idwed_vars)

    cat("== INDICATORS SPATIALISATION -- OK ==\n")
    cat("==================================================\n")
    gc()

    return(spat_ind)
}

#' rpc
#'
#' This function calculates the global epidemic risk indicator (`Rpheno`), measuring the probability of successful development of a bark beetle generation based on the indicators `awakening_doy`, `swarming_doy`, and `maturing_doy`, and the hydraulic stress endured by spruce trees throughout the year.
#'
#' @param spat_ind The raster stack returned by the kpi() function.
#' @return A probability raster (`Rpheno`), ranging from 0 to 1, indicating the global epidemic risk of bark beetle generation development and attack on water-stress spruce forests.
#' @examples
#' \dontrun{
#'  rpheno <- rpc(spat_ind)
#' }
#' @export
rpc <- function(spat_ind) {
    cat("===== Rpheno CALCULATION =====\n")

    prob_awakening <- (-0.0028 * spat_ind$awakening_doy) + 1.0282
    prob_swarming <- (-0.0028 * spat_ind$swarming_doy) + 1.0311
    prob_maturing <- (-0.0032 * spat_ind$maturing_doy) + 1.1699
    prob_hsi <- spat_ind$hsi

    rpheno <- terra::rast(spat_ind$awakening_doy)
    terra::values(rpheno) <- round(terra::values(prob_awakening) * terra::values(prob_swarming) * terra::values(prob_maturing) * terra::values(prob_hsi), 2)
    rpheno <- terra::clamp(rpheno, lower = 0, upper = 1, values = TRUE)
    names(rpheno) <- "Rpheno"

    cat("== Rpheno CALCULATION -- OK ==\n")
    cat("==============================\n")
    gc()

    return(rpheno)
}

#' pipeline
#'
#' Runs the entire B2SPM pipeline, from computing under-phloem temperatures to spatializing the attack risk.
#'
#' @param topography The raster stack returned by the topo_comp() function.
#' @param drias_table The DRIAS table processed either by the drias_reader() or the drias_fetcher() function.
#' @param return_tables Whether the intermediate phenological tables (`awakening_table`, `swarming_table` and `maturing_table`) should be returned as ouputs of the pipeline or not.
#' @param precision The level of computation precision of the `pcc()` function. See the `precision` paramater in the `ppc()` function documentation for more details.
#' @return A list of the raster stack containing the spatialized phenological indicators (`awakening_doy`, `swarming_doy`, `maturing_doy`), the attack risk (`Rpheno`), and the maximum number of generations per year (`max_gen`), and the intermediate phenological tables if `return_tables = TRUE`.
#' @examples
#' \dontrun{
#'  results <- pipeline(topography, drias_table)
#' }
#' @export
pipeline <- function(topography, drias_table, return_tables = FALSE, precision = 1) {
    cat("\n")
    cat("+--------------------------------------------------------------------------------------------------+\n")
    cat("|-------------------------------- B2SPM PIPELINE INITIALISATION... --------------------------------|\n")
    cat("+--------------------------------------------------------------------------------------------------+\n")
    cat("\n")
    gc()

    start_time <- Sys.time()

    if (!inherits(topography, "SpatRaster") && any(names(topography) != c("spruce_forests", "alt", "slope", "aspect", "tpi"))) {
        stop("!! ERROR - 'topography' is invalid. It must necessarily be the result of the topo_comp() function.")
    }
    if (!inherits(drias_table, "data.frame") && any(names(drias_table) != c("id", "X93", "Y93", "date", "doy", "tmin", "tmax", "tmean", "tot_pr", "spec_hum", "vis_solrad", "ir_solrad", "wind", "pet"))) {
        stop("!! ERROR - 'drias_table' is invalid. It must necessarily be the result of the drias_reader() or drias_fetcher() function.")
    }
    if (!is.logical(return_tables)) {
        stop("!! ERROR - 'return_tables' must be TRUE or FALSE.")
    }
    if (precision != 1 && precision != 2 && precision != 3) {
        stop("!! ERROR - 'precision' must be set to 1, 2 or 3.")
    }

    cat("Checking available RAM...\n")
    if (round(ps::ps_system_memory()$total / 1024^3, 0) < 4) {
        cat("<!> Warning - The required 4GB of RAM is not available for pipeline to launch.\n")
        clear_mem <- readline("This will clear the current R environment, excepted 'drias_table' and 'topography' objects. Do you accept? (y/n): ")
        while (toupper(clear_mem) != "y" && toupper(clear_mem) != "n") {
            cat("Incorrect input, you must answer by yes (y) or no (n). Try again.\n")
            clear_mem <- readline("(y/n): ")
        }

        if (clear_mem == "y") {
            cat("The current R environment will be cleared.\n")
            env <- ls(envir = .GlobalEnv)
            env <- setdiff(env, c("drias_table", "topography"))
            rm(env)
            gc()
            cat("The current R environment has been correctly cleared.\n")

            if (round(ps::ps_system_memory()$total / 1024^3, 0) < 4) {
                stop("!! ERROR - The required 4GB of RAM is still not available for pipeline to launch, even after a complete R environment clearing. The code will abort.")
            } else {
                cat("The required 4GB of RAM is now available for pipeline to launch.")
            }
            
        } else {
            stop("The pipeline couldn't be launch due to insufficient RAM. Please clear some space first.\n")
        }
    } else {
        cat("Available RAM -- OK\n")
        cat("\n")
    }

    cat("Checking data compatibility...\n")
    if (terra::ext(topography$alt) < terra::ext(terra::vect(drias_table, geom = c("X93", "Y93"), crs = terra::crs(topography))) + 250) {
        stop("!! ERROR - The input DEM in the `topo_comp()` function is smaller than the ROI + the 250m required buffer.\n")
    } else {
        cat("Data compatibility -- OK\n")
        cat("\n")
    }

    drias_table <- phloem_rm(drias_table)
    cat("\n")

    drias_table <- ppc(drias_table, topography, precision = precision)
    cat("\n")

    awakening_table <- awakening(drias_table, topography)
    cat("\n")
    swarming_table <- swarming(drias_table, awakening_table, topography)
    cat("\n")
    maturing_table <- maturing(drias_table, swarming_table, topography)
    cat("\n")

    hsi_table <- hsc(drias_table)
    cat("\n")

    spat_ind <- kpi(awakening_table, swarming_table, maturing_table, hsi_table, topography)
    cat("\n")

    rpheno <- rpc(spat_ind)
    cat("\n")

    if (return_tables == FALSE) {
        results <- list(spat_ind, rpheno)
        names(results) <- c("spat_ind", "rpheno")
    } else {
        results <- list(awakening_table, swarming_table, maturing_table, hsi_table, spat_ind, rpheno)
        names(results) <- c("awakening_table", "swarming_table", "maturing_table", "hsi_table", "spat_ind", "rpheno")
    }

    end_time <- Sys.time()
    process_time <- end_time - start_time
    cat(paste0("Processing time: ", process_time, "\n"))

    cat("\n")
    cat("+--------------------------------------------------------------------------------------------------+\n")
    cat("|-------------------------------------- B2SPM PIPELINE -- OK --------------------------------------|\n")
    cat("+--------------------------------------------------------------------------------------------------+\n")
    cat("\n")
    gc()

    return(results)
}