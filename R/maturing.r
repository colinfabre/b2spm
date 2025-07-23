#' maturing
#'
#' This function calculates the maturing date (`maturing`) of bark beetle larvae after swarming (`swarming`), weighted by the topographic modulators CSIadj and SVI.
#'
#' @param drias_table The DRIAS table processed by the phloem_rm() function.
#' @param swarming_table The data.frame returned by the swarming() function.
#' @param topography The raster stack returned by the topo_comp() function.
#' @return A data.frame with the columns: `id`, `X93`, `Y93`, `maturing_doy`, `CSI_adj` and `SVI`.
#' @examples
#' \dontrun{
#'  maturing_table <- maturing(drias_data, swarming_table, topography)
#' }
#' @export
maturing <- function(drias_table, swarming_table, topography) {
    cat("Calculating maturing doy...\n")

    stations <- unique(drias_table$id)

    results <- lapply(stations, function(station) {
        station_data <- drias_table[drias_table$id == station, ]
        swarming_day <- swarming_table$swarming_doy[swarming_table$id == station]
        station_data <- station_data[station_data$doy >= swarming_day, ]
        station_data$dev_cumsum <- cumsum(station_data$tmean_phloem)

        maturing_day <- ifelse(length(station_data$doy[station_data$dev_cumsum >= 145]) > 0, min(station_data$doy[station_data$dev_cumsum >= 145]), 0)

        if (swarming_day != 0) {
            window <- swarming_day - 50
            start_day <- max(1, window)
            topomod_data <- station_data[station_data$doy >= start_day & station_data$doy <= swarming_day, ]
            station_point <- terra::vect(topomod_data, geom = c("X93", "Y93"), crs = "EPSG:2154")[1]
            topomod_data$alt <- terra::extract(topography$alt, station_point)[, 2]
            topomod_data$aspect <- terra::extract(topography$aspect, station_point)[, 2]
            topomod_data$tpi <- terra::extract(topography$tpi, station_point)[, 2]

            highest_spruces <- max(topography$alt[topography$spruce_forests == 1])

            csi_table <- stats::aggregate(cbind(vis_solrad, ir_solrad) ~ id, data = topomod_data, sum, na.rm = TRUE)
            csi_table$csi <- csi_table$vis_solrad + csi_table$ir_solrad
            csi_table$csi_norm <- (csi_table$csi - 100) / 29900
            topo_table <- stats::aggregate(cbind(alt, aspect) ~ id, data = topomod_data, stats::median, na.rm = TRUE)
            csi_table <- merge(csi_table, topo_table, by = "id")
            csi_table$spruces_alt_norm <- (((highest_spruces - csi_table$alt) / highest_spruces) + 34) / (-34)
            csi_table$csi_adj_norm <- pmin(pmax(round(mean(csi_table$csi_norm + (1 + cos(csi_table$aspect - pi )) + csi_table$spruces_alt_norm), 3), 0), 1)

            svi_table <- stats::aggregate(cbind(tot_pr, tmax) ~ id, data = topomod_data, sum, na.rm = TRUE)
            svi_table$t_pr <- svi_table$tmax / svi_table$tot_pr
            svi_table$t_pr_norm <- (svi_table$t_pr + 2.5) / (400 - 2.5)
            topo_table <- stats::aggregate(cbind(alt, tpi) ~ id, data = topomod_data, stats::median, na.rm = TRUE)
            svi_table <- merge(svi_table, topo_table, by = "id")
            svi_table$spruces_alt_norm <- (((highest_spruces - csi_table$alt) / highest_spruces) + 34) / (-34)
            svi_table$tpi_norm <- (1 + (svi_table$tpi / 100) + 1) / 2
            svi_table$svi_norm <- pmin(pmax(round(mean(svi_table$t_pr_norm + svi_table$spruces_alt_norm + svi_table$tpi_norm), 3), 0), 1)

            if (csi_table$csi_adj_norm + svi_table$svi_norm <= 1) {
                maturing_day_adj <- ceiling(maturing_day + 100 - 45 * (csi_table$csi_adj_norm + svi_table$svi_norm))
            } else {
                maturing_day_adj <- ceiling(maturing_day + 65 - 10 * (csi_table$csi_adj_norm + svi_table$svi_norm))
            }
            if (maturing_day_adj - swarming_day < 45) {
                maturing_day_adj <- swarming_day + 45
            }
        } else {
            station_point <- terra::vect(station_data, geom = c("X93", "Y93"), crs = "EPSG:2154")[1]
            station_data$alt <- terra::extract(topography$alt, station_point)[, 2]
            station_data$aspect <- terra::extract(topography$aspect, station_point)[, 2]
            station_data$tpi <- terra::extract(topography$tpi, station_point)[, 2]

            highest_spruces <- max(topography$alt[topography$spruce_forests == 1])

            csi_table <- stats::aggregate(cbind(vis_solrad, ir_solrad) ~ id, data = station_data, sum, na.rm = TRUE)
            csi_table$csi <- csi_table$vis_solrad + csi_table$ir_solrad
            csi_table$csi_norm <- (csi_table$csi - 100) / 29900
            topo_table <- stats::aggregate(cbind(alt, aspect) ~ id, data = station_data, stats::median, na.rm = TRUE)
            csi_table <- merge(csi_table, topo_table, by = "id")
            csi_table$spruces_alt_norm <- (((highest_spruces - csi_table$alt) / highest_spruces) + 34) / (-34)
            csi_table$csi_adj_norm <- pmin(pmax(round(mean(csi_table$csi_norm + (1 + cos(csi_table$aspect - pi )) + csi_table$spruces_alt_norm), 3), 0), 1)

            svi_table <- stats::aggregate(cbind(tot_pr, tmax) ~ id, data = station_data, sum, na.rm = TRUE)
            svi_table$t_pr <- svi_table$tmax / svi_table$tot_pr
            svi_table$t_pr_norm <- (svi_table$t_pr + 2.5) / (400 - 2.5)
            topo_table <- stats::aggregate(cbind(alt, tpi) ~ id, data = station_data, stats::median, na.rm = TRUE)
            svi_table <- merge(svi_table, topo_table, by = "id")
            svi_table$spruces_alt_norm <- (((highest_spruces - csi_table$alt) / highest_spruces) + 34) / (-34)
            svi_table$tpi_norm <- (1 + (svi_table$tpi / 100) + 1) / 2
            svi_table$svi_norm <- pmin(pmax(round(mean(svi_table$t_pr_norm + svi_table$spruces_alt_norm + svi_table$tpi_norm), 3), 0), 1)

            maturing_day_adj <- 0
            cat(paste0("<!> Warning - maturing() returned a 0-doy for station ", station, " meaning bark-beetle pupae won't mature this year.\n"))
            cat("<!> Warning - It's caused by an impossible swarming day.\n")
        }

        if (maturing_day_adj > 365) {
            maturing_day_adj <- 0
            cat(paste0("<!> Warning - maturing() returned a 0-doy for station ", station, " meaning bark-beetle pupae won't mature this year.\n"))
            cat("<!> Warning - It's caused by inadapted climate or site conditions.\n")
        }

        coords <- station_data[1, c("X93", "Y93")]
        maturing_table <- data.frame(id = station, X93 = coords$X93, Y93 = coords$Y93, maturing_doy = maturing_day_adj, maturing_doy_raw = maturing_day, csi_adj_norm = csi_table$csi_adj_norm, svi_norm = svi_table$svi_norm)
        names(maturing_table) <- c("id", "X93", "Y93", "maturing_doy", "maturing_doy_raw", "csi_adj", "svi")
        return(maturing_table)
    })
    
    cat("Calculating maturing doy -- OK\n")
    gc()
    
    return(do.call(rbind, results))
}