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
    cat("Calculating the awakening doy...\n")
    
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
        
        highest_spruces <- max(topography$alt[topography$spruce_forests == 1])

        csi_table <- stats::aggregate(cbind(vis_solrad, ir_solrad) ~ id, data = station_data, sum, na.rm = TRUE)
        csi_table$csi <- csi_table$vis_solrad + csi_table$ir_solrad
        csi_table$csi_norm <- (csi_table$csi - 100) / 29900
        topo_table <- stats::aggregate(cbind(alt, aspect) ~ id, data = station_data, stats::median, na.rm = TRUE)
        csi_table <- merge(csi_table, topo_table, by = "id")
        csi_table$spruces_alt_norm <- (((highest_spruces - csi_table$alt) / highest_spruces) + 34) / (-34)
        csi_table$csi_adj_norm <- pmin(pmax(round(mean(csi_table$csi_norm + (1 + cos(csi_table$aspect - pi )) + csi_table$spruces_alt_norm), 3), 0), 1)

        awakening_day_adj <- ceiling(awakening_day + 45 * (1 - csi_table$csi_adj_norm)^0.5)
        if (awakening_day_adj > 365) {
            awakening_day_adj <- 0
            cat(paste0("<!> Warning - awakening() returned a 0-doy for the station ", station, ", meaning adult bark-beetles won't awake this year.\n"))
            cat("<!> Warning - It's caused by inadapted climate or site conditions.\n")
        }
        coords <- station_data[1, c("X93", "Y93")]
        awakening_table <- data.frame(id = station, X93 = coords$X93, Y93 = coords$Y93, awakening_doy = awakening_day_adj, awakening_day_raw = awakening_day, csi_adj_norm = csi_table$csi_adj_norm)
        names(awakening_table) <- c("id", "X93", "Y93", "awakening_doy", "awakening_doy_raw", "csi_adj")
        return(awakening_table)
    })

    cat("Calculating awakening doy -- OK\n")
    gc()
    
    return(do.call(rbind, results))
}