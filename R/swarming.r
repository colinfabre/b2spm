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
    cat("Calculating the swarming doy...\n")

    stations <- unique(drias_table$id)

    results <- lapply(stations, function(station) {
        station_data <- drias_table[drias_table$id == station, ]
        awakening_day <- awakening_table$awakening_doy[awakening_table$id == station]

        swarming_day <- suppressWarnings(max(min(station_data$doy[station_data$doy > awakening_day & station_data$tmean >= 15 & station_data$tmax <= 30 & station_data$tot_pr == 0 & station_data$wind < 3]), awakening_day + 1))
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

    cat("Calculating swarming doy -- OK\n")
    gc()
    
    return(do.call(rbind, results))
}