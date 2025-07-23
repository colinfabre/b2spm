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
    cat("Reading the input DRIAS table path...\n")
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

    cat("Reading DRIAS table -- OK\n")
    gc()

    return(drias_table)
}