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
    cat("Spatializing phenological indicators and hydraulic stress index...\n")

    data <- Reduce(function(x, y) merge(x, y, by = c("id", "X93", "Y93"), all = TRUE), list(awakening_table, swarming_table, maturing_table, hsi_table))
    data <- data[, c(1:4, 7:8, 12)]
    data <- data["awakening_doy" != 0, ]
    data_sf <- sf::st_as_sf(data, coords = c("X93", "Y93"), crs = 2154)
    vars <- c("awakening_doy", "swarming_doy", "maturing_doy", "hsi")

    grid <- terra::rast(terra::ext(topography), crs = terra::crs(topography), resolution = 250, vals = 0)
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

    cat("Spatializing indicators -- OK\n")
    gc()

    return(spat_ind)
}