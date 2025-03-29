#' topo_comp
#'
#' This function loads the spruce forest mask and computes topographic data for the study area using the provided DEM.
#'
#' @param dem A raster (preferably a SpatRaster in EPSG:2154) representing elevation, covering the entire study area plus a 5-pixel buffer, at a 25m-pixel spatial resolution. It has to be located within the French alpine arc (see README.md for further information).
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

            slope <- terra::terrain(dem, v = "slope", unit = "degrees")
            aspect <- terra::terrain(dem, v = "aspect", unit = "radians")
            tpi <- terra::terrain(dem, v = "TPI")
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
#' @param drias_txt_path The path to the text file containing DRIAS climatic data. The data must be comma-separated, the date must be in DD/MM/YYYY format, and it must include tmin (K), tmax (K), tmean (K), tot_pr (kg/m2/s), spec_hum (kg/kg), vis_solrad (W/m2), ir_solrad (W/m2), and wind (m/s).
#' @param smoothing Should the data be averaged around the central year (must provide an odd number of years). Simulated climate data are usually averaged with a 10-year range on either side of the central year to analyze.
#' @return A data.frame containing the columns: `id`, `X93`, `Y93`, `date`, `doy`, `tmin`, `tmax`, `tmean`, `tot_pr`, `spec_hum`, `vis_solrad`, `ir_solrad`, `wind`.
#' @examples
#' \dontrun{
#'  drias_table <- drias_reader("Chablais_2030.txt")
#' }
#' @export
drias_reader <- function(drias_txt_path, smoothing = FALSE) {
    if (!is.logical(smoothing)) {
        stop("!! ERROR - 'smoothing' must be TRUE or FALSE.")
    }

    drias_table <- utils::read.table(drias_txt_path, sep = ",", row.names = NULL)
    names(drias_table) <- c("id", "X93", "Y93", "date", "tmin", "tmax", "tmean", "tot_pr", "spec_hum", "vis_solrad", "ir_solrad", "wind")

    cat("===== DAY OF YEAR =====\n")
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

            mean_values <- stats::aggregate(temp[, 6:13], by = list(doy = temp$doy), FUN = mean)
            yta_mask <- as.numeric(substr(drias_table$date, 7, 10)) == yta
            doy_indices <- match(drias_table$doy[yta_mask], mean_values$doy)
            drias_table[yta_mask, 6:13] <- mean_values[doy_indices, -1]
        }
    } else {drias_table <- temp}

    drias_table <- cbind(drias_table[, 1:4], drias_table[, "doy"], drias_table[, 5:12])
    names(drias_table) <- c("id", "X93", "Y93", "date", "doy", "tmin", "tmax", "tmean", "tot_pr", "spec_hum", "vis_solrad", "ir_solrad", "wind")

    drias_table$tmin <- drias_table$tmin - 273.15
    drias_table$tmax <- drias_table$tmax - 273.15
    drias_table$tmean <- drias_table$tmean - 273.15
    drias_table$tot_pr <- drias_table$tot_pr * 86400

    cat("== DAY OF YEAR -- OK ==\n")
    cat("=======================\n")
    gc()

    return(drias_table)
}

#' phloem_rm
#'
#' This function uses a constrained non-linear radiative model calibrated for spruce to calculate the temperatures beneath the phloem regulating bark beetle development.
#'
#' @param drias_table The DRIAS table processed by the drias_reader() function.
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
    
    cat("== RADIATIVE MODEL FOR UNDER-PHLOEM TEMPERATURE CALCULATION -- OK ==\n")
    cat("====================================================================\n")
    gc()

    drias_table <- cbind(drias_table[, 1:5], drias_table[, "tmin"], drias_table[, "tmin_phloem"], drias_table[, "tmax"], drias_table[, "tmax_phloem"], drias_table[, "tmean"], drias_table[, "tmean_phloem"], drias_table[, 12:16])
    names(drias_table) <- c("id", "X93", "Y93", "date", "doy", "tmin", "tmin_phloem", "tmax", "tmax_phloem", "tmean", "tmean_phloem", "tot_pr", "spec_hum", "vis_solrad", "ir_solrad", "wind")

    return(drias_table)
}

#' awakening
#'
#' This function calculates the awakening day (`awakening`) of adult bark beetles, weighted by the topographic modulator CSI_adj.
#'
#' @param drias_table The DRIAS table processed by the phloem_rm() function.
#' @param topography The raster stack returned by the data_fetcher() function.
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
        station_point <- terra::vect(station_data, geom = c("X93", "Y93"), crs = "EPSG:27572")
        terra::crs(station_point) <- "EPSG:27572"
        station_point <- terra::project(station_point, terra::crs(topography))
        station_data$alt <- terra::extract(topography$alt, station_point)[, 2]
        station_data$aspect <- terra::extract(topography$aspect, station_point)[, 2]

        csi_table <- stats::aggregate(cbind(vis_solrad, ir_solrad) ~ id, data = station_data, sum, na.rm = TRUE)
        csi_table$csi <- csi_table$vis_solrad + csi_table$ir_solrad
        topo_table <- stats::aggregate(cbind(alt, aspect) ~ id, data = station_data, stats::median, na.rm = TRUE)
        csi_table <- merge(csi_table, topo_table, by = "id")
        csi_table$csi_adj <- pmax(round(csi_table$csi * cos(csi_table$aspect * pi / 180) * (1 - (csi_table$alt / 3500)), 2), 0)
        csi_table$csi_adj_log <- log1p(csi_table$csi_adj)

        awakening_day <- ceiling(awakening_day * exp(0.2 * csi_table$csi_adj_log))
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
#' @param topography The raster stack returned by the data_fetcher() function.
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

        swarming_day <- max(station_data$doy[station_data$doy > awakening_day & station_data$tmean >= 16.11 & station_data$tmax <= 31.29], awakening_day + 1)
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
#' @param topography The raster stack returned by the data_fetcher() function.
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
            station_point <- terra::vect(topomod_data, geom = c("X93", "Y93"), crs = "EPSG:27572")
            terra::crs(station_point) <- "EPSG:27572"
            station_point <- terra::project(station_point, terra::crs(topography))
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
            station_point <- terra::vect(station_data, geom = c("X93", "Y93"), crs = "EPSG:27572")
            terra::crs(station_point) <- "EPSG:27572"
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

#' kpi
#'
#' This function spatializes phenological indicators (`awakening_doy`, `swarming_doy`, `maturing_doy`) using IDW or ordinary kriging with an exponential model, depending on the number of DRIAS points present in the study area.
#'
#' @param awakening_table The data.frame returned by the awakening() function.
#' @param swarming_table The data.frame returned by the swarming() function.
#' @param maturing_table The data.frame returned by the maturing() function.
#' @param topography The raster stack returned by the data_fetcher() function.
#' @return A raster stack containing the interpolated phenological indicators (`awakening_doy`, `swarming_doy`, `maturing_doy`), restricted to pure spruce forests.
#' @examples
#' \dontrun{
#'  pheno_ind <- kpi(awakening_table, swarming_table, maturing_table, topography)
#' }
#' @export
kpi <- function(awakening_table, swarming_table, maturing_table, topography) {
    cat("===== PHENOLOGICAL INDICATORS SPATIALISATION =====\n")

    pheno_data <- Reduce(function(x, y) merge(x, y, by = c("id", "X93", "Y93"), all = TRUE), list(awakening_table, swarming_table, maturing_table))
    pheno_data <- pheno_data[, c(1:4, 6:7)]
    pheno_data <- pheno_data["awakening_doy" != 0, ]
    pheno_sf <- sf::st_as_sf(pheno_data, coords = c("X93", "Y93"), crs = 27572)
    vars <- c("awakening_doy", "swarming_doy", "maturing_doy")

    grid <- topography$spruce_forests
    terra::values(grid) <- NA
    grid_df <- as.data.frame(grid, xy = TRUE, na.rm = FALSE)
    grid_sf <- sf::st_as_sf(grid_df, coords = c("x", "y"), crs = 27572)
    pheno_ind <- terra::rast()

    if (length(pheno_data$id) < 50) {
        cat("<!> Warning - The number of DRIAS points is not enough to fit a semi-variogram and proceed to a kriging spatialization.\n")
        cat("<!> Warning - An IDW algorithm is used instead.\n")
        idw_spationer <- function(var) {
            formula <- stats::as.formula(paste(var, "~1"))

            idw_result <- gstat::idw(formula = formula, locations = pheno_sf, newdata = grid_sf, idp = 0.5, nmax = 8, debug.level = -1)
            idw_df <- cbind(sf::st_coordinates(idw_result), idw_result$var1.pred)
            colnames(idw_df) <- c("x", "y", var)

            idwed_ind <- terra::rast(idw_df, type = "xyz", crs = terra::crs(grid))
            idwed_ind <- terra::resample(idwed_ind, grid, method = "near")
            idwed_ind[is.na(topography$spruce_forests)] <- NA
            idwed_ind[idwed_ind < 0] <- 0
            terra::values(idwed_ind) <- ceiling(terra::values(idwed_ind))
            names(idwed_ind) <- var

            return(idwed_ind)
        }

        idwed_vars <- lapply(vars, idw_spationer)
        pheno_ind <- do.call(c, idwed_vars)
    }

    if (length(pheno_data$id) >= 50) {
        cat("<!> Warning - The number of DRIAS points is enough to fit a semi-variogram and proceed to a kriging spatialization.\n")
        krig_spationer <- function(var) {
            formula <- stats::as.formula(paste(var, "~1"))
            vgm_model <- gstat::variogram(formula, data = pheno_sf)
            vgm_fit <- gstat::fit.variogram(vgm_model, gstat::vgm("Exp", range = 500, nugget = NA, psill = NA), debug.level = -1)

            krig_result <- gstat::krige(formula = formula, locations = pheno_sf, newdata = grid_sf, model = vgm_fit, nmax = 8, debug.level = -1)
            krig_df <- cbind(sf::st_coordinates(krig_result), krig_result$var1.pred)
            colnames(krig_df) <- c("x", "y", var)

            kriged_ind <- terra::rast(krig_df, type = "xyz", crs = terra::crs(grid))
            kriged_ind <- terra::resample(kriged_ind, grid, method = "near")
            kriged_ind[is.na(topography$spruce_forests)] <- NA
            kriged_ind[kriged_ind < 0] <- 0
            terra::values(kriged_ind) <- ceiling(terra::values(kriged_ind))
            names(kriged_ind) <- var

            return(kriged_ind)
        }

        kriged_vars <- lapply(vars, krig_spationer)
        pheno_ind <- do.call(c, kriged_vars)
    }

    cat("== PHENOLOGICAL INDICATORS SPATIALISATION -- OK ==\n")
    cat("==================================================\n")
    gc()

    return(pheno_ind)
}

#' rpc
#'
#' This function calculates the global epidemic risk indicator (`Rpheno`), measuring the probability of successful development of a bark beetle generation based on the indicators `awakening_doy`, `swarming_doy`, and `maturing_doy`, and thus the risk of attack.
#'
#' @param pheno_ind The raster stack returned by the kpi() function.
#' @return A probability raster (`Rpheno`), ranging from 0 to 1, indicating the global epidemic risk of bark beetle generation development and thus the risk of attack.
#' @examples
#' \dontrun{
#'  rpheno <- rpc(pheno_ind)
#' }
#' @export
rpc <- function(pheno_ind) {
    cat("===== Rpheno CALCULATION =====\n")

    prob_awakening <- (-0.0028 * pheno_ind$awakening_doy) + 1.0282
    prob_swarming <- (-0.0028 * pheno_ind$swarming_doy) + 1.0311
    prob_maturing <- (-0.0032 * pheno_ind$maturing_doy) + 1.1699

    rpheno <- round(prob_awakening * prob_swarming * prob_maturing, 2)
    rpheno[rpheno < 0] <- 0
    names(rpheno) <- "Rpheno"

    cat("== Rpheno CALCULATION -- OK ==\n")
    cat("==============================\n")
    gc()

    return(rpheno)
}

#' pipeline
#'
#' Runs the entire B2SPM pipeline, from reading DRIAS data to spatializing the attack risk and number of attacks per year.
#'
#' @param drias_txt_path Path to the DRIAS file containing daily meteorological data.
#' @param dem A raster (preferably a SpatRaster in EPSG:2154) representing elevation, covering the entire study area plus a 5-pixel buffer.
#' @return A raster stack containing the spatialized phenological indicators (`awakening_doy`, `swarming_doy`, `maturing_doy`), the attack risk (`Rpheno`), and the maximum number of generations per year (`max_gen`).
#' @examples
#' \dontrun{
#'  dem <- terra::rast("dem_roi.tif")
#'  results <- pipeline("drias.txt", bbox)
#' }
#' @export
pipeline <- function(drias_txt_path, dem) {
    cat("\n")
    cat("+--------------------------------------------------------------------------------------------------+\n")
    cat("|-------------------------------- B2SPM PIPELINE INITIALISATION... --------------------------------|\n")
    cat("+--------------------------------------------------------------------------------------------------+\n")
    cat("\n")

    topography <- topo_comp(dem)
    cat("\n")

    drias_table <- drias_reader(drias_txt_path)
    cat("\n")

    drias_table <- phloem_rm(drias_table)
    cat("\n")

    awakening_table <- awakening(drias_table, topography)
    cat("\n")
    swarming_table <- swarming(drias_table, awakening_table, topography)
    cat("\n")
    maturing_table <- maturing(drias_table, swarming_table, topography)
    cat("\n")

    pheno_ind <- kpi(awakening_table, swarming_table, maturing_table, topography)
    cat("\n")

    rpheno <- rpc(pheno_ind)
    cat("\n")

    results <- c(pheno_ind, rpheno)

    cat("\n")
    cat("+--------------------------------------------------------------------------------------------------+\n")
    cat("|-------------------------------------- B2SPM PIPELINE -- OK --------------------------------------|\n")
    cat("+--------------------------------------------------------------------------------------------------+\n")
    cat("\n")

    return(results)
}