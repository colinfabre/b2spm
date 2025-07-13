#' pipeline
#'
#' Runs the entire B2SPM pipeline, from computing topography cast shadowing to spatializing the attack risk.
#'
#' @param topography The raster stack returned by the topo_comp() function.
#' @param drias_table The DRIAS table processed either by the drias_reader() or the drias_fetcher() function.
#' @param return_tables Whether the intermediate phenological tables (`awakening_table`, `swarming_table` and `maturing_table`) should be returned as ouputs of the pipeline or not.
#' @return A list of the raster stack containing the spatialized phenological indicators (`awakening_doy`, `swarming_doy`, `maturing_doy`), the attack risk (`Rpheno`), and the maximum number of generations per year (`max_gen`), and the intermediate phenological tables if `return_tables = TRUE`.
#' @examples
#' \dontrun{
#'  results <- pipeline(topography, drias_table)
#' }
#' @export
pipeline <- function(topography, drias_table, return_tables = FALSE) {
    cat("\n")
    cat("+--------------------------------------------------------------------------------------------------+\n")
    cat("|----------------------------------------- B2SPM PIPELINE -----------------------------------------|\n")
    cat("+--------------------------------------------------------------------------------------------------+\n")
    cat("\n")
    gc()

    Sys.sleep(2)
    cat("Pipeline initialisation...\n")
    cat("\n")
    start_time <- Sys.time()
    Sys.sleep(2)

    cat("Checking available RAM...\n")
    if (round(ps::ps_system_memory()$total / 1024^3, 0) < 4) {
        cat("<!> Warning - The required 4GB of RAM is not available for pipeline to launch.\n")
        clear_mem <- readline("This will clear the current R environment, excepted 'drias_table' and 'topography' objects. Do you accept? (y/n): ")
        while (tolower(clear_mem) != "y" && tolower(clear_mem) != "n") {
            cat("Incorrect input, you must answer by yes (y) or no (n). Try again.\n")
            clear_mem <- readline("(y/n): ")
        }

        if (clear_mem == "y") {
            cat("The current R environment will be cleared.\n")
            conf <- readline("Are you sure? (y/n): ")
            if (tolower(conf) == "y") {
                for (i in range(5)) {
                    Sys.sleep(1)
                    cat(paste0("Cleaning of the environment in ", 6 - i, "\n"))
                }

                env <- ls(envir = .GlobalEnv)
                env <- setdiff(env, c("drias_table", "topography"))
                rm(env)
                gc()
                cat("The current R environment has been correctly cleared.\n")

                if (round(ps::ps_system_memory()$total / 1024^3, 0) < 4) {
                    stop("!! ERROR - The required 4GB of RAM is still not available for pipeline to launch, even after a complete R environment clearing. The code will abort.")
                } else {
                cat("The required 4GB of RAM is now available for pipeline to launch.\n")
                }
            }
        } else {
            stop("The pipeline couldn't be launch due to insufficient RAM. Please clear some space first.\n")
        }
    } else {
        cat("Available RAM -- OK\n")
        cat("\n")
    }

    Sys.sleep(2)

    cat("Checking data compatibility...\n")
    if (!inherits(topography, "SpatRaster") && any(names(topography) != c("spruce_forests", "alt", "slope", "aspect", "tpi"))) {
        stop("!! ERROR - 'topography' is invalid. It must necessarily be the result of the topo_comp() function.")
    }
    if (!inherits(drias_table, "data.frame") && any(names(drias_table) != c("id", "X93", "Y93", "date", "doy", "tmin", "tmax", "tmean", "tot_pr", "spec_hum", "vis_solrad", "ir_solrad", "wind", "pet"))) {
        stop("!! ERROR - 'drias_table' is invalid. It must necessarily be the result of the drias_reader() or drias_fetcher() function.")
    }
    if (!is.logical(return_tables)) {
        stop("!! ERROR - 'return_tables' must be TRUE or FALSE.")
    }
    if (terra::ext(topography$alt) < terra::ext(terra::vect(drias_table, geom = c("X93", "Y93"), crs = terra::crs(topography))) + 250) {
        stop("!! ERROR - The input DEM in the `topo_comp()` function is smaller than the ROI + the 250m required buffer.\n")
    } else {
        cat("Data compatibility -- OK\n")
        cat("\n")
    }

    Sys.sleep(2)
    cat("Launching the pipeline...\n")
    cat("\n")
    Sys.sleep(2)

    drias_table <- ppc(drias_table, topography, return_tables)
    cat("\n")

    drias_table <- phloem_rm(drias_table)
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
    cat("Pipeline -- OK")
    process_time <- round(end_time - start_time, 0) - (5 * 2 + 5) ### number of Sys.sleep(2)
    cat(paste0("Processing time: ", process_time, "mins\n"))

    cat("\n")
    cat("+--------------------------------------------------------------------------------------------------+\n")
    cat("|----------------------------------------- B2SPM PIPELINE -----------------------------------------|\n")
    cat("+--------------------------------------------------------------------------------------------------+\n")
    cat("\n")
    gc()

    return(results)
}