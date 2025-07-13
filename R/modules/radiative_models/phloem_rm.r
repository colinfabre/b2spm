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
    cat("Calculating under-phloem temperatures using a calibrated spruce-wood radiative model...\n")

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

    cat("Calculating under-phloem temperatures -- OK\n")
    gc()

    return(drias_table)
}