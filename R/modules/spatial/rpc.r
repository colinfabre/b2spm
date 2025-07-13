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
    cat("Calculating Rpheno global epidemic risk indicator...\n")

    prob_awakening <- (-8 * 10^-5 * (spat_ind$awakening_doy)^2) + (0.0175 * spat_ind$awakening_doy) + 0.0498
    prob_swarming <- (-8 * 10^-5 * (spat_ind$swarming_doy)^2) + (0.018 * spat_ind$swarming_doy) - 0.0033
    prob_maturing <- (-7 * 10^-5 * (spat_ind$maturing_doy)^2) + (0.022 * spat_ind$maturing_doy) - 0.7852
    prob_hsi <- spat_ind$hsi

    rpheno <- terra::rast(spat_ind$awakening_doy)
    terra::values(rpheno) <- round(terra::values(prob_awakening) * terra::values(prob_swarming) * terra::values(prob_maturing) * terra::values(prob_hsi), 2)
    rpheno <- terra::clamp(rpheno, lower = 0, upper = 1, values = TRUE)
    names(rpheno) <- "Rpheno"

    cat("Calculating Rpheno -- OK\n")
    gc()

    return(rpheno)
}