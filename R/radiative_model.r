jakoby_consts <- function() {
    k1 <- 0.01
    k2 <- 0.5
    k3 <- 0.2

    return(list(k1 = k1, k2 = k2, k3 = k3))
}



phloem_calc <- function(climate_stack) {

    jakoby_consts <- jakoby_consts()
    k1 <- jakoby_consts$k1
    k2 <- jakoby_consts$k2
    k3 <- jakoby_consts$k3
    vis_rad_max <- vis_rad_calc(vis_rad = climate_stack$vis_rad)
    wind_vars <- wind_calc(climate_stack$wind)
    wind_day <- wind_vars$wind_day
    wind_night <- wind_vars$wind_night
    rh_vars <- rh_calc(climate_stack$spec_hum, climate_stack$t_mean)
    rh_day <- rh_vars$rh_day
    rh_night <- rh_vars$rh_night

    delta_t_max <- k1 * vis_rad_max - k2 * wind_day - k3 * (rh_day - 50)
    delta_t_min <- -k2 * wind_night - k3 * (rh_night - 50)

    t_max_phloem <- pmin(climate_stack$t_max + delta_t_max, climate_stack$t_max + 10)
    t_min_phloem <- pmax(climate_stack$t_min + delta_t_min, climate_stack$t_min - 2)
    t_mean_phloem <- (t_max_phloem + t_min_phloem) / 2

    return(list(t_max_phloem = t_max_phloem, t_min_phloem = t_min_phloem, t_mean_phloem = t_mean_phloem))
}