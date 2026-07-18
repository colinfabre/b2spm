rh_calc <- function(spec_hum, t_mean) {
    rh <- (spec_hum / (0.622 + 0.378 * spec_hum)) * exp(17.27 * t_mean / (t_mean + 237.3))

    rh_day <- rh * 0.8
    rh_night <- rh * 1.2

    return(rh, rh_day, rh_night)
}



wind_calc <- function(wind) {
    wind_day <- wind * 1.2
    wind_night <- wind * 0.8

    return(wind_day, wind_night)
}



vis_rad_calc <- function(vis_rad) {
    vis_rad_max <- vis_rad * 1.5

    return(vis_rad_max)
}