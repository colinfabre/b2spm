library(terra)



drias_download <- function(roi, output) {
    if (!dir.exists(output)) {
        dir.create(output, recursive = TRUE)
    }

    if (!inherits(roi, "SpatVector")) {
        roi <- terra::vect(roi)
    }

    drias_urls <- c(
        "https://object.files.data.gouv.fr/meteofrance-drias/DRIAS2020_NEW/CNRM-CERFACS-CNRM-CM5/CNRM-ALADIN63/rcp85/day/prtotAdjust/prtotAdjust_France_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CNRM-ALADIN63_v2_MF-ADAMONT-SAFRAN-1980-2011_day_20060101-21001231.nc",
        "https://object.files.data.gouv.fr/meteofrance-drias/DRIAS2020_NEW/CNRM-CERFACS-CNRM-CM5/CNRM-ALADIN63/rcp85/day/tasAdjust/tasAdjust_France_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CNRM-ALADIN63_v2_MF-ADAMONT-SAFRAN-1980-2011_day_20060101-21001231.nc",
        "https://object.files.data.gouv.fr/meteofrance-drias/DRIAS2020_NEW/CNRM-CERFACS-CNRM-CM5/CNRM-ALADIN63/rcp85/day/tasmaxAdjust/tasmaxAdjust_France_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CNRM-ALADIN63_v2_MF-ADAMONT-SAFRAN-1980-2011_day_20060101-21001231.nc",
        "https://object.files.data.gouv.fr/meteofrance-drias/DRIAS2020_NEW/CNRM-CERFACS-CNRM-CM5/CNRM-ALADIN63/rcp85/day/tasminAdjust/tasminAdjust_France_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CNRM-ALADIN63_v2_MF-ADAMONT-SAFRAN-1980-2011_day_20060101-21001231.nc",
        "https://object.files.data.gouv.fr/meteofrance-drias/DRIAS2020_NEW/CNRM-CERFACS-CNRM-CM5/CNRM-ALADIN63/rcp85/day/hussAdjust/hussAdjust_France_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CNRM-ALADIN63_v2_MF-ADAMONT-SAFRAN-1980-2011_day_20060101-21001231.nc",
        "https://object.files.data.gouv.fr/meteofrance-drias/DRIAS2020_NEW/CNRM-CERFACS-CNRM-CM5/CNRM-ALADIN63/rcp85/day/sfcWindAdjust/sfcWindAdjust_France_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CNRM-ALADIN63_v2_MF-ADAMONT-SAFRAN-1980-2011_day_20060101-21001231.nc",
        "https://object.files.data.gouv.fr/meteofrance-drias/DRIAS2020_NEW/CNRM-CERFACS-CNRM-CM5/CNRM-ALADIN63/rcp85/day/rldsAdjust/rldsAdjust_France_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CNRM-ALADIN63_v2_MF-ADAMONT-SAFRAN-1980-2011_day_20060101-21001231.nc",
        "https://object.files.data.gouv.fr/meteofrance-drias/DRIAS2020_NEW/CNRM-CERFACS-CNRM-CM5/CNRM-ALADIN63/rcp85/day/rsdsAdjust/rsdsAdjust_France_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CNRM-ALADIN63_v2_MF-ADAMONT-SAFRAN-1980-2011_day_20060101-21001231.nc",
        "https://object.files.data.gouv.fr/meteofrance-drias/DRIAS2020_NEW/CNRM-CERFACS-CNRM-CM5/CNRM-ALADIN63/rcp85/day/evspsblpotAdjust/evspsblpotAdjust_France_CNRM-CERFACS-CNRM-CM5_rcp85_r1i1p1_CNRM-ALADIN63_v2_MF-ADAMONT-SAFRAN-1980-2011_day_20060101-21001231_Hg0175.nc"
    )

    params <- c("drias_urls", "output", "process_nc")
    iterover <- seq_len(drias_urls)
    iterwith <- function(i) {
        process_nc(drias_urls[i], roi, output)
    }
    parallelism(params, iterover, iterwith)
}