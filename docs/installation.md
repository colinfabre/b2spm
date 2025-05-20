---
layout: default
title: Installation & Dependencies
permalink: /installation/
---

# Installation

B2SPM is compatible with Windows, Mac, and Linux.
You can install B2SPM from GitHub using:

```r
# install.packages("devtools")
devtools::install_github("colinfabre/b2spm")
```

# Dependencies

## Internal Functionality

B2SPM relies on several other packages to function; make sure to download them beforehand to ensure proper operation:

```r
install.packages("utils", "stats", "terra", "sf", "gstat", "ps")
```

There are no specific package versions required.

## External Data

B2SPM requires a DEM (Digital Elevation Model) covering at least the entire study area where the modeling will be conducted. The resolution of this DEM should not exceed 25m (the final analysis grid is 250m, and a resolution that is too fine (e.g., 5m) significantly increases computation times without improving model accuracy).

The phenological model relies on climate data from the [DRIAS portal](https://www.drias-climat.fr), which provides regionalized climate projections for mainland France. The corrected DRIAS-2020 data are available for several different models, but the most accurate for France is MétéoFrance's global ARPEGE-Climat model, forced by the regional climate model (RCM) and used in the ALADIN63 simulation. The preferred scenario is RCP8.5, as it aligns with the current greenhouse gas emission curve.

<p align="center">
  <img src="/figures/drias_model.png" width="1000" height="450">
</p>

The user can select the desired study year and extract all DRIAS points contained within their study area, adding a one-row and one-column margin on all sides to ensure proper spatialization of phenological indicators.

<p align="center">
  <img src="/figures/drias_time.png" width="1000" height="285">
  <img src="/figures/drias_points.png" width="1000" height="855">
</p>

The climate parameters required to run B2SPM are:

- Daily minimum temperature (K)
- Daily maximum temperature (K)
- Daily average temperature (K)
- Total daily precipitation (kg/m²/s)
- Specific humidity (kg/kg)
- Incident visible solar radiation (W/m²)
- Incident infrared solar radiation (W/m²)
- Wind speed (m/s)
- Potential evapotranspiration (kg/m²/s)

<p align="center">
  <img src="/figures/drias_params.png" width="1000" height="580">
</p>

Data downloaded from DRIAS-2020 must be provided as `.txt` files (native format) following this structure:

- Columns:
  - POSITION: DRIAS point identifier
  - DATE: measurement date (format DD/MM/YYYY)
  - Climate parameters (one column per selected parameter)

- POSITION Format:
  - Point identifier
  - Lambert93 coordinates (LambertX, LambertY)

- Separator: fields must be comma-separated (,)

Example of an expected DRIAS file structure:

```r
point_id, X_LambII, Y_LambII,date,tmin,tmax,tmoy,tot_pr,spec_hum,vis_solrad,ir_solrad,wind,pet
ID001,908000,2129000,01/01/2030,-2.3,5.7,1.5,0.0,0.0025,150.0,220.0,2.1,0.0
ID001,908000,2129000,02/01/2030,-1.8,6.2,2.1,0.0001,0.0024,152.5,225.0,1.8,0.00005
```

<p align="center">
  <img src="/figures/drias_format.png" width="1000" height="725">
</p>