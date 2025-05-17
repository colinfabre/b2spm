---
layout: default
title: Main Functions
permalink: /functions/
---

## Main Functions

`topo_comp(dem)`: Extracts spruce forest areas and topographic data from the study area.

`drias_reader(drias_txt_path)` or `drias_fetcher(topography, year)`: Reads downloaded DRIAS data or fetches DRIAS database corresponding to the specified year and study area.

**`pipeline(topography, drias_table)`: Runs the complete B2SPM pipeline on the study area from the radiative model computation to the risk index calculation.**