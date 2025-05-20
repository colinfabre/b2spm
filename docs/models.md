---
layout: default
title: Models
permalink: /models/
---

## Radiative Model for Sub-Phloem Temperature

The bark beetle hibernates and then digs its egg-laying galleries where larvae mature under the bark of spruce trees. Calculating sub-phloem temperatures is therefore crucial to understanding and accurately modeling the insect's phenological cycle. This temperature is influenced by several factors related to forest structure and environmental conditions:

- Solar radiation (visible and infrared): the bark absorbs part of the radiation, locally increasing the temperature. However, in a dense spruce forest, shading reduces this effect.

- Wind: stronger winds promote heat dissipation from the trunk, reducing heating.

- Specific humidity and precipitation: high humidity and regular precipitation increase the thermal inertia of the bark, maintaining a more stable temperature (buffer effect).

- Diurnal thermal gradient: the thermal inertia of the phloem reduces daily temperature amplitude, with limited elevation of $$T_{max_{phloem}}$$ and moderate cooling of $$T_{min_{phloem}}$$.

A constrained non-linear radiative model accounting for the impact of these environmental variables has been directly integrated into B2SPM to model the evolution of sub-phloem temperature, which dictates bark beetle phenology. It is structured as follows:

### 1. Heat Conduction Equation for Spruce Bark and Phloem

$$\frac{\partial t}{\partial T_{phloem}} = \alpha \left( \frac{\partial^2 T_{phloem}}{\partial x^2} \right) + S(x,t)$$

where:
- $$T_{phloem}(x,t)$$ temperature at depth $$x$$ in the bark/phloem at time $$t$$
- $$\alpha \approx 1.2 \times 10^{-7}{m^2/s}$$ thermal diffusivity of wood $$\frac{k}{\rho C_p}$$, as a function of conductivity $$k \approx 0.15{W/m\cdot K}$$, density $$\rho \approx 600{kg/m^3}$$ (for decaying wood) and heat capacity $$C_p \approx 1600{J/kg \cdot K}$$
- $$S(x,t)$$ source term representing solar radiation absorption in the phloem

Assuming a steady-state equilibrium $$(\frac{\partial T}{\partial t} = 0)$$, and an exponential attenuation of absorbed heat in the bark, we obtain:

$$T_{phloem} = T_{ext} + \beta_{vis.solrad} \times e^{-\lambda x} \times vis.solrad + \beta_{ir.solrad} \times e^{-\lambda x} \times {ir.solrad}$$

where:
- $$\lambda \approx 0.2{mm^{-1}}$$ radiation attenuation coefficient in the bark (related to bark thickness and cellular structure)

### 2. Modeling Cooling by Convection and Evapotranspiration

Newton's Law of Cooling: $$\frac{dT_{phloem}}{dt} = -h \times (T_{phloem} - T_{ext})$$

where:
- $$h$$ convection coefficient (dependent on wind and humidity)

At steady-state:

$$T_{phloem} = T_{ext} + \frac{S}{h}$$

where:
- $$S$$ absorbed heat

But we know that $$h$$ depends on wind and humidity:

$$h = h_{0} \times \left(1 + \gamma_{wind} \times e^{\delta \times {wind}} \right) \times \left(1 + \gamma_{spec.hum} \times {spec.hum}^b \right)$$

where:
- $$h_0 \approx 10\text{W/m}^2 \cdot \text{K}$$ base convection coefficient
- $$\gamma_{wind} \approx 0.1$$ wind adjustment coefficient
- $$\gamma_{spec.hum} \approx 0.3$$ humidity adjustment coefficient
- $$b \approx 1.3$$ humidity adjustment exponent

### 3. Final Radiative Model

Combining (1) and (2), we obtain:

$$T_{phloem} = T_{ext} + \frac{\beta_{vis.solrad} \times e^{-\lambda x} \times {vis.solrad} + \beta_{ir.solrad} \times e^{-\lambda x} \times {ir.solrad}}{h_0 \times \left(1 + \gamma_{wind} \times e^{\delta \times {wind}}\right) \times \left(1 + \gamma_{spec.hum} \times {spec.hum}^b\right)}$$

with the physical constraint:

$$0 \leq T_{phloem} - T_{ext} \leq 2.5^\circ C$$

## Calculation of the Photoperiod Moderated by Cast Shadows

Swarming and laying of bark-beetles are controlled by atmospherical conditions (such as temperature and precipitation), and the cessation of these activities is mainly dictated by two conditions:

- a minimum threshold photoperiod (10h)
- an average daily air temperature threshold (15°C)

When the photoperiod threshold, and to a lesser extent the temperature threshold, are reached, the bark beetle's phenological activity decreases and it prepares for hibernation. Diapause is finally triggered after 5 days of maintaining these conditions. This marks the end of bark beetle attacks on spruce trees for that year.

The photoperiod depends on the day of the year ($$doy$$) and the latitude. By simplifying the Fourier series representing the path of the sun, we can thus deduce the exposure time of a given point. But the surrounding topography induces cast shadows on this point, seeing its theorical solar exposure being reduced. Therefore, the theorical photoperiod is corrected by integrating topographic shading over daytime hours using the input DEM.

### 1. Theorical Solar Exposure

This method uses an astronomical solar geometry model (simplified Fourier's Series) to estimate the theorical solar exposure at a given point over daytime hours.

$$\gamma = \frac{2\pi}{365}(doy - 1 + \frac{h}{24})$$

where:
- $$\gamma$$ day angle (ie fraction fo the year in radians in the Fourier's Series, modulating the trigonometrical components of the declination as a function of the Earth's position on its orbit)
- $$doy$$ day of the year (between 1 and 366)
- $$h = 12$$

$$\delta = 0.006918 - 0.399912\cos(\gamma) + 0.070257\sin(\gamma) - 0.006758\cos(2\gamma) + 0.000907\sin(2\gamma) - 0.002697\cos(3\gamma) + 0.00148\sin(3\gamma)$$

where:
- $$\delta$$ solar declination angle

$$pp(doy) = \frac{24}{\pi} \times \arccos(-\tan(\phi) \times \tan(\delta))$$

where:
- $$\phi$$ station latitude (radians)
- $$pp(doy)$$ theorical photoperiod (hours) for the given doy

### 2. Cast Shadowing

$$\gamma = \frac{2\pi}{365}(doy - 1 + \frac{h}{24})$$

where:
- $$\gamma$$ day angle (ie fraction fo the year in radians in the Fourier's Series, modulating the trigonometrical components of the declination as a function of the Earth's position on its orbit)
- $$doy$$ day of the year (between 1 and 366)
- $$h$$ calculation hour in 06.00 (sunrise at 45-47°N), 09.00, 12.00 (solar apex), 15.00 and 18.00 (sunset at 45-47°N)

$$\delta = 0.006918 - 0.399912\cos(\gamma) + 0.070257\sin(\gamma) - 0.006758\cos(2\gamma) + 0.000907\sin(2\gamma) - 0.002697\cos(3\gamma) + 0.00148\sin(3\gamma)$$

where:
- $$\delta$$ solar declination angle (radians)

$$angle(h) = \pi \times \frac{h - 12}{12}$$

where:
- $$angle(h)$$ solar angle from North (radians) as a function of the calculation hour $$h$$

$$\alpha_{s} = \arcsin[\sin(\phi) \times \sin(\delta) + \cos(\phi) \times \cos(\delta) \times \cos(angle_{hour})]$$

where:
- $$\alpha_{s}$$ solar altitude angle (radians); only positive solar altitude angles are considered for illumination.

$$\theta_{s} = [atan2(-cos(\delta) \times sin(angle_{hour}), sin(\delta) \times cos(\phi) - cos(\delta) \times sin(\phi) \times cos(angle_{hour})) \times \frac{180}{\pi}] \bmod 360$$

where:
- $$\theta_{s}$$ solar azimuth angle (radians) within [0; 360]

For each calculation hour $$h$$ and its corresponding solar altitude angle $$\alpha_{s}$$ and azimuth angle $$\theta_{s}$$, the hillshade is computed using `terra::shade()` and slope and aspect derived from the input DEM.

$$cshd(doy) = \frac{1}{5} \sum_{h=9}^{18} \mathbb{I}_{\text{sunlit}}(\alpha_s(h), \theta_s(h))$$

where:
- $$cshd(doy)$$ clear-sky shading coefficient between 0 (complete shadowing) and 1 (complete clearing)

### 3. Real Photoperiod

$$pp_{cshd}(doy) = pp(doy) \times cshd(doy)$$

where:
- $$pp_{cshd}(doy)$$ effective photoperiod (in hours) for the given doy, modulated by the surrounding topography

As this whole process can be very computationally expensive (around 20 hours for the 808 DRIAS points of the alpine arc for one year), there are 3 precision levels:
- `precision = 1` means the photoperiod is computed as a constant through a 15-days window.
- `precision = 2` means the cast shadows are daily computed at 12.00 (minimum shadowing).
- `precision = 3` means the cast shadows are daily computed throughout the course of the sun (at 06.00, 09.00, 12.00, 15.00 and 18.00); in that case, expect the computation time to be multiplied by 2.5 from `precision = 2`.