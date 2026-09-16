# Appliances Energy Prediction Data

Energy consumption and environmental measurements collected at 10-minute
intervals from a low-energy residential building. The data include
appliance energy use together with indoor and outdoor temperature,
humidity, and weather measurements.

## Usage

``` r
energydata
```

## Format

A data frame with 19,735 observations and 29 variables:

- date:

  Date and time of the observation.

- Appliances:

  Appliance energy use, in Wh.

- lights:

  Energy use of light fixtures, in Wh.

- T1:

  Temperature in the kitchen area.

- RH_1:

  Relative humidity in the kitchen area.

- T2:

  Temperature in the living room area.

- RH_2:

  Relative humidity in the living room area.

- T3:

  Temperature in the laundry room area.

- RH_3:

  Relative humidity in the laundry room area.

- T4:

  Temperature in the office room.

- RH_4:

  Relative humidity in the office room.

- T5:

  Temperature in the bathroom.

- RH_5:

  Relative humidity in the bathroom.

- T6:

  Outdoor temperature measured on the north side of the building.

- RH_6:

  Outdoor relative humidity measured on the north side of the building.

- T7:

  Temperature in the ironing room.

- RH_7:

  Relative humidity in the ironing room.

- T8:

  Temperature in teenager room 2.

- RH_8:

  Relative humidity in teenager room 2.

- T9:

  Temperature in the parents' room.

- RH_9:

  Relative humidity in the parents' room.

- T_out:

  Outdoor temperature from the weather station.

- Press_mm_hg:

  Atmospheric pressure in mm Hg.

- RH_out:

  Outdoor relative humidity.

- Windspeed:

  Wind speed.

- Visibility:

  Visibility.

- Tdewpoint:

  Dew-point temperature.

- rv1:

  Random variable 1.

- rv2:

  Random variable 2.

## Source

UCI Machine Learning Repository, Appliances Energy Prediction Dataset:
<https://archive.ics.uci.edu/dataset/374/appliances+energy+prediction>

## References

Candanedo, L. M., Feldheim, V., and Deramaix, D. (2017). Data driven
prediction models of energy use of appliances in a low-energy house.
*Energy and Buildings*, 140, 81–97.
