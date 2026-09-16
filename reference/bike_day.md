# Bike Sharing Daily Data

Daily bike-sharing counts together with weather and calendar information
from the Capital Bikeshare system. The data contain daily observations
for 2011 and 2012.

## Usage

``` r
bike_day
```

## Format

A data frame with 731 observations and 16 variables:

- instant:

  Record index.

- dteday:

  Calendar date.

- season:

  Season indicator.

- yr:

  Year indicator, with 0 corresponding to 2011 and 1 to 2012.

- mnth:

  Month, from 1 to 12.

- holiday:

  Indicator for whether the day is a holiday.

- weekday:

  Day-of-week indicator.

- workingday:

  Indicator for whether the day is a working day.

- weathersit:

  Weather situation category.

- temp:

  Normalized temperature.

- atemp:

  Normalized apparent temperature.

- hum:

  Normalized humidity.

- windspeed:

  Normalized wind speed.

- casual:

  Number of casual bike users.

- registered:

  Number of registered bike users.

- cnt:

  Total number of bike rentals.

## Source

UCI Machine Learning Repository, Bike Sharing Dataset:
<https://archive.ics.uci.edu/dataset/275/bike+sharing+dataset>

## References

Fanaee-T, H. and Gama, J. (2013). Event labeling combining ensemble
detectors and background knowledge.
