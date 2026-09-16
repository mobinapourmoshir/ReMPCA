\name{energydata}
\alias{energydata}
\docType{data}

\title{Appliances Energy Prediction Data}

\description{
Energy consumption and environmental measurements collected at
10-minute intervals from a low-energy residential building.
The data include appliance energy use together with indoor and
outdoor temperature, humidity, and weather measurements.
}

\usage{
energydata
}

\format{
A data frame with 19,735 observations and 29 variables:
\describe{
  \item{date}{Date and time of the observation.}
  \item{Appliances}{Appliance energy use, in Wh.}
  \item{lights}{Energy use of light fixtures, in Wh.}
  \item{T1}{Temperature in the kitchen area.}
  \item{RH_1}{Relative humidity in the kitchen area.}
  \item{T2}{Temperature in the living room area.}
  \item{RH_2}{Relative humidity in the living room area.}
  \item{T3}{Temperature in the laundry room area.}
  \item{RH_3}{Relative humidity in the laundry room area.}
  \item{T4}{Temperature in the office room.}
  \item{RH_4}{Relative humidity in the office room.}
  \item{T5}{Temperature in the bathroom.}
  \item{RH_5}{Relative humidity in the bathroom.}
  \item{T6}{Outdoor temperature measured on the north side of the building.}
  \item{RH_6}{Outdoor relative humidity measured on the north side of the building.}
  \item{T7}{Temperature in the ironing room.}
  \item{RH_7}{Relative humidity in the ironing room.}
  \item{T8}{Temperature in teenager room 2.}
  \item{RH_8}{Relative humidity in teenager room 2.}
  \item{T9}{Temperature in the parents' room.}
  \item{RH_9}{Relative humidity in the parents' room.}
  \item{T_out}{Outdoor temperature from the weather station.}
  \item{Press_mm_hg}{Atmospheric pressure in mm Hg.}
  \item{RH_out}{Outdoor relative humidity.}
  \item{Windspeed}{Wind speed.}
  \item{Visibility}{Visibility.}
  \item{Tdewpoint}{Dew-point temperature.}
  \item{rv1}{Random variable 1.}
  \item{rv2}{Random variable 2.}
}
}

\source{
UCI Machine Learning Repository, Appliances Energy Prediction Dataset:
\url{https://archive.ics.uci.edu/dataset/374/appliances+energy+prediction}
}

\references{
Candanedo, L. M., Feldheim, V., and Deramaix, D. (2017).
Data driven prediction models of energy use of appliances in a
low-energy house.
\emph{Energy and Buildings}, 140, 81--97.
}

\keyword{datasets}
