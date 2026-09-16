\name{bike_day}
\alias{bike_day}
\docType{data}

\title{Bike Sharing Daily Data}

\description{
Daily bike-sharing counts together with weather and calendar information
from the Capital Bikeshare system. The data contain daily observations
for 2011 and 2012.
}

\usage{
bike_day
}

\format{
A data frame with 731 observations and 16 variables:
\describe{
  \item{instant}{Record index.}
  \item{dteday}{Calendar date.}
  \item{season}{Season indicator.}
  \item{yr}{Year indicator, with 0 corresponding to 2011 and 1 to 2012.}
  \item{mnth}{Month, from 1 to 12.}
  \item{holiday}{Indicator for whether the day is a holiday.}
  \item{weekday}{Day-of-week indicator.}
  \item{workingday}{Indicator for whether the day is a working day.}
  \item{weathersit}{Weather situation category.}
  \item{temp}{Normalized temperature.}
  \item{atemp}{Normalized apparent temperature.}
  \item{hum}{Normalized humidity.}
  \item{windspeed}{Normalized wind speed.}
  \item{casual}{Number of casual bike users.}
  \item{registered}{Number of registered bike users.}
  \item{cnt}{Total number of bike rentals.}
}
}

\source{
UCI Machine Learning Repository, Bike Sharing Dataset:
\url{https://archive.ics.uci.edu/dataset/275/bike+sharing+dataset}
}

\references{
Fanaee-T, H. and Gama, J. (2013).
Event labeling combining ensemble detectors and background knowledge.
}

\keyword{datasets}
