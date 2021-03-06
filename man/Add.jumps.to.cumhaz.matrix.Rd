% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/average_cumulative_hazard.R
\name{Add.jumps.to.cumhaz.matrix}
\alias{Add.jumps.to.cumhaz.matrix}
\title{Function to extend length of a Kaplan-Meier curve}
\usage{
Add.jumps.to.cumhaz.matrix(cumhaz.matrix, maxevents)
}
\arguments{
\item{cumhaz.matrix}{a data.frame with information on Kaplan-Meier (Nelson-Aalen/Breslow) values at each event-time point, along point-wise confidence interval (CI) values and treatment/control indicator. By convention such information must be arranged in columns with name 'km' (survival probability), 'time' (event-time), 'km_l' (lower CI value for 'km'), 'km_u' (upper CI value for 'km'), 'trt' (binary indicator for a treatment/control grouping, can be NA).}

\item{maxevents}{number of (group-specific) event-times (points) used to extend curve length.}
}
\description{
Function to extend length of a Kaplan-Meier curve
}
\details{
The last available KM value is stretched as a constant between the last and future event-time. The event-time line is "evenly" stretched between the last.time.point + mean(diff(time)) and last.time.point + (mean(diff(time))*number.of.point.to.add). This function can be useful when having several simulated curves that have different length due to randomness. Then, the function normalizes all simulated curves to a common 'maxevents' length (e.g. the max length). This can be later used to average across curves, for example.
Basic information for 'cumhaz.matrix' can be, for instance, computed via 'survfit' from package 'survival'.
}
