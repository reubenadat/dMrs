% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dMrs.R
\name{refData_match}
\alias{refData_match}
\title{refData_match}
\usage{
refData_match(wDAT, rDAT)
}
\arguments{
\item{wDAT}{A working dataset data.frame containing required 
columns \code{age}, \code{time}, \code{delta}, 
\code{datediag_yr}, \code{dateEvent_yr}, and \code{sex} 
corresponding to age, observed time, event status, 
diagnosis year, observed year, and sex (coded 'female'/'male'), 
respectively.}

\item{rDAT}{A reference dataset data.frame containing
required columns \code{Year}, \code{Age}, \code{qx},
and \code{sex} corresponding to reference year, 
reference age, event hazard, and sex, respectively.}
}
\description{
This function takes as input
a working dataset of interest and a 
reference dataset.
}
