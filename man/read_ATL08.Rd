% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/read_ATL08.R
\name{read_ATL08}
\alias{read_ATL08}
\title{Read ATL08 h5 file as a data frame}
\usage{
read_ATL08(
  atl08_h5,
  beam = c("gt1l", "gt1r", "gt2l", "gt2r", "gt3l", "gt3r"),
  beam_strength = c("weak", "strong"),
  lat_range,
  odir,
  atl03_pc = FALSE
)
}
\arguments{
\item{atl08_h5}{Path to h5 file}

\item{beam}{Character vector indicating beams to process}

\item{beam_strength}{Character vector indicating the strength of beams to process}

\item{lat_range}{Numeric vector. Lower and upper latitude to return}

\item{odir}{Character. Output directory}

\item{atl03_pc}{Logical indicating if the classification of ATL03 photons should be returned}
}
\description{
Read ATL08 h5 file as a data frame
}
