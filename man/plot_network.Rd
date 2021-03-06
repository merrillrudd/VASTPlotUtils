% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_network.R
\name{plot_network}
\alias{plot_network}
\title{Plot stream network and optional observations}
\usage{
plot_network(Network_sz_LL, Data = NULL, byYear = FALSE,
  byValue = FALSE, value_label = NULL, FilePath = NULL,
  FileName = "Network", arrows = FALSE, root = FALSE,
  obs_color = NULL)
}
\arguments{
\item{Data}{data frame with observations}

\item{byYear}{plot observations by year, default = FALSE}

\item{byValue}{plot observations representing size of point by their value}

\item{value_label}{label for the value observed}

\item{FilePath}{path to figure directory}

\item{FileName}{label for plot}

\item{arrows}{default = FALSE do not plot segments}

\item{root}{default TRUE to show root nodes, FALSE in case there are other root nodes that are not meaningful.}

\item{obs_color}{option for hard-coding the colors for figures, used to match colors for plotting only one category when other figures are multiple categories (optional)}
}
\value{
Figure plotting stream network and observations
}
\description{
\code{plot_network} generates figure showing network nodes, directional arrows, and observations
}
