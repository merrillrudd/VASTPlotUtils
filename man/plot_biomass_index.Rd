% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_biomass_index.R
\name{plot_biomass_index}
\alias{plot_biomass_index}
\title{{Plot index of abundance}}
\usage{
plot_biomass_index(TmbData, Sdreport, DirName = paste0(getwd(), "/"),
  PlotName = "Index", interval_width = 1, strata_names = NULL,
  category_names = NULL, use_biascorr = TRUE, plot_legend = TRUE,
  total_area_km2 = NULL, plot_log = FALSE, width = 4, height = 4,
  create_covariance_table = FALSE, ...)
}
\arguments{
\item{TmbData}{Formatted data inputs, from `VAST::make_data(...)`}

\item{Sdreport}{output from fit$parameter_estimates$SD}

\item{DirName}{Directory for saving plot and table}

\item{PlotName}{Name for plot}

\item{interval_width}{width for confidence intervals}

\item{strata_names}{names for spatial strata}

\item{category_names}{names for categories or a label if no categories are used}

\item{use_biascorr}{Boolean, whether to use bias-corrected estimates if available}

\item{plot_legend}{Add legend for labelling colors}

\item{total_area_km2}{Total area for calculating a design-based estimator using one design-stratum (only recommended for model exploration)}

\item{plot_log}{Boolean, whether to plot y-axis in log-scale}

\item{width}{plot width in inches}

\item{height}{plot height in inches}

\item{...}{Other inputs to `par()`}
}
\value{
Return Tagged list of output
\describe{
  \item{Table}{table of index estimates by stratum and year, e.g., for including in an assessment model}
}
}
\description{
{ plots an index proportion to population abundance}
}
