#' @title
#' Plot Pearson residuals on map
#'
#' @description
#' \code{plot_residuals} shows average Pearson residual for every knot for encounter probability and positive catch rate components
#'
#' @param fit model fit from fit_model
#' @param Data input data in format of Data_Geostat
#' @param Network_sz_LL stream network info
#' @param category_names category names
#' @param FilePath path to save figures
#' @return A tagged list of Pearson residuals
#' \describe{
#'   \item{Q1_xt}{Matrix of average residuals for encounter/non-encounter component by site \code{x} and year \code{t}}
#'   \item{Q2_xt}{Matrix of average residuals for positive-catch-rate component by site \code{x} and year \code{t}}
#' }

#' @export
plot_residuals = function( fit, Data, Network_sz_LL, category_names, FilePath ){

  ##################
  # Basic inputs
  ##################

  ## raw residuals
  pred_dens <- fit$Report$D_gcy
  obs_dens <- array(NA, dim=dim(pred_dens))

  years <- unique(Data$Year[order(Data$Year)])
  n_t <- length(years)
  n_c <- dim(pred_dens)[2]
  n_g <- dim(pred_dens)[1]
  for(t in 1:n_t){
    for(c in 1:n_c){
      sub <- Data %>% filter(Year==years[t]) %>% filter(CategoryNum == c)
      obs_dens[sub$Knot,c,t] <- sub$Catch_KG/sub$AreaSwept_km2
    }
  }

  resid_dens <- obs_dens - pred_dens
  resid_dens_list <- lapply(1:n_c, function(x){
    obs_dens_x <- data.frame(obs_dens[,x,])
    colnames(obs_dens_x) <- years
    obs_dens_x <- obs_dens_x %>% mutate('child_s'=1:nrow(obs_dens_x)) %>% gather(key = "Year", value = "Observed", -child_s)

    pred_dens_x <- data.frame(pred_dens[,x,])
    colnames(pred_dens_x) <- years
    pred_dens_x <- pred_dens_x %>% mutate('child_s'=1:nrow(pred_dens_x)) %>% gather(key = "Year", value = "Predicted", -child_s)

    op <- full_join(obs_dens_x, pred_dens_x)

    resid_dens_x <- op %>% mutate("RawResidual"=Observed - Predicted) %>% mutate("PearsonResidual" = (Observed - Predicted)/sqrt(Predicted))
    resid_dens_xll <- full_join(resid_dens_x, Network_sz_LL) %>% mutate("Category" = category_names[x])
    return(resid_dens_xll)
  })
  resid_dens <- do.call(rbind, resid_dens_list)
   
   for(c in 1:n_c){
    presid <- ggplot(resid_dens %>% filter(Category == category_names[c]) %>% filter(is.na(PearsonResidual)==FALSE)) +
            geom_point(data = Network_sz_LL, aes(x = Lon, y = Lat), color = "gray") +
            geom_point(aes(x = Lon, y = Lat, color = PearsonResidual), cex=2) +
            facet_wrap(.~Year) +
            xlab("Longitude") +  ylab("Latitude") +
            ggtitle(paste0(category_names[c], " Pearson residuals")) +
            mytheme()
    ggsave(file.path(FilePath, paste0("Pearson_residuals_", category_names[c], ".png")), presid, width = 8, height=6)
   }

  pfit <- ggplot(resid_dens %>% filter(is.na(PearsonResidual)==FALSE)) +
        geom_point(aes(x = Predicted, y = PearsonResidual)) +
        facet_wrap(.~Category, scales="free") +
        mytheme()
  ggsave(file.path(FilePath, paste0("Pearson_residuals_vs_predicted.png")), pfit)

# return(resid_dens)
}
