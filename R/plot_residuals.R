#' @title
#' Plot Pearson residuals on map
#'
#' @description
#' \code{plot_residuals} shows average Pearson residual for every knot for encounter probability and positive catch rate components
#'
#' @param ObsModel
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
plot_residuals = function( ObsModel, fit, Data, Network_sz_LL, category_names, FilePath = NULL ){

  ##################
  # Basic inputs
  ##################
  Report <- fit$Report
  Data_Geostat <- Data

## Pearson resids for detection and catch rate
  D_i <- Report$R1_i*Report$R2_i
  PR1_i <- PR2_i <- rep(NA, length(D_i))
  for(i in 1:length(D_i)){
    
    ## bernoulli for presence
      mui <- Report$R1_i[i]
      obs <- as.numeric(Data_Geostat$Catch_KG[i]>0)
      PR1_i[i] <- (obs-mui)/sqrt(mui*(1-mui)/1)
    
    ## NA for 0 observations
    obs <- Data_Geostat$Catch_KG[i]
    if(obs>0){
      ## make sure to use the right variance as this depends on gear type
      # gr <- as.numeric(Data_Geostat$Gear[i])
      if(ObsModel == 1) PR2_i[i] <- (log(obs)-log(Report$R2_i[i])+Report$SigmaM[1]^2/2)/Report$SigmaM[1]
      if(ObsModel %in% c(7,11)) PR2_i[i] <- (log(obs)-log(Report$R2_i[i]))/log(Report$R2_i[i])
    }
  }
  df <- cbind(Data_Geostat, PR1=PR1_i, PR2=PR2_i, positive=ifelse(Data_Geostat$Catch_KG>0,1,0))
  xlim <- range(df$Lon); ylim <- range(df$Lat)
  # tmp <- 1:3
  # if(results$model!='combined') tmp <- 1
  # for(gr in tmp){
    # gt <- levels(Data_Geostat$Gear)[1]

    # g <- ggplot(subset(df, positive==1), aes(Lon, Lat, size=abs(PR2), color=PR2>0))+
    #   geom_point(alpha=.25) + facet_wrap('Year') + xlim(xlim) + ylim(ylim)+
    #   scale_size('Pearson Resid', range=c(0,3))  + theme_bw()
  # g <- plot_network(Network_sz_LL = Network_sz_LL, Data = df, plot_type = 2, byYear=TRUE, byValue=TRUE)
  p1 <- ggplot(df %>% filter(positive == 1)) +
    geom_point(data = Network_sz_LL, aes(x = Lon, y = Lat), col = "gray", alpha = 0.6) +
    geom_point(aes(x = Lon, y = Lat, col = PR1>0, size = abs(PR1)), alpha = 0.7) +
    facet_wrap(Year~Category) +
    xlim(xlim) + ylim(ylim) +
    scale_color_brewer(palette = "Set1") +
    scale_size("Pearson Residual", range = c(0,max(df$PR1, na.rm = TRUE))) +
    xlab("Longitude") + ylab("Latitude") +
    ggtitle("First component") + 
    theme_bw(base_size = 14)
  if(is.null(FilePath)==FALSE){
    ggsave(file.path(FilePath, "Pearson_resid_firstcomponent.png"), p1, width = 12, height = 10)
  } else {
    print(p1)
  }

  p2 <- ggplot(df %>% filter(positive == 1)) +
    geom_point(data = Network_sz_LL, aes(x = Lon, y = Lat), col = "gray", alpha = 0.6) +
    geom_point(aes(x = Lon, y = Lat, col = PR2>0, size = abs(PR2)), alpha = 0.7) +
    facet_wrap(Year~Category) +
    xlim(xlim) + ylim(ylim) +
    scale_color_brewer(palette = "Set1") +
    scale_size("Pearson Residual", range = c(0,max(df$PR2, na.rm = TRUE))) +
    xlab("Longitude") + ylab("Latitude") +
    ggtitle("Second component") + 
    theme_bw(base_size = 14)
  if(is.null(FilePath)==FALSE){
    ggsave(file.path(FilePath, "Pearson_resid_secondcomponent.png"), p2, width = 12, height = 10)
  } else {
    print(p2)
  }

  # if(is.null(FilePath)==FALSE) ggsave(file.path(FilePath, 'Pearson_resid_catchrate.png'), plot=g,
  #          width=7, height=5)
  # if(is.null(FilePath)){
  #   dev.new()
  #   print(g)
  # }
  #   # g <- ggplot(subset(df), aes(Lon, Lat, size=abs(PR1), color=PR1>0))+
  #   #   geom_point(alpha=.25) + facet_wrap('Year') + xlim(xlim) + ylim(ylim)+
  #   #   scale_size('Pearson Resid', range=c(0,3))  + theme_bw()
  # g <- plot_network(Network_sz_LL = Network_sz_LL, Data=df, plot_type=1, byYear=TRUE, byValue = TRUE)
  # if(is.null(FilePath)==FALSE)  ggsave(file.path(FilePath, 'Pearson_resid_encounter.png'), plot=g,
  #          width=7, height=5)
  # if(is.null(FilePath)){
  #   dev.new()
  #   print(g)
  # }

sub <- subset(df, positive==1)
sresid <- sum(sub$PR2, na.rm=TRUE )
sresid2 <- sum(sub$PR1, na.rm=TRUE)
sdf <- data.frame("P1" = sum(sub$PR1, na.rm=TRUE), "P2"=sum(sub$PR2, na.rm=TRUE))
if(is.null(FilePath)==FALSE) write.csv(sdf, file=file.path(FilePath, "Sum_resid.csv"))


# if(type == "density"){
#   pred_dens <- fit$Report$D_gcy  
#   ## density
#   obs_dens <- array(NA, dim=dim(pred_dens))

#   years <- unique(Data$Year[order(Data$Year)])
#   n_t <- length(years)
#   n_c <- dim(pred_dens)[2]
#   n_g <- dim(pred_dens)[1]
#   for(t in 1:n_t){
#     for(c in 1:n_c){
#       sub <- Data %>% filter(Year==years[t]) %>% filter(CategoryNum == c)
#       obs_dens[sub$Knot,c,t] <- sub$Catch_KG/sub$AreaSwept_km2
#     }
#   }
# }


#   resid_dens <- obs_dens - pred_dens
#   resid_dens_list <- lapply(1:n_c, function(x){
#     obs_dens_x <- data.frame(obs_dens[,x,])
#     colnames(obs_dens_x) <- years
#     obs_dens_x <- obs_dens_x %>% mutate('child_s'=1:nrow(obs_dens_x)) %>% gather(key = "Year", value = "Observed", -child_s)

#     pred_dens_x <- data.frame(pred_dens[,x,])
#     colnames(pred_dens_x) <- years
#     pred_dens_x <- pred_dens_x %>% mutate('child_s'=1:nrow(pred_dens_x)) %>% gather(key = "Year", value = "Predicted", -child_s)

#     op <- full_join(obs_dens_x, pred_dens_x)

#     resid_dens_x <- op 
#     resid <- unlist(sapply(1:nrow(resid_dens_x), function(y){
#       obs <- resid_dens_x$Observed[y]
#       if(is.na(obs)==FALSE){
#         if(obs > 0){
#           out <- (log(resid_dens_x$Observed[y]) - resid_dens_x$Predicted[y])/sqrt(resid_dens_x$Predicted[y])
#         } else{ out <- NA}
#       } else{ out <- NA}
#     })) 
#     resid_dens_x$PearsonResidual <- resid
#     resid_dens_xll <- full_join(resid_dens_x, Network_sz_LL) %>% mutate("Category" = category_names[x])
#     return(resid_dens_xll)
#   })
#   resid_dens <- do.call(rbind, resid_dens_list)
   
#    for(c in 1:n_c){
#     presid <- ggplot(resid_dens %>% filter(Category == category_names[c]) %>% filter(is.na(PearsonResidual)==FALSE)) +
#             geom_point(data = Network_sz_LL, aes(x = Lon, y = Lat), color = "gray") +
#             geom_point(aes(x = Lon, y = Lat, color = PearsonResidual), cex=2) +
#             facet_wrap(.~Year) +
#             xlab("Longitude") +  ylab("Latitude") +
#             ggtitle(paste0(category_names[c], " Pearson residuals")) +
#             mytheme()
#     ggsave(file.path(FilePath, paste0("Pearson_residuals_", category_names[c], ".png")), presid, width = 8, height=6)
#    }

#   pfit <- ggplot(resid_dens %>% filter(is.na(PearsonResidual)==FALSE)) +
#         geom_point(aes(x = Predicted, y = PearsonResidual)) +
#         facet_wrap(.~Category, scales="free") +
#         mytheme()
#   ggsave(file.path(FilePath, paste0("Pearson_residuals_vs_predicted.png")), pfit)

# return(resid_dens)
}
