#' @title{Plot standard maps}
#'
#' @description {plots a standard set of diagnostic maps}
#'
#' @param plot_set integer-vector defining plots to create
#' \describe{
#'   \item{plot_set=1}{Probability of encounter/non-encounter}
#'   \item{plot_set=2}{Log-expected positive catch rate}
#'   \item{plot_set=3}{Log-predicted density (product of encounter probability and positive catch rates)}
#'   \item{plot_set=4}{Log-positive catch rates (rescaled)}
#'   \item{plot_set=5}{Log-predicted density (rescaled)}
#'   \item{plot_set=6}{Spatio-temporal variation in encounter probability}
#'   \item{plot_set=7}{Spatio-temporal variation in log-positive catch rates}
#'   \item{plot_set=8}{Linear predictor for encounter probability}
#'   \item{plot_set=9}{Linear predictor for positive catch rates}
#'   \item{plot_set=10}{Coefficient of variation for predicted density (available only if \code{Data_Fn(...,Options=c('SD_site_logdensity'=1,...))}}
#'   \item{plot_set=11}{Covariates that are included in the model}
#'   \item{plot_set=12}{Total biomass across all categories (only useful in a multivariate model)}
#'   \item{plot_set=13}{Covariate effects on encounter probability}
#'   \item{plot_set=14}{Covariate effects on positive catch rates}
#'   \item{plot_set=15}{Individual covariate effects on encounter probability}
#'   \item{plot_set=16}{Individual covariate effects on positive catch rates}
#' }
#' @param Report tagged list of outputs from TMB model via \code{Obj$report()}
#' @param TmbData optional tagged list of outputs from fit$data_list or VAST::make_data
#' @param spatial_list required for Method == Stream_network, optional for other spatial models, tagged list of outputs from \code{make_spatial_info}
#' @param Sdreport Standard deviation outputs from TMB model via \code{sdreport(Obj)}
#' @param Panel Whether to plot years for a given category (\code{Panel="Category"}) or categories for a given year ((\code{Panel="Year"})  in each panel figure
#' @param Ylim ylimits for each panel
#' @param Xlim xlimits for each panel
#' @param FileName Directory (absolute path) and base for filenames of plots
#' @param year_labels Year names for labeling panels
#' @param years_to_plot integer vector, specifying positions of \code{year_labelss} for plotting (used to avoid plotting years with no data, etc.)
#' @param category_names character vector specifying names for different categories (only used for R package \code{VAST})
#' @param covar_names character vector specifying covariate names for labeling figures
#' @param legend Boolean whether to plot legend or not
#' @param format image format
#' @param ... arguments passed to \code{par}
#'
#' @return Mat_xt a matrix (rows: modeled knots; column: modeled year) for plotted output of last element of \code{plot_set}
#'

#' @export
plot_maps <-
function(plot_set=3, Report, Sdreport=NULL, Xlim, Ylim,
         TmbData=NULL, spatial_list=NULL, Panel="Category",
         FileName=paste0(getwd(),"/"), year_labels=NULL, years_to_plot=NULL, Format="png",
         category_names=NULL, covar_names=NULL,
         legend=TRUE, ...){

  # local functions
  logsum = function(vec){ max(vec) + log(sum(exp(vec-max(vec)))) }

  # Fill in missing inputs
  if( "D_xt" %in% names(Report)){
    # SpatialDeltaGLMM
    if( is.null(year_labels) ) year_labels = 1:ncol(Report$D_xt)
    if( is.null(years_to_plot) ) years_to_plot = 1:ncol(Report$D_xt)
    category_names = "singlespecies"
    Ncategories = length(category_names)
    Nyears = dim(Report$D_xt)[2]
  }
  if( "D_xct" %in% names(Report)){
    # VAST Version < 2.0.0
    if( is.null(year_labels) ) year_labels = 1:dim(Report$D_xct)[3]
    if( is.null(years_to_plot) ) years_to_plot = 1:dim(Report$D_xct)[3]
    if( is.null(category_names) ) category_names = 1:dim(Report$D_xct)[2]
    Ncategories = dim(Report$D_xct)[2]
    Nyears = dim(Report$D_xct)[3]
  }
  if( "D_xcy" %in% names(Report)){
    # VAST Version >= 2.0.0
    if( is.null(year_labels) ) year_labels = 1:dim(Report$D_xcy)[3]
    if( is.null(years_to_plot) ) years_to_plot = 1:dim(Report$D_xcy)[3]
    if( is.null(category_names) ) category_names = 1:dim(Report$D_xcy)[2]
    Ncategories = dim(Report$D_xcy)[2]
    Nyears = dim(Report$D_xcy)[3]
  }
  if( "D_gcy" %in% names(Report)){
    # VAST Version >= 8.0.0
    if( is.null(year_labels) ) year_labels = 1:dim(Report$D_gcy)[3]
    if( is.null(years_to_plot) ) years_to_plot = 1:dim(Report$D_gcy)[3]
    if( is.null(category_names) ) category_names = 1:dim(Report$D_gcy)[2]
    Ncategories = dim(Report$D_gcy)[2]
    Nyears = dim(Report$D_gcy)[3]
  }
  if("dhat_ktp" %in% names(Report)){
    # MIST Version <= 14
    if( is.null(year_labels) ) year_labels = 1:dim(Report$dhat_ktp)[2]
    if( is.null(years_to_plot) ) years_to_plot = 1:dim(Report$dhat_ktp)[2]
    if( is.null(category_names) ) category_names = 1:dim(Report$dhat_ktp)[3]
    Ncategories = dim(Report$dhat_ktp)[3]
    Nyears = dim(Report$dhat_ktp)[2]
  }
  if("dpred_ktp" %in% names(Report)){
    # MIST Version >= 15
    if( is.null(year_labels) ) year_labels = 1:dim(Report$dpred_ktp)[2]
    if( is.null(years_to_plot) ) years_to_plot = 1:dim(Report$dpred_ktp)[2]
    if( is.null(category_names) ) category_names = 1:dim(Report$dpred_ktp)[3]
    Ncategories = dim(Report$dpred_ktp)[3]
    Nyears = dim(Report$dpred_ktp)[2]
  }

  # Errors
  if( Nyears != length(year_labels) ){
    stop("Problem with `year_labels`")
  }
  if( Ncategories != length(category_names) ){
    stop("Problem with `category_names`")
  }

  # Extract elements
  plot_codes <- c("Pres", "Pos", "Dens", "Pos_Rescaled", "Dens_Rescaled", "Eps_Pres", "Eps_Pos", "LinPred_Pres", "LinPred_Pos", "Dens_CV", "Covariates", "Total_dens", "Cov_effects_Pres", "Cov_effects_Pos", "Indiv_cov_effects_Pres", "Indiv_cov_effects_Pos")
  plot_names <- c("presence_absence", "log-positive catch rates", "log-density", "positive catch rates", "density", "epsilon for presence_absence", "epsilon for positive catch rates", "encounter probability linear predictor", "positive catch rates linear predictor", "density CV", "covariates", "total density", "encounter probability covariate effects", "catch rates covariate effects", "encounter probability individual covariate effects", "catch rates individual covariate effects")
  if( is.null(textmargin)){
    textmargin <- c("Probability of encounter", "Density, ln(kg. per square km.)", "Density, ln(kg. per square km.)", "", "", "", "", "", "", "CV of density (dimensionless)", "Covariate value", "Density, ln(kg. per square km.)", "", "", "", "")
  }
  # Select locations to plot
  # if( Nknots<Inf ){
  #   NN_plot = stats::kmeans(x=PlotDF[,c("Lon","Lat")], centers=Nknots, iter.max=50, nstart=2, trace=0)
  #   Match = match( 1:Nknots, NN_plot$cluster)
  #   PlotDF = PlotDF[Match,]
  #   message( "Restricted plotting locations to ", Nknots, " locations" )
  # }

  # Loop through plots
  for(plot_num in plot_set){

    # Extract matrix to plot
    if(plot_num==1){
      # Presence/absence ("Pres")
      if("D_xt"%in%names(Report)) Array_xct = Report$R1_xt
      if("D_xct"%in%names(Report)) Array_xct = Report$R1_xct
      if("D_xcy"%in%names(Report)) Array_xct = Report$R1_xcy
      if("D_gcy"%in%names(Report)) Array_xct = Report$R1_gcy
      # if("D_gcy"%in%names(Report)) stop("`plot_maps` not implemented for requested plot_num")
      if(any(c("dhat_ktp","dpred_ktp")%in%names(Report))) stop("Not implemented for SpatialVAM")
      message( "plot_num=1 doesn't work well when using ObsModel[2]==1" )
    }
    if(plot_num==2){
      # Positive values ("Pos")
      if("D_xt"%in%names(Report)) Array_xct = log(Report$R2_xt)
      if("D_xct"%in%names(Report)) Array_xct = log(Report$R2_xct)
      if("D_xcy"%in%names(Report)) Array_xct = log(Report$R2_xcy)
      if("D_gcy"%in%names(Report)) Array_xct = log(Report$R2_gcy)
      # if("D_gcy"%in%names(Report)) stop("`plot_maps` not implemented for requested plot_num")
      if(any(c("dhat_ktp","dpred_ktp")%in%names(Report)))  stop("Not implemented for SpatialVAM")
    }
    if(plot_num==3){
      # Density ("Dens")
      if("D_xt"%in%names(Report)) Array_xct = log(Report$D_xt)
      if("D_xct"%in%names(Report)) Array_xct = log(Report$D_xct)
      if("D_xcy"%in%names(Report)) Array_xct = log(Report$D_xcy)
      if("D_gcy"%in%names(Report)) Array_xct = log(Report$D_gcy)
      if("dhat_ktp" %in% names(Report)) Array_xct = aperm(Report$dhat_ktp[,,cI],c(1,3,2))
      if("dpred_ktp" %in% names(Report)) Array_xct = aperm(Report$dpred_ktp[,,cI],c(1,3,2))
    }
    if(plot_num==4){
      # Positive values rescaled ("Pos_Rescaled")
      if("D_xt"%in%names(Report)) Array_xct = log(Report$R2_xt+quantile(Report$R2_xt,0.25))
      if("D_xct"%in%names(Report)) Array_xct = log(Report$R2_xct+quantile(Report$R2_xct,0.25))
      if("D_xcy"%in%names(Report)) Array_xct = log(Report$R2_xcy+quantile(Report$R2_xcy,0.25))
      if("D_gcy"%in%names(Report)) Array_xct = log(Report$R2_gcy+quantile(Report$R2_gcy,0.25))
      # if("D_gcy"%in%names(Report)) stop("`plot_maps` not implemented for requested plot_num")
      if(any(c("dhat_ktp","dpred_ktp")%in%names(Report)))  stop("Not implemented for SpatialVAM")
    }
    if(plot_num==5){
      # Density rescaled ("Dens_Rescaled")
      if("D_xt"%in%names(Report)) Array_xct = log(Report$D_xt+quantile(Report$D_xt,0.25))
      if("D_xct"%in%names(Report)) Array_xct = log(Report$D_xct+quantile(Report$D_xct,0.25))
      if("D_xcy"%in%names(Report)) Array_xct = log(Report$D_xcy+quantile(Report$D_xcy,0.25))
      if("D_gcy"%in%names(Report)) Array_xct = log(Report$D_gcy+quantile(Report$D_gcy,0.25))
      # if("D_gcy"%in%names(Report)) stop("`plot_maps` not implemented for requested plot_num")
      if(any(c("dhat_ktp","dpred_ktp")%in%names(Report)))  stop("Not implemented for SpatialVAM")
    }
    if(plot_num==6){
      # Epsilon for presence/absence ("Eps_Pres")
      if("D_xt"%in%names(Report)) Array_xct = Report$Epsilon1_st
      if("D_xct"%in%names(Report)) Array_xct = Report$Epsilon1_sct
      if("D_xcy"%in%names(Report)) Array_xct = Report$Epsilon1_sct
      if("D_gcy"%in%names(Report)) Array_xct = Report$Epsilon1_gct
      if(any(c("dhat_ktp","dpred_ktp")%in%names(Report)))  stop("Not implemented for SpatialVAM")
    }
    if(plot_num==7){
      # Epsilon for positive values ("Eps_Pos")
      if("D_xt"%in%names(Report)) Array_xct = Report$Epsilon2_st
      if("D_xct"%in%names(Report)) Array_xct = Report$Epsilon2_sct
      if("D_xcy"%in%names(Report)) Array_xct = Report$Epsilon2_sct
      if("D_gcy"%in%names(Report)) Array_xct = Report$Epsilon2_gct
      if(any(c("dhat_ktp","dpred_ktp")%in%names(Report)))  stop("Not implemented for SpatialVAM")
    }
    if(plot_num==8){
      # Linear predictor for probability of encounter
      if("D_xt"%in%names(Report)) Array_xct = Report$P1_xt
      if("D_xct"%in%names(Report)) Array_xct = Report$P1_xct
      if("D_xcy"%in%names(Report)) Array_xct = Report$P1_xcy
      if("D_gcy"%in%names(Report)) stop("`plot_maps` not implemented for requested plot_num")
      if(any(c("dhat_ktp","dpred_ktp")%in%names(Report)))  stop("Not implemented for SpatialVAM")
    }
    if(plot_num==9){
      # Linear predictor for positive catch rates
      if("D_xt"%in%names(Report)) Array_xct = Report$P2_xt
      if("D_xct"%in%names(Report)) Array_xct = Report$P2_xct
      if("D_xcy"%in%names(Report)) Array_xct = Report$P2_xcy
      if("D_gcy"%in%names(Report)) stop("`plot_maps` not implemented for requested plot_num")
      if(any(c("dhat_ktp","dpred_ktp")%in%names(Report)))  stop("Not implemented for SpatialVAM")
    }
    if(plot_num==10){
      # Density ("Dens") CV             # Index_xtl
      if( is.null(Sdreport) ) stop("Must supply 'Sdreport' if 'plot_num=10'")
      if("D_xt"%in%names(Report)){
        if( !("log(Index_xtl)" %in% rownames(TMB::summary.sdreport(Sdreport))) ) stop("Please re-run with Options('SD_site_logdensity'=1,...) to use 'plot_num=10' in 'SpatialDeltaGLMM'")
        Array_xct = array( TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))=="log(Index_xtl)"),], dim=c(dim(Report$D_xt),ncol(Report$Index_tl),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) )[,,1,'Std. Error']
      }
      if("D_xct"%in%names(Report)){
        if( !("log(Index_xctl)" %in% rownames(TMB::summary.sdreport(Sdreport))) ) stop("Please re-run with Options('SD_site_logdensity'=1,...) to use 'plot_num=10' in 'VAST'")
        Array_xct = array( TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))=="log(Index_xctl)"),], dim=c(dim(Report$D_xct),dim(Report$Index_ctl)[3],2), dimnames=list(NULL,NULL,NULL,NULL,c('Estimate','Std. Error')) )[,,,1,'Std. Error']
      }
      if("D_xcy"%in%names(Report)){
        if( !("log(Index_xcyl)" %in% rownames(TMB::summary.sdreport(Sdreport))) ) stop("Please re-run with Options('SD_site_logdensity'=1,...) to use 'plot_num=10' in 'VAST'")
        Array_xct = array( TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))=="log(Index_xcyl)"),], dim=c(dim(Report$D_xcy),dim(Report$Index_cyl)[3],2), dimnames=list(NULL,NULL,NULL,NULL,c('Estimate','Std. Error')) )[,,,1,'Std. Error']
      }
      if("D_gcy"%in%names(Report)){
        if( !("log(Index_gcyl)" %in% rownames(TMB::summary.sdreport(Sdreport))) ) stop("Please re-run with Options('SD_site_logdensity'=1,...) to use 'plot_num=10' in 'VAST'")
        Array_xct = array( TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))=="log(Index_gcyl)"),], dim=c(dim(Report$D_gcy),dim(Report$Index_cyl)[3],2), dimnames=list(NULL,NULL,NULL,NULL,c('Estimate','Std. Error')) )[,,,1,'Std. Error']
      }
      if(any(c("dhat_ktp","dpred_ktp")%in%names(Report))) stop("'plot_num=10' not implemented for 'SpatialVAM'")
      # Convert to CV
      Array_xct = sqrt( exp(Array_xct^2) - 1 )
      if("D_gcy"%in%names(Report)) stop("`plot_maps` not implemented for requested plot_num")
    }
    if(plot_num==11){
      if(is.null(TmbData)) stop( "Must provide `TmbData` to plot covariates" )
      if(!("X_gtp" %in% names(TmbData))) stop( "Can only plot covariates for VAST version >= 2.0.0" )
      Array_xct = aperm( TmbData$X_gtp, perm=c(1,3,2) )
      category_names = 1:dim(Array_xct)[2]
    }
    if(plot_num==12){
      # Total density ("Dens")
      if("D_xt"%in%names(Report)) Array_xct = log(Report$D_xt)
      if("D_xct"%in%names(Report)) Array_xct = log(apply(Report$D_xct,FUN=sum,MARGIN=c(1,3)))
      if("D_xcy"%in%names(Report)) Array_xct = log(apply(Report$D_xcy,FUN=sum,MARGIN=c(1,3)))
      if("D_gcy"%in%names(Report)) Array_xct = log(apply(Report$D_gcy,FUN=sum,MARGIN=c(1,3)))
      if("dhat_ktp" %in% names(Report)) Array_xct = apply(aperm(Report$dhat_ktp,c(1,3,2)),FUN=logsum,MARGIN=c(1,3))
      if("dpred_ktp" %in% names(Report)) Array_xct = apply(aperm(Report$dpred_ktp,c(1,3,2)),FUN=logsum,MARGIN=c(1,3))
    }
    if(plot_num==13){
      # Covariate effects for probability of encounter
      if("D_xt"%in%names(Report)) stop()
      if("D_xct"%in%names(Report)) stop()
      if("D_xcy"%in%names(Report)) Array_xct = Report$eta1_xct
      if("D_gcy"%in%names(Report)) Array_xct = Report$eta1_gct
      if("dhat_ktp" %in% names(Report)) stop()
      if("dpred_ktp" %in% names(Report)) stop()
    }
    if(plot_num==14){
      # Covariate effects for positive catch rates
      if("D_xt"%in%names(Report)) stop()
      if("D_xct"%in%names(Report)) stop()
      if("D_xcy"%in%names(Report)) Array_xct = Report$eta2_xct
      if("D_gcy"%in%names(Report)) Array_xct = Report$eta2_gct
      if("dhat_ktp" %in% names(Report)) stop()
      if("dpred_ktp" %in% names(Report)) stop()
    }
    if(plot_num==15){
      if(is.null(TmbData)) stop("Must provide `TmbData` to plot covariates.")
      if(!("X_gtp" %in% names(TmbData))) stop( "Can only plot covariates for VAST version >= 2.0.0")

        X <- TmbData$X_gtp
        # gamma1_ctp <- summary(Sdreport)[which(grepl("gamma1", rownames(summary(Sdreport)))),"Estimate"]
        gamma1_ctp <- Sdreport$par.fixed[grepl("gamma1_ctp", names(Sdreport$par.fixed))]

        n_p <- dim(X)[3]
        category_names <- covar_names
        Array_xct <- array(NA, dim=c(dim(X)[1], n_p, dim(X)[2]))
        for(i in 1:n_p){
          eta <- gamma1_ctp[i] * X[,,i]
          Array_xct[,i,] <- eta
        }
    }
    if(plot_num==16){
      if(is.null(TmbData)) stop("Must provide `TmbData` to plot covariates.")
      if(!("X_gtp" %in% names(TmbData))) stop( "Can only plot covariates for VAST version >= 2.0.0")

        X <- TmbData$X_gtp
        # gamma1_ctp <- summary(Sdreport)[which(grepl("gamma1", rownames(summary(Sdreport)))),"Estimate"]
        gamma2_ctp <- Sdreport$par.fixed[grepl('gamma2_ctp', names(Sdreport$par.fixed))]

        n_p <- dim(X)[3]
        category_names <- covar_names
        Array_xct <- array(NA, dim=c(dim(X)[1], n_p, dim(X)[2]))
        for(i in 1:n_p){
          eta <- gamma2_ctp[i] * X[,,i]
          Array_xct[,i,] <- eta
        }
    }

          ## indicates Method == "Stream_network"
    if(all(is.null(TmbData$parent_s))==FALSE){
      require(ggplot2)
      if(all(is.null(spatial_list))) stop("add spatial_list (output from FishStatsUtils::make_spatial_info) to use ggplot2 plots and/or plot stream network.")

            mytheme <- function (base_size = 14, base_family = "") 
            {
                theme_grey(base_size = base_size, base_family = base_family) %+replace%
                theme(axis.title.x = element_text(margin = margin(10,0,0,0)),
                      #axis.title.x = element_text(vjust = -1.5),
                      #axis.title.y = element_text(margin = margin(0,20,0,0)),
                      #axis.title.y = element_text(vjust = -0.1),
                      axis.text = element_text(size = rel(0.8)),
                      axis.ticks = element_line(colour = "black"), 
                      legend.key = element_rect(colour = "grey80"),
                      panel.background = element_rect(fill = "white", colour = NA),
                      panel.border = element_rect(fill = NA, colour = "grey50"),
                      panel.grid.major = element_line(colour = "grey90", size = 0.2),
                      panel.grid.minor = element_line(colour = "grey98", size = 0.5),
                      strip.background = element_rect(fill = "grey80", colour = "grey50", size = 0.2))
            }
    }


    # Plot for each category
    if( tolower(Panel)=="category" ){
      if(length(dim(Array_xct))==2) Nplot = 1
      if(length(dim(Array_xct))==3) Nplot = dim(Array_xct)[2]
      for( cI in 1:Nplot){
        if(length(dim(Array_xct))==2) Return = Mat_xt = Array_xct
        if(length(dim(Array_xct))==3) Return = Mat_xt = Array_xct[,cI,]

        if(all(is.null(TmbData$parent_s))){
          # Do plot
          if( is.null(mfrow)) mfrow = c(ceiling(sqrt(length(years_to_plot))), ceiling(length(years_to_plot)/ceiling(sqrt(length(years_to_plot)))))
          if(add==FALSE) par( mfrow=mfrow )
          PlotMap_Fn( MappingDetails=MappingDetails, Mat=Mat_xt[,years_to_plot,drop=FALSE], PlotDF=PlotDF, MapSizeRatio=MapSizeRatio, Xlim=Xlim, Ylim=Ylim, FileName=paste0(FileName,plot_codes[plot_num],ifelse(Nplot>1,paste0("--",category_names[cI]),"")), year_labels=year_labels[years_to_plot], Rescale=Rescale, Rotate=Rotate, Format=Format, Res=Res, zone=zone, Cex=Cex, textmargin=textmargin[plot_num], add=add, pch=pch, Legend=Legend, mfrow=mfrow, plot_legend_fig=plot_legend_fig, ...)
        } else {
          ## matrix is number of nodes by number of years
          if(is.null(dim(Mat_xt)))  n_t = 1 else n_t <- dim(Mat_xt)[2]
          if(n_t != length(year_labels)) stop("number of years in density array does not match Data_Geostat years")
          if(n_t > 1) {
            xct <- lapply(1:n_t, function(x){
              out <- data.frame('value'=Mat_xt[,x], 'year'=year_labels[x], spatial_list$loc_x, 'category'=category_names[cI])
              return(out)
            })
            xct <- do.call(rbind, xct)
          } else xct <- data.frame('value'=Mat_xt, 'year'=year_labels, spatial_list$loc_x, "category"=category_names[1])

          p <- ggplot(xct) +
              geom_point(aes(x = E_km, y = N_km, color = value), cex=Cex) +#, ...) +
              scale_color_distiller(palette = "Spectral") +
              scale_x_continuous(breaks=quantile(xct$E_km, prob=c(0.1,0.5,0.9)), labels=round(quantile(xct$E_km, prob=c(0.1,0.5,0.9)),0)) +
              # guides(color=guide_legend(title=plot_codes[plot_num])) +
              facet_wrap(~year) + 
              mytheme() +
              xlab("Eastings") + ylab("Northings")
          if(Nplot!=1) p <- p + ggtitle(paste(category_names[cI], plot_names[plot_num]))
          if(Nplot==1) p <- p + ggtitle(paste(plot_names[plot_num]))

          if(!is.null(FileName)){
            if(Nplot!=1) ggsave(file.path(FileName, paste0(plot_names[plot_num], "_", cI, "_byCat.png")), p, width=8,height=8)
            if(Nplot==1) ggsave(file.path(FileName, paste0(plot_names[plot_num], "_byCat.png")), p, width=8,height=8)
          }          
        }
      }
    }
    # Plot for each year
    if( tolower(Panel)=="year" ){
      Nplot = length(years_to_plot)
      for( tI in 1:Nplot){
        if(length(dim(Array_xct))==2) Mat_xc = Array_xct[,years_to_plot[tI],drop=TRUE]
        if(length(dim(Array_xct))==3) Mat_xc = Array_xct[,,years_to_plot[tI],drop=TRUE]
        if(is.null(dim(Mat_xc)) & is.vector(Mat_xc)){
          Ncategories = 1
        } else { Ncategories = dim(Mat_xc)[2] }
        Return = Mat_xc = array( as.vector(Mat_xc), dim=c(dim(Array_xct)[1],Ncategories)) # Reformat to make sure it has same format for everything

        if(all(is.null(TmbData$parent_s))){
          # Do plot
          if( is.null(mfrow)) mfrow = c(ceiling(sqrt(length(category_names))), ceiling(length(category_names)/ceiling(sqrt(length(category_names)))))
          if(add==FALSE) par( mfrow=mfrow )
          PlotMap_Fn( MappingDetails=MappingDetails, Mat=Mat_xc, PlotDF=PlotDF, MapSizeRatio=MapSizeRatio, Xlim=Xlim, Ylim=Ylim, FileName=paste0(FileName,plot_codes[plot_num],ifelse(Nplot>1,paste0("--",year_labels[years_to_plot][tI]),"")), year_labels=category_names, Rescale=Rescale, Rotate=Rotate, Format=Format, Res=Res, zone=zone, Cex=Cex, textmargin=textmargin[plot_num], add=add, pch=pch, Legend=Legend, mfrow=mfrow, plot_legend_fig=plot_legend_fig, ...)
        } else { 
         ## matrix is number of nodes by number of years
          n_c <- dim(Mat_xc)[2]
          if(n_c > 1) {
            xct <- lapply(1:n_c, function(x){
              out <- data.frame('value'=Mat_xc[,x], 'year'=years_to_plot[tI], spatial_list$loc_x, 'category'=category_names[x])
              return(out)
            })
            xct <- do.call(rbind, xct)
          } else xct <- data.frame('value'=Mat_xc, 'year'=years_to_plot[tI], spatial_list$loc_x, "category"='total')

          p <- ggplot(xct) +
              geom_point(aes(x = E_km, y = N_km, color = value), cex=2) +
              scale_color_distiller(palette = "Spectral") +
              scale_x_continuous(breaks=quantile(xct$E_km, prob=c(0.1,0.5,0.9)), labels=round(quantile(xct$E_km, prob=c(0.1,0.5,0.9)),0)) +
              # guides(color=guide_legend(title=plot_codes[plot_num])) +
              facet_wrap(~category) + 
              mytheme() +
              xlab("Eastings") + ylab("Northings")

          if(n_c == 1) width = 6; height = 5
          if(n_c > 1) width = 10; height = 5
          if(Nplot!=1) p <- p + ggtitle(paste(year_labels[tI], plot_names[plot_num]))
          if(Nplot==1) p <- p + ggtitle(paste(plot_names[plot_num]))

          if(!is.null(FileName)){
            if(Nplot!=1) ggsave(file.path(FileName, paste0(plot_names[plot_num], "_", tI, "_byYear.png")), p, width=width, height=height)
            if(Nplot==1) ggsave(file.path(FileName, paste0(plot_names[plot_num], "_byYear.png")), p, width = width, height = height)
          }          
        }
      }
    }
  }

  return( invisible(Return) )
}