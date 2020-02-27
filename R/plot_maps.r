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
#' @param fit tagged list of outputs from TMB model via \code{Obj$report()}
#' @param TmbData optional tagged list of outputs from fit$data_list or VAST::make_data
#' @param spatial_list required for Method == Stream_network, optional for other spatial models, tagged list of outputs from \code{make_spatial_info}
#' @param Sdreport Standard deviation outputs from TMB model via \code{sdreport(Obj)}
#' @param Panel Whether to plot years for a given category (\code{Panel="Category"}) or categories for a given year ((\code{Panel="Year"})  in each panel figure
#' @param Ylim ylimits for each panel
#' @param Xlim xlimits for each panel
#' @param Zlim value scale limits
#' @param DirName Directory (absolute path)
#' @param PlotName plot names are automatically generated but option to add a modifier
#' @param category_names character vector specifying names for different categories (only used for R package \code{VAST})
#' @param covar_names character vector specifying covariate names for labeling figures
#' @param legend Boolean whether to plot legend or not
#' @param textmargin option to include y-axis text
#' @param option to add arrows between network nodes
#' @param cex size of points for map
#' @param ... arguments passed to \code{par}
#'
#' @return Mat_xt a matrix (rows: modeled knots; column: modeled year) for plotted output of last element of \code{plot_set}
#'

#' @export
plot_maps <-
function(plot_set=3, fit, Sdreport=NULL, Xlim=NULL, Ylim=NULL, Zlim = NULL, 
         TmbData=NULL, spatial_list=NULL, Panel="Category",
         DirName=NULL, PlotName=NULL,
         category_names=NULL, covar_names=NULL,
         legend=TRUE, textmargin=NULL, arrows=FALSE, cex=0.5,...){

  # local functions
  logsum = function(vec){ max(vec) + log(sum(exp(vec-max(vec)))) }

  year_labels = fit$year_labels
  years_to_plot = fit$years_to_plot
  Report <- fit$Report
  # Fill in missing inputs
  if( "D_xt" %in% names(Report)){
    # SpatialDeltaGLMM
    category_names = "singlespecies"
    Ncategories = length(category_names)
    Nyears = dim(Report$D_xt)[2]
  }
  if( "D_xct" %in% names(Report)){
    # VAST Version < 2.0.0
    if( is.null(category_names) ) category_names = 1:dim(Report$D_xct)[2]
    Ncategories = dim(Report$D_xct)[2]
    Nyears = dim(Report$D_xct)[3]
  }
  if( "D_xcy" %in% names(Report)){
    # VAST Version >= 2.0.0
    if( is.null(category_names) ) category_names = 1:dim(Report$D_xcy)[2]
    Ncategories = dim(Report$D_xcy)[2]
    Nyears = dim(Report$D_xcy)[3]
  }
  if( "D_gcy" %in% names(Report)){
    # VAST Version >= 8.0.0
    if( is.null(category_names) ) category_names = 1:dim(Report$D_gcy)[2]
    Ncategories = dim(Report$D_gcy)[2]
    Nyears = dim(Report$D_gcy)[3]
  }
  if("dhat_ktp" %in% names(Report)){
    # MIST Version <= 14
    if( is.null(category_names) ) category_names = 1:dim(Report$dhat_ktp)[3]
    Ncategories = dim(Report$dhat_ktp)[3]
    Nyears = dim(Report$dhat_ktp)[2]
  }
  if("dpred_ktp" %in% names(Report)){
    # MIST Version >= 15
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

    inp_Zlim <- Zlim
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
      # if(is.null(Zlim)) inp_Zlim <- c(0,1)
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
      # if("D_xt"%in%names(Report)) Array_xct = log(Report$D_xt+quantile(Report$D_xt,0.25))
      # if("D_xct"%in%names(Report)) Array_xct = log(Report$D_xct+quantile(Report$D_xct,0.25))
      # if("D_xcy"%in%names(Report)) Array_xct = log(Report$D_xcy+quantile(Report$D_xcy,0.25))
      # if("D_gcy"%in%names(Report)) Array_xct = log(Report$D_gcy+quantile(Report$D_gcy,0.25))
      if("D_xt"%in%names(Report)) Array_xct = Report$D_xt
      if("D_xct"%in%names(Report)) Array_xct = Report$D_xct
      if("D_xcy"%in%names(Report)) Array_xct = Report$D_xcy
      if("D_gcy"%in%names(Report)) Array_xct = Report$D_gcy

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

      # require(ggplot2)
      if(all(is.null(spatial_list))) stop("add spatial_list (output from FishStatsUtils::make_spatial_info) to use ggplot2 plots.")


    # Plot for each category
    if( tolower(Panel)=="category" ){
      if(length(dim(Array_xct))==2) Nplot = 1
      if(length(dim(Array_xct))==3) Nplot = dim(Array_xct)[2]
      for( cI in 1:Nplot){
        if(length(dim(Array_xct))==2) Return = Mat_xt = Array_xct
        if(length(dim(Array_xct))==3) Return = Mat_xt = Array_xct[,cI,]

          ## matrix is number of nodes by number of years
          if(is.null(dim(Mat_xt)))  n_t = 1 else n_t <- dim(Mat_xt)[2]
          if(n_t != length(year_labels)) stop("number of years in density array does not match Data_Geostat years")
          if(n_t > 1) {
            xct <- lapply(1:n_t, function(x){
              out <- data.frame('value'=Mat_xt[,x], 'year'=year_labels[x], spatial_list$loc_g, 'category'=category_names[cI])
              return(out)
            })
            xct <- do.call(rbind, xct)
          } else xct <- data.frame('value'=Mat_xt, 'year'=year_labels, spatial_list$loc_g, "category"=category_names[1])

          if(all(is.null(Xlim))) Xlim = c(min(xct$E_km),max(xct$E_km))
          if(all(is.null(Ylim))) Ylim = c(min(xct$N_km),max(xct$N_km))
          p <- ggplot(xct)
          if(arrows == TRUE){
            Network_sz_EN <- data.frame('parent_s'=fit$data_list$parent_s, 'child_s'=fit$data_list$child_s, fit$spatial_list$loc_g)
            l2 <- lapply(1:nrow(Network_sz_EN), function(x){
              parent <- Network_sz_EN$parent_s[x]
              find <- Network_sz_EN %>% filter(child_s == parent)
              if(nrow(find)>0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=find$E_km, 'N2'=find$N_km)
              if(nrow(find)==0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=NA, 'N2'=NA)
              # if(nrow(find)>0) out <- cbind.data.frame(Network_sz_EN[x,], 'long2'=find$long, 'lat2'=find$lat)
              # if(nrow(find)==0) out <- cbind.data.frame(Network_sz_EN[x,], 'long2'=NA, 'lat2'=NA)
              return(out)
            })
            l2 <- do.call(rbind, l2)
            p <- p + geom_segment(data=l2, aes(x = E_km,y = N_km, xend = E2, yend = N2), arrow=arrow(length=unit(0.2,"cm")), col="gray")
          }
          if(is.null(Zlim)) inp_Zlim = quantile(xct$value, prob = c(0,1), na.rm=TRUE)
          p <- p +
              geom_point(aes(x = E_km, y = N_km, color = value), cex = cex) +
              scale_color_distiller(palette = "Spectral", limits = inp_Zlim) +
              coord_cartesian(xlim = Xlim, ylim = Ylim) +
              scale_x_continuous(breaks=quantile(xct$E_km, prob=c(0.1,0.5,0.9)), labels=round(quantile(xct$E_km, prob=c(0.1,0.5,0.9)),0)) +
              # guides(color=guide_legend(title=plot_codes[plot_num])) +
              facet_wrap(~year) + 
              mytheme() +
              xlab("Eastings") + ylab("Northings")
          if(Nplot!=1) p <- p + ggtitle(paste(category_names[cI], plot_names[plot_num], PlotName))
          if(Nplot==1) p <- p + ggtitle(paste(plot_names[plot_num], PlotName))

          if(!is.null(DirName)){
            if(Nplot!=1) ggsave(file.path(DirName, paste0(plot_names[plot_num], PlotName, "_", cI, "_byCat.png")), p, width=8,height=8)
            if(Nplot==1) ggsave(file.path(DirName, paste0(plot_names[plot_num], PlotName, "_byCat.png")), p, width=8,height=8)
          }      
          if(is.null(DirName)){
            dev.new()
            print(p)
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

         ## matrix is number of nodes by number of years
          n_c <- dim(Mat_xc)[2]
          if(n_c > 1) {
            xct <- lapply(1:n_c, function(x){
              out <- data.frame('value'=Mat_xc[,x], 'year'=year_labels[years_to_plot[tI]], spatial_list$loc_g, 'category'=category_names[x])
              return(out)
            })
            xct <- do.call(rbind, xct)
          } else xct <- data.frame('value'=Mat_xc, 'year'=year_labels[years_to_plot[tI]], spatial_list$loc_g, "category"='total')
          if(all(is.null(Xlim))) Xlim = c(min(xct$E_km),max(xct$E_km))
          if(all(is.null(Ylim))) Ylim = c(min(xct$N_km),max(xct$N_km))
          p <- ggplot(xct)
          if(arrows == TRUE){
            Network_sz_EN <- data.frame('parent_s'=fit$data_list$parent_s, 'child_s'=fit$data_list$child_s, fit$spatial_list$loc_g)
            l2 <- lapply(1:nrow(Network_sz_EN), function(x){
              parent <- Network_sz_EN$parent_s[x]
              find <- Network_sz_EN %>% filter(child_s == parent)
              if(nrow(find)>0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=find$E_km, 'N2'=find$N_km)
              if(nrow(find)==0) out <- cbind.data.frame(Network_sz_EN[x,], 'E2'=NA, 'N2'=NA)
              # if(nrow(find)>0) out <- cbind.data.frame(Network_sz_EN[x,], 'long2'=find$long, 'lat2'=find$lat)
              # if(nrow(find)==0) out <- cbind.data.frame(Network_sz_EN[x,], 'long2'=NA, 'lat2'=NA)
              return(out)
            })
            l2 <- do.call(rbind, l2)
            p <- p + geom_segment(data=l2, aes(x = E_km,y = N_km, xend = E2, yend = N2), arrow=arrow(length=unit(0.2,"cm")), col="gray")
          }
          if(is.null(Zlim)) inp_Zlim = quantile(xct$value, prob = c(0,1), na.rm=TRUE)
          p <- p +
              geom_point(aes(x = E_km, y = N_km, color = value), cex=cex) +
              scale_color_distiller(palette = "Spectral", limits = inp_Zlim) +
              coord_cartesian(xlim = Xlim, ylim = Ylim) +
              scale_x_continuous(breaks=quantile(xct$E_km, prob=c(0.1,0.5,0.9)), labels=round(quantile(xct$E_km, prob=c(0.1,0.5,0.9)),0)) +
              # guides(color=guide_legend(title=plot_codes[plot_num])) +
              facet_wrap(~category) + 
              mytheme() +
              xlab("Eastings") + ylab("Northings")

          if(n_c == 1) width = 6; height = 5
          if(n_c > 1) width = 10; height = 5
          if(Nplot!=1) p <- p + ggtitle(paste(year_labels[years_to_plot[tI]], plot_names[plot_num], PlotName))
          if(Nplot==1) p <- p + ggtitle(paste(plot_names[plot_num], PlotName))

          if(!is.null(DirName)){
            if(Nplot!=1) ggsave(file.path(DirName, paste0(plot_names[plot_num], PlotName, "_", tI, "_byYear.png")), p, width=width, height=height)
            if(Nplot==1) ggsave(file.path(DirName, paste0(plot_names[plot_num], PlotName, "_byYear.png")), p, width = width, height = height)
          }       
          if(is.null(DirName)){
            dev.new()
            print(p) 
          }  
      }
    }
  }

  return( invisible(Return) )
}
