
#' @title
#' Plot shifts in distribution and area occupied
#'
#' @description
#' \code{plot_range_index} plots center-of-gravity, kernel-area occupied, and effective-area occupied
#'
#' @param Sdreport TMB output from `TMB::sdreport(Obj)`
#' @param Report Reporting output from `Obj$report()`
#' @param TmbData Formatted data inputs, from `VAST::Data_Fn(...)`
#' @param Year_Set Year names for plotting
#' @param PlotDir Directory for plots
#' @param FileName_COG Full filename (including directory) for center-of-gravity plot
#' @param FileName_Area Full filename (including directory) for center-of-gravity plot
#' @param FileName_EffArea Full filename (including directory) for center-of-gravity plot
#' @param Znames Names for center-of-gravity summary statistics
#' @inheritParams plot_biomass_index
#'
#' @return Return Tagged list of output
#' \describe{
#'   \item{COG_Table}{table of center-of-gravity ("COG") estimates by year}
#'   \item{KernelArea_Table}{table of kernel-area approximation to area occupied, estimated by year}
#'   \item{EffectiveArea_Table}{table of effective-area approximation to area occupied, estimated by year, recommended over \code{KernelArea_Table}}
#' }
#'

#' @export
plot_range_index = function( Sdreport, Report, TmbData, Year_Set=NULL, Years2Include=NULL, strata_names=NULL,
  PlotDir=paste0(getwd(),"/"), FileName_COG=paste0(PlotDir,"/center_of_gravity.png"),
  FileName_Area=paste0(PlotDir,"/Area.png"), FileName_EffArea=paste0(PlotDir,"/Effective_Area.png"), Znames=rep("",ncol(TmbData$Z_xm)), use_biascorr=TRUE,
  category_names=NULL, interval_width=1, ...){

  # Informative errors
  if(is.null(Sdreport)) stop("Sdreport is NULL; please provide Sdreport")

  # Which parameters
  if( "ln_Index_tl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # SpatialDeltaGLMM
    CogName = "mean_Z_tm"
    EffectiveName = "effective_area_tl"
    TmbData[['n_c']] = 1
  }
  if( "ln_Index_ctl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # VAST Version < 2.0.0
    CogName = "mean_Z_ctm"
    EffectiveName = "effective_area_ctl"
  }
  if( "ln_Index_cyl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # VAST Version >= 2.0.0
    CogName = "mean_Z_cym"
    EffectiveName = "effective_area_cyl"
    TmbData[["n_t"]] = nrow(TmbData[["t_yz"]])
  }

  # Default inputs
  if( is.null(Year_Set)) Year_Set = 1:TmbData$n_t
  if( is.null(Years2Include) ) Years2Include = 1:TmbData$n_t
  if( is.null(strata_names) ) strata_names = 1:TmbData$n_l
  if( is.null(category_names) ) category_names = 1:TmbData$n_c
  Return = list( "Year_Set"=Year_Set )

  # Plot distribution shift and kernal-area approximation to area occupied if necessary outputs are available
  if( !any(c("mean_Z_tm","mean_Z_ctm","mean_Z_cym") %in% names(Report)) ){
    message( "To plot range-shifts and kernal-approximation to area occupied, please re-run with Options['Calculate_Range']=1" )
  }else{
    message( "Plotting center-of-gravity..." )

    # Extract index  (using bias-correctino if available and requested)
    SD = TMB::summary.sdreport(Sdreport)
    SD_mean_Z_ctm = array( NA, dim=c(unlist(TmbData[c('n_c','n_t','n_m')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) )
    if( use_biascorr==TRUE && "unbiased"%in%names(Sdreport) ){
      SD_mean_Z_ctm[] = SD[which(rownames(SD)==CogName),c('Est. (bias.correct)','Std. Error')]
    }
    if( !any(is.na(SD_mean_Z_ctm)) ){
      message("Using bias-corrected estimates for center of gravity...")
    }else{
      message("Not using bias-corrected estimates for center of gravity...")
      SD_mean_Z_ctm[] = SD[which(rownames(SD)==CogName),c('Estimate','Std. Error')]
    }

    # Plot center of gravity
    png( file=FileName_COG, width=6.5, height=TmbData$n_c*2, res=200, units="in")
      par( mar=c(2,2,1,0), mgp=c(1.75,0.25,0), tck=-0.02, oma=c(1,1,0,1.5), mfrow=c(TmbData$n_c,dim(SD_mean_Z_ctm)[[3]]), ... )  #
      for( cI in 1:TmbData$n_c ){
      for( mI in 1:dim(SD_mean_Z_ctm)[[3]]){
        Ybounds = (SD_mean_Z_ctm[cI,,mI,'Estimate']%o%rep(interval_width,2) + SD_mean_Z_ctm[cI,,mI,'Std. Error']%o%c(-interval_width,interval_width))
        Ylim = range(Ybounds,na.rm=TRUE)
        plot_lines(x=Year_Set, y=SD_mean_Z_ctm[cI,,mI,'Estimate'], ybounds=Ybounds, col_bounds=rgb(1,0,0,0.2), fn=plot, type="l", lwd=2, col="red", bounds_type="shading", ylim=Ylim, xlab="", ylab="", main="")
        if( cI==1 ) mtext(side=3, text=Znames[mI], outer=FALSE )
        if( mI==dim(SD_mean_Z_ctm)[[3]] & TmbData$n_c>1 ) mtext(side=4, text=category_names[cI], outer=FALSE, line=0.5)
      }}
      mtext( side=1:2, text=c("Year","Location"), outer=TRUE, line=c(0,0) )
    dev.off()

    # Write to file
    COG_Table = NULL
    for( cI in 1:TmbData$n_c ){
    for( mI in 1:dim(SD_mean_Z_ctm)[[3]]){
      Tmp = cbind("m"=mI, "Year"=Year_Set, "COG_hat"=SD_mean_Z_ctm[cI,,mI,'Estimate'], "SE"=SD_mean_Z_ctm[cI,,mI,'Std. Error'])
      if( TmbData$n_c>1 ) Tmp = cbind( "Category"=category_names[cI], Tmp)
      COG_Table = rbind(COG_Table, Tmp)
    }}

    # Plot area
    #KernelArea_Table = cbind("Year"=Year_Set, "KernelArea"=SD_log_area_Z_tmm[,2,1,1], "SE"=SD_log_area_Z_tmm[,2,1,2])
    #png( file=FileName_Area, width=4, height=4, res=200, units="in")
    #  par( mfrow=c(1,1), mar=c(3,3,2,0), mgp=c(1.75,0.25,0), tck=-0.02, oma=c(0,0,0,0))
    #  plot_lines( x=Year_Set, y=SD_log_area_Z_tmm[,2,1,1], ybounds=SD_log_area_Z_tmm[,2,1,1]%o%rep(1,2)+SD_log_area_Z_tmm[,2,1,2]%o%c(-1,1), fn=plot, bounds_type="shading", col_bounds=rgb(1,0,0,0.2), col="red", lwd=2, xlab="Year", ylab="ln(km^2)", type="l", main="Kernel approximation to area occupied")
    #dev.off()

    # Return stuff
    Return = c(Return, list("SD_mean_Z_ctm"=SD_mean_Z_ctm, "COG_Table"=COG_Table))
  }

  #
  # Only run if necessary outputs are available
  if( !any(c("effective_area_tl","effective_area_ctl","effective_area_cyl") %in% names(Report)) ){
    message( "To plot effective area occupied, please re-run with Options['Calculate_effective_area']=1" )
  }else{
    message( "Plotting effective area occupied..." )

    # Extract estimates
    SD = TMB::summary.sdreport(Sdreport)
    SD_effective_area_ctl = SD_log_effective_area_ctl = array( NA, dim=c(unlist(TmbData[c('n_c','n_t','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) )
    # Effective area
    if( use_biascorr==TRUE && "unbiased"%in%names(Sdreport) ){
      SD_effective_area_ctl[] = SD[which(rownames(SD)==EffectiveName),c('Est. (bias.correct)','Std. Error')]
    }
    if( !any(is.na(SD_effective_area_ctl)) ){
      message("Using bias-corrected estimates for effective area occupied (natural scale)...")
    }else{
      message("Not using bias-corrected estimates for effective area occupied (natural scale)...")
      SD_effective_area_ctl[] = SD[which(rownames(SD)==EffectiveName),c('Estimate','Std. Error')]
    }
    # Log-Effective area
    if( use_biascorr==TRUE && "unbiased"%in%names(Sdreport) ){
      SD_log_effective_area_ctl[] = SD[which(rownames(SD)==paste0("log_",EffectiveName)),c('Est. (bias.correct)','Std. Error')]
    }
    if( !any(is.na(SD_log_effective_area_ctl)) ){
      message("Using bias-corrected estimates for effective area occupied (log scale)...")
    }else{
      message("Not using bias-corrected estimates for effective area occupied (log scale)...")
      SD_log_effective_area_ctl[] = SD[which(rownames(SD)==paste0("log_",EffectiveName)),c('Estimate','Std. Error')]
    }

    # Plot area
    plot_index( Index_ctl=array(SD_log_effective_area_ctl[,,,'Estimate'],dim(SD_log_effective_area_ctl)[1:3]),
      sd_Index_ctl=array(SD_log_effective_area_ctl[,,,'Std. Error'],dim(SD_log_effective_area_ctl)[1:3]),
      Year_Set=Year_Set, Years2Include=Years2Include, strata_names=strata_names, category_names=category_names,
      DirName="", PlotName=FileName_EffArea, scale="uniform",
      interval_width=interval_width, xlab="Year", ylab="Effective area occupied, ln(km^2)", Yrange=c(NA,NA),
      width=ceiling(TmbData$n_c/ceiling(sqrt(TmbData$n_c)))*4, height=ceiling(sqrt(TmbData$n_c))*4 )
    #png( file=FileName_EffArea, width=ceiling(TmbData$n_c/ceiling(sqrt(TmbData$n_c)))*2.5, height=ceiling(sqrt(TmbData$n_c))*2.5, res=200, units="in")
    #  par( mfrow=c(1,1), mar=c(2,2,1,0), mgp=c(1.75,0.25,0), tck=-0.02, oma=c(1,1,1,0), mfrow=c(ceiling(sqrt(TmbData$n_c)),ceiling(TmbData$n_c/ceiling(sqrt(TmbData$n_c)))))
    #  for( cI in 1:TmbData$n_c ){
    #    Ybounds = SD_log_effective_area_ctl[cI,,1,1]%o%rep(interval_width,2) + SD_log_effective_area_ctl[cI,,1,2]%o%c(-interval_width,interval_width)
    #    plot_lines( x=Year_Set, y=SD_log_effective_area_ctl[cI,,1,1], ybounds=Ybounds, ylim=range(Ybounds), fn=plot, bounds_type="shading", col_bounds=rgb(1,0,0,0.2), col="red", lwd=2, xlab="", ylab="", type="l", main=category_names[cI])
    #  }
    #  mtext( side=1:3, text=c("Year","ln(km^2)","Effective area occupied"), outer=TRUE, line=c(0,0,0) )
    #dev.off()

    # Write to file
    EffectiveArea_Table = NULL
    for( cI in 1:TmbData$n_c ){
      Tmp = cbind("Year"=Year_Set, "EffectiveArea"=SD_log_effective_area_ctl[cI,,1,1], "SE"=SD_log_effective_area_ctl[cI,,1,2])
      if( TmbData$n_c>1 ) Tmp = cbind( "Category"=category_names[cI], Tmp)
      EffectiveArea_Table = rbind(EffectiveArea_Table, Tmp)
    }

    # Return stuff
    Return = c(Return, list("SD_effective_area_ctl"=SD_effective_area_ctl, "SD_log_effective_area_ctl"=SD_log_effective_area_ctl, "EffectiveArea_Table"=EffectiveArea_Table))
  }

  # Return list of stuff
  return( invisible(Return) )
}

