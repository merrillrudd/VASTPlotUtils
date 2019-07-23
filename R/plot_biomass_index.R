
#' @title{Plot index of abundance}
#'
#' @description{ plots an index proportion to population abundance}
#'
#' @param TmbData Formatted data inputs, from `VAST::make_data(...)`
#' @param Sdreport output from fit$parameter_estimates$SD
#' @param DirName Directory for saving plot and table
#' @param PlotName Name for plot
#' @param interval_width width for confidence intervals
#' @param strata_names names for spatial strata
#' @param category_names names for categories
#' @param use_biascorr Boolean, whether to use bias-corrected estimates if available
#' @param plot_legend Add legend for labelling colors
#' @param total_area_km2 Total area for calculating a design-based estimator using one design-stratum (only recommended for model exploration)
#' @param plot_log Boolean, whether to plot y-axis in log-scale
#' @param width plot width in inches
#' @param height plot height in inches
#' @param ... Other inputs to `par()`
#' @inheritParams plot_maps
#'
#' @return Return Tagged list of output
#' \describe{
#'   \item{Table}{table of index estimates by stratum and year, e.g., for including in an assessment model}
#' }
#'

#' @export
plot_biomass_index <-
function( TmbData, Sdreport, Year_Set=NULL, Years2Include=NULL, DirName=paste0(getwd(),"/"), PlotName="Index", interval_width=1,
  strata_names=NULL, category_names=NULL, use_biascorr=TRUE, plot_legend=TRUE, total_area_km2=NULL, plot_log=FALSE, width=4, height=4,
  create_covariance_table=FALSE, ... ){

  require(ggplot2)

  # Informative errors
  if(is.null(Sdreport)) stop("Sdreport is NULL; please provide Sdreport")
  if(!is.null(category_names) && length(category_names)!=TmbData$n_c ) stop("`category_names` must have same length as `TmbData$n_c`")
  if(!is.null(strata_names) && length(strata_names)!=TmbData$n_l ) stop("`strata_names` must have same length as `TmbData$n_l`")

  # Which parameters
  if( "ln_Index_tl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # SpatialDeltaGLMM
    ParName = "Index_tl"
    TmbData[['n_c']] = 1
  }
  if( "ln_Index_ctl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # VAST Version < 2.0.0
    ParName = "Index_ctl"
  }
  if( "ln_Index_cyl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # VAST Version >= 2.0.0
    ParName = "Index_cyl"
    TmbData[["n_t"]] = nrow(TmbData[["t_yz"]])
  }
  if( "Index_tp" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    # SpatialVAM
    ParName = "Index_tp"
    TmbData[["n_l"]] = 1
    TmbData[["n_c"]] = TmbData[["n_p"]]
  }

  # Add t_iz if missing (e.g., from earlier version of VAST, or SpatialDeltaGLMM)
  if( !("t_iz" %in% names(TmbData)) ){
    TmbData$t_iz = matrix( TmbData$t_i, ncol=1 )
  }

  # Add in t_yz if missing (e.g., from earlier version of VAST, or SpatialDeltaGLMM)
  if( !("t_yz" %in% names(TmbData)) ){
    TmbData$t_yz = matrix(1:TmbData$n_t - 1, ncol=1)
  }

  # Fill in missing
  Year_Set = 1:TmbData$n_t
  Years2Include = 1:TmbData$n_t
  if( is.null(strata_names) ) strata_names = 1:TmbData$n_l
  if( is.null(category_names) ) category_names = 1:TmbData$n_c

  # Logical check
  if( "unbiased"%in%names(Sdreport) ){
    if( all(is.na(Sdreport$unbiased$value)) ){
      stop("You appear to be using bias-correction, but all values are NA. Please report problem to package author.")
    }
  }

  # Defaults
  if( "treat_nonencounter_as_zero" %in% names(TmbData$Options_list$Options) ){
    treat_missing_as_zero = TmbData$Options_list$Options["treat_nonencounter_as_zero"]
  }else{
    treat_missing_as_zero = FALSE
  }

  # Objects
  SD = TMB::summary.sdreport(Sdreport)
  if( !"report" %in% names(as.list(args(TMB:::as.list.sdreport))) ){
    warning( "package `TMB` should be updated to easily access standard errors")
  }
  SD_stderr = TMB:::as.list.sdreport( Sdreport, what="Std. Error", report=TRUE )
  SD_estimate = TMB:::as.list.sdreport( Sdreport, what="Estimate", report=TRUE )
  if( use_biascorr==TRUE && "unbiased"%in%names(Sdreport) ){
    SD_estimate_biascorrect = TMB:::as.list.sdreport( Sdreport, what="Est. (bias.correct)", report=TRUE )
  }
  if( any(is.na(SD_estimate)) | any(is.na(SD_stderr)) ){
    stop( "Problem: Standard errors contain NAs")
  }

  # Extract index (using bias-correctino if available and requested)
  if( ParName %in% c("Index_tl","Index_ctl","Index_cyl")){
    Index_ctl = log_Index_ctl = array( NA, dim=c(unlist(TmbData[c('n_c','n_t','n_l')]),2), dimnames=list(category_names,Year_Set,strata_names,c('Estimate','Std. Error')) )
    # Index
    if( use_biascorr==TRUE && "unbiased"%in%names(Sdreport) ){
      Index_ctl[] = SD[which(rownames(SD)==ParName),c('Est. (bias.correct)','Std. Error')]
    }
    if( !any(is.na(Index_ctl)) ){
      message("Using bias-corrected estimates for abundance index (natural-scale)...")
    }else{
      message("Not using bias-corrected estimates for abundance index (natural-scale)...")
      Index_ctl[] = SD[which(rownames(SD)==ParName),c('Estimate','Std. Error')]
    }
    # Log-Index
    if( use_biascorr==TRUE && "unbiased"%in%names(Sdreport) ){
      log_Index_ctl[] = SD[which(rownames(SD)==paste0("ln_",ParName)),c('Est. (bias.correct)','Std. Error')]
    }
    if( !any(is.na(log_Index_ctl)) ){
      message("Using bias-corrected estimates for abundance index (log-scale)...")
    }else{
      message("Not using bias-corrected estimates for abundance index (log-scale)...")
      log_Index_ctl[] = SD[which(rownames(SD)==paste0("ln_",ParName)),c('Estimate','Std. Error')]
    }
  }
  if( ParName %in% c("Index_tp")){
    if( use_biascorr==TRUE && "unbiased"%in%names(Sdreport) ){
      Index_ctl = aperm( array( c(Sdreport$unbiased$value[which(names(Sdreport$value)==ParName)],TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))==ParName),'Std. Error']), dim=c(unlist(TmbData[c('n_t','n_c','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) ), perm=c(2,1,3))
      if( "ln_Index_tp" %in% rownames(TMB::summary.sdreport(Sdreport))){
        log_Index_ctl = aperm( array( c(Sdreport$unbiased$value[which(names(Sdreport$value)==paste0("ln_",ParName))],TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))==paste0("ln_",ParName)),'Std. Error']), dim=c(unlist(TmbData[c('n_t','n_c','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) ), perm=c(2,1,3))
      }else{
        log_Index_ctl = log( Index_ctl )
        log_Index_ctl[,,,'Std. Error'] = log_Index_ctl[,,,'Std. Error'] / log_Index_ctl[,,,'Estimate']
        warning( "Using kludge for log-standard errors of index, to be replaced in later versions of 'MIST'" )
      }
    }else{
      Index_ctl = aperm( array( TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))==ParName),], dim=c(unlist(TmbData[c('n_t','n_c','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) ), perm=c(2,1,3,4))
      if( "ln_Index_tp" %in% rownames(TMB::summary.sdreport(Sdreport))){
        log_Index_ctl = aperm( array( TMB::summary.sdreport(Sdreport)[which(rownames(TMB::summary.sdreport(Sdreport))==paste0("ln_",ParName)),], dim=c(unlist(TmbData[c('n_t','n_c','n_l')]),2), dimnames=list(NULL,NULL,NULL,c('Estimate','Std. Error')) ), perm=c(2,1,3,4))
      }else{
        log_Index_ctl = log( Index_ctl )
        log_Index_ctl[,,,'Std. Error'] = log_Index_ctl[,,,'Std. Error'] / log_Index_ctl[,,,'Estimate']
        warning( "Using kludge for log-standard errors of index, to be replaced in later versions of 'MIST'" )
      }
    }
  }

  # Extract biomass ratio Bratio_cty if available (only available if >= V5.3.0 and using spatial Gompertz model features)
  if( "Bratio_cyl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    Bratio_ctl = array( NA, dim=c(unlist(TmbData[c('n_c','n_t','n_l')]),2), dimnames=list(category_names,Year_Set,strata_names,c('Estimate','Std. Error')) )
    if( use_biascorr==TRUE && "unbiased"%in%names(Sdreport) ){
      Bratio_ctl[] = SD[which(rownames(SD)=="Bratio_cyl"),c('Est. (bias.correct)','Std. Error')]
    }
    if( !any(is.na(Bratio_ctl)) ){
      message("Using bias-corrected estimates for biomass ratio (natural-scale)...")
    }else{
      message("Not using bias-corrected estimates for biomass ratio (natural-scale)...")
      Bratio_ctl[] = SD[which(rownames(SD)=="Bratio_cyl"),c('Estimate','Std. Error')]
    }
  }else{
    Bratio_ctl = NULL
  }
  if( "ln_Bratio_cyl" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    log_Bratio_ctl = array( NA, dim=c(unlist(TmbData[c('n_c','n_t','n_l')]),2), dimnames=list(category_names,Year_Set,strata_names,c('Estimate','Std. Error')) )
    if( use_biascorr==TRUE && "unbiased"%in%names(Sdreport) ){
      log_Bratio_ctl[] = SD[which(rownames(SD)=="ln_Bratio_cyl"),c('Est. (bias.correct)','Std. Error')]
    }
    if( !any(is.na(log_Bratio_ctl)) ){
      message("Using bias-corrected estimates for biomass ratio (log-scale)...")
    }else{
      message("Not using bias-corrected estimates for biomass ratio (log-scale)...")
      log_Bratio_ctl[] = SD[which(rownames(SD)=="ln_Bratio_cyl"),c('Estimate','Std. Error')]
    }
  }else{
    log_Bratio_ctl = NULL
  }

  # Extract Fratio
  if( "Fratio_ct" %in% rownames(TMB::summary.sdreport(Sdreport)) ){
    Fratio_ct = array( NA, dim=c(unlist(TmbData[c('n_c','n_t')]),2), dimnames=list(category_names,Year_Set,c('Estimate','Std. Error')) )
    Fratio_ct[] = SD[which(rownames(SD)=="Fratio_ct"),c('Estimate','Std. Error')]
    #Fratio_ct = abind::abind( SD_estimate$Fratio, SD_stderr$Fratio, along=3 )
    #dimnames(Fratio_ct) = list(category_names,Year_Set,c('Estimate','Std. Error'))
  }else{
    Fratio_ct = NULL
  }

  # Calculate design-based
  if( !is.null(total_area_km2) & TmbData$n_c==1 ){
    message( "Calculating naive design-based index -- do not use this, its intended only for comparison purposes" )
    Calc_design = TRUE
    Design_t = tapply( TmbData$b_i/TmbData$a_i, INDEX=TmbData$t_i, FUN=mean ) * total_area_km2 / 1000 # Convert to tonnes
    Design_t = cbind( "Estimate"=Design_t, "Std. Error"=sqrt(tapply(TmbData$b_i/TmbData$a_i,INDEX=TmbData$t_i,FUN=var)/tapply(TmbData$b_i/TmbData$a_i,INDEX=TmbData$t_i,FUN=length))*total_area_km2/1000)
    Design_t = cbind( Design_t, "CV"=Design_t[,'Std. Error'] / Design_t[,'Estimate'] )
  }else{
    Calc_design = FALSE
  }

  # Fix at zeros any years-category combinations with no data
  if( treat_missing_as_zero==TRUE ){
    # Determine year-category pairs with no data
    Num_ct = tapply( TmbData$b_i, INDEX=list(factor(TmbData$c_i,levels=1:TmbData$n_c-1),factor(TmbData$t_i[,1],levels=1:TmbData$n_t-1)), FUN=function(vec){sum(!is.na(vec))} )
    Num_ct = ifelse( is.na(Num_ct), 0, Num_ct )
    # Replace values with 0 (estimate) and NA (standard error)
    Index_ctl[,,,'Estimate'] = ifelse(Num_ct%o%rep(1,TmbData$n_l)==0, 0, Index_ctl[,,,'Estimate'])
    Index_ctl[,,,'Std. Error'] = ifelse(Num_ct%o%rep(1,TmbData$n_l)==0, NA, Index_ctl[,,,'Std. Error'])
    log_Index_ctl[,,,'Estimate'] = ifelse(Num_ct%o%rep(1,TmbData$n_l)==0, -Inf, log_Index_ctl[,,,'Estimate'])
    log_Index_ctl[,,,'Std. Error'] = ifelse(Num_ct%o%rep(1,TmbData$n_l)==0, NA, log_Index_ctl[,,,'Std. Error'])
  }

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

  # Plot biomass and Bratio
  if(all(is.integer(TmbData$b_i))){
    Plot_suffix = "Count"
  } else { Plot_suffix = "Biomass"}
  if( !is.null(Bratio_ctl) ) Plot_suffix = c( Plot_suffix, "Bratio" )

  for( plotI in 1:length(Plot_suffix) ){


      if( Plot_suffix[plotI]=="Count" ){ Array_ctl = Index_ctl * 1000; log_Array_ctl = log_Index_ctl }
      if( Plot_suffix[plotI]=="Biomass" ){ Array_ctl = Index_ctl; log_Array_ctl = log_Index_ctl }
      if( Plot_suffix[plotI]=="Bratio" ){ Array_ctl = Bratio_ctl; log_Array_ctl = log_Bratio_ctl }

      byCat <- lapply(1:TmbData$n_c, function(x){
        byStrat <- lapply(1:TmbData$n_l, function(y){
          est <- Array_ctl[x,Years2Include,y,1]
          ybounds = (Array_ctl[x,Years2Include,y,'Estimate']%o%c(1,1))*exp(log_Array_ctl[x,Years2Include,y,'Std. Error']%o%c(-interval_width,interval_width))
          df <- data.frame("Year"=names(est), "Category"=category_names[x], "Stratum"=strata_names[y], "Estimate"=est, "Ybound_low"=ybounds[,1], "Ybound_high"=ybounds[,2])
          return(df)
        })
        byStrat <- do.call(rbind, byStrat)
      })
      byCat <- do.call(rbind, byCat)
      byCat$Year <- as.numeric(byCat$Year)

      Ylim = lapply(1:TmbData$n_c, function(x){
        c(0, max(Array_ctl[x,Years2Include,,'Estimate']%o%c(1,1) * exp(log_Array_ctl[x,Years2Include,,'Std. Error']%o%c(-interval_width,interval_width)),na.rm=TRUE)*1.05 )
      })
      blank_info <- data.frame(Category = c(sapply(1:TmbData$n_c, function(x) rep(category_names[x], 2))), x = 0, y = c(sapply(1:TmbData$n_c, function(x) Ylim[[x]])))
      
      name <- switch(Plot_suffix[plotI], "Biomass"="Abundance (metric tonnes)", "Count"="Abundance (individuals)", "Bratio"="Biomass ratio")
      byCat$Stratum <- factor(byCat$Stratum)
      p <- ggplot(byCat) +
        geom_segment(aes(x = Year, y = Ybound_low, xend = Year, yend = Ybound_high, color=Stratum)) +
        geom_line(aes(x = Year, y=Estimate, color = Stratum)) +
        geom_point(aes(x = Year, y=Estimate, color = Stratum), cex=2) +
        geom_blank(data=blank_info, aes(x=x, y=y)) +
        facet_grid(Category~., scales="free_y") +
        scale_x_continuous(breaks = seq(min(byCat$Year),max(byCat$Year),by=5), labels=Year_Set[Years2Include][seq(min(byCat$Year),max(byCat$Year),by=5)]) +
        expand_limits(y = 0) +
        scale_y_continuous(expand = c(0,0)) +
        # coord_cartesian(ylim = Ylim) +
        scale_color_brewer(palette = "Set1") +
        ylab(name) +
        mytheme()
      if(length(unique(strata_names))==1) p <- p + guides(color = FALSE)
      if(!is.null(DirName)) ggsave(paste0(DirName,"/",PlotName,"-",Plot_suffix[plotI],".png"), p)
    }



  # Write to file
  Table = NULL
  for( cI in 1:TmbData$n_c ){
    Tmp = data.frame( "Year"=Year_Set, "Unit"=1, "Fleet"=rep(strata_names,each=TmbData$n_t), "Estimate_metric_tons"=as.vector(Index_ctl[cI,,,'Estimate']), "SD_log"=as.vector(log_Index_ctl[cI,,,'Std. Error']), "SD_mt"=as.vector(Index_ctl[cI,,,'Std. Error']) )
    if( TmbData$n_c>1 ) Tmp = cbind( "Category"=category_names[cI], Tmp)
    Table = rbind( Table, Tmp )
  }
  if(!is.null(total_area_km2)) Table = cbind(Table, "Naive_design-based_index"=Design_t)
  write.csv( Table, file=paste0(DirName,"/Table_for_SS3.csv"), row.names=FALSE)

  # Return stuff
  Return = list( "Table"=Table, "log_Index_ctl"=log_Index_ctl, "Index_ctl"=Index_ctl, "Ylim"=Ylim )

  # Extract and save covariance
  if( "cov"%in%names(Sdreport) & create_covariance_table==TRUE ){
    DF = expand.grid( "Category"=1:TmbData$n_c, "Year"=1:TmbData$n_t, "Stratum"=1:TmbData$n_l )
    Which = which( names(Sdreport$value)==ParName )
    Cov = Sdreport$cov[Which,Which]
    Corr = cov2cor(Cov) - diag(nrow(Cov))
    rowcolDF = cbind( "RowNum"=row(Corr)[lower.tri(Corr,diag=TRUE)], "ColNum"=col(Corr)[lower.tri(Corr,diag=TRUE)] )
    Table = cbind( DF[rowcolDF[,'ColNum'],], DF[rowcolDF[,'RowNum'],] )
    colnames(Table) = paste0(colnames(Table), rep(c(1,2),each=3))
    Table = cbind( Table, "Correlation"=cov2cor(Cov)[lower.tri(Corr,diag=TRUE)], "Covariance"=Cov[lower.tri(Corr,diag=TRUE)] )
    Table = cbind( Table, "Index1"=Index_ctl[as.matrix(cbind(DF[rowcolDF[,'ColNum'],],1))], "Index2"=Index_ctl[as.matrix(cbind(DF[rowcolDF[,'RowNum'],],1))] )
    WhichZero = which( (Table[,'Index1']*Table[,'Index2']) == 0 )
    Table[WhichZero,c('Correlation','Covariance')] = 0
    Return = c( Return, "Table_of_estimted_covariance"=Table )
  }
  if( !is.null(Bratio_ctl)) Return = c( Return, list("Bratio_ctl"=Bratio_ctl) )
  if( !is.null(log_Bratio_ctl)) Return = c( Return, list("log_Bratio_ctl"=log_Bratio_ctl) )
  if( !is.null(Fratio_ct)) Return = c( Return, list("Fratio_ct"=Fratio_ct) )

  return( invisible(Return) )
}
