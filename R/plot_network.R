#' @title
#' Plot stream network and optional observations
#' 
#' @description
#' \code{plot_network} generates figure showing network nodes, directional arrows, and observations
#' 
#' @param Network_sz_LL
#' @param Data data frame with observations
#' @param plot_type 0 = data (Catch_KG), 1 = PR1 residuals for first linear predictor, or 2 PR2 residuals for second linear predictor
#' @param byYear plot observations by year, default = FALSE
#' @param byValue plot observations representing size of point by their value
#' @param value_label label for the value observed
#' @param root default TRUE to show root nodes, FALSE in case there are other root nodes that are not meaningful.
#' @param arrows default = FALSE do not plot segments
#' @param FilePath path to figure directory
#' @param FileName label for plot
#' @param obs_color option for hard-coding the colors for figures, used to match colors for plotting only one category when other figures are multiple categories (optional)
#' @return Figure plotting stream network and observations
#' @export
plot_network <- function(Network_sz_LL, Data=NULL, plot_type = 0, byYear = FALSE, byValue=FALSE, value_label=NULL, FilePath=NULL, FileName="Network", arrows = FALSE, root = FALSE, obs_color=NULL){

  if(byYear==FALSE){
    ### across years
    aa <- ggplot(Network_sz_LL) +
      mytheme()
    
    ## add roots underneath points
    if(root == TRUE) aa <- aa + geom_point(data = Network_sz_LL %>% filter(parent_s==0), aes(x = Lon, y = Lat), color="goldenrod", cex=5, alpha=0.5)
    
    ## add points and complete figure
    aa <- aa +
      geom_point(data = Network_sz_LL, aes(x = Lon, y = Lat), color = "gray", alpha=0.6) +
      xlab("Longitude") + ylab("Latitude")  


    ## option to add arrows
    if(arrows == TRUE){
      l2 <- lapply(1:nrow(Network_sz_LL), function(x){
        parent <- Network_sz_LL$parent_s[x]
        find <- Network_sz_LL %>% filter(child_s == parent)
        if(nrow(find)>0) out <- cbind.data.frame(Network_sz_LL[x,], 'Lon2'=find$Lon, 'Lat2'=find$Lat)
        if(nrow(find)==0) out <- cbind.data.frame(Network_sz_LL[x,], 'Lon2'=NA, 'Lat2'=NA)
        # if(nrow(find)>0) out <- cbind.data.frame(Network_sz_LL[x,], 'long2'=find$long, 'lat2'=find$lat)
        # if(nrow(find)==0) out <- cbind.data.frame(Network_sz_LL[x,], 'long2'=NA, 'lat2'=NA)
        return(out)
      })
      l2 <- do.call(rbind, l2)
      aa <- aa + geom_segment(data=l2, aes(x = Lon,y = Lat, xend = Lon2, yend = Lat2), arrow=arrow(length=unit(0.2,"cm")), col="gray")
    }

    ## option to add observations
    if(all(is.null(Data))==FALSE & byYear==FALSE){
      if(byValue==FALSE){
        aa <- aa + 
          geom_point(data = Data, aes(x = Lon, y = Lat, color=Category), cex=2) + 
          scale_color_brewer(palette = "Set1")      
      }
      if(byValue==TRUE){
        if(all(is.null(value_label))) stop("please include label for value type")
        aa <- aa + 
          scale_color_brewer(palette = "Set1") +
          guides(size=guide_legend(title=value_label))           
        if(plot_type == 0) aa <- aa +  geom_point(data = Data, aes(x = Lon, y = Lat, color=Category, size=Catch_KG), alpha=0.6)
        if(plot_type == 1) aa <- aa +  geom_point(data = Data, aes(x = Lon, y = Lat, color=PR1>0, size=abs(PR1)), alpha=0.25) + scale_size("Pearson residual", range = c(0,3))
        if(plot_type == 2) aa <- aa +  geom_point(data = Data %>% filter(positive == 1), aes(x = Lon, y = Lat, color=PR2>0, size=abs(PR2)), alpha=0.25) + scale_size("Pearson residual", range = c(0,3))
      } 

    }
    if(all(is.null(Data)) | (length(unique(Data$Category))==1 & byValue==FALSE)){
      aa <- aa + guides(fill = FALSE)
      width <- 9
      height <- 8
    } else {
      width = 10
      height = 9
    }
    if(all(is.null(FilePath))==FALSE) ggsave(file.path(FilePath, paste0(FileName, ".png")), aa, width = width, height = height)
    if(all(is.null(FilePath))) print(aa)
    return(aa)
  }
    
  if(byYear == TRUE){
    
    if(all(is.null(Data))) stop("Error: must include observations in data frame 'Data' to plot network by year.")
    
    ### by year
    years <- unique(Data$Year)[order(unique(Data$Year))]	
      
    Network_sz_LL_wYear <- lapply(1:length(years), function(x){
      out <- cbind.data.frame(Network_sz_LL, "Year"=years[x])
      return(out)
    })
    Network_sz_LL_wYear <- do.call(rbind, Network_sz_LL_wYear)	
    
    bb <- ggplot(Network_sz_LL_wYear) +
      facet_wrap(~Year) +
      mytheme() 	
    
    ## add roots underneath points
    if(root == TRUE) bb <- bb + geom_point(data = Network_sz_LL %>% filter(parent_s==0), aes(x = Lon, y = Lat), color="goldenrod", cex=5)	
    
    ## add points and complete figure
    bb <- bb +
      geom_point(data = Network_sz_LL, aes(x = Lon, y = Lat), color = "gray", cex=0.5, alpha=0.6) +
      xlab("Longitude") + ylab("Latitude")	
    
    ## option to add observations
    ## option to add observations
    bb <- bb + scale_x_continuous(breaks=round(quantile(Data$Lon,prob=c(0.2,0.8)),1), labels=round(quantile(Data$Lon,prob=c(0.2,0.8)),1)) 
    if(byValue==FALSE){
      if(is.null(obs_color)){
        bb <- bb + 
          geom_point(data = Data, aes(x = Lon, y = Lat, color=Category), cex=1.8, alpha=0.6) + 
          scale_color_brewer(palette = "Set1")
      }
      if(all(is.null(obs_color))==FALSE){
        if(length(obs_color)!=length(unique(Data$Category))) stop("input observation colors must match number of categories in data to plot")
        bb <- bb + 
          geom_point(data = Data, aes(x = Lon, y = Lat, color=Category), cex=1.8, alpha=0.6, fill=obs_color)
      }      
    }
    if(byValue==TRUE){
      if(is.null(obs_color)){
        bb <- bb + 
          # geom_point(data = Data, aes(x = Lon, y = Lat, fill=Category, size=Catch_KG), alpha=0.6) + 
          scale_color_brewer(palette = "Set1")
        if(plot_type == 0) bb <- bb +  geom_point(data = Data, aes(x = Lon, y = Lat, color=Category, size=Catch_KG), alpha=0.6)
        if(plot_type == 1) bb <- bb +  geom_point(data = Data, aes(x = Lon, y = Lat, color=PR1>0, size=abs(PR1)), alpha=0.25) + scale_size("Pearson residual", range = c(0,3))
        if(plot_type == 2) bb <- bb +  geom_point(data = Data %>% filter(positive == 1), aes(x = Lon, y = Lat, color=PR2>0, size=abs(PR2)), alpha=0.25) + scale_size("Pearson residual", range = c(0,3))
      }
      if(all(is.null(obs_color))==FALSE){
        if(length(obs_color)!=length(unique(Data$Category))) stop("input observation colors must match number of categories in data to plot")
        bb <- bb + 
          # geom_point(data = Data, aes(x = Lon, y = Lat, fill=Category, size=Catch_KG), alpha=0.6, fill=obs_color) +
          guides(size=guide_legend(title=value_label))    
        if(plot_type == 0) bb <- bb +  geom_point(data = Data, aes(x = Lon, y = Lat, color=Category, size=Catch_KG))
        if(plot_type == 1) bb <- bb +  geom_point(data = Data, aes(x = Lon, y = Lat, color=PR1>0, size=abs(PR1))) + scale_size("Pearson residual", range = c(0,3))
        if(plot_type == 2) bb <- bb +  geom_point(data = Data %>% filter(positive == 1), aes(x = Lon, y = Lat, color=PR2>0, size=abs(PR2))) + scale_size("Pearson residual", range = c(0,3))
       
      }        
    }

    if(length(unique(Data$Category))==1 & byValue==FALSE){
      bb <- bb + guides(color = FALSE)
      width <- 10
      height = 10
    } else {
      width = 12
      height = 10
    }	
    
    if(all(is.null(FilePath))==FALSE) ggsave(file.path(FilePath, paste0(FileName, ".png")), bb, width = width, height = height)
    if(all(is.null(FilePath))) print(bb)
    return(bb)
  }
  
}