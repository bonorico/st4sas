## average_cumulative_hazard.R contains modules to avarage over several simulated Kaplan-Meier / cumulative hazard curves.
## Copyright (C) 2019 Federico Bonofiglio

## This file is part of st4sas.

## st4sas is free software: you can redistribute it and/or modify
## it under the terms of the GNU General Public License as published by
## the Free Software Foundation, either version 3 of the License, or
## (at your option) any later version.

## st4sas is distributed in the hope that it will be useful,
## but WITHOUT ANY WARRANTY; without even the implied warranty of
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
## GNU General Public License for more details.

## You should have received a copy of the GNU General Public License
## along with st4sas.  If not, see <https://www.gnu.org/licenses/>.






#' @name count.arm.events
#' @title Time-events (points) count of a (stratified) Kaplan-Meier (Nelson-Aalen/Breslow) estimate. 
#'
#' @param cumhaz.matrix a data.frame with information on Kaplan-Meier (Nelson-Aalen/Breslow) values at each event-time point, along point-wise confidence interval (CI) values and treatment/control indicator. By convention such information must be arranged in columns with name 'km' (survival probability), 'time' (event-time), 'km_l' (lower CI value for 'km'), 'km_u' (upper CI value for 'km'), 'trt' (binary indicator for a treatment/control grouping, can be NA).
#'
#' @details basic information can be, for instance, computed via 'survfit' from package 'survival'.


count.arm.events <- function( cumhaz.matrix ){
    
    if ( all(!is.na(cumhaz.matrix$trt)) ) {
        
        nevents <- as.vector( table( cumhaz.matrix$trt ) )  # trt 0 , trt 1  

        data.frame( N.trt0 = nevents[1], N.trt1 = nevents[2] )

    } else {

        dim(cumhaz.matrix)[1] # unstratified length

    }
    
    
}




                                        # extend curve/time length to pregiven size

extend.time.jumps <- function( time, maxn ){

                                        # set an average jump length
    mean.dt <- mean( diff(time), na.rm = T ) 
    
    time.length <- length(time)
    last.time.point <- time[time.length]
    
    if ( maxn >=  time.length )
        N.to.add <- maxn - time.length
    else
        stop( "maxn is lesser than cumhaz length" )

    next.last.time <- last.time.point + mean.dt # last time point + mean(dt)
    new.last.timepoint <- last.time.point + ( mean.dt*N.to.add ) # accrue last time point by N times the mean dt
    
    added.time <- seq( next.last.time, new.last.timepoint, length =  N.to.add) # no overlapping with last time point

    extended.time.line <- c( time, added.time )

    return( extended.time.line )
    
} 

                                        #

extend.cumhaz.values <- function( cumhaz, maxn){

    curve.length <- length(cumhaz)
    lastvalue <- cumhaz[curve.length]
    
    if ( maxn >=  curve.length )
        N.to.add <- maxn - curve.length
    else
        stop( "maxn is lesser than cumhaz length" )
    
    add.values <- rep( lastvalue, N.to.add )

    extended.cumhaz <- c( cumhaz, add.values )  # add values on right tail

    return( extended.cumhaz )
    
}



#' @name Add.jumps.to.cumhaz.matrix
#' @title Function to extend length of a Kaplan-Meier curve
#'
#' @param cumhaz.matrix a data.frame with information on Kaplan-Meier (Nelson-Aalen/Breslow) values at each event-time point, along point-wise confidence interval (CI) values and treatment/control indicator. By convention such information must be arranged in columns with name 'km' (survival probability), 'time' (event-time), 'km_l' (lower CI value for 'km'), 'km_u' (upper CI value for 'km'), 'trt' (binary indicator for a treatment/control grouping, can be NA).
#'
#' @param maxevents number of (group-specific) event-times (points) used to extend curve length. 
#' @details The last available KM value is stretched as a constant between the last and future event-time. The event-time line is "evenly" stretched between the last.time.point + mean(diff(time)) and last.time.point + (mean(diff(time))*number.of.point.to.add). This function can be useful when having several simulated curves that have different length due to randomness. Then, the function normalizes all simulated curves to a common 'maxevents' length (e.g. the max length). This can be later used to average across curves, for example.
#' Basic information for 'cumhaz.matrix' can be, for instance, computed via 'survfit' from package 'survival'.

Add.jumps.to.cumhaz.matrix <- function( cumhaz.matrix, maxevents ){

    
    K <- length( maxevents )
    trt.level <- unique(cumhaz.matrix$trt) # 0,1 OR NA. problem if both
    if ( K != length( trt.level )  )
        stop( "maxevents strata do not match trt strata" )

    if ( length(trt.level) == 1) {
        if ( is.na(trt.level) ){
            trt.level <- 1
            cumhaz.matrix$trt <- is.na(cumhaz.matrix$trt)
        } }
    
                                        # (split-modify-recompose dataset)
    
    out <- do.call( "rbind",  # recompose
                   lapply( 1:K,
                          function(i){
                              splitdata <- cumhaz.matrix[cumhaz.matrix$trt == trt.level[i], ]

                              ## here you need follow some columns numbering/labelling convention ...
                              
                              data.frame( km = extend.cumhaz.values( splitdata$km, maxevents[i]),
                                         time = extend.time.jumps( splitdata$time, maxevents[i] ),
                                         km_l = extend.cumhaz.values( splitdata$km_l, maxevents[i]),
                                         km_u = extend.cumhaz.values( splitdata$km_u, maxevents[i]),
                                         trt = trt.level[i]  )
                          }         )                               
                   )
    
    return(out)
    
}



##

                                        #
compute.average.curve <- function( extended.km.matrix, extended.time.matrix, extended.km_l.matrix, extended.km_u.matrix, extended.trt ){
    
    mean.km <- apply( extended.km.matrix, 1, mean, na.rm = T ) # MC average cumhaz
    mean.time <- apply( extended.time.matrix, 1, mean, na.rm = T )  # MC average time
    mean.km_l <- apply( extended.km_l.matrix, 1, mean, na.rm = T )
    mean.km_u <- apply( extended.km_u.matrix, 1, mean, na.rm = T ) 
    averaged.frame <- data.frame( km = mean.km, time = mean.time, km_l = mean.km_l, km_u = mean.km_u, trt = extended.trt )

    return(averaged.frame)

}
                                        #



#' @name Average.MC.cumhaz
#' @title Average several Kaplan-Meier (or Nelson-Aalen/Breslow) curves
#'
#' @param MClist a list of 'cumhaz.matrix' objects (see coun.arm.events), as for instance obtained via a simulation/resampling procedure, and used to produce a single averaged Kaplan-Meier (or Nelson-Aalen/Breslow) curve
#'
#' @examples
#' 
#'help("st4sas")
#' 


Average.MC.cumhaz <- function( MClist ){ # MClist is a list of MC breslow curves
    
                                        # step 0 : count arm-specific events
    nevents <- do.call( "rbind", lapply( MClist, function(x) count.arm.events(x) ) )     
    
                                        # step 1 : determine longest curve in MClist, fix maxlength

    maxNevents <- apply( nevents, 2, max , na.rm = T )
    
                                        # step 2 : extend shorter curves to maxlenght and fill empty elements with maxcumhaz, and stack into matrix 
    
    extended.MClist <- lapply( MClist, function(x) Add.jumps.to.cumhaz.matrix( x, maxNevents ) )
    extended.MClist <- extended.MClist[sapply(extended.MClist, function(x)!is.null(x))]  # drp eventual NULL vals

    extended.trt <- extended.MClist[[1]]$trt  # all frames should have equal trt tabulation (NOT TRUE, SOME ARE MISSING) ???? TODO(me) : check this claim ... update: WHAT do u mean ?? (29/06/17)

    extended.km <- lapply(extended.MClist, function(x) x$km) 
    extended.time <- lapply(extended.MClist, function(x) x$time) 
    extended.km_l <- lapply(extended.MClist, function(x) x$km_l) 
    extended.km_u <- lapply(extended.MClist, function(x) x$km_u)                                                    
                                        # step 3 : compute rowwise averages (average cumhaz)

    averaged.frame <- compute.average.curve( do.call("cbind",extended.km), do.call("cbind",extended.time), do.call("cbind",extended.km_l), do.call("cbind",extended.km_u), extended.trt ) 
                                        # step 4 (see below)

    H <- length(extended.MClist)
    sampleNr <- rep(paste("sample", 1:H), rep(length(extended.km[[1]]), H))

    raw.samples <- data.frame(km = unlist( extended.km), time = unlist( extended.time), km_l = unlist( extended.km_l), km_u = unlist( extended.km_u), trt = extended.trt, sampleNr =  sampleNr)

    out <- list( "mc.average" = averaged.frame, # MC average (extended event time- to be chipped)
                "raw.extended" =  raw.samples  )                                                                                                            
    return( out )
    
}






