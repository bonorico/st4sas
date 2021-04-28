## coxoptim.R contains modules to maximize a Cox partial likelihood having numerator 'y' in the form of a sufficient statistic.
## Copyright (C) 2018 Federico Bonofiglio

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




#' @name coxlik
#' @title Cox partial likelihood maximization using sufficient statistics formulation.
#'
#' @description Performs usual Cox partial likelihood optimization (using package 'maxLik') but uses sufficient statistics where possible (e.g. in the log likelihood numerator) instead of raw data, that is instead used in the covariates-based risk-set (log likelihood denominator). This setting can be useful in simulaiton settings where the numerator can be kept fixed and denominator raw data is simulated instead. For instance one might recover summary information on number-of-events in a treatment arm, or number-of-events times blood pressure values, but not raw data on treatment assignement and blood pressure values. Under such incomplete data situation one might still run a Cox regression noting that the summary information is a sufficient statistic in the Cox likelihood and trying to simulate raw covariate information in the denominator, akin to [1] or [2]. 
#' @param y (vector valued) sufficient statistic. It automatically enters the Cox log likelihood nominator
#'
#' @param x raw covariates data (one observation per individual) in matrix form 
#'
#' @param str categorical vector indicating several strata to be used in a stratified Cox refression.
#'
#'
#' @references [1] Borgan, Ø. and Keogh, R. (2015). Nested case–control studies: should one break the matching? Lifetime Data Analysis, 21(4):517–541.
#' 
#' [2] Keogh, R. H., Seaman, S. R., Bartlett, J. W., and Wood, A. M. (2018). Multiple imputation of missing data in nested case-control and case-cohort studies. Biometrics.



coxlik <- function(y, x, str=1){

    maxLik(lCox,  grad=grCox, hess=hessCox,
           start=rep(0.0, length(y)),
           y=y, x=x , str=str )
    
}



#' @name CollectLoglikResults
#' @title collects results from an optimized log likelihood.
#'
#' @description extracts basic optimization output and computes meaningful summaries like MLE estimates, its covariance, normally approximates 95\% confidence intervals, log likelihood deviance, and the AIC.
#' 
#' @param obj an output from the 'maxLik' function 
#'
#'
#'



CollectLoglikResults <- function(obj, summary=TRUE){

    
    if(class(obj[[1]])[1]=="maxLik")
        
        res <- do.call("rbind",

                       lapply(obj , function(x){

                           beta <- if(!is.null(x$est)){
                                       x$est
                                   }else{
                                       NA }  # fix here must return NA of correct length
                           N <- length(beta)
                           max <- ifelse(!is.null(x$max),
                                         x$max, NA )
                           
                           hess <-
                               try( -solve(x$hess), T  ) 

                           vbeta <- if( class(hess)=="try-error"){
                                        rep(NA,N)
                                    }else{
                                        diag(hess) 
                                    }

                           aic <- 2*length(beta) - 2*max

                           low <- beta - 1.96*sqrt(vbeta)
                           up <- beta + 1.96*sqrt(vbeta)

                           data.frame(beta=t(beta),
                                      vbeta=
                                          t( vbeta  ),
                                      low = t(low),
                                      up=t(up),
                                      max=t(max),
                                      aic = t(aic)
                                      )
                       }
                       )    )

    else  # LM
        
        res <- do.call("rbind",
                       lapply(obj , function(x){
                           
                           vbeta <- ifelse(x$vbeta<0, NA, x$vbeta)
                           aic <- ifelse(is.nan(x$aic), NA, x$aic)
                           low <- ifelse(is.nan(x$low), NA, x$low)
                           up <- ifelse(is.nan(x$up), NA, x$up)
                           max <- ifelse(is.nan(x$max), NA, x$max)
                           data.frame(beta= t(x$beta),
                                      vbeta= t( vbeta  ),                   
                                      low = t(low),
                                      up=t(up),
                                      max=t(max),
                                      aic = t(aic))
                       }
                       )  )  
    
    sumres <- NULL
    if(summary){


        MCerr <- apply(res, 2, sd, na.rm=TRUE)
        
        
        MCmean <- apply(res, 2, mean, na.rm=TRUE)

        sumres <- rbind(MCerr, #  MCmode,
                        MCmean )
        rownames(sumres) <- c("var", #  "mode",
                              "mean")
        
    }

    output <- list("mc"=as.data.frame(res),
                   "mcsum"=as.data.frame(sumres))
    output
}
























