# FA fit for rxCovCor object
factanal.rxCovCor <- function( x, factors = 2, start = NULL, 
                               rotation = "varimax", control = NULL, ... )
{ 
    # check to see if the rxCovCor object is a covariance or correlation
    if( !inherits( x, "rxCovCor" ) ){
        stop( "'x' must be of class 'rxCovCor'" )
    }

    n.obs <- x$valid.obs
    fit   <- factanal(  scores   = "none", 
                        covmat   = x$CovCor, 
                        n.obs    = n.obs,
                        factors  = factors, 
                        start    = start, 
                        rotation = rotation, 
                        control  = control, 
                        ... )
    res   <- list( "fit" = fit, "scores" = NULL )
    class( res ) <- "rsrFactanal"
    return( res )
}

# FA fit for xdf object
rsrFactanal <- function( formula,
                         factors,
                         inFile, 
                         outFile   = NULL,  
                         append    = "none",
                         overwrite = FALSE, 
                         scores    = c( "none", "regression", "Bartlett" ),
                         start     = NULL, 
                         rotation  = "varimax",
                         blocksPerRead  = rxGetOption( "blocksPerRead" ),
                         reportProgress = rxGetOption( "reportProgress" ),
                         ... )
{
    scores <- match.arg( scores )
    if( reportProgress > 0 ){
        cat( "Calculating Covariance Matrix:\n" )
    }
    CovMat <- rxCovCor( formula        = formula, 
                        data           = inFile,
                        blocksPerRead  = blocksPerRead,
                        reportProgress = reportProgress )
    if( reportProgress > 0 ){
        cat( "Fitting the model:\n" )
    }
    res    <- factanal.rxCovCor( x        = CovMat, 
                                 factors  = factors ,
                                 start    = start, 
                                 rotation = rotation,
                                 ... )
    if( scores != "none" ){
        if( reportProgress > 0 ){
            cat( "Calculating Factor Scores:\n" )
        }
        if( is.null( outFile ) ) { 
            outFile <- inFile
            append  <- TRUE 
        }
        res$scores <- outFile
        rsrFactanalPredictXdf( res$fit, 
                               scores         = scores,
                               inFile         = inFile, 
                               outFile        = outFile, 
                               append         = append,
                               overwrite      = overwrite, 
                               blocksPerRead  = blocksPerRead,
                               reportProgress = reportProgress,
                               ... )
    } 
    return( res )
}

# print 
print.rsrFactanal <- function( x ){
    if( inherits( x, "rsrFactanal" ) ){
        print( x$fit )
    } else {
        print( x )
    }
}
# summary
summary.rsrFactanal <- function( x ){
    if( inherits( x, "rsrFactanal" ) ){
        summary( x$fit )
    } else {
        summary( x )
    }
}
# scores
scores.rsrFactanal <- function( x ){
    
    if( inherits( x, "rsrFactanal" ) ){
        x$scores
    } else {
        scores( x )
    }
}
# loadings
loadings.rsrFactanal<- function( x ){

    if( inherits( x, "rsrFactanal" ) ){
        x$fit$loadings
    } else {
        loadings( x )
    }
}
# plot 
scoreplot <- function (x, ...) {UseMethod( "scoreplot", x )}

scoreplot.rsrFactanal   <- 
function( x, choices = 1L:2L, label = NULL, label.nm, alpha = 1,
          numBreaksX = 100, numBreaksY = 100,
          main = NULL, xlab, ylab, xlim, ylim, ... ){
    if( inherits( x, "rsrFactanal" ) ){
        nm <- dimnames( x$fit$loadings )[[2]][choices]
        if( is.numeric( choices ) ){
            choice.nm  <- nm[choices]
            loading.sc <- x$fit$loadings[,choices]
        } else {
            if( !all( choices %in%  nm) ){
                stop( "choices you specified are not incorrect" ) 
            } else {
                choice.nm  <- choices
                loading.sc <- x$fit$loadings[,choices]
            }
        }
        if( is.null( x$scores ) ){ stop( "object does not have scores" ) }
        binned.count <- rsrBinXdf( x$scores, 
                                    label = label, label.nm =label.nm, 
                                    choices = choice.nm,
                                    numBreaksX = numBreaksX,
                                    numBreaksY = numBreaksY )
        op <- par( pty = "s" )
        if( !is.null( main ) )
            op <- c(op, par( mar = par( "mar" ) + c( 0, 0, 1, 0 ) ) )

        unsigned.range <- function( x ){
            c( -abs( min( x, na.rm = TRUE ) ), abs( max( x, na.rm = TRUE ) ) )
        }
        rangx1 <- unsigned.range( binned.count$x.breaks ) 
        rangx2 <- unsigned.range( binned.count$y.breaks )

        if( missing( xlim ) && missing( ylim ) ){
            xlim <- ylim <- rangx1 <- rangx2 <- range( rangx1, rangx2 )
        } else if( missing( xlim ) ){ 
            xlim <- rangx1
        } else if( missing( ylim ) ){ 
            ylim <- rangx2
        }
        if( missing( xlab ) ){
            xlab = choice.nm[1]
        } 
        if( missing( ylab ) ){
            ylab= choice.nm[2]
        }
        df.scores <- rsrGridMap( binned.count, alpha = alpha, 
                                xlab = xlab, ylab = ylab,
                                xlim = xlim, ylim = ylim, ... )
    } else {
        stop("x must be of class rsrFactanal")
    }
    invisible( df.scores )
}
rsrFactorPlot<-function( x, factorName = "Factor1", ... ){
    if( inherits( x, "rsrFactanal" ) ){
        if( is.null( x$scores ) ){
            stpt("Factor score is not calculated.")
        }
        idx <- 1:x$fit$n.obs
        rxLinePlot( formula( paste( factorName, " ~ idx" ) ), 
                    x$scores, ... )
    } else {
        stop( "only works on rsrFactanal object" )
    }
}

# internal function for data %*% loading operation
rsrFactanalPredictXdf <-
function( factanal,                  # an object of class factanal
          scores, 
          inFile,                    # Rest of args same as rxDataStepXdf
          outFile,
          rowSelection   = NULL,
          transforms     = NULL,
          append         = "none",
          overwrite      = FALSE,
          removeMissings = FALSE,
          computeLowHigh = TRUE,
          blocksPerRead  = rxGetOption( "blocksPerRead" ),
          reportProgress = rxGetOption( "reportProgress" ) )
{
    # Check princomp input object type
    if( !inherits( factanal, "factanal" ) ){
        stop( "'factanal' must be of class 'factanal'" )
    }
    # Prepare expression-based arguments
    rowSelection <- substitute( rowSelection )
    transforms   <- substitute( transforms )
    # Using rxDataStepXdf nomenclature, the variables required for generating
    # the scores are specified the transformVars argument.
    transformVars <- rownames( factanal$loadings )
    if( !is.null( transformVars ) ){
        # Get variable names from file
        possibleNames <- rxGetVarNames( inFile )
        # Append new variables names from transforms
        if( !is.null( transforms ) && 
            is.expression( transforms ) &&
            length( transforms ) == 1L && 
            length( transforms[[1L]] ) > 0L &&
            deparse( transforms[[1L]][[1L]] ) == "list" ){
            possibleNames <- c( possibleNames, names( transforms )[[1L]][-1L] )
        }
        if(!all( transformVars %in% possibleNames ) ){
            stop( "'inFile' does not have named columns matching one or more of the original columns" )
        }
    }
    newTransformFunc <- function( factanal , scores)
    {
        # Create a closure variable containing the number of columns
        p  <- NROW( factanal$loadings )
        # Create function with one argument dataList, which is a list object
        function( dataList )
        {
            nm      <- rownames( factanal$loadings )
            xMatrix <- do.call( cbind, dataList[nm] )
            # Ensure xMatrix has the appropriate number of columns
            if( NCOL( xMatrix ) != p )
                stop( "'inFile' and 'transforms' do not produce the appropriate columns" )
                
            if( scores != "none" ){
                Lambda <- factanal$loadings
                zz <- scale( xMatrix, TRUE, TRUE )
                switch(scores,
                       regression = {
                           sc <- zz %*% solve( factanal$correlation, Lambda )
                           if( !is.null( Phi <- attr( Lambda, "covariance" ) ) ){
                               sc <- sc %*% Phi
                           }
                       },
                       Bartlett = {
                            d   <- 1 / factanal$uniquenesses
                            tmp <- t( Lambda * d )
                            sc  <- t( solve( tmp %*% factanal$loadings, tmp %*% t( zz ) ) )
                       })
                rownames( sc ) <- rownames( xMatrix )
                colnames( sc ) <- colnames( Lambda )
                score <- sc
            }
            c( dataList, as.data.frame( score ) )
        }
    }
    # Create a transformFunc that contains princomp input in its closure
    transformFunc <- newTransformFunc( factanal, scores )

    # Call rxDataStepXdf, do.call is used here for safe processing of
    # rowSelection and transforms arguments
    do.call( rxDataStepXdf,
             list( inFile         = inFile,
                   outFile        = outFile,
                   rowSelection   = rowSelection,
                   transforms     = transforms,
                   transformFunc  = transformFunc,
                   transformVars  = transformVars,
                   append         = append,
                   overwrite      = overwrite,
                   removeMissings = removeMissings,
                   computeLowHigh = computeLowHigh,
                   blocksPerRead  = blocksPerRead,
                   reportProgress = reportProgress ) )
}
