
# PCA fit for rxCovCor object
princomp.rxCovCor <- function( x, ... ){
    # check to see if the rxCovCor object is a covariance or correlation
    if( !inherits( x, "rxCovCor" ) ){
        princomp( x )
    } else {
        if( x$params$Type == "cov" ){
            corFlag <- FALSE
        } else if ( x$params$Type == "cor" ){
            corFlag <- TRUE
        } else {
            stop( gettextf( "x must be rxCovCor matrix with type cov or cor" ) )
        }
        fit <- princomp( cor    = corFlag, 
                         scores = FALSE, 
                         covmat = list( cov    = x$CovCor, 
                                        n.obs  = x$valid.obs, 
                                        center = x$Means ), 
                         subset = subset, ... )
        fit$center   <- x$Means
        res          <- list( "fit"    = fit, 
                              "scores" = NULL )
        class( res ) <- "rsrPrincomp"
        return( res )
    }
}

# PCA fit for xdf object
rsrPrincomp <- function( formula,
                        inFile, 
                        score          = FALSE, 
                        outFile        = NULL,  
                        append         = "none",
                        overwrite      = FALSE, 
                        blocksPerRead  = rxGetOption( "blocksPerRead" ),
                        reportProgress = rxGetOption( "reportProgress" ),
                        ... )
{
    if( reportProgress > 0 ){
        cat( "Calculating Covariance Matrix:\n" )
    }
    CovMat <- rxCovCor( formula = formula, data = inFile,   
                        blocksPerRead  = blocksPerRead,
                        reportProgress = reportProgress )
    if( reportProgress > 0 ){
        cat( "Fitting the model:\n" )
    }
    res    <- princomp.rxCovCor( x = CovMat, ... )
    if( score == TRUE ){
        if( reportProgress > 0 ){
            cat( "Calculating Scores:\n" )
        }
        if( is.null( outFile ) ) { 
            outFile <- inFile
            append  <- TRUE
        }
        res$scores <- outFile
        rsrPrincompPredictXdf( res$fit, 
                               inFile         = inFile, 
                               outFile        = outFile, 
                               append         = append,
                               overwrite      = overwrite,
                               blocksPerRead  = blocksPerRead,
                               reportProgress = reportProgress
        )
    }
    return( res )
}

# print 
print.rsrPrincomp   <- function( x ){
    if( inherits( x, "rsrPrincomp" ) ){
        print( x$fit )
    } else {
        print( x )
    }
}

# summary
summary.rsrPrincomp <- function( x ){
    if( inherits( x, "rsrPrincomp" ) ){
        summary( x$fit )
    } else {
        summary( x )
    }
}
scores <- function (x, ...){ 
    UseMethod( "scores", x ) # x is object
}

# default method
scores.default <- function( x ){
    if( !is.null( x$scores ) ){
        return( x$scores )
    } else {
        return( x )
    }
}
# scores
scores.rsrPrincomp  <- function( x ){
    if( inherits( x, "rsrPrincomp" ) ){
        x$scores
    } else {
        scores( x ) 
    }
}

# loadings
loadings <- function (x, ...) {UseMethod( "loadings", x )}

loadings.rsrPrincomp<- function( x ){
    if( inherits( x, "rsrPrincomp" ) ){
        x$fit$loadings
    } else {
        loadings( x ) 
    }
}

rsrComponentPlot<-function( x, componentName = "Comp.1", ... ){ #,,... ){
    if( inherits( x, "rsrPrincomp" ) ){
        if( is.null( x$scores ) ){
            stpt("Principal Component score is not calculated.")
        }
        possibleNames <- rxGetVarNames( x$scores )
        if( !all( componentName %in% possibleNames ) ){
            stop( "x does not have component with specified name" )
        }
        idx <- 1:x$fit$n.obs
        rxLinePlot( formula( paste( componentName, " ~ idx" ) ), x$scores, ... )
    } else {
        stop( "only works on rsrPrincomp object" )
    }
}
# plot 
biplot.rsrPrincomp    <- function( x, choices = 1L:2L, label = NULL, label.nm = NULL, 
                                   numBreaksX = 100, numBreaksY = 100, 
                                   expand = 1, scale = 1, pc.biplot = FALSE, 
                                   arrow.len = 0.1, alpha = 1,
                                   main = NULL, xlab, ylab, xlim, ylim, ... ){
    if( inherits( x, "rsrPrincomp" ) ){
        nm <- dimnames( x$fit$loadings )[[2]]
        if( is.numeric( choices ) ){
            choice.nm  <- nm[choices]
            loading.sc <- x$fit$loadings[,choices]
        } else {
            if( !all( choices %in%  nm ) ){
                stop( "choices you specified are not incorrect" ) 
            } else {
                choice.nm  <- choices
                loading.sc <- x$fit$loadings[,choices]
            }
        }
        if( is.null( x$scores ) ){ stop( "object does not have scores" ) }
        lam <- x$fit$sdev[choices]
        if( is.null( n <- x$fit$n.obs ) ) n <- 1
        lam <- lam * sqrt(n)
        if( scale < 0 || scale > 1 ) warning( "'scale' is outside [0, 1]" )
        if( scale != 0 ) lam <- lam^scale else lam <- 1
        if( pc.biplot )  lam <- lam / sqrt( n )
        cat("Binning the scores \n")
        binned.count <- rsrBinXdf( x$scores, scale = lam,
                                    label = label, label.nm = label.nm, 
                                    choices = choice.nm,
                                    numBreaksX = numBreaksX,
                                    numBreaksY = numBreaksY )
        op <- par( pty = "s" )
        if( !is.null( main ) )
            op <- c(op, par( mar = par( "mar" ) + c( 0, 0, 1, 0 ) ) )

        l.sc.scaled <- loading.sc * rep( lam, each = dim( loading.sc )[1],
                                         byrow = FALSE )
        unsigned.range <- function( x ){
            c( -abs( min( x, na.rm = TRUE ) ), abs( max( x, na.rm = TRUE ) ) )
        }
        rangx1 <- unsigned.range( binned.count$x.breaks ) 
        rangx2 <- unsigned.range( binned.count$y.breaks )
        rangy1 <- unsigned.range( l.sc.scaled[, 1L] )
        rangy2 <- unsigned.range( l.sc.scaled[, 2L] )
    
        if( missing( xlim ) && missing( ylim ) ){
            xlim <- rangx1 
            ylim <- rangx2 #<- range( rangx1, rangx2 )
        } else if( missing( xlim ) ){ 
            xlim <- rangx1
        } else if( missing( ylim ) ){ 
            ylim <- rangx2
        }
        ratio <- max( rangy1 / rangx1, rangy2 / rangx2 ) / expand
        if( missing( xlab ) ){
            xlab = choice.nm[1]
        } 
        if( missing( ylab ) ){
            ylab= choice.nm[2]
        }
        cat("Plotting grid \n")
        df.scores <- rsrGridMap( binned.count, alpha = alpha, 
                                xlab = xlab, ylab = ylab,
                                xlim = xlim, ylim = ylim, ... )
        cat("Plotting lines \n")
        par( new = TRUE )
        plot( l.sc.scaled, axes = FALSE, type = "n", 
                xlim = xlim * ratio, ylim = ylim * ratio,
                xlab = "", ylab = "", ... )
        axis( 3, ... )
        axis( 4, ... )
        arrows( 0, 0, l.sc.scaled[,1] * 0.8, l.sc.scaled[,2] * 0.8, 
                length = arrow.len )
        text( l.sc.scaled, labels = dimnames( l.sc.scaled )[[1]] )
    } else {
        biplot( x )
    }
    invisible( df.scores )
}
plot.rsrPrincomp <- function( x, main = deparse( substitute( x ) ), ...){
  screeplot( x$fit, main = main, ... )
}

# prediction
predict.rsrPrincomp <- function( object, 
                                newdata, 
                                outFile   = NULL,
                                append    = "none",
                                overwrite = FALSE,
                                ... )
{
    if( inherits( object, "rsrPrincomp" ) ){
        if( missing( newdata ) ) return( object$scores )
        if( is.null( outFile ) ){ 
            outFile = newdata
            append  = TRUE 
        }
        rsrPrincompPredictXdf(  object$fit, 
                                inFile    = newdata, 
                                outFile   = outFile, 
                                append    = append,
                                overwrite = overwrite,
                                ...)
        return( outFile )
    } else {
        predict( object, newdata )
    }
}

# internal function for data %*% loading operation
rsrPrincompPredictXdf <-
function( princomp,                  # an object of class princomp
          inFile,                    # Rest of args same as rxDataStepXdf
          outFile,
          rowSelection   = NULL,
          transforms     = NULL,
          append         = "none",
          overwrite      = FALSE,
          removeMissings = FALSE,
          computeLowHigh = TRUE,
          blocksPerRead  = rxGetOption("blocksPerRead"),
          reportProgress = rxGetOption("reportProgress"))
{
    # Check princomp input object type
    if( !inherits( princomp, "princomp" ) ){
        stop( "'princomp' must be of class 'princomp'" )
    }
    # Prepare expression-based arguments
    rowSelection <- substitute( rowSelection )
    transforms   <- substitute( transforms )

    # Using rxDataStepXdf nomenclature, the variables required for generating
    # the scores are specified the transformVars argument.
    transformVars <- rownames( princomp$loadings )
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
    newTransformFunc <- function( princomp )
    {
        # Create a closure variable containing the number of columns
        p  <- NCOL( princomp$loadings )
        # Create function with one argument dataList, which is a list object
        function( dataList )
        {
            nm      <- rownames( princomp$loadings )
            xMatrix <- do.call( cbind, dataList[nm] )
            # Ensure xMatrix has the appropriate number of columns
            if( NCOL( xMatrix ) != p )
                stop( "'inFile' and 'transforms' do not produce the appropriate columns" )
            scores  <- scale( xMatrix, princomp$center, princomp$scale ) %*% princomp$loadings
            colnames( scores ) <- colnames( princomp$loadings )
            c( dataList, as.data.frame( scores ) )
        }
    }
    # Create a transformFunc that contains princomp input in its closure
    transformFunc <- newTransformFunc( princomp )

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
