
# PCA fit for rxCovCor object
princomp.rxCovCor <- function( CovCorMat, ... ){
    # check to see if the rxCovCor object is a covariance or correlation
    if( !inherits( CovCorMat, "rxCovCor" ) ){
        princomp( CovCorMat )
    } else {
        if( CovCorMat$params$Type == "cov" ){
            corFlag <- FALSE
        } else if ( CovCorMat$params$Type == "cor" ){
            corFlag <- TRUE
        } else {
            stop( gettextf( "CovCorMat must be rxCovCor matrix with type cov or cor" ) )
        }
        fit <- princomp( cor    = corFlag, 
                         scores = FALSE, 
                         covmat = list( cov    = CovCorMat$CovCor, 
                                        n.obs  = CovCorMat$valid.obs, 
                                        center = CovCorMat$Means ), 
                         subset = subset, ... )
        fit$center   <- CovCorMat$Means
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
    CovMat <- rxCovCor( formula = formula, data = inFile,   
                        blocksPerRead  = blocksPerRead,
                        reportProgress = reportProgress )
    res    <- princomp.rxCovCor( CovCorMat = CovMat, ... )
    if( score == TRUE ){
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
scores <- function (x, ...)  {
  if( is.null( attr( x, "class" ) ) ){
    return( x )
  }
  else  UseMethod( "scores", x ) # x is object
}
# Step 3
scores.default <- function( x ){
    if(!is.null(x$scores)){
        return(x$scores)
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

loadings <- function (x, ...)  {
  if( is.null( attr( x, "class" ) ) ){
    return( x )
  }
  else  UseMethod( "loadings", x ) # x is object
}
# loadings
loadings.rsrPrincomp<- function( x ){
    if( inherits( x, "rsrPrincomp" ) ){
        x$fit$loadings
    } else {
        loadings( x ) 
    }
}

rsrPlotScores<-function( x, componentName,xlab= "Index", main="Principal Component Score", ... ){
    if( inherits( x, "rsrPrincomp" ) ){
        if( is.null( x$scores ) ){
            stpt("Principal Component score is not calculated.")
        }
        idx <- 1:x$fit$n.obs
        rxLinePlot( formula( paste( componentName, " ~ idx" ) ), x$scores, 
             xlab = xlab, main = main, ... )
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
    
       # if( missing( xlim ) && missing( ylim ) ){
       #     xlim <- rangx1 
       #     ylim <- rangx2 #<- range( rangx1, rangx2 )
       # } else if( missing( xlim ) ){ 
       #     xlim <- rangx1
       # } else if( missing( ylim ) ){ 
       #     ylim <- rangx2
       # }
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
        #par( new = TRUE )
        #plot( l.sc.scaled, axes = FALSE, type = "n", 
        #        xlim = xlim * ratio, ylim = ylim * ratio,
        #        xlab = "", ylab = "", ... )
        #axis( 3, ... )
        #axis( 4, ... )
        #arrows( 0, 0, l.sc.scaled[,1] * 0.8, l.sc.scaled[,2] * 0.8, 
        #        length = arrow.len )
        #text( l.sc.scaled, labels = dimnames( l.sc.scaled )[[1]] )
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
        rsrPrincompPredictXdf(   object$fit, 
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

rsrLoadingPlot<- function( loading, main = "Loading matrix", rangex, rangey, xbin, ybin, 
                        zrange=c( -1, 1 ), zbin = 11, digits=2, cex.col=0.7, cex.var=0.5,
                        xlabels = NULL, ylabels = NULL, color = rev( brewer.pal( zbin, "RdBu" ) )
                        ){
                            
    z.breaks = seq( zrange[1], zrange[2], length.out = zbin + 1 )
    z = array( as.double( cut( loading, breaks = z.breaks, labels = 1:zbin ) ), dim( loading ) )
    d = dim( loading )
    y <- 1:( d[1] + 1 )
    x <- 1:( d[2] + 1 )
    if( length( color ) < zbin ){ stop("color must be same length as zbin") }
    pal <- color
    layout( matrix( c( 2, 1 ), 1, 2, byrow = FALSE ), c( 10.5, 1.5 ) )
    par( mar = c( 0.5, 0.1, 2, 0.1 ), pty = "m" )
    plot( c( 0, 1 ), c( min( z.breaks ), max( z.breaks ) ), type = "n",
        bty = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n" )
    for( i in 2:( length( z.breaks ) ) ) {
        rect( xleft = 0.5, ybottom = z.breaks[i - 1], 
              xright = 1, ytop = z.breaks[i], col = pal[i - 1] )
        text( x = 0.45, y = z.breaks[i - 1], 
              labels = format(round(z.breaks[i - 1], digits)), 
              cex = cex.col, adj = 1, xpd = TRUE)
    }
    rect( xleft = 0.5, ybottom = z.breaks[length( z.breaks )], 
          xright = 1, ytop = z.breaks[length( z.breaks )], col = pal[length( pal )])
    text( x = 0.45, y = z.breaks[length( z.breaks )], 
          labels = format(round(z.breaks[length( z.breaks )], digits ) ), 
          cex = cex.col, adj = 1, xpd = TRUE )
    par( mar=c( 1, 4, 4, 1 ) )
    plot( range( x ),range( y ), axes = FALSE, type = "n",
          xlim = range( x ), ylim = range( y ),
          xaxs = "i", yaxs = "i", xlab = "", ylab = "" )
    for(i in 1:d[2]){
        for(j in d[1]:1){
            rect( x[i], y[j], x[i+1], y[j+1], col = pal[z[j,i]], 
                  border = pal[z[j,i]] )
        }
}
    if( is.null( xlabels ) ){ xlabels <- 1:d[2] }
    if( is.null( ylabels ) ){ ylabels <- 1:d[1] }
    axis( 2, at = ( ( 1:d[1] ) + 0.5 ), labels = rev( ylabels ), 
             las = 2, tick = FALSE, line = -1 )
    axis( 3, at = ( ( 1:d[2] )+ 0.5 ), labels = xlabels ,
             las = 3, tick = FALSE, line = -1 )
    title( main = main )
    invisible( z )
}