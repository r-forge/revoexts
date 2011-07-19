rsrBinXdf <- 
function( inFile, choices = 1L:2L, scale = c( 1, 1 ), 
          label = NULL, label.nm = NULL,
          numBreaksX = 100, numBreaksY = 100, ...)
{
    # parameter checks
    if( length( choices ) != 2L ){ stop( "length of choices must be 2" ) }
    if( is.numeric( choices ) ){
        nm   <- rxGetVarNames( inFile )[choices]
    }
    if( is.character( choices ) ){
        if( !all( choices %in% rxGetVarNames( inFile ) ) ){
            stop( gettextf( "The choices you specified does not exist in '%s'.", 
                            deparse( substitute( inFile ) ) ), 
                  domain = NA ) 
        }
        nm <- choices
    }
    
    form <- paste( "~", paste( nm,collapse = " + " ) )
    rangeData <- rxSummary( formula( form ), inFile, reportProgress = 0 )[["sDataFrame"]][c( "Name", "Min", "Max" )]
    rownames( rangeData ) <- rangeData[["Name"]]
    rangeData <- rangeData[-1L]
    x <- nm[1]
    y <- nm[2]
    breaksX <- seq( rangeData[x, "Min"]/scale[1], 
                    rangeData[x, "Max"]/scale[1], length = numBreaksX + 1 )
    breaksY <- seq( rangeData[y, "Min"]/scale[2], 
                    rangeData[y, "Max"]/scale[2], length = numBreaksY + 1 )
    numClass <- 1
    ClassLabels <- NULL

    if( !is.null( label ) ){
        if( rxGetInfoXdf( inFile )$numRows != rxGetInfoXdf( label )$numRows ){
            stop( "Label must have the same length as the data.", domain = NA )
        }
        # this only works for categorical variable!!!!
        if( is.null( label.nm ) ){ 
            ClassLabels <- rxGetVarInfoXdf( label )[[1]]$levels
        } else {
            labelSummary <- rxSummary( formula( paste( "~F(",label.nm,")" ) ),
                                       label )
            ClassLabels  <- levels( labelSummary$categorical[[1]][,1]) 
        }
        numClass  <- length( ClassLabels )
        labSource <- RxXdfData( label )
        rxOpen( labSource )
    }

    # Calculate the number of observations within each "tile"
    if( numClass == 1 ){
        dimension <- c( length( breaksX ) - 1, length( breaksY ) - 1 ) 
    } else {
        dimension <- c( length( breaksX ) - 1, length( breaksY ) - 1, numClass ) 
    }
    counts    <- array( 0, dimension )
    simSource <- RxXdfData( inFile )
    cat("Opening data \n")
    rxOpen( simSource )
    while( NROW( chunk <- rxReadNext( simSource ) ) > 0 ){
        if( !is.null( label ) ){ labchunk <- rxReadNext( labSource ) 
            for( i in 1:numClass ){
                # note right cut off must be inf because of precision
                counts[,,i] <-counts[,,i] +
                 unclass( table( cut( 
                  chunk[[x]][labchunk[[label.nm]] == ClassLabels[i]]/scale[1], 
                  breaks = c( breaksX[-length( breaksX )], Inf ), 
                  right = FALSE ),
                                 cut( 
                  chunk[[y]][labchunk[[label.nm]] == ClassLabels[i]]/scale[2], 
                  breaks = c( breaksY[-length( breaksY )], Inf ), 
                  right = FALSE ) ) )
            }
        } else {
            # note right cut off must be inf because of precision
            counts[,] <-counts[,] +
                 unclass( table( cut( chunk[[x]] / scale[1], 
                                      breaks = c( breaksX[-length( breaksX )], 
                                                  Inf ), 
                                      right = FALSE ),
                                 cut( chunk[[y]] / scale[2], 
                                      breaks = c( breaksY[-length( breaksY )], 
                                                  Inf ), 
                                      right = FALSE ) ) )
        }
} 
    cat("Closing data \n")
    rxClose( simSource )
    if( !is.null( label ) ){ 
        rxClose( labSource )
        dimnames( counts ) <- list( paste( "[", breaksX[- length( breaksX )], 
                                           ",", breaksX[- 1], "]" ),
                                    paste( "[", breaksY[- length( breaksY )], 
                                           ",", breaksY[- 1], "]" ), 
                                    ClassLabels )
    } else {
        dimnames( counts ) <- list( paste( "[", breaksX[- length( breaksX )], 
                                           ",", breaksX[- 1], "]" ),
                                    paste( "[", breaksY[- length( breaksY )], 
                                           ",", breaksY[- 1], "]" ) )
    }
    result <- list ( counts = counts, x.breaks = breaksX, y.breaks = breaksY, 
                     labels = ClassLabels, varname = nm )
    class( result ) <- "rsrBinnedCount"
    return( result )
}
