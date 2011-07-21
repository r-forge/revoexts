
rsrLoadingPlot<- function( loading, main = "Loading matrix", ylabels, xlabels,
                           mar = c( 0, 5, 4, 0 ), ... )
{
    par( mar = mar, oma =c( 0, 0, 1, 0 ) )
    image( t( loading[dim( loading )[1]:1,] ),   
           col = rev( brewer.pal( 11, "RdBu" ) ) , axes = FALSE, ... )
    if( missing( ylabels ) && !is.na( dimnames( loading )[[1]] ) ){
        ylabels = dimnames( loading )[[1]]
    }
    if( missing( xlabels ) && !is.na( dimnames( loading )[[2]])){
        xlabels = dimnames( loading )[[2]]
    }
    if( !missing( ylabels ) && !is.null( ylabels ) ){
        if( dim( loading )[1]!= length( ylabels ) ){ stop("ylabel must be as long as 1st dimention of loading") }
        axis( 2, at = seq( 0 , 1, length.out= dim(loading)[1]), labels = rev( ylabels ), 
             las = 2, tick = FALSE, line = -1 )
    }
    if( !missing( xlabels ) && !is.null( xlabels ) ){
        if( dim( loading )[2]!= length( xlabels ) ){ stop("xlabel must be as long as 2st dimention of loading") }
        axis( 3, at = seq( 0 , 1, length.out= dim(loading)[2]), labels = xlabels ,
             las = 3, tick = FALSE, line = -1 )
    }
    title( main = main, outer = TRUE )
    #imageplot( matrix = loading, main = main, ... )
}

pairs.loadings<- function( loading, include=1:3, col = "black", labels = NULL, compnames = NULL, ... )
{
    clen = length( include )
    ldim = dim( loading )
    dnam = dimnames( loading )
    if( length( col ) == 1 ){ col = rep( col, ldim[1] ) }
    
    if( is.null( labels ) ){
        if(!is.null( dnam[[1]] ) ) { 
            labels = dnam[[1]] 
        } else {
            labels = 1:ldim[1]
        }
    }
    if( is.null( compnames ) ){
        if(!is.null( dnam[[2]] ) ){
            compnames = dnam[[2]][include]
        } else {
            compnames = include
        }
    }
    if( length( dnam[[1]] ) != length( col ) ){ 
        stop("col must be either one color or as long as the number of rows")
    }
    par( mfrow = c( clen, clen ) )
    par( mar = c( 2, 2, 1, 1 ), mgp = c( 2, 0.7, 0 ), tck = -0.01 )
     for( j in 1:clen ){
        for( k in 1:clen ){
            if( k != j ){
                plot( loading[,include[k]], loading[,include[j]], type="n", 
                       xlab="", ylab="", ... )
                abline( 0, 0, col = "grey89" )
                abline( v = 0,col = "grey89" )
                for( i in 1:ldim[1] ){
                    text(loading[i,include[k]], loading[i,include[j]], labels[i])#, col = col[i] )
                }
                #dev.off()
            } else {
                plot.new()#(0,0,type="n")
                text( 0.5, 0.5, compnames[k] )
            }
        }
    }
}

