
color.gradient<-
function( n, hue, alpha = 1 ) 
{
    hsv( h     = hue, 
         s     = seq( from = 0, to = 1, length.out = n ), 
         v     = 1, 
         alpha = alpha )
}

# The col that one can specify is one that corresponds to hue for full saturation and value. 
# One may specify other colors but it will come out as its full saturation and value color.
# List of the colors that comes out correctly are 
# "blue", "blue1", "chartreuse", "chartreuse1", "cyan", "cyan1", "darkorange", 
# "darkorange1", "deepskyblue", "deepskyblue1", 
# "gold", "gold1", "green", "green1", "magenta", "magenta1", "orange", 
# "orange1", "orangered", "orangered1", "red", "red1", "springgreen", 
# "springgreen1", "turquoise1", "yellow", "yellow1"
# by setting parameters for jitter(), "x.factor", "y.factor", "x.amount", and 
# "y.amount", to have no jittering, this function becomes an approximate 
# version of heatmap. 

rsrGridMap <- 
function( count, alpha = 1, main = NULL, xlab, ylab, col,
            mar = c( 3.5, 3.5, 2, 2 ), mgp = c( 2, 0.7, 0 ), tck = -0.01,
            x.factor = 1, y.factor = 1, x.amount = NULL, y.amount = NULL, xlim = NULL, ylim = NULL, ... )
{
    par ( mar = mar, mgp = mgp, tck = tck )
    if( inherits( count, "rsrBinnedCount" ) ){
        count.array <- count$counts 
        dm  <- dim( count.array )
        ldm <- length ( dm ) 
        x.center <- ( count$x.breaks[-1] + 
                      count$x.breaks[-length( count$x.breaks )] ) / 2
        y.center <- ( count$y.breaks[-1] + 
                      count$y.breaks[-length( count$y.breaks )] ) / 2
        x.side   <- min( count$x.breaks[-1] - 
                         count$x.breaks[-length( count$x.breaks )] )
        y.side   <- min( count$y.breaks[-1] - 
                         count$y.breaks[-length( count$y.breaks )] ) 
        if( is.null( count$labels ) ){
            dimnames( count.array ) <- list( x.center, y.center )
        } else {
            dimnames( count.array ) <- list( x.center, y.center, 
                                     as.character(1:length( count$labels ) ) )
        }
        if( missing( xlab ) ) { xlab = count$varname[1]}
        if( missing( ylab ) ) { ylab = count$varname[2]}
        maxcnt <- max( as.vector( count.array ) )
        dfz <- as.data.frame( as.table( count.array ) )
        if ( ldm == 2 ){
            dfz <- dfz[dfz[,3]!=0,]
            # random ordering to avoid one color always being on the top
            dfz <- dfz[sample(dim(dfz)[1]),] 
            x   <- as.double( levels( dfz[,1] ) )[dfz[,1]]
            y   <- as.double( levels( dfz[,2] ) )[dfz[,2]]
            cls <- rep( 1, length( x ) )
            cnt <- as.double( dfz[,3] )
        } else if ( ldm == 3 ){
            dfz <- dfz[dfz[,4]!=0,]
            # random ordering to avoid one color always being on the top
            dfz <- dfz[sample(dim(dfz)[1]),] 
            x   <- as.double( levels( dfz[,1] ) )[dfz[,1]]
            y   <- as.double( levels( dfz[,2] ) )[dfz[,2]]
            cls <- as.double( dfz[,3] )
            cnt <- as.double( dfz[,4] )
        }
        ucl     <- length( unique( cls ) )
        pal     <- matrix( NA, maxcnt, ucl ) # color pallet
        if( missing( col ) ){
            hu  <- ( 1:ucl ) / ucl
        } else {
            if( length( col ) == ucl ){ 
                hu <- rgb2hsv( col2rgb( col ) )[1,]
            } else if( length( col ) == 1 ){
                hu <- rep( col, ucl )    
            } else {
                stop("col must be specified for each subset")
            }
        }
        for( i in 1:ucl ){
            pal[,i] <- color.gradient( maxcnt + 1, hue = hu[i], alpha = alpha )[-1]
        }
        if( ucl == 1){
            jx <- x
            jy <- y
        } else {
            jx <- jitter( x, factor = x.factor, amount = x.amount )
            jy <- jitter( y, factor = y.factor, amount = y.amount )
        }
        if( missing( xlim ) ){ 
            xlim <- range(jx) + x.side*c(-1,1)
        }
        if( missing( ylim ) ){ 
            ylim <- range(jy) + y.side*c(-1,1)
        }
        par( bg="gray89" )
        plot( jx, jy, type = "n", pch = 19, 
              col = pal[as.matrix( cbind( cnt, cls ) )], 
              xlab= xlab, ylab=ylab, xlim=xlim, ylim=ylim )
        par( bg="white" )
        symbols( jx, jy, add=T, inches =F,
                 rectangles = cbind( rep( x.side,length(jx)), 
                                     rep(y.side,length(jy) ) ),
                 fg = pal[as.matrix( cbind( cnt, cls ) )],
                 bg = pal[as.matrix( cbind( cnt, cls ) )] )
        invisible( cbind( x, y, cls, cnt ) )
    } else {
        stop( "count must be an object with class rsrBinnedCount" )
    }
}
