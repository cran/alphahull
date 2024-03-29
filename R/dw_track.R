dw_track <-
function (x, y = NULL, eps, nps = 20000, sc = 100) 
{
    X <- xy.coords(x, y)
    dat <- data.frame(location_long = X$x, location_lat = X$y)
    colnames(dat) <- c("location_long", "location_lat")
    np <- dim(dat)[1]
    seg <- spatstat.geom::psp(dat[-np, 1], dat[-np, 2], dat[-1, 1], dat[-1, 
        2], window = spatstat.geom::owin(range(dat[, 1]), range(dat[, 2])))
    pseg <- spatstat.random::runifpointOnLines(nps, seg)
    xah <- c(dat[, 1], pseg$x)
    yah <- c(dat[, 2], pseg$y)
    xah <- jitter(xah) * sc
    yah <- jitter(yah) * sc
    dw <- dw(cbind(xah, yah), eps = eps * sc)
    aux <- dw
    na <- dim(aux)[1]
    if (na >= 1) {
        dfr <- list()
        for (j in 1:na) {
            c <- aux[j, 1:2]
            r <- aux[j, 3]
            v <- aux[j, 4:5]
            theta <- aux[j, 6]
            angles <- anglesArc(v, theta)
            seqang <- seq(angles[1], angles[2], length = 100)
            x <- c[1] + r * cos(seqang)
            y <- c[2] + r * sin(seqang)
            daux <- data.frame(x = x/sc, y = y/sc)
            dfr[[j]] <- geom_path(data = daux, aes(x = x, y = y), 
                col = 4)
        }
        return(dfr)
    }
    else {
        warning("Please, choose a different value of eps")
    }
}
