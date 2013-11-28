inahull <-
function (ahull.obj, p) 
{
    compl <- ahull.obj$complement
    halfpl <- which(compl[, "r"] < 0)
    n.halfpl <- length(halfpl)
    ball <- which(compl[, "r"] > 0)
    n.ball <- length(ball)
    in.compl <- FALSE
    if (n.halfpl >= 1) {
        h <- 1
        while ((h <= n.halfpl) & in.compl == FALSE) {
            sig = compl[halfpl[h], 3]
            a = compl[halfpl[h], 1]
            b = compl[halfpl[h], 2]
            if (sig <= -3) {
                if (p[1] > a) {
                  if (sig == -3) {
                    in.compl <- TRUE
                  }
                }
                else if (p[1] < a) {
                  if (sig == -4) {
                    in.compl <- TRUE
                  }
                }
            }
            else {
                if (p[2] > a + b * p[1]) {
                  if (sig == -1) {
                    in.compl <- TRUE
                  }
                }
                else if (p[2] < a + b * p[1]) {
                  if (sig == -2) {
                    in.compl <- TRUE
                  }
                }
            }
            h <- h + 1
        }
    }
    if (in.compl == FALSE) {
        k <- 1
        while ((k <= n.ball) & in.compl == FALSE) {
            r = compl[ball[k], 3]
            c1 = compl[ball[k], 1]
            c2 = compl[ball[k], 2]
            d <- sqrt((p[1] - c1)^2 + (p[2] - c2)^2)
            if (d < r) {
                in.compl <- TRUE
            }
            k <- k + 1
        }
    }
    return(in.ahull = !in.compl)
}
