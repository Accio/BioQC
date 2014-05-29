## rankSumTestWithCorrelation function from the limma package (version 3.18.13)

## authors: Gordon Smyth and Di Wu, following Zar, JD (1999) Biostatistical Analysis 4th Edition
## used under GPL(>=2) license. The function has been slihtly modified to allow reporting results

rankSumTestWithCorrelation <- function (index, statistics, correlation = 0, df = Inf) 
{
    n <- length(statistics)
    r <- rank(statistics)
    r1 <- r[index]
    n1 <- length(r1)
    n2 <- n - n1
    U <- n1 * n2 + n1 * (n1 + 1)/2 - sum(r1)
    mu <- n1 * n2/2
    if (correlation == 0 || n1 == 1) {
        sigma2 <- n1 * n2 * (n + 1)/12
    }
    else {
        sigma2 <- asin(1) * n1 * n2 + asin(0.5) * n1 * n2 * (n2 - 
            1) + asin(correlation/2) * n1 * (n1 - 1) * n2 * (n2 - 
            1) + asin((correlation + 1)/2) * n1 * (n1 - 1) * 
            n2
        sigma2 <- sigma2/2/pi
    }
    TIES <- (length(r) != length(unique(r)))
    if (TIES) {
        NTIES <- table(r)
        prod <- sum(NTIES * (NTIES + 1) * (NTIES - 1))
        denom <- (n * (n + 1) * (n - 1))
        adjustment <- prod/denom
        sigma2 <- sigma2 * (1 - adjustment)
    }
    zlowertail <- (U + 0.5 - mu)/sqrt(sigma2)
    zuppertail <- (U - 0.5 - mu)/sqrt(sigma2)
    less <-pt(zuppertail, df = df, lower.tail = FALSE)
    greater <- pt(zlowertail, df = df)
    res <- c(U=U,
             mu=mu,
             n1=n1,
             n2=n2,
             sigma2=sigma2,
             r1sum=sum(r1),
             zlt=zlowertail,
             zut=zuppertail,
             less = less,
             greater = greater)
    return(res)
}
