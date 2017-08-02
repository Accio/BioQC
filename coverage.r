cov <- covr::package_coverage(quiet = F)
as.data.frame(cov)
missing <- zero_coverage(cov)
cdf <- as.data.frame(cov)
ccdf <- cdf
ccdf[, "value"] <- (cdf[, "value"] != 0)
um <- rle(cdf[, 2])
uml <- um[[1]]
pef1 <- c(um[[2]][1], mean(ccdf[1:um[[1]][1], "value"]))
ma <- matrix(pef1, nrow = 76, ncol = 2)
k <- 1
for (j in 1:76) {
  ma[j, ] <-  c(um[[2]][j], mean(ccdf[k:(um[[1]][j]+k-1), "value"]))
  k <- k + um[[1]][j]
}
b <- print(cov, group = "function", by = "expression")
print(cov, group = "filename", by = "expression")
print(cov, group = "function", by = "line")
print(cov, group = "filename", by = "line")
print(cov, by = "expression")
print(cov)


