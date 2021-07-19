#' Prettify default signature names
#' @param names Character strings, signature names
#' @param includeNamespace Logical, whether the namespace of the signatures should be included
#' @return Character strings, pretty signature names
#' @examples
#' sig <- readCurrentSignatures()
#' prettyNames <- prettySigNames(names(sig))
#' @export
prettySigNames <- function(names, includeNamespace=TRUE) {
  res <- gsub("_[0-9\\.]*_[0-9]$", "", names)
  res <- gsub("_", " ", res)
  res <- gsub("\\s+", " ", res)
  nss <- c("NR", "Roche", "Linsley", "Fantom5")
  for (ns in nss) {
    isNS <- grepl(ns, res)
    fmt <- paste0(ns, ".*$")
    if(includeNamespace) {
      res[isNS] <- gsub(fmt, paste0("[", ns, "]"), res[isNS])
    } else {
      res[isNS] <- gsub(fmt, "", res[isNS])
    }
  }
  return(res)
}
