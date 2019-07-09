appendOneGmtList <- function(gmtList, newGmtList) {
  res <- gmtList
  res@.Data <- c(res@.Data, newGmtList@.Data)
  names(res) <- c(gsName(gmtList), gsName(newGmtList))
  return(res)
}

#' Append a GmtList object to another one
#' 
#' @param gmtList A \code{GmtList} object
#' @param newGmtList Another \code{GmtList} object to be appended
#' @param ... Further \code{GmtList} object to be appended
#' 
#' @return A new \code{GmtList} list, with all elements in the input appended in the given order
#' 
#' @examples 
#' test_gmt_file<- system.file("extdata/test.gmt", package="BioQC")
#' testGmtList1 <- readGmt(test_gmt_file, namespace="test1")
#' testGmtList2 <- readGmt(test_gmt_file, namespace="test2")
#' testGmtList3 <- readGmt(test_gmt_file, namespace="test3")
#' testGmtAppended <- appendGmtList(testGmtList1, testGmtList2, testGmtList3)
#' @export
appendGmtList <- function(gmtList, newGmtList, ...) {
  res <- appendOneGmtList(gmtList, newGmtList)
  newList <- list(...)
  for (ngs in newList) {
    res <- appendOneGmtList(res, ngs)
  }
  return(res)
}
