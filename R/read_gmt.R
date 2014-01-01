read_gmt <- function(filename) {
  stopifnot(file.exists(filename))
  filename <- path.expand(filename)
  .Call("read_gmt", filename)
}
