## --- Trace

Trace = function(A) {
  if (is.matrix(A)) {
    if (nrow(A) != ncol(A)) {
      stop('Matrix should be square in Trace()')
    }
  }
  return(sum(diag(A)))
}