Warning <- function(...)warning(..., call.=FALSE)
Stop <- function(...)stop(..., call.=FALSE)


recycle <- function(x, n){
  rep(x, ceiling(n / length(x)))[1:n]
}
