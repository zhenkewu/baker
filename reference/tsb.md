# generate stick-breaking prior (truncated) from a vector of random probabilities

generate stick-breaking prior (truncated) from a vector of random
probabilities

## Usage

``` r
tsb(u)
```

## Arguments

- u:

  a vector of probabilities, with the last element 1.

## Value

a vector of the same length as u; sum to 1.

## Examples

``` r
oldpar <- graphics::par(mfrow=c(3,3),oma=c(0,1,5,0),
   mar=c(1,2,1,1))
for (iter in 1:9){
 u   <- c(rbeta(9,1,0.8),1)
 res <- tsb(u)
 barplot(res,ylim=c(0,1),main=paste0("Random Sample #", iter),ylab="Probability")
}
graphics::mtext("Truncated Stick-Breaking Dist. (10 segments)",3,
     outer=TRUE,cex=1.5,line=1.5)

par(oldpar)
```
