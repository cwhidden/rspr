# rspR 
A couple R wrappers to rspr.
rspr() is a simple tree versus tree comparison. If exact == TRUE the exact distance is reported.
If not, then the approximation is returned.
```R
rspr(tree1, tree2, exact = TRUE) 
```
rspr.matrix() computes the matrix of distances in a list of trees (of class "multiPhylo").

```R
rspr.matrix(treelist, type = c("restricted", "onerow", "full"), maxdist = 1) 
```
The function will
+ Compute censored SPR distances, i.e., maxdist if dij <= maxdist and -1 if dij> maxdist (type = "restricted"), or
+ compute the set of SPR distances from the first tree (type = "onerow") or
+ compute the full matrix (type = "full"). Fun for small data sets, may be impractical for large ones.

See [rspR_examples](https://github.com/maxbiostat/rspr/blob/master/R/rspR_examples.r) for usage examples.
