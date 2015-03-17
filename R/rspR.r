############
### Wrappers for the great rpsr (https://github.com/cwhidden/rspr) by Chris Whidden (FHCRC)
### Note that the functions here assume rspr is installed in the system, ie., at /usr/bin/rspr
### Copyleft (or the one to blame): Carvalho, LMF (2015)
### Last update: github is your friend.
############
rspr <- function(t1, t2, exact = TRUE){ 
  ## Simple two-trees comparison
  ## takes two trees and returns the exact/or approximate rSPR distance between them 
  ## note that this assumes rspr is installed in the system, ie., at /usr/bin/rspr
  tts <- list(t1, t2)
  class(tts) <- "multiPhylo"
  raw <- system2(input = paste(write.tree(tts)), command = "rspr", stdout = TRUE)
  if(exact){    
    d <- strsplit(raw[grep("total exact", raw)], "=")[[1]][2]
  }else{
    d <- strsplit(raw[grep("approx drSPR", raw)], "=")[[1]][2]
  }
  return(as.numeric(d))
}
rspr.matrix <- function(trees, type = c("restricted", "onerow", "full"), maxdist = 1){ 
  require(ape)
  ## takes a list of trees and returns a matrix with either
  ## (0) censored distances (maxdist if dij = maxdist and -1 if dij> maxdist) or;
  ## (1) distances from every tree to the first one in the list;
  ## (2) the actual spr distances.
  type <- match.arg(type)
  if(class(trees)!="multiPhylo") stop("Input is not of class multiPhylo")   
  root.test <- unlist(lapply(trees, is.rooted)) # testing if trees are rooted
  tmp <- tempfile("rsprmatrix", fileext = ".csv") # Unfortunately, we'll have to resort to this hack for the time being.
  if(type == "restricted"){
    if(any(root.test)){ # If any of the trees in the set is rooted, uses rooted (default) version
      dists <- system2(input = paste(write.tree(trees)), command = "rspr",
                       args = paste("-pairwise_max", maxdist), stdout = tmp)
    }else{
      dists <- system2(input = paste(write.tree(trees)), command = "rspr",
                       args = paste("-unrooted -pairwise_max", maxdist),  stdout = tmp)
    }    
  }else{
    if(type == "onerow"){
      if(any(root.test)){ 
        dists <- system2(input = paste(write.tree(trees)), command = "rspr",  args = "-pairwise 0 1", stdout = tmp)
      }else{
        dists <- system2(input = paste(write.tree(trees)), command = "rspr",  args = "-unrooted -pairwise 0 1", stdout = tmp)
      }
    }else{
      if(any(root.test)){
        dists <- system2(input = paste(write.tree(trees)), command = "rspr",  args = "-pairwise", stdout = tmp)
      }else{
        dists <- system2(input = paste(write.tree(trees)), command = "rspr",  args = "-unrooted -pairwise", stdout = tmp)
      }
    }
  }
  m <- as.matrix(read.table(tmp, sep = ",", header = FALSE, fill = TRUE))
  if(any(dim(m)>1))  m[lower.tri(m)] <- t(m)[lower.tri(t(m))]
  rownames(m) <- colnames(m) <- NULL
  return(m)
}
