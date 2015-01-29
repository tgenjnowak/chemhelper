
diff.mat2letter.not <- function(g){
    
  n <- nrow(g)
  
  # Re-arrange data into an "edge list" for use in igraph (i.e. which groups are "connected") - Solution from "David Eisenstat" ()
  same <- which(g)
  g2 <- data.frame(N1=((same-1) %% n) + 1, N2=((same-1) %/% n) + 1)
  g2 <- g2[order(g2[[1]]),] # Get rid of loops and ensure right naming of vertices
  g3 <- simplify(graph.data.frame(g2,directed = FALSE))
  
  
  
  # Calcuate the maximal cliques - these are groupings where every node is connected to all others
  cliq <- maximal.cliques(g3) # Solution from "majom" ()
  
  # Reorder by level order - Solution from "MrFlick" ()
  ml<-max(sapply(cliq, length))
  reord <- do.call(order, data.frame(
    do.call(rbind, 
            lapply(cliq, function(x) c(sort(x), rep.int(0, ml-length(x))))
    )
  ))
  cliq <- cliq[reord]
  cliq
  
  # Generate labels to  factor levels
  lab.txt <- vector(mode="list", n) # empty list
  lab <- letters[seq(cliq)] # clique labels
  for(i in seq(cliq)){ # loop to concatenate clique labels
    for(j in cliq[[i]]){
      lab.txt[[j]] <- paste0(lab.txt[[j]], lab[i])
    }
  }
  
  
  
  
  return(unlist(lab.txt))
  
}



edges2diff.mat <- function(lower,upper){
  
  n <- length(lower)
  g <- outer(lower, upper,"-")
  g <- !(g<0)
  g <- g + t(g) # not necessary, but make matrix symmetric
  g <- g!=1
  rownames(g) <- 1:n # change row names
  colnames(g) <- 1:n # change column names
  
  return(g)
}

