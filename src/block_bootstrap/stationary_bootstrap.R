

library(dplyr)
library(readxl)
data <- read_xlsx("Google Drive/Mit drev/Speciale/Block bootstrap/data.xlsx", col_names = F)
colnames(data) <- "v1"


stationary_bootstrap <- function(data, B, w){
  t <- dim(data)[1]
  k <- dim(data)[2]
  
  #Define the probability of a new block
  p <- 1/w
  #Set up the bsdata and indices
  indices = matrix(0, nrow = t, ncol = B)
  #Initial positions
  indices[1,] <- ceiling(t*matrix(runif(B),1,B))
  #Set up the random numbers
  select <- matrix(runif(t*B),t,B) < p
  select_indices <- which(select, arr.ind = T)
  indices[select_indices] <- ceiling(runif(sum(select), max = t))
  
  for(i in 2:t){
    # Determine whether we stay (rand>p) or move to a new starting value
    # (rand<p)
    indices[i,!select[i,]] <- indices[i-1,!select[i,]] + 1
  }
  indices[indices>t] <- indices[indices>t] - t
  
  bsdata <- matrix(sapply(indices, function(y) data[y,]), ncol= ncol(indices))
  bsdata <- data.frame(bsdata)
  
  indices <- data.frame(indices)
  
  # Output BSDATA and INDICES
  return(list(bsdata = bsdata, indices = indices))
}

test <- stationary_bootstrap(data,100,5)




