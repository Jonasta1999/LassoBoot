

library(dplyr)
library(readxl)
data <- read_xlsx("Block bootstrap/data.xlsx", col_names = F)
colnames(data) <- "v1"


block_bootstrap2 <- function(data, B, w){
  t <- dim(data)[1]
  k <- dim(data)[2]
  
  #Compute the number of blocks nedded
  s <- ceiling(t/w)
  #Generate the starting points
  bs <- floor(matrix(runif(s*B)*t, nrow = s, ncol = B)) + 1
  indices <- matrix(0, nrow = s*w, ncol = B)
  index <- 1
  #Adder is a variable that needs to be added each loop
  adder <- matrix(rep(0:(w-1), each = B), ncol = B, byrow = T)
  
  for(i in seq(1, t, w)){
    indices[i:(i+w-1),] <- matrix(rep(bs[index,], each = w), ncol = B) + adder
    index <- index + 1
  }
  
  indices <- indices[1:t,]
  indices[indices > t] <- indices[indices > t] -t
  
  bsdata <- matrix(sapply(indices, function(y) data[y,]), ncol = ifelse(!is.null(ncol(indices)), ncol(indices), 1))
  bsdata <- matrix(bsdata, nrow = t) 
  bsdata <- bsdata %>% t() %>% matrix(nrow = t) %>% data.frame()
  
  indices <- data.frame(indices)
  
  # Output BSDATA and INDICES
  return(list(bsdata = bsdata, indices = indices))
}




