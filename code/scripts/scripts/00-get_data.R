# retrieve dataset
source(args[[1]])
x <- fun()

# number of features & observations
dim(x)



dimnames(x) <- list(
  paste0("gene", seq_len(nrow(x))),
  paste0("cell", seq_len(ncol(x))))

saveRDS(x, args[[2]])