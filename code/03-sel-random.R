fun <- \(x){
    s <- x[[1]]
    o <- sample(seq_len(nrow(s)), size = round(nrow(s)*0.1))
}