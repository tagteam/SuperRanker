res <- NULL
for(i in 1:1500) {
  filename <- paste('SuperFinal-', i, '.RData', sep="")
  if(file.exists(filename)) {
    load(file=filename)
    res <- c(res, list(out))
  } else {
    cat(filename, " doesn't exist\n")
  }
}

save(res, file="joinedOutput.RData")