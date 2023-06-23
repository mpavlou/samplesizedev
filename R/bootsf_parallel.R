
bootsf_parallel <-function(data, nboot=200){
  cal <-NULL
  cores <- parallel::detectCores()
  cl    <- parallel::makeCluster(cores[1]-1)

  doParallel::registerDoParallel(cl)

  a<-   foreach::foreach(j = 1:nboot, packages="speedglm") %dopar% {
    bs <- sample(nrow(data), replace=T)
    databs=data[bs,]
    xvarsbs=databs[,-1];ybs<-databs[,1]
    fitbs <- speedglm::speedglm(ybs~xvarsbs, family=binomial())

    eta_est    <- as.matrix(cbind(1,data[,-1]))%*%coef(fitbs)
    fitcal     <- speedglm::speedglm(data[,1]~eta_est, family=binomial())
    cal[j]     <- as.vector(coef(fitcal)[2])
  }

  parallel::stopCluster(cl)

  cs <- matrix(unlist(a), byrow=TRUE, nrow=nboot)
  return(median(cs,na.rm=TRUE))
}
