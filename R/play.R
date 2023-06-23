es    <- NULL
es_sd <- NULL

i<-0

for (n.predictors in seq(4, 30, 4)) {

  i<-i+1

n_actual <- samplesizedev(outcome="Binary", S = 0.9, p = 0.2, c = 0.75, n.predictors = n.predictors,  nsim = 1000, parallel = TRUE)$actual

expected <- expected_cs_mape_binary(n = n_actual, p = 0.2, c = 0.75, n.predictors = n.predictors, nsim = 1000, parallel = TRUE)

es[i]    <- expected[2]
es_sd[i] <- expected[3]

print(expected)

}

es    <- unlist(es)
es_sd <- unlist(es_sd)
n.predictors <- seq(4, 30, 2)

wide <- data.frame(es, es_sd, n.predictors)


ggplot(wide, aes(x = n.predictors, y = es)) + geom_line(size=1) +
  geom_errorbar(aes(ymin = es - 1.96*es_sd, ymax=es + 1.96*es_sd), width=.2,position=position_dodge(.2))+
  ylab("Estimated CS")+ xlab("Number of predictors") +
  geom_hline(yintercept = 0.9, linetype="dashed",size=1) +   geom_line(color="blue")

ggplot(wide, aes(x = n.predictors, y = es_sd)) + geom_line(size=1) +
  ylab("SD(CS)")+ xlab("Number of predictors") + geom_line(color="blue") + scale_x_continuous( breaks=n.predictors)


  ggplot(wide, aes(x = n.predictors, y = es)) + geom_line(size=1) +
           geom_errorbar(aes(ymin = es - 0.68*es_sd, ymax=es + 0.68*es_sd), width=.2,position=position_dodge(.2))+
           ylab("Estimated CS")+ xlab("Number of predictors") +
           geom_hline(yintercept = 0.9, linetype="dashed",size=1) +   geom_line(color="blue") +
           geom_hline(yintercept = 1, linetype="dashed",size=1, col="red") +   geom_line(color="blue") + scale_x_continuous( breaks=n.predictors)


  es          <- NULL
  es_sd       <- NULL
  es_rmsd     <- NULL

  es_lsf      <- NULL
  es_lsf_sd   <- NULL
  es_lsf_rmsd <- NULL



  for (n.predictors in seq(5, 30, 5)) {

    i<-i+1

    n_actual <- samplesizedev(outcome="Binary", S = 0.9, p = 0.2, c = 0.75, n.predictors = n.predictors,  nsim = 1000, parallel = TRUE)$actual

    system.time(expected <- expected_cs_mape_binary(n = n_actual, p = 0.2, c = 0.75, n.predictors = n.predictors, nsim = 1000, parallel = TRUE))

    system.time (expected_lsf <- expected_cs_mape_binary(n = n_actual, p = 0.2, c = 0.75, n.predictors = n.predictors, nsim = 1000, parallel = TRUE, method = "LSF"))

    print(expected)

    print(expected_lsf)


    es[i]      <- expected[2]
    es_sd[i]   <- expected[3]
    es_rmsd[i] <- expected[4]

    es_lsf[i]      <- expected_lsf[2]
    es_lsf_sd[i]   <- expected_lsf[3]
    es_lsf_rmsd[i] <- expected_lsf[4]

  }

  es          <- unlist(es)
  es_sd       <- unlist(es_sd)
  es_rmsd     <- unlist(es_rmsd)

  es_lsf      <- unlist(es_lsf)
  es_lsf_sd   <- unlist(es_lsf_sd)
  es_lsf_rmsd <- unlist(es_lsf_rmsd)


 long <- data.frame(rbind(cbind(es, es_sd, es_rmsd), cbind(es_lsf, es_lsf_sd, es_lsf_rmsd)))
 long$method <- c( rep(0,6), rep(1,6) )
 long$method = factor(long$method)
 long$predictors = rep(seq(5, 30, 5),2)
 levels(long$method) = c("MLE", "LSF")

 names(long)=c("cs", "cs_sd", "cs_rmsd", "method", "predictors")


 ggplot(long, aes(x=predictors, y=cs_rmsd, group=method, col=method)) + geom_line(size=1) +
   ylab("Root-Mean Square Deviation from target value of CS")+ xlab("Number of predictors") +
   geom_vline(xintercept = 10, linetype="dashed",size=1) + scale_x_continuous( breaks=long$predictors)



 ggplot(long, aes(x=as.factor(predictors), y=sqrt(cs_sd), group=method, col=method)) + geom_line(size=1) +
   ylab("Standard Deviation of CS")+ xlab("Number of predictors") +
   geom_vline(xintercept = 3, linetype="dashed",size=1) + + scale_y_continuous(limits = c(0.8, 1.3))


 ggplot(long, aes(x=as.factor(predictors), y=cs, group=method, col=method)) + geom_line(size=1) +
   ylab("Mean CS")+ xlab("Number of predictors") +
 #  geom_vline(xintercept = 3, linetype="dashed",size=1) +
   geom_errorbar(aes(ymin = cs - 0.68*cs_sd, ymax = cs + 0.68*cs_sd), width=.2,position=position_dodge(.2))

 #
 # ggplot(long, aes(x=as.factor(predictors), y=sqrt(mape_lsf), group=method, col=method)) + geom_line(size=1) +
 #   ylab("MAPE")+ xlab("Number of predictors") +
 #   geom_vline(xintercept = 3, linetype="dashed",size=1)

 save.image("n_predictors.Rdata")

 load("n_predictors.Rdata")


 ###############################

 n.predictors <- 20


 es          <- NULL
 es_sd       <- NULL
 es_rmsd     <- NULL

 es_lsf      <- NULL
 es_lsf_sd   <- NULL
 es_lsf_rmsd <- NULL

 mape        <- NULL
 mape_lsf    <- NULL



 for (n_actual in c(1580/4, 1580/2, 1580*3/4, 1580)) {

   i<-i+1

   system.time(expected <- expected_cs_mape_binary(n = n_actual, p = 0.2, c = 0.75, n.predictors = n.predictors, nsim = 1000, parallel = TRUE))

   system.time (expected_lsf <- expected_cs_mape_binary(n = n_actual, p = 0.2, c = 0.75, n.predictors = n.predictors, nsim = 1000, parallel = TRUE, method = "LSF"))

   print(expected)

   print(expected_lsf)


   es[i]      <- expected[2]
   es_sd[i]   <- expected[3]
   es_rmsd[i] <- expected[4]

   es_lsf[i]      <- expected_lsf[2]
   es_lsf_sd[i]   <- expected_lsf[3]
   es_lsf_rmsd[i] <- expected_lsf[4]

   mape[i]        <- expected[5]
   mape_lsf[i]    <- expected_lsf[5]

 }

 es          <- unlist(es)
 es_sd       <- unlist(es_sd)
 es_rmsd     <- unlist(es_rmsd)
 mape        <- unlist(mape)

 es_lsf      <- unlist(es_lsf)
 es_lsf_sd   <- unlist(es_lsf_sd)
 es_lsf_rmsd <- unlist(es_lsf_rmsd)
 mape_lsf    <- unlist(mape_lsf)



 long <- data.frame(rbind(cbind(es, es_sd, es_rmsd, mape), cbind(es_lsf, es_lsf_sd, es_lsf_rmsd, mape_lsf)))
 long$method <- c( rep(0,4), rep(1,4) )
 long$method = factor(long$method)
 long$n = rep(c(1580/4, 1580/2, 1580*3/4, 1580), 2)
 levels(long$method) = c("MLE", "LSF")

 names(long)=c("cs", "cs_sd", "cs_rmsd", "mape", "method", "n")


 ggplot(long, aes(x=n, y=cs_rmsd, group=method, col=method)) + geom_line(size=1) +
   ylab("Root-Mean Square Deviation from target value of CS")+ xlab("Number of predictors") +
   geom_hline(yintercept = long$cs_rmsd[4], linetype="dashed",size=1) + scale_x_continuous( breaks=long$n)

 ggplot(long, aes(x=n, y=mape, group=method, col=method)) + geom_line(size=1) +
   ylab("Root-Mean Square Deviation from target value of CS")+ xlab("Number of predictors") +
   geom_hline(yintercept = long$mape[4], linetype="dashed",size=1) + scale_x_continuous( breaks=long$n)

 ggplot(long, aes(x=n, y=sqrt(cs_sd), group=method, col=method)) + geom_line(size=1) +
   ylab("Standard Deviation of CS")+ xlab("Number of predictors") +
   geom_vline(xintercept = 3, linetype="dashed",size=1) + scale_y_continuous(limits = c(0.8, 1.3))

 ggplot(long, aes(x=n, y=cs, group=method, col=method)) + geom_line(size=1) +
 ylab("Mean CS")+ xlab("Number of predictors") +
   #  geom_vline(xintercept = 3, linetype="dashed",size=1) +
 geom_errorbar(aes(ymin = cs - 0.68*cs_sd, ymax = cs + 0.68*cs_sd), width=.2,position=position_dodge(.2))









