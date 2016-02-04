#source("http://bioconductor.org/biocLite.R") 
#biocLite("LBE") 
#install.packages("PhViD")
library(PhViD)
library(plyr)
result <- read.csv("/home/xuefliang/PhViD/trt.csv",stringsAsFactors = T)
#table(result$trt,result$pt)
result$newvar <- rep(1,length(result$trt))
#aggregate(result$newvar,by = list(result$trt,result$pt),FUN=sum,na.rm=T)
signaldat <- ddply(result,.(trt,pt),summarize,sum=sum(newvar))
#参数数据框要求：第一列药品 第二列 副反应 第三列药品产生了几个副反应
#MARGIN.THRES：边际频数小于MARGIN.THRES的值将被忽略
signaldat <- as.PhViD(signaldat,MARGIN.THRES = 1)

#对象PhViD L数据框记录药品和反应 N合计的报告数
#对象PhViD data n11药品产生了几个反应 n1. 记录药品总反应发生数 n.1记录各种特定的副反应的发生数，
#比如休克1个，发热3个

#ROR(signaldat)
#n11代表公式中的a，n10代表公式中的b，n10 <- n1. - n11表示某药品总发生数-某药品的特定反应数
#n01代表公式中的c，n.1 - n11表示特定的反应数-某药品的特定反应数，n00代表公式中的d

DATABASE <- signaldat

ror <- function (DATABASE, OR0 = 1, MIN.n11 = 1, DECISION = 1, DECISION.THRES = 0.05, 
                 RANKSTAT = 1) 
{
  require("LBE")
  if (RANKSTAT == 2 & DECISION == 1) 
    stop("The FDR can't be used as decision rule with this ranking statistic")
  DATA <- DATABASE$data
  N <- DATABASE$N
  L <- DATABASE$L
  n11 <- DATA[, 1]
  n1. <- DATA[, 2]
  n.1 <- DATA[, 3]
  n10 <- n1. - n11
  n01 <- n.1 - n11
  n00 <- N - (n11 + n10 + n01)
  E <- n1. * n.1/N
  #MIN.n11生成可能信息的配对的最小数量，默认为1
  if (MIN.n11 > 1) {
    E <- E[n11 >= MIN.n11]
    n1. <- n1.[n11 >= MIN.n11]
    n.1 <- n.1[n11 >= MIN.n11]
    n10 <- n10[n11 >= MIN.n11]
    n01 <- n01[n11 >= MIN.n11]
    n00 <- n00[n11 >= MIN.n11]
    LL <- data.frame(drugs = L[, 1], events = L[, 2], n11)
    LL1 <- LL[, 1][n11 >= MIN.n11]
    LL2 <- LL[, 2][n11 >= MIN.n11]
    rm(list = "L")
    L <- data.frame(LL1, LL2)
    n11 <- n11[n11 >= MIN.n11]
  }
  Nb.Cell <- length(n11)
  #计算ROR
  logROR <- log(n11 * n00/(n10 * n01)) 
  var.logROR <- 1/n11 + 1/n10 + 1/n01 + 1/n00
  #计算95%CI下限
  LCI <- exp(logROR-1.96*sqrt(var.logROR)) 
  #p  Probability的缩写，表示概率函数.pnorm(0)计算标准正态分布从负无穷大到0的概率，即计算0左边的曲线下面积
  #d  Density的缩写，表示密度函数。标准正态分布x=0对应的值可以用dnorm(0)计算，即计算对应的y值。
  #q  Quantile的缩写，表示分位函数。如果知道标准正态分布从负无穷大到x的概率是0.9678，想要知道这个x的值，可以通过qnorm(0.9678)计算。
  #r Random的缩写，表示随机函数。用于随机生成符合正太分布的数值，如果想随机生成10个符合标准正态分布的函数，可以用rnorm(10)来获得
  #计算mean=log(OR0)，sd=sqrt(var.logROR)正态分布的曲线下面积
  #1-pnorm计算右侧概率
  pval.logOR.uni <- 1 - pnorm(logROR, log(OR0), sqrt(var.logROR))
  #OR0值即odds ratio，也称优势比、比数比，默认为1
  #petit_rankstat与OR0进行比较后的统计量
  petit_rankstat <- (logROR - log(OR0))/sqrt(var.logROR)
  pval.uni <- pval.logOR.uni
  pval.uni[pval.uni > 1] <- 1
  pval.uni[pval.uni < 0] <- 0
  PVAL.UNI <- pval.uni
  #不输出图
  LBE.res <- LBE(2 * apply(cbind(pval.uni, 1 - pval.uni), 
                           1, min), plot.type = "none")
  pi.c <- LBE.res$pi0
  fdr <- pi.c * sort(pval.uni[pval.uni <= 0.5])/(c(1:sum(pval.uni <= 
                                                           0.5))/Nb.Cell)
  fdr <- c(fdr, pi.c/(2 * ((sum(pval.uni <= 0.5) + 1):Nb.Cell)/Nb.Cell) + 
             1 - sum(pval.uni <= 0.5)/((sum(pval.uni <= 0.5) + 1):Nb.Cell))
  FDR <- apply(cbind(fdr, 1), 1, min)
  if (RANKSTAT == 2) {
    FDR <- rep(NaN, length(n11))
  }
  LB <- qnorm(0.025, logROR, sqrt(var.logROR))
  #RANKSTAT 相同配对的排名规则： 1 按P值，默认值 2 按ROR值95%CI下限
  if (RANKSTAT == 1) 
    RankStat <- PVAL.UNI
  if (RANKSTAT == 2) 
    RankStat <- LB
  #DECISION 信号生成的决策依据：1 FDR（错误发现率） 默认值，2 信号的数量 3 秩统计
  if (DECISION == 1 & RANKSTAT == 1) 
    Nb.signaux <- sum(FDR <= DECISION.THRES)
  if (DECISION == 2) 
    Nb.signaux <- min(DECISION.THRES, Nb.Cell)
  if (DECISION == 3) {
    if (RANKSTAT == 1) 
      Nb.signaux <- sum(RankStat <= DECISION.THRES, na.rm = TRUE)
    if (RANKSTAT == 2) 
      Nb.signaux <- sum(RankStat >= DECISION.THRES, na.rm = TRUE)
  }
  RES <- vector(mode = "list")
  RES$INPUT.PARAM <- data.frame(OR0, MIN.n11, DECISION, DECISION.THRES, 
                                RANKSTAT)
  RES$ALLSIGNALS <- data.frame(L[, 1][order(petit_rankstat, decreasing = TRUE)],
                               L[, 2][order(petit_rankstat, decreasing = TRUE)], 
                               n11[order(petit_rankstat, decreasing = TRUE)], 
                               E[order(petit_rankstat, decreasing = TRUE)], 
                               RankStat[order(petit_rankstat, decreasing = TRUE)], 
                               exp(logROR)[order(petit_rankstat,decreasing = TRUE)], 
                               n1.[order(petit_rankstat, decreasing = TRUE)],
                               n.1[order(petit_rankstat, decreasing = TRUE)], FDR)
  colnames(RES$ALLSIGNALS) <- c("drug code", "event effect", 
                                "count", "expected count", "p-value", "ROR", "drug margin", 
                                "event margin", "FDR")
  if (RANKSTAT == 2) {
    colnames(RES$ALLSIGNALS)[5] <- "LB95(log(ROR))"
  }
  RES$SIGNALS <- RES$ALLSIGNALS[1:Nb.signaux, ]
  RES$NB.SIGNALS <- Nb.signaux
  RES
}

#PRR()

#FDR的含义是阳性检验结果中判断错误的比例。
#FDR是假设检验错误率的控制指标，可根据需要灵活选取。
#FDR是筛选出的差异变量的评价指标。
#FDR控制时决定一个显著性水平的界值，从而使FDR被限制在某一固定水平上。通常采用线性向上的控制方法，假设变量间独立。
#FDR估计是指在某一“有统计学意义”阈值情况下所有选入的有差异变量的假发现率的估计值。
#$F(p)=p_{0}F_{0}(p)+p_{1}F_{1}(p)$
#$p_{0}$为真实无效假设所占总检验次数的比列，无差异变量占总变量个数的比例。
#$F_{0}(p)$为无效假设下p的分布函数
#$p_{1}$为实际有差异变量在所有变量中所占的比例
#$F_{1}(p)$为备则假设成立下p值的分布函数
#给出的FDR的估计，称其为q值。
#ror(signaldat)


bcpnn <- function (DATABASE, RR0 = 1, MIN.n11 = 1, DECISION = 1, DECISION.THRES = 0.05, 
          RANKSTAT = 1, MC = FALSE, NB.MC = 10000) 
{
  DATA <- DATABASE$data
  N <- DATABASE$N
  L <- DATABASE$L
  #n11代表公式中的a
  n11 <- DATA[, 1]
  #n1.代表a+b
  n1. <- DATA[, 2]
  #n.1代表a+c
  n.1 <- DATA[, 3]
  #n10代表公式中的b
  n10 <- n1. - n11
  #n01代表公式中的c
  n01 <- n.1 - n11
  #n00代表公式中的d
  n00 <- N - (n11 + n10 + n01)
  E <- n1. * n.1/N
  if (MIN.n11 > 1) {
    E <- E[n11 >= MIN.n11]
    n1. <- n1.[n11 >= MIN.n11]
    n.1 <- n.1[n11 >= MIN.n11]
    n10 <- n10[n11 >= MIN.n11]
    n01 <- n01[n11 >= MIN.n11]
    n00 <- n00[n11 >= MIN.n11]
    L <- L[n11 >= MIN.n11, ]
    n11 <- n11[n11 >= MIN.n11]
  }
  Nb.Cell <- length(n11)
  if (MC == FALSE) {
    post.H0 <- matrix(nrow = Nb.Cell, ncol = length(RR0))
    p1 <- 1 + n1.
    p2 <- 1 + N - n1.
    q1 <- 1 + n.1
    q2 <- 1 + N - n.1
    r1 <- 1 + n11
    r2b <- N - n11 - 1 + (2 + N)^2/(q1 * p1)
    EICb <- log(2)^(-1) * (digamma(r1) - digamma(r1 + r2b) - 
                             (digamma(p1) - digamma(p1 + p2) + digamma(q1) - digamma(q1 + 
                                                                                       q2)))
    VICb <- log(2)^(-2) * (trigamma(r1) - trigamma(r1 + r2b) + 
                             (trigamma(p1) - trigamma(p1 + p2) + trigamma(q1) - 
                                trigamma(q1 + q2)))
    post.H0 <- pnorm(log(RR0), EICb, sqrt(VICb))
    LB <- qnorm(0.025, EICb, sqrt(VICb))
  }
  if (MC == TRUE) {
    #MCMCpack包提供了进行贝叶斯推断和贝叶斯计算的工具（特别是MCMC）。 MCMCpack包的设计思想是针对特定的模型运用MCMC方法。
    require(MCMCpack)
    n1. <- n11 + n10
    n.1 <- n11 + n01
    Nb_Obs <- length(n11)
    q1. <- (n1. + 0.5)/(N + 1)
    q.1 <- (n.1 + 0.5)/(N + 1)
    q.0 <- (N - n.1 + 0.5)/(N + 1)
    q0. <- (N - n1. + 0.5)/(N + 1)
    a.. <- 0.5/(q1. * q.1)
    a11 <- q1. * q.1 * a..
    a10 <- q1. * q.0 * a..
    a01 <- q0. * q.1 * a..
    a00 <- q0. * q.0 * a..
    g11 <- a11 + n11
    g10 <- a10 + n10
    g01 <- a01 + n01
    g00 <- a00 + n00
    g1. <- g11 + g10
    g.1 <- g11 + g01
    post.H0 <- vector(length = length(n11))
    LB <- vector(length = length(n11))
    quantile <- vector("numeric", length = length(n11))
    for (m in 1:length(n11)) {
      p <- rdirichlet(NB.MC, c(g11[m], g10[m], g01[m], 
                               g00[m]))
      p11 <- p[, 1]
      p1. <- p11 + p[, 2]
      p.1 <- p11 + p[, 3]
      IC_monte <- log(p11/(p1. * p.1))
      temp <- IC_monte < log(RR0)
      post.H0[m] <- sum(temp)/NB.MC
      LB[m] <- sort(IC_monte)[round(NB.MC * 0.025)]
    }
    rm(p11, p1., p.1, IC_monte, temp)
    gc()
  }
  if (RANKSTAT == 1) 
    RankStat <- post.H0
  if (RANKSTAT == 2) 
    RankStat <- LB
  if (RANKSTAT == 1) {
    FDR <- (cumsum(post.H0[order(post.H0)])/(1:length(post.H0)))
    FNR <- rev(cumsum((1 - post.H0)[order(1 - post.H0)]))/(Nb.Cell - 
                                                             1:length(post.H0))
    Se <- cumsum((1 - post.H0)[order(post.H0)])/(sum(1 - 
                                                       post.H0))
    Sp <- rev(cumsum(post.H0[order(1 - post.H0)]))/(Nb.Cell - 
                                                      sum(1 - post.H0))
  }
  if (RANKSTAT == 2) {
    FDR <- (cumsum(post.H0[order(LB, decreasing = TRUE)])/(1:length(post.H0)))
    FNR <- rev(cumsum((1 - post.H0)[order(1 - LB, decreasing = TRUE)]))/(Nb.Cell - 
                                                                           1:length(post.H0))
    Se <- cumsum((1 - post.H0)[order(LB, decreasing = TRUE)])/(sum(1 - 
                                                                     post.H0))
    Sp <- rev(cumsum(post.H0[order(1 - LB, decreasing = TRUE)]))/(Nb.Cell - 
                                                                    sum(1 - post.H0))
  }
  if (DECISION == 1) 
    Nb.signaux <- sum(FDR <= DECISION.THRES)
  if (DECISION == 2) 
    Nb.signaux <- min(DECISION.THRES, Nb.Cell)
  if (DECISION == 3) {
    if (RANKSTAT == 1) 
      Nb.signaux <- sum(RankStat <= DECISION.THRES, na.rm = TRUE)
    if (RANKSTAT == 2) 
      Nb.signaux <- sum(RankStat >= DECISION.THRES, na.rm = TRUE)
  }
  RES <- vector(mode = "list")
  RES$INPUT.PARAM <- data.frame(RR0, MIN.n11, DECISION, DECISION.THRES, 
                                RANKSTAT)
  if (RANKSTAT == 1) {
    RES$ALLSIGNALS <- data.frame(L[, 1][order(RankStat)], 
                                 L[, 2][order(RankStat)], n11[order(RankStat)], E[order(RankStat)], 
                                 RankStat[order(RankStat)], (n11/E)[order(RankStat)], 
                                 n1.[order(RankStat)], n.1[order(RankStat)], FDR, 
                                 FNR, Se, Sp)
    colnames(RES$ALLSIGNALS) <- c("drug code", "event effect", 
                                  "count", "expected count", "post.H0", "n11/E", "drug margin", 
                                  "event margin", "FDR", "FNR", "Se", "Sp")
  }
  if (RANKSTAT == 2) {
    RES$ALLSIGNALS <- data.frame(L[, 1][order(RankStat, decreasing = TRUE)], 
                                 L[, 2][order(RankStat, decreasing = TRUE)], n11[order(RankStat, 
                                                                                       decreasing = TRUE)], E[order(RankStat, decreasing = TRUE)], 
                                 RankStat[order(RankStat, decreasing = TRUE)], (n11/E)[order(RankStat, 
                                                                                             decreasing = TRUE)], n1.[order(RankStat, decreasing = TRUE)], 
                                 n.1[order(RankStat, decreasing = TRUE)], FDR, FNR, 
                                 Se, Sp, post.H0[order(RankStat, decreasing = TRUE)])
    colnames(RES$ALLSIGNALS) <- c("drug code", "event effect", 
                                  "count", "expected count", "Q_0.025(log(IC))", "n11/E", 
                                  "drug margin", "event margin", "FDR", "FNR", "Se", 
                                  "Sp", "postH0")
  }
  RES$SIGNALS <- RES$ALLSIGNALS[1:Nb.signaux, ]
  RES$NB.SIGNALS <- Nb.signaux
  RES
}


bcpnn(signaldat, RR0 = 1, MIN.n11 = 3, DECISION = 3, DECISION.THRES = 0.05, 
      RANKSTAT = 2, MC = FALSE, NB.MC = 10000)