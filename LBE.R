library(LBE)
data(golub.pval)
#golub.pval <- head(golub.pval,n = 30)
#LBE零假设所占比例的估计程序，错误发现率和q值
res=LBE(golub.pval)

#$\pi_{0}(a)=(1/m)*\sum_{i=1}^{m}[-ln(1-\pi_{0})^{a}]/\Gamma(a+1)$
function (pval, a = NA, l = 0.05, ci.level = 0.95, qvalues = TRUE, 
          plot.type = "main", FDR.level = 0.05, n.significant = NA) 
#pval:p值向量
#a:$[-ln(1-\pi_{0})^{a}]$,如果a==NA（默认），a值自动计算，标准差渐近线的上限，pi0值小于为1的阈值。a>=1,a值直接应用在$[-ln(1-pi)^{a}]$公式中。a<1,对p值进行转化。

#l：在a=NA时，阈值，渐近标准差的上限
#ci.level:pi0的置信区间
#qvalues：q值，是否估计q值和fdr值。qvalue=F，仅估计pi0的零假设比例。
#plot.type：plot.type = "none"，不输出图。plot.type = "main"，估计在纵轴为q值-横轴为p值的图。
#plot.type = "multiple"，输出1、p值的直方图。2、估计的q值-p值。
#3、针对每个q值为cutoff的显著性检验的数量。4、针对每个显著性检验的错误发现率(fdr)
#FDR.level：fdr的控制标准。
#n.significant：如果有值，fdr由“n.significant”中定义的拒绝域中的最小的p值进行估计。
data(hedenfalk.pval)
res=LBE(hedenfalk.pval)