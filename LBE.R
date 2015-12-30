library(LBE)
data(golub.pval)
#golub.pval <- head(golub.pval,n = 30)
#LBE零假设所占比例的估计程序，错误发现率和q值
res=LBE(golub.pval,plot.type = "none")

#$\pi_{0}(a)=(1/m)*\sum_{i=1}^{m}[-ln(1-\pi_{0})^{a}]/\Gamma(a+1)$
function (pval, a = NA, l = 0.05, ci.level = 0.95, qvalues = TRUE, 
          plot.type = "main", FDR.level = 0.05, n.significant = NA) 
#pval:p值向量
#a:$[-ln(1-\pi_{0})^{a}]$,如果a==NA（默认），a值自动计算，为渐近标准差的上界，pi0为判断为真的零假设的比例，其值小于1。。
#a>=1,a值直接应用在$[-ln(1-pi)^{a}]$公式中。a<1,对p值进行转化。
#l：在a=NA时，阈值，渐近标准差的上界
#ci.level:pi0的置信区间的水平
#qvalues：估计q值或fdr值的逻辑值。qvalue=F，判断为真的零假设的比例。
#plot.type：plot.type = "none"，不输出图。plot.type = "main"，估计在纵轴为q值-横轴为p值的图。
#plot.type = "multiple"，输出1、p值的直方图。2、估计的q值-p值。
#3、针对每个q值为cutoff的显著性检验的数量。4、针对每个显著性检验的错误发现率(fdr)
#FDR.level：fdr的控制标准,仅在n.significant = NA时使用。
#n.significant：如果指定值，由“n.significant”的最小的p值估计fdr的拒绝域中。
 
#结果
#call 计算公式回显
#FDR n.significant == NA时为fdr的控制线，n.significant != NA时为fdr的估计值
#pi0 pi0的估计值，为真实无效假设所占总检验次数的比列。
#pi0.ci pi0的置信区间
#ci.level pi0的置信区间的水平
#a 这个$[-ln(1-\pi_{0})^{a}]$公式中用到的值
#l 渐近标准差的上界
#qvalues 估计出的一组q值向量
#pvalues 输入的p值向量
#significant 是否拒绝零假设的指标
#n.significant 拒绝零假设的数目
data(hedenfalk.pval)
res=LBE(hedenfalk.pval)