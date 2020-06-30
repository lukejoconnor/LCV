# LCV example script written by Katie Siewert

#Start with data munged using the ldsc package
trait1File=".../my_trait1_summary_statistics_munged.sumstats.gz"

trait2File=".../my_trait2_summary_statistics_munged.sumstats.gz"

#Load trait 1 data and calculate Zs
d1 = na.omit(read.table(gzfile(trait1File),header=TRUE,sep="\t",stringsAsFactors = FALSE))

#Load trait 2 data and calculate Zs
d2 = na.omit(read.table(gzfile(trait2File),header=TRUE,sep="\t",stringsAsFactors = FALSE))

#Load LD scores
ldscoresFile="".../unannotated_LDscores.l2.ldsc""
d3=read.table(ldscoresFile,header=TRUE,sep=' ',stringsAsFactors=FALSE)

#Merge
m = merge(d3,d1,by="SNP")
data = merge(m,d2,by="SNP")

#Sort by position 
data = data[order(data[,"CHR"],data[,"BP"]),]

#Flip sign of one z-score if opposite alleles-shouldn't occur with UKB data
#If not using munged data, will have to check that alleles match-not just whether they're opposite A1/A2
mismatch = which(data$A1.x!=data$A1.y,arr.ind=TRUE)
data[mismatch,]$Z.y = data[mismatch,]$Z.y*-1
data[mismatch,]$A1.y = data[mismatch,]$A1.x
data[mismatch,]$A2.y = data[mismatch,]$A2.x


#Run LCV-need to setwd to directory containing LCV package
setwd(".../LCV/")
source(".../LCV/RunLCV.R")

LCV = RunLCV(data$L2,data$Z.x,data$Z.y)
sprintf("Estimated posterior gcp=%.2f(%.2f), log10(p)=%.1f; estimated rho=%.2f(%.2f)",LCV$gcp.pm, LCV$gcp.pse, log(LCV$pval.gcpzero.2tailed)/log(10), LCV$rho.est, LCV$rho.err)

