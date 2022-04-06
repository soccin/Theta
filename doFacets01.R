args=commandArgs(trailing=T)

if(len(args)<1) {
    cat("\n\n   usage: doFacets01.R COUNTFILE\n\n")
    quit()
}

pngCairo<-function(filename,width=14,height=8.5,pointsize=12,res=150) {

    png(filename,type="cairo",units="in",
        width=width,height=height,pointsize=pointsize,res=res)

}

library(readr)
library(yaml)

SEED=20211010
set.seed(20211010)

#countFile="Counts/counts_ASHPC_0001_Pa_P_2___ASHPC_0001_Pa_R_1_.txt.gz"
countFile=args[1]
samplePair=gsub("_.txt.*","",gsub("counts_","",basename(countFile)))
tumor=gsub("___.*","",samplePair)
normal=gsub(".*___","",samplePair)
oBase=cc("facets",samplePair)

params=list()
params$args=list(
    SEED=SEED,
    countFile=countFile,
    tumor=tumor,
    normal=normal)

oDir=file.path("out",samplePair)
dir.create(oDir, showWarnings = FALSE, recursive = TRUE)

library(facets)
rcmat=readSnpMatrix(countFile)

stats=list()
stats$normalDepth=as.list(quantile(rcmat$NOR.DP))
stats$tumorDepth=as.list(quantile(rcmat$TUM.DP))

xx = preProcSample(rcmat,snp.nbhd=1000)
oo=procSample(xx,cval=500)
fit=emcncf(oo)

md=c("loglik","purity","ploidy","dipLogR","emflags")
res=list()
res$fit=fit[md]

pngCairo(file.path(oDir,paste0(oBase,"_CopyNum.png")),height=11,width=8.5)
plotSample(x=oo,emfit=fit)
dev.off()

pngCairo(file.path(oDir,paste0(oBase,"_diLog.png")),height=8.5,width=8.5)
logRlogORspider(oo$out, oo$dipLogR)
dev.off()

write_csv(fit$cncf,file.path(oDir,paste0(oBase,"_CNCF.csv")))

params$stats=stats
params$results=res

write_yaml(params,file.path(oDir,paste0(oBase,".yaml")))


