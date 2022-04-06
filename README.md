# Project Theta

Code to run FACETS on WGS samples

FACETS (//github.com/mskcc/facets/) version Oct 8, 2021 (commit 3058bba) was used to process Whole Genome Samples (WGS). To generate the snp count files the included snp-pileup command was run using a cleaned version of dbSnp version 137 (filter to contain only SNV and with duplicates removed). The resulting tumor/normal pileups were then processed using the standard facets method as outlined in the manual. Specifically:
```
mat = readSnpMatrix(countFile)
xx  = preProcSample(mat,snp.nbhd=1000)
oo  = procSample(xx,cval=500)
fit = emcncf(oo)
```

