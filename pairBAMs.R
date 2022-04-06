require(tidyverse)
require(fs)

args=commandArgs(trailing=T)
bams=scan(args[1],"")

manifest=tibble(BAM=bams) %>%
    mutate(Sample=basename(BAM) %>% gsub(".bam","",.)) %>%
    mutate(Patient=str_extract(Sample,"^([^_]+)_([^_]+)")) %>%
    mutate(Type=ifelse(grepl("_R$|_R_",Sample),"Normal","Tumor"))

normals=filter(manifest,Type=="Normal")

pairs=filter(manifest,Type=="Tumor") %>%
    left_join(normals,by="Patient") %>%
    filter(!is.na(Type.y)) %>%
    select(Normal=BAM.y,Tumor=BAM.x)

write_tsv(pairs,"pairs",col_names=F)

unpairedNormals=normals %>% filter(!BAM %in% pairs$Normal)
unpairedTumors=filter(manifest,Type=="Tumor") %>% filter(!BAM %in% pairs$Tumor)

unpaired=bind_rows(unpairedNormals,unpairedTumors) %>% arrange(Patient,Type)

write_tsv(unpaired,"unpaired")
