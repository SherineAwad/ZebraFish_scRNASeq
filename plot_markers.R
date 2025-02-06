library(Signac)
library(Seurat)
library(GenomeInfoDb)
library(BSgenome.Drerio.UCSC.danRer11)
library(future)
library(ggplot2)
library(patchwork)
library(Matrix)
library(dplyr)
set.seed(1234)
library(AcidGenomes)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)
myRDS <- args[1]
split_string <- strsplit(myRDS, '_')[[1]]
mysample <- split_string[1]
print(mysample) 

suffix=""
if (substring_exists <- grepl("subset", myRDS) )
{
suffix= "S"
}
print(suffix) 
myObject <- readRDS(myRDS)

DefaultAssay(myObject) <- "RNA"

gene_list <- rownames(myObject)
head(gene_list) 

genes1 = c("notch1a", "notch1b","notch3")
genes2 = c("her4.2","her4.2","her4.3","her4.4")
genes3 = c("dla","dlb","dlc","dld")
genes4 = c("dll4","jag1b","jag2b","dtx2")
genes5 = c("dtx4a","dtx4b","mib1","mib2")
genes6 = c("numb","numbl","clcf1","crlf1a")
genes7 = c("lepr","lepa","lepb" ,"mmp9","ascl1a","insm1a")
genes8 = c("her12","sox2","ascl1a","ascl1b","gli1","sfrp2","fgf19", "stat3") 
genes9 = c("rho","nrl","otx2","crx","guca1b","rom1a","rom1b")
genes10 = c("opn1mw1", "opn1mw2","opn1mw3","opn1mw4")
genes11 = c("sebox","bhlhe23","cabp5a","cabp5b","vsx2") 
genes12 = c("gad1a","gad1b","gad2","slc6a9","tfap2b","prox1a")
genes19 = c("pax6a","calb2a","calb2b","pcp4a","pcp4b","elavl3","isl1")
genes13 = c("lhx1a","cbln4","calb1","nefla","neflb","nefma","nefmb") 
genes14 = c("nefla", "neflb", "nefma","nefmb")
genes15 = c("ptprc", "csf2rb", "csf1ra")
genes16 = c("pax2a","pax2b","igf2a", "igf2b","gfap")
genes17 = c("tie1","mbpa", "rpe65a", "chata", "th", "acta2")
genes18 = c("vcx2", "ccnd1", "atoh7", "neurod1", "pax6")
genes19 = c("pax6a","calb2a","calb2b","pcp4a","pcp4b","elavl3","isl1")
genes20 = c("ebf3a","rbfox3a","rbfox3b","isl1","isl2a")
genes21 = c("isl2b","pou4f1","pou4f2","pou4f3","rbpms")
genes22 = c("opn1sw1", "opn1sw2","arr3a","thrb")
genes23 = c("prkca","pcp4a","pcp4b","isl1") 
genes24 = c("her4.1", "mpeg1.1", "csf1ra")
genes25 = c("kcnj8","rgs5a", "rgs5b")
genes26 = c("Onecut1","Onecut2","nr2e3","sncga","sncgb","thy1")
figure_name <- ""
figure_name <- paste(mysample, "_markers.pdf", sep=suffix)
pdf(file =figure_name)
FeaturePlot(myObject , features =genes1, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5, order =TRUE)
FeaturePlot(myObject , features =genes2, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5, order =TRUE)
FeaturePlot(myObject , features =genes3, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes4, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes5, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes6, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes7, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes8, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes9, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes10, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes11, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes12, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes13, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes14, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes15, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes16, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes17, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes18, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes19, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes20, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes21, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes22, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes23, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes24, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes25, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)
FeaturePlot(myObject , features =genes26, reduction = "umap", cols = c("lightgrey", "red"), pt.size = 0.5,order =TRUE)


dev.off()



