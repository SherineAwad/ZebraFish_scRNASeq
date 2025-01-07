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

split_string <- strsplit(myRDS, "_")[[1]]
mysample <- split_string[1]

myObject <- readRDS(myRDS)



DefaultAssay(myObject) <- "RNA"

myObject%>%dplyr::glimpse()
print("===============================")
print(myRDS)
table(Idents(myObject)) 
prop.table(table(Idents(myObject)))
table(myObject@meta.data[,'sample'])
print("===============================")



table(myObject@meta.data$sample, myObject@active.ident)
#Rod   BC  RGC   AC Progenitors   HC   MG Cones


df <- data.frame(CellType = c("Rod","BC","RGC","AC","Progenitors","HC","MG","Cones"), 
S1=c(2224,657,190,452,1738,757,4944,414), 
S2=c(2673,1201,537,1871,1966,1160,5805,565)) 
library(tidyr)
df_long <- df %>% gather(key = "Condition", value = "Value", -CellType)
df_percent <- df_long %>%
  group_by(Condition) %>%
  mutate(Percent = Value / sum(Value) * 100)


stallion = c("Rod"="#D51F26","BC"="#003782","RGC"="#208A42","AC"="#820078","Progenitors"="#F47D2B", "HC"="#FFA500","MG"="#8A9FD1","Cones" ="#E6C122")
pdf(file = "Cellratio.pdf", width=4, height=4, onefile=FALSE)
ggplot(df_percent, aes(x = Condition, y = Percent, fill = CellType)) +
  geom_bar(stat = "identity", position = "fill") +  # position = "fill" makes it a percent stacked barplot
  scale_y_continuous(labels = scales::percent) +    # Show y-axis as percentage
  labs(y = "Cell Ratio", x = "Sample", fill = "Cell Type") +
  theme_minimal() +
  ggtitle("Celltype Ratio") +
   scale_fill_manual(values = stallion) +
   theme(axis.text.x = element_text(angle = 45, hjust = 1),
   axis.line = element_line(size = 1.5))


dev.off()




