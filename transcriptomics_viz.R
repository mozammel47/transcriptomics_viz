###Importing library
library(ggplot2)
library(plotly)
library(DESeq2)
library(tibble)
library(DT)
library(reshape2)
library(dplyr)


###We start by running DESeq2 to find differentially expressed genes.
##Loading data
countData <- as.matrix(read.csv(file="rna_count_data.txt",sep="\t",row.names="locus_tag"))

conditions <- c(rep("control",3), rep("treatment", 3))

mockRna.dse <- DESeqDataSetFromMatrix(countData,
            colData = as.data.frame(conditions), design = ~ conditions)

colData(mockRna.dse)$conditions <- factor(colData(mockRna.dse)$conditions,
                                          levels=c("control", "treatment"))

## Get a DESeqDataSet object
mockRna.dse <- DESeq(mockRna.dse)

##Extracting table

resultsNames(mockRna.dse) # lists the coefficients
res <- results(mockRna.dse, name="conditions_treatment_vs_control") 

resLFC <- lfcShrink(mockRna.dse, coef="conditions_treatment_vs_control", type="apeglm")

write.table(resLFC, file = "rna_fold-change-table.txt", sep = "\t")


############  Table  ####################

lfc_table= read.table(file="rna_fold-change-table.txt",sep = "\t", header = TRUE,
                     row.names = 1) %>% rownames_to_column(var = "locus_tag")

fc_table_sig = lfc_table %>% filter(pvalue<0.05 & log2FoldChange<=-1 | log2FoldChange>=1)

fc_table_sig = (format(fc_table_sig, digits=3))

fc_table_sig = datatable(fc_table_sig)

fc_table_sig = datatable(fc_table_sig, options = list(
  columnDefs = list(list(className = 'dt-center', targets="_all"))))


#############  VolcanoPlot #################
lfc_table$diffexpressed <- "NO"
lfc_table$diffexpressed[lfc_table$log2FoldChange > 1 & lfc_table$pvalue < 0.05] <- "UP"
lfc_table$diffexpressed[lfc_table$log2FoldChange < -1 & lfc_table$pvalue < 0.05] <- "DOWN"

volcano_plot = ggplotly(lfc_table %>% ggplot(aes(x=log2FoldChange, y=-log10(pvalue), z= locus_tag, col=diffexpressed))+
           geom_point(alpha =0.7, size=2)+ 
           scale_color_manual(values=c("darkorchid", "grey50", "salmon"))+
           labs(x="log2 fold change", y="-log10 p-value",title = "")+ 
           theme_my_custom()+
           scale_y_continuous(limits = c(0,50),breaks=c(0,10,20,30,40,50))+
           scale_x_continuous(limits = c(-5,5),breaks=c(-4,-2,0,2,4))+ 
           geom_vline(xintercept=c(-1, 1), alpha = 0.9, color="black", linetype = 2)+
           geom_hline(yintercept=-log10(0.05), alpha = 0.9, color="black", linetype = 2)
)  
######### PCA Plot##################

pca_plot = plotPCA(rlog(mockRna.dse), intgroup="conditions")+
  labs(x="PC1 (76.0 %)", y="PC2 (9.0 %)",title = "")+ 
  scale_color_manual(values = c("salmon", "darkorchid"), labels = c("Control", "Treatment"))+
  scale_x_continuous(breaks=c(-5,0,5))+scale_y_continuous(breaks=c(-3,0,3))+
  theme_my_custom()+
  theme(
    legend.position="bottom", legend.title=element_blank(),
    panel.grid = element_blank(),axis.line = element_blank())


#############  Heatmap  ###############

# Transform count data using the variance stablilizing transform
deseq2VST <- vst(mockRna.dse)

# Convert the DESeq transformed object to a data frame
deseq2VST <- assay(deseq2VST)
deseq2VST <- as.data.frame(deseq2VST)
deseq2VST$Gene <- rownames(deseq2VST)

#Keep only the significantly differentiated genes where the fold-change was at least 3
resOrdered1 <- res[complete.cases(res),]  #remove any rows with NA
resOrdered1 <- resOrdered1[order(resOrdered1$padj),]

sigGenes <- rownames(resOrdered1[resOrdered1$padj <= 10e-10 & abs(resOrdered1$log2FoldChange) > 2,])
deseq2VST <- deseq2VST[deseq2VST$Gene %in% sigGenes,]

# Now overwrite our original data frame with the long format
deseq2VST <- melt(deseq2VST, id.vars=c("Gene"))

# Make a heatmap

heatmap = ggplot(deseq2VST, aes(x=variable, y=Gene, fill=log2(value))) + 
           geom_tile(color="gray10") + 
           scale_fill_distiller(palette = "PRGn", "log2 fold change")+ theme_bw()+
           theme(axis.text.x=element_text(angle=65, hjust=1, vjust = 1),
                 axis.title.x = element_blank(),
                 axis.text = element_text(colour = "black"))

##### Building custom theme #####################
#Define theme_my_custom() function
theme_my_custom <- function(){ 
  font <- "Arial"   #assign font family up front
  theme_bw() %+replace%    #replace elements we want to change
    theme(
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black",fill = NA, size = 1),      
      plot.title = element_text(family = font, size = 14, hjust = 0, vjust = 2,color = "black"),              
      plot.subtitle = element_text(family = font, size = 12, color = "black"),
      plot.caption = element_text(family = font, size = 12, hjust = 1, color = "black"),
      axis.title = element_text(family = font, size = 12, color = "black"),
      axis.text = element_text(family = font, size = 12, color = "black"),
      axis.text.x = element_text(color = "black", margin=margin(5, b = 10)),
      legend.title = element_blank(),              
      legend.text = element_text(family = font, size = 11, color = "black"),
      strip.text = element_text(family = font, size = 12, color = "black", 
                                margin = margin(0.2,0,0.2,0, "cm"))     
    )
}
