# Packages
library( "DESeq2" )
library(ggplot2)

# Read data
countData <- read.csv( '/Users/innaa/Documents/PEGBM_ns/DESeq2 inputs/counts_Villi_EVT_A.csv', header = TRUE, sep = ",")
head(countData)
metaData <- read.csv( '/Users/innaa/Documents/PEGBM_ns/DESeq2 inputs/meta_data_Villi_EVT_A.csv', header = TRUE, sep = ",")
metaData

# prog non prog 
dds <- DESeqDataSetFromMatrix(countData=countData, 
                              colData=metaData, 
                              design=~SegmentDisplayName  , tidy = TRUE)



dds <- DESeq(dds)
resultsNames(dds)

res1<-results(dds, name = resultsNames(dds)[2]) 
summary(res1)
res1 <- res1[order(res1$pvalue),]
head(res1)
res_tab = as.data.frame(res1)
# Export
write.csv(as.data.frame(res1), 
          file=paste("/Users/innaa/Documents/PEGBM_ns/DESeq2 inputs/Results/", 
                     resultsNames(dds)[2], "_villi_united.csv", sep =""))



# Get normalized counts
normalized_counts <- counts(dds, normalized=TRUE)


gene_of_interest <- "FUT8"

# Check distribution of expression values
summary(normalized_counts[gene_of_interest,])

# Create violin plot with log2 transformation
plot_data <- data.frame(
  Expression = normalized_counts[gene_of_interest,],
  Group = metaData$SegmentDisplayName
)

# Perform Kruskal-Wallis test
kw_test <- kruskal.test(Expression ~ Group, data = plot_data)
pval <- kw_test$p.value
pval_text <- ifelse(pval < 0.001, 
                    paste("p < 0.001"),
                    paste("p =", sprintf("%.3f", pval)))

# Plot with log2 transformation
ggplot(plot_data, aes(x=Group, y=log2(Expression + 1), fill=Group)) +
  geom_violin(trim=FALSE) +
  geom_boxplot(width=0.1, fill="white") +
  geom_jitter(width=0.1, size=1, alpha=0.5) +
  theme_bw() +
  labs(title=paste("Expression of", gene_of_interest),
       y="log2(Normalized Expression + 1)",
       x="Diagnostic Group") +
  theme(
    plot.title = element_text(hjust = 0.5, size=14),
    axis.text = element_text(size=12),
    axis.title = element_text(size=12)
  ) +
  annotate("text", 
           x = mean(1:length(unique(plot_data$Group))), 
           y = max(log2(plot_data$Expression + 1)) * 1.1,
           label = pval_text,
           size = 4)


