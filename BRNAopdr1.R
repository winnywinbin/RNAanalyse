
library(limma)
library(edgeR)
library(Glimma)
library(gplots)
#library(org.Mm.eg.db)
library(RColorBrewer)
# De data wordt ingeladen als seqdata.
seqdata <- read.delim("/Users/winnythoen/Desktop/BioInformatica/Jaar\ 3/BRNA/chronische_inflammatie_rawcounts.csv",
                      stringsAsFactors = FALSE, sep = ",")
# De sample info wordt ingeladen als sampleinfo
sampleinfo <- read.delim("/Users/winnythoen/Desktop/BioInformatica/Jaar\ 3/BRNA/chronische_inflammatie_metadata.csv",
                         sep = ",")

# Bekijk de eerste 6 regels.
head(seqdata)
# Hoeveel rows en columns zijn er?
dim(seqdata)
# 64102     9
sampleinfo
# Countdata is een nieuw data object, wat alle counts van de 8 samples.
countdata <- seqdata[,-(1)]
head(countdata)
rownames(countdata) <- seqdata[,1]
colnames(countdata)

# Filtering to remove lowly expressed genes
# We gebruiken threshold 0.5 bij minstens 2 samples.

myCPM <- cpm(countdata)
head(myCPM)
thresh <- myCPM > 0.5
head(thresh)
# Summery hoe vaak TRUE voorkomt per row
table(rowSums(thresh))
# 14197 row hebben alles TRUE
# WE willen alles houden wat ten minste 2 TRUEs in een row
keep <- rowSums(thresh) >= 2
# Nu de subset aanpassen om alleen de hoge te houden
counts.keep <- countdata[keep,]
summary(keep)
dim(counts.keep)
# Als de teller kleiner is dan 10-15 wordt dat als laag beschouwd
# Wat aangeeft dat het geassocieerde gen niet in dat sample tot expressie wordt gebracht.
plot(myCPM[,1], countdata[,1])
# nu wordt het verkleint om te kijken van dichtbij wat er daadwerkelijk gebeurd
# bij kleinere waardes
plot(myCPM[,1],countdata[,1],ylim=c(0,50),xlim=c(0,3))
# We voegen een lijn toe om goed te kunnen kijken waar de grens ligt.
abline(v = 0.5, col = "blue")
abline(h = 10, col = "blue")

# We maken een DGEList -> hierin bewaren we data.
# We gebruiken de counts.keep data , dit is een subset met alle hoge data.
y <- DGEList(counts.keep)
# Nu hebben een een goede lijst zonder de lage expressie genen,
# en de counts zijn opgeslagen in een DGEList.
# Nu gaan we kijken naar verschillende plots om te checken of 
# de data echt van goede kwaliteit is en wat we verwachten.
# Eerst kijken hoeveel reads we hebben in y
y$samples
# We can also plot the library sizes as a barplot to see whether there are any major discrepancies between the samples more easily.
barplot(y$samples$lib.size, names=colnames(y), las=2)
title("Barplot of library sizes")
# Count data is not normally distributed, 
# so if we want to examine the distributions of the raw counts we need to log the counts. 
# Next we’ll use box plots to check the distribution of the read counts on the log2 scale. 
# We can use the cpm function to get log2 counts per million, which are corrected for the different library sizes. 
# The cpm function also adds a small offset to avoid taking log of zero.
# Get log2 counts per million
logcounts <- cpm(y, log = TRUE)
# Check distributions of samples using boxplots
boxplot(logcounts, xlab="", ylab="Log2 counts per million", las=2)
# horizontale lijn, dat corresponds to the median  logCPM
abline(h=median(logcounts), col="blue")
abline(h=mean(logcounts), col="red", pch=18)
title("Unnormalised logCPM")

# Je ziet dat de samples erg dichtbij de horizontale lijn liggen
# Dat is een goed teken!!
# Als er iets verder vanaf zat zouden we dat verder moeten onderzoeken
plotMDS(y)
title("MDS plot ")
# We specify the option to let us plot two plots side-by-sde
par(mfrow=c(1,1))
# Let's set up colour schemes for CellType
# How many cell types and in what order are they stored?
levels(sampleinfo$celltype)
col.cell <- c("orange", "purple", "red", "blue")[sampleinfo$celltype]
data.frame(sampleinfo$celltype, col.cell)
# Redo de MDS with colouring cells
plotMDS(y, col=col.cell)
legend("bottomleft", fill = c("orange", "purple", "red", "blue"), 
       legend = unique(sampleinfo$celltype))
title("Cell type")
levels(sampleinfo$dex)
col.dex <- c("pink", "green")[sampleinfo$dex]
plotMDS(y, col=col.dex)
legend("bottomleft", fill = c("pink", "green"),
       legend = unique(sampleinfo$dex),
       cex = 0.8)
title("DEX")

# Nu gaan we proberen het 1 te maken
col.cell <- c("dark green", "purple", "dark red", "light blue")[sampleinfo$celltype]
shps <- c(15,16)[sampleinfo$dex]
plotMDS(y, dim=c(1,2), col=col.cell, pch = shps, cex = 2)
legend("bottomright", fill = c("dark green", "purple", "dark red", "light blue"),
       legend = unique(sampleinfo$celltype))
legend("topleft", legend = unique(sampleinfo$dex),
       pch = shps)
title("Cell type & Dex")

#Another alternative is to generate an interactive MDS plot 
#using the Glimma package. 
#This allows the user to interactively explore the different dimensions.
labels <- paste(sampleinfo$id, sampleinfo$celltype, sampleinfo$dex)
group <- paste(sampleinfo$celltype, sampleinfo$dex, sep = ".")
group <- factor(group)
glMDSPlot(y, labels = labels, groups = group, folder = "mds")

# Normaliseren
y <- calcNormFactors(y)

par(mfrow=c(1,2))
plotMD(logcounts, column = 1)
abline(h=0, col="grey")
plotMD(logcounts, column = 4)
abline(h=0, col="grey")

par(mfrow=c(1,2))
plotMD(y, column = 1)
abline(h=0, col="grey")
plotMD(y, column = 4)
abline(h=0, col="grey")

par(mfrow=c(1,2))
plotMD(logcounts, column = 4)
abline(h=0, col="grey")
plotMD(y, column = 4)
abline(h=0, col="grey")

# Clustering
var_genes <- apply(logcounts, 1, var)
select_var <- names(sort(var_genes, decreasing = TRUE))[1:500]
highly_variable_lcpm <- logcounts[select_var,]
mypallet <- brewer.pal(11, "RdYlBu")
morecols <- colorRampPalette(mypallet)
col.cell <- c("green", "blue", "grey", "pink")[sampleinfo$celltype]
heatmap.2(highly_variable_lcpm, col = rev(morecols(50)), trace = "none", 
          main = "Top 500 most variable genes across samples.", ColSideColors = col.cell,
          scale = "row")

# Create the design matrix
group

design <- model.matrix(~sampleinfo$celltype + sampleinfo$dex)
colnames(design) <- c("Intercept", "N061011", "N080611", "N61311", "dextreated")
design

# Voom transformation the data
#The voom transformation uses the experiment design matrix, and produces an EList object.
par(mfrow=c(1,1))
v <- voom(y, design, plot = TRUE)

par(mfrow=c(1,2))
boxplot(logcounts, xlab = "", ylab="Log2 counts per million",
        las=2, main="Unnormalised logCPM")
abline(h=median(logcounts), col="blue")
boxplot(v$E, xlab="", ylab="Log2 counts per million",
        las=2, main="Voom transformed logCPM")
abline(h=median(v$E), col="blue")

# Fit the linear Model
fit <-lmFit(v)
names(fit)
head(coef(fit))

cont.matrix <- makeContrasts(N08 = N080611, N6 = N61311, N06 = N061011, levels = design)
fit2 <- contrasts.fit(fit, cont.matrix)
fit2 <- eBayes(fit2)
summa.fit <- decideTests(fit2)
vennDiagram(summa.fit)
summary(summa.fit)

# Top 10 genen van elke sample weergegeven gesorteerd op p-value
topTable(fit2, coef = "N08", sort.by = "p")
topTable(fit2, coef = "N6", sort.by = "p")
topTable(fit2, coef = "N06", sort.by = "p")

# Vervolgens gaan we annotatie bijvoegen.
annotation <- read.delim("/Users/winnythoen/Desktop/BioInformatica/Jaar\ 3/BRNA/annotables_unique_grch37.csv",
                         sep = ",")
annotable <- annotation[annotation$ensgene %in% rownames(fit),] 
fit$genes <- annotable
names(fit$genes)
ann <- annotable[,c(1,3,9)]
fit2$genes <- ann

lim_N08 <- topTable(fit2, coef = "N08", sort.by = "p", n="Inf")
lim_N6 <- topTable(fit2, coef = "N6", sort.by = "p", n="Inf")
lim_N06 <- topTable(fit2, coef = "N06", sort.by = "p", n="Inf")
write.csv(lim_N08, file = "N08Results.csv", row.names = FALSE)
write.csv(lim_N6, file = "N6Results.csv", row.names = FALSE)
write.csv(lim_N06, file = "N06Results.csv", row.names = FALSE)


# We want to highlight the significant genes. We can get this from decideTests.
par(mfrow=c(1,2))
plotMD(fit2, coef = 1, status = summa.fit[,"N08"], values = c(-1,1))
plotMD(fit2, coef = 1, status = summa.fit[,"N6"], values = c(-1,1))
plotMD(fit2, coef = 1, status = summa.fit[,"N06"], values = c(-1,1))

# Visualiseren genen set/ Enrichment plots

names(fit2)
head(fit2$coefficients)
head(fit2$t)
par(mfrow=c(1,1))



volcanoplot(fit2, coef = 1, highlight = 100, names = fit2$genes$symbol)