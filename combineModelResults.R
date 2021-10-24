library(ggfortify)
library(caret)
library(pROC)
library(missForest)
library(randomForest)
library(Boruta)

#This script make prediction for each gene use model trained in trainModel.R
set.seed(8)
#minium number of mutations (in each of class) required for gene-specific model, genes with mutations less than this number will be predicted by non-gene specific model
min_gene_num <- 45
#weight ratio of target gene to non-target gene in gene specific model
targetGeneWeight <- 3.0

PID_SNV_BOX_file <- "./test/Rapid.HGMD_OMIM_PID.4varcards.tsv.anno.allChrs.4SNVBox.out.forR.CLIN_SIG.uniq4R.test.txt"
nonPID_SNV_BOX_file <- "./test/gnomad.benignMuts.4SNVBox.out.forR.CLIN_SIG.uniq4R.highAC.test.txt"
#read Varcards annotated data
PID_VarCards_file <- "./test/Rapid.HGMD_OMIM_PID.4varcards.tsv.anno.allChrs.CLIN_SIG.uniq.test.txt"
nonPID_VarCards_file <- "./test/gnomad.benignMuts.varcardsAnno.noFileName.CLIN_SIG.uniq.highAC.test.txt"

args <- commandArgs(T)
if(length(args) < 1){
        print("Usage: Rscript combineModelResults.R <output directory>")
        stop()
}
output_dir <- args[1]

#output folder
BorutaOut <- paste0(output_dir, "/BorutaOut")
#output folder for model object files
ModelObjs <- paste0(output_dir, "/ModelObjs")

#get classes balanced training dataset
source("readSNVdata.R")
source("balanceTwoGroups.R")
SNVDataObj <- readData(PID_SNV_BOX_file, nonPID_SNV_BOX_file, PID_VarCards_file, nonPID_VarCards_file, paste0(BorutaOut, "/all_data"))
balancedDataObj <- getIndex4BalancedData(SNVDataObj)
df_SNV_Box <- balancedDataObj$df_SNV_Box
varCards_data_withNA <- balancedDataObj$varCards_data_withNA
geneSymbols <- balancedDataObj$geneSymbols
dataset <- balancedDataObj$dataset
labels <- dataset$Group

#Only train gene specfic model for genes with >= min_gene_num PID and nonPID mutations
geneFreq.pid <- as.data.frame(table(geneSymbols[labels=="PID"]))
genesWithEnoughMuts_PID <- geneFreq.pid[geneFreq.pid$Freq > min_gene_num, 1]

geneFreq.nonPID <- as.data.frame(table(geneSymbols[labels=="Non_PID"]))
genesWithEnoughMuts_nonPID <- geneFreq.nonPID[geneFreq.nonPID$Freq > min_gene_num, 1]

genesWithEnoughMuts <- intersect(genesWithEnoughMuts_nonPID, genesWithEnoughMuts_PID)
genesWithEnoughMutIndices <- geneSymbols %in% genesWithEnoughMuts
otherGeneIndices <- !genesWithEnoughMutIndices

#new gene symbol list with genes without enough mutations renamed to "Other"
newGeneSymbols <- geneSymbols
newGeneSymbols[otherGeneIndices] <- "Other"


write.table(cbind(dataset, newGeneSymbols, as.character(labels)), paste0(BorutaOut, "/data4training.txt"), sep="\t", quote=F)

write.table(cbind(dataset, newGeneSymbols, geneSymbols, as.character(labels)), paste0(BorutaOut, "/AllDataUseForModel.tmpOriGeneName.txt"), sep="\t", quote=F)

#tmpIndices <- 1:length(newGeneSymbols)%%7==1;dataset=dataset[tmpIndices, ];newGeneSymbols=newGeneSymbols[tmpIndices]

trainModelByGene <- function(dataset, newGeneSymbols, useAllGenes=F){
	geneModels <- list();
	pred4AllGenes <- data.frame()
	allGenes <- unique(newGeneSymbols)
	if(useAllGenes){
		allGenes <- c("all_genes")
	}
	#loop through all genes
	for(gene in allGenes){
		print(gene)
		remainingGenesData <- dataset[newGeneSymbols!=gene,]
		geneData <- dataset[newGeneSymbols==gene,]
		nonOtherData <- dataset[newGeneSymbols!="Other",]
		selectedData <- NULL
		weights <- NULL
		needed_index_cutoff <-NULL
		#training data for model 'all_genes', uniform weights
		if(gene == "all_genes"){
			selectedData <- dataset
			geneData <- dataset
			weights <- rep(1.0, dim(selectedData)[1])
			needed_index_cutoff <- 1
		#training data for model 'Other', uniform weights. Use all genes for training
		}else if(gene == "Other"){
			selectedData <- rbind(remainingGenesData, geneData)
			weights <- rep(1.0, dim(selectedData)[1])
			needed_index_cutoff <- dim(remainingGenesData)[1]+1
		}else{
			#print("gene specific")
			selectedData <- rbind(remainingGenesData, geneData)
			weights <- rep(1.0, dim(selectedData)[1])
			needed_index_cutoff <- dim(remainingGenesData)[1]+1
			weights[needed_index_cutoff:length(weights)] <- targetGeneWeight
		}
		fit.rf.pred <- readRDS(file = paste0(ModelObjs, "/", gene, ".pred.rds"))
		targetPred <- fit.rf.pred[fit.rf.pred$rowIndex>=needed_index_cutoff,]
		fit.rf <- NULL
		pred4AllGenes <- rbind(pred4AllGenes, targetPred)	
	}
	colnames(pred4AllGenes)[3] <- "obs"
	
	list(pred_df=pred4AllGenes)
}

geneSpecificModelData <- trainModelByGene(dataset, newGeneSymbols)
pred4AllGenes <- geneSpecificModelData$pred_df

all_gene_ModelData <- trainModelByGene(dataset, newGeneSymbols, useAllGenes=T)
pred_dfWithOldMethodsUseAllGenes <- all_gene_ModelData$pred_df




pdf(paste0(BorutaOut, "/AUC.pdf"), width=5, height=5)
par(mar=c(3,3,2,1))
sz_col1 <- rgb(241/255,88/255,84/255)
sz_col2 <- rgb(93/255,165/255,218/255)
#plot(roc(pred4AllGenes$obs, pred4AllGenes$PID), col="red")
rr_geneSpecific <- roc(pred4AllGenes$obs, pred4AllGenes$PID)
plot(rr_geneSpecific$specificities, rr_geneSpecific$sensitivities, col=sz_col1, xlim=c(1,0), type="l", xaxt="n", yaxt="n", lwd=2, xlab="", ylab="")
lines(c(1,0), c(0,1), col="grey", lwd=1.6)
axis(1, mgp=c(1.5, .55, 0))
axis(2, mgp=c(1.5, .55, 0))
mtext("Specificity", side=1, line=1.55)
mtext("Sensitivity", side=2, line=1.55)
#plot(roc(pred_dfOnlySNVBox$obs, pred_dfOnlySNVBox$PID), col="green", add=T)
rr_NonGeneSpecific <- roc(pred_dfWithOldMethodsUseAllGenes$obs, pred_dfWithOldMethodsUseAllGenes$PID)
plot(rr_NonGeneSpecific, col=sz_col2, add=T)
#plot(roc(labels, varCards_pred_data[,1]))
methodNames <- sub("_score", "", colnames(varCards_data_withNA))

uniq_col <- c(rgb(93/255,165/255,218/255), rgb(96/255,189/255,104/255), rgb(241/255,124/255,176/255), rgb(77/255,77/255,77/255), rgb(178/255,145/255,47/255), rgb(178/255,118/255,178/255), rgb(241/255,88/255,84/255))


uniq_col_num <- 7
allColors <- rep(uniq_col,10)
auc_oldMethods <- c(1.0)
for(i in 2:dim(varCards_data_withNA)[2]){
	
	rr_oldMethods <- roc(labels, varCards_data_withNA[,i])
	auc_oldMethods <- c(auc_oldMethods, rr_oldMethods$auc)
	#print(paste(colnames(varCards_data_withNA)[i], rr_oldMethods$auc))
	plot(rr_oldMethods, add=T, col=allColors[i], lty=as.integer(i/uniq_col_num)+2)
	print(paste(methodNames[i], rr_oldMethods$auc))
}
dev.off()
pdf(paste0(BorutaOut, "/AUC.legend.pdf"), width=3,heigh=5)
par(mar=c(0,0,0,0))
plot(c(-10,-10), xlim=c(0,0.52), ylim=c(26,1), xaxt="n", yaxt="n", xlab="", ylab="")
x1 <- 0
x2 <- 0.15
x3 <- 0.18
ypos <- 1
lines(c(x1,x2), c(ypos,ypos), col=sz_col1, lwd=2)
text(x3, ypos, paste("Gene specific (", signif(rr_geneSpecific$auc,2), ")", sep=""), adj=0)
ypos <- 2
lines(c(x1,x2), c(ypos,ypos), col=sz_col2, lwd=2)
text(x3, ypos, paste("Non gene specific (", signif(rr_NonGeneSpecific$auc,2), ")", sep=""), adj=0)

auc_order <- order(auc_oldMethods, decreasing=T)
ypos <- 2
for(i in auc_order){
	if(i== 1){
		next
	}
	ypos <- ypos +1

	lines(c(x1,x2), c(ypos,ypos), col=allColors[i], lty=as.integer(i/uniq_col_num)+2, lwd=2)
	text(x3, ypos, paste(methodNames[i], " (", signif(auc_oldMethods[i],2), ")", sep=""), adj=0)
}
dev.off()

pdf(paste0(BorutaOut, "/featurePlot.pdf"), width=13,heigh=13)
#featurePlot(x=df2, y=labels, plot="box")
scales <- list(x=list(relation="free", cex=1.3), y=list(relation="free", cex=1.3))
featurePlot(x=df_SNV_Box, y=labels, plot="density", 
scales=scales, plot.points = FALSE, layout = c(5, 5), 
auto.key = list(columns = 2), adjust = 1.5, lwd=3)
dev.off();
