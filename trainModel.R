library(ggfortify)
library(caret)
library(pROC)
library(missForest)
library(randomForest)
library(Boruta)

#This script trains a model for a single gene (or all_genes for non-gene-specific model), then output the model to an rds file
set.seed(8)
#name of target gene
args <- commandArgs(T)
if(length(args) < 2){
	print("Usage: Rscript trainModel.R <target gene name> <output directory>")
	stop()
}
targetGene <- args[1]
output_dir <- args[2]

#!!Please change the input training data before running the script
#SNVbox annotated data for PID and nonPID mutations
PID_SNV_BOX_file <- "./test/Rapid.HGMD_OMIM_PID.4varcards.tsv.anno.allChrs.4SNVBox.out.forR.CLIN_SIG.uniq4R.test.txt"
nonPID_SNV_BOX_file <- "./test/gnomad.benignMuts.4SNVBox.out.forR.CLIN_SIG.uniq4R.highAC.test.txt"
#read Varcards annotated data
PID_VarCards_file <- "./test/Rapid.HGMD_OMIM_PID.4varcards.tsv.anno.allChrs.CLIN_SIG.uniq.test.txt"
nonPID_VarCards_file <- "./test/gnomad.benignMuts.varcardsAnno.noFileName.CLIN_SIG.uniq.highAC.test.txt"


#minium number of mutations (in each of class) required for gene-specific model, genes with mutations less than this number will be predicted by non-gene specific model
min_gene_num <- 45
#weight ratio of target gene to non-target gene in gene specific model
targetGeneWeight <- 3.0

#output folder
BorutaOut <- paste0(output_dir, "/BorutaOut")
#output folder for model object files
ModelObjs <- paste0(output_dir, "/ModelObjs")

dir.create(file.path(BorutaOut), showWarnings = FALSE)
dir.create(file.path(ModelObjs), showWarnings = FALSE)

#get classes balanced training dataset
source("readSNVdata.R")
source("balanceTwoGroups.R")
#read the training data files
SNVDataObj <- readData(PID_SNV_BOX_file, nonPID_SNV_BOX_file, PID_VarCards_file, nonPID_VarCards_file, paste0(BorutaOut, "/all_data"))
#balance the two classes in the training dataset
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

#new gene symbol list, genes without enough mutations renamed to "Other"
newGeneSymbols <- geneSymbols
newGeneSymbols[otherGeneIndices] <- "Other"


write.table(cbind(dataset, newGeneSymbols, as.character(labels)), paste0(BorutaOut, "/data4training.txt"), sep="\t", quote=F)


dir.create(file.path(BorutaOut), showWarnings = FALSE)
dir.create(file.path(ModelObjs), showWarnings = FALSE)

#tmpIndices <- 1:length(newGeneSymbols)%%7==1;dataset=dataset[tmpIndices, ];newGeneSymbols=newGeneSymbols[tmpIndices]


trainModelByGene <- function(targetGene, dataset, newGeneSymbols, useAllGenes=F){
	#targetGene <- "BTK";useAllGenes=F	
	#A list, gene symbol as name, model as value
	geneModels <- list();
	
	#model prediction results for all genes, result for a single gene if the script is used to predict one gene
	pred4AllGenes <- data.frame()
	
	allGenes <- targetGene
	if(useAllGenes){
		allGenes <- c("all_genes")
	}
	for(gene in allGenes){
		#gene <- allGenes
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

		
		#Use Boruta to select important features
		boruta_output <- Boruta(Group ~ ., data=na.omit(geneData), doTrace=2, maxRuns = 50)  
		boruta_signif <- names(boruta_output$finalDecision[boruta_output$finalDecision %in% c("Confirmed", "Tentative")]) 
		DataFeaturesSelected <- selectedData[,c(boruta_signif, "Group")]
		
		#Output and visualize Boruta results
		borutaPrefix <- paste0(BorutaOut, "/ML.results.", gene, ".performance")
		pdf(file=paste0(borutaPrefix, ".boruta_output.pdf"), width=5, height=15)
		par(mar=c(4,11,1,1))
		plot(boruta_output, cex.axis=0.8, las=2, ylab="", xlab="", main="", horizontal=T)  
		mtext(paste("Importance of SNV features\n in", gene), side = 1, line = 2.7, cex=0.9)
		dev.off()
		write.table(boruta_signif, file=paste0(borutaPrefix, ".v_imp.sorted.2txt"),sep="\t",quote=F, row.names=F)
		write.table(DataFeaturesSelected, file=paste0(borutaPrefix, ".DataFeaturesSelected.txt"),sep="\t")
		write.table(boruta_output$ImpHistory, file=paste0(borutaPrefix, ".ImpHistory.txt"),sep="\t",quote=F)
		print("Dimension of feature selected training dataset:")
		print(dim(DataFeaturesSelected))

		#Select the best mtry parameters
		bestmtry <- tuneRF(DataFeaturesSelected[,1:(dim(DataFeaturesSelected)[2]-1)], DataFeaturesSelected[,c("Group")], stepFactor=2.0, improve=1e-5, ntree=500)
		bestmtry_val <- bestmtry[order(bestmtry[,2]),1][1]
		max.mtry <-20
		if(bestmtry_val > max.mtry){
			bestmtry_val <- max.mtry
		}
		control <- trainControl(method="cv", 
				        summaryFunction=twoClassSummary, 
				        classProbs=T,
				        savePredictions=T)
		metric <- "Accuracy"
		tunegrid <- expand.grid(.mtry=c(bestmtry_val))
		
		#Model training
		fit.rf <- train(Group~., data=DataFeaturesSelected, weights=weights,method="cforest", metric=metric, trControl=control, tuneGrid=tunegrid)
		#geneModels[[gene]] <- fit.rf
		targetPred <- fit.rf$pred[fit.rf$pred$rowIndex>=needed_index_cutoff, ]
		pred4AllGenes <- rbind(pred4AllGenes, targetPred)
		#pred_out <- predict(fit.rf, newdata=geneData, type="prob")
		#GeneObs <- geneData[,"Group"]
		#pred4AllGenes <- rbind(pred4AllGenes, cbind(pred_out, GeneObs) )
		saveRDS(fit.rf$pred, file = paste0(ModelObjs, "/", gene, ".pred.rds"))	
		saveRDS(fit.rf$finalModel, file = paste0(ModelObjs, "/", gene, ".finalModel.rds"))
	}
	colnames(pred4AllGenes)[3] <- "obs"
	list(pred_df=pred4AllGenes)
}

geneSpecificModelData <- trainModelByGene(targetGene, dataset, newGeneSymbols)

