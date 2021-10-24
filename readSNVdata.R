library(Hmisc)
set.seed(8)
#this function read the training data file, and perform preprocessing
readData <- function(PID_SNV_BOX_file, nonPID_SNV_BOX_file, PID_VarCards_file, nonPID_VarCards_file, out_dir){
	borutaPrefix <- out_dir
	#read the SNVbox annotated data for PID and nonPID mutations
	PID_SNV_BOX_raw <- read.csv(PID_SNV_BOX_file, sep="", stringsAsFactors=F, row.names=1)
	nonPID_SNV_BOX_raw <- read.csv(nonPID_SNV_BOX_file, sep="", stringsAsFactors = F, row.names=1)
	SNV_BOX_raw <- as.data.frame(rbind(PID_SNV_BOX_raw, nonPID_SNV_BOX_raw))

	#SNVbox data filtering
	SNV_BOX_raw2 <- SNV_BOX_raw[,-c(1,2,88)]
	non_na_feature_index <- apply(SNV_BOX_raw2,2,function(x){sum(is.na(x))})==0
	row_indices <- SNV_BOX_raw$Group != "Non_PID"
	df_SNV_Box <- SNV_BOX_raw2[row_indices,non_na_feature_index]

	#read Varcards annotated data
	PID_VarCards <- read.csv(PID_VarCards_file, sep="\t", stringsAsFactors = F, row.names=1, na.strings="-")
	nonPID_VarCards <- read.csv(nonPID_VarCards_file, sep="\t", stringsAsFactors = F, row.names=1, na.strings="-")

	#data preprocessing, use only the information about prediction results of existing tools
	both_VarCards <- rbind(PID_VarCards, nonPID_VarCards)
	colnames(both_VarCards)[colnames(both_VarCards) == "REVEL"] <- "REVEL_score"
	all_score_indices <- grepl("score", colnames(both_VarCards))
	revel_index <- colnames(both_VarCards) == "REVEL"
	geneScoreIndices <- grepl("Gene_score_", colnames(both_VarCards))
	rc_varCards_indices <- (all_score_indices | revel_index) & !geneScoreIndices & !colnames(both_VarCards) == "Damaging_score"

	#Rename the group labels
	#Missing value imputation
	group_labels <- as.factor(ifelse(SNV_BOX_raw$Group == "PID", "PID", "Non_PID")[row_indices])
	varCards_data_withNA <- both_VarCards[rownames(df_SNV_Box), rc_varCards_indices]
	geneSymbols <- both_VarCards[rownames(df_SNV_Box), "Gene_symbol"]
	moreDetailedGroup <- both_VarCards[rownames(df_SNV_Box), "Group"]

	impute_df <- function(x){
		for(i in 1:dim(x)[2]){
			x[,i] <- impute(x[,i])
		}
		x
	}
	varCards_pred_data <- impute_df(varCards_data_withNA)

	#Data for class balancing
	dataset <- cbind(df_SNV_Box, varCards_pred_data, group_labels)
	colnames(dataset)[dim(dataset)[2]] <- "Group"

	write.table(cbind(dataset,geneSymbols), file=paste0(borutaPrefix, ".unfiltered_unbalanced.all_SNVs.dataset.txt"), sep="\t", quote=F)
	list(df_SNV_Box=df_SNV_Box,
	     varCards_data_withNA=varCards_data_withNA, 
	     geneSymbols=geneSymbols, 
	     dataset=dataset
	     )
}
