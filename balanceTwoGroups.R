#This function take output from readSNVdata.R, balance two classes, to make the nonPID/PID mutation ratio roughly equals to ctrlCaseRatio
getIndex4BalancedData <- function(dataList, ctrlCaseRatio = 1.2){
	set.seed(8)
	df_SNV_Box <- dataList$df_SNV_Box
	varCards_data_withNA <- dataList$varCards_data_withNA
	geneSymbols <- dataList$geneSymbols
	dataset <- dataList$dataset
	    
	labels <- dataset$Group
	geneGroup <- cbind(geneSymbols, as.character(labels))
	geneSet <- unique(geneGroup[,1])

	all.selected.indices <- c()
	
	for(gene in geneSet){
		row.indices <- 1:dim(geneGroup)[1]
		#get the original indices of nonPID and PID mutations
		control.bool.indices <- geneGroup[,1]==gene & geneGroup[,2]=="Non_PID"
		pid.bool.indices <- geneGroup[,1]==gene & geneGroup[,2]=="PID"

		pid.count <- sum(pid.bool.indices)
		control.count <- sum(control.bool.indices)
		#number of nonPID mutations > ctrlCaseRatio * PID_mutations
		if(control.count > pid.count * ctrlCaseRatio){
			ctrl.index.selected <- sample(row.indices[control.bool.indices], pid.count*ctrlCaseRatio)
		#number of nonPID mutations <= ctrlCaseRatio * PID_mutations
		}else{
			ctrl.index.selected <- row.indices[control.bool.indices]
		}
		#keep all PID mutations.
		#all.selected.indices contain indices for multiple genes
		all.selected.indices <- c(all.selected.indices, ctrl.index.selected, row.indices[pid.bool.indices])
	}
	selectedIndices <- sort(all.selected.indices)
	list(df_SNV_Box=df_SNV_Box[selectedIndices,],
	     varCards_data_withNA=varCards_data_withNA[selectedIndices,],
	     geneSymbols=geneSymbols[selectedIndices],
	     dataset=dataset[selectedIndices,]
	)

}





