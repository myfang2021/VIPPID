# Model training for VIPPID v1.0

## About
VIPPID is a tool to predict PID (primary immunodeficiency) mutations. It is applicable to missense single nucleotide substitutions, and it uses information from [Varcards] (http://159.226.67.237/sun/varcards/about/download) Varcards and SNVBox [SNVBox] (https://karchinlab.org/apps/appSnvBox.html) for training. It trains gene specific model for genes with more than certain number (45 by default) of mutations, for other genes, non-gene-specific model will be trained.


## Dependencies
#### Tested environment:
* R version 3.6.3 (2020-02-29)
* Platform: x86_64-pc-linux-gnu (64-bit)
* Running under: Ubuntu 18.04.2 LTS
### Required packages:
* Boruta (7.0.0)
* pROC (1.17.0.1)
* ggfortify (0.4.11)
* Hmisc (4.4-2)
* missForest (1.4)
* caret (6.0-86)
* randomForest (4.6-14)

## Usage
### Prepare the input data
#### SNVBox annotation data files
The first type of input files are for the prediction information from SNVbox, after getting prediction output from SNVBox, you have to insert one column to the beginning, in the format of Chr_pos_refBase_altBase (e.g. chr10_14977469_C_T), and append a column to the end of file to indicate the class of the mutation (i.e. PID or benign), and add a header line like in the test data.
You have to generate one SNVBox data file for the PID mutations, and another file for the non-PID mutations.
#### Varcards annotation data files
The second type of input files are for the prediction information from Varcards.
Use VarcardsPreprocess/preprocessVarcardAnno4R.sh to preprocess the data from Varcards.
Usage:
```sh preprocessVarcardAnno4R.sh Varcard_input.txt <mutation class>```

### Train the model
The training takes two steps:
**Step 1**: train gene-specific model for each gene with enough mutations, then train non-gene-specific model. 
To train gene-specific model, run
```Rscript trainModel.R <gene symbol> <output directory>```
To train non-gene-specific model, run both:

```Rscript trainModel.R Other <output directory>```

```Rscript trainModel.R all_genes <output directory>```

### Example
```Rscript trainModel.R FAS ./test/```

or

```
for i in FAS PRF1 Other all_genes;do 
  /usr/bin/Rscript trainModel.R $i ./test/;
done
```
**Please note that you have to change the input files in trainModel.R before running the script.**
**Step 2**: make prediction, collect prediction results and draw ROC graph.
run:

```Rscript combineModelResults.R <output directory>```
### Example

```Rscript combineModelResults.R ./test/```

Then the output will be within folders BorutaOut and ModelObjs in the specified output directory.

