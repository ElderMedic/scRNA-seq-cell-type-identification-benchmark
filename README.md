# scRNA-seq Cell Type Identification Benchmark
## Github repo for thesis project 
> Discussion on Cell Type Identification Methods for Single-cell RNA Sequencing Data 

> 有关单细胞RNA-seq的细胞标识方法探讨

## Abstract
>In recent years, single-cell RNA sequencing technology has developed rapidly and is also widely used in biological research. Modeling and identifying the cell types based on the characteristics of single-cell RNA-seq data is one of the important challenges in the field of bioinformatics. In this regard, supervised classification models allow for the automatic and efficient labeling of cells, and as a result a number of methods based on supervised classification have been developed. 

>This paper investigates the currently published methods for cell type identification, designs and constructs test experiments, systematically analyzes and compares them from various aspects, and aims to provide researchers with some references or suggestions. First, we selected 10 popular automatic cell identification tools according to their underlying algorithm. Then, we used human pancreatic and mouse Tabula Muris single-cell RNA-seq data to comprehensively evaluate the above methods from multiple aspects such as accuracy, stability, generalization ability, operating efficiency, and other factors (including cell number, feature selection, and batch effects). We find that most of the tested methods can complete the task well, but each has its own advantages and disadvantages. To summarize, we recommend SVM rejection, SingleCellNet. SingleR and CHETAH. We also discussed future directions that may address the problem of cell clustering and classification, and shed some light on subsequent methodological studies and data analysis. 

>Keywords: scRNA-seq; cell type identification; classifier; machine learning 

## Description
- All python and R script used to run and benchmark tools in the study is presented in "Script" folder.
- For reproduction: run Cross\_Validation for each scRNA-seq dataset first, if you need feature selection, run "rank\_gene\_dropout.py" now or never, then run the exact "run_XXX.X" script you want, use "evaluate.R" to generate comprehensive evaluation result for the classifier finally.
- Test results and figure raw data is stored in "Results_preliminary" and "scRNA-seq Benchmark datasets" folder.
- For stability and generalization ability tests, all required input file are already in place.
- Due to the lack of computational resource, some experiments planned in advance have never been actually tested, we are glad if anyone can make them happen. 

### How to Use
- Generate training and test cell indices: 
	Cross-validation R script → .RData file
- Input: Expression matrix and true cell-type labels .csv file, RData file above
- Feature selection: run gene dropout rate ranking Python script to get gene indices
- Run each classifier wrapper in R/Python with default parameter
- Evaluate the classifier prediction by metrics calculation R script
- Save recorded data and draw figures and tables
