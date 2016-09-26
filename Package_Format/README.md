#Package Format 

For users who prefer a package implementation, it can be found here. 

To estimate DNAm age, use predictor(dataset).

If you do not wish to normalize your dataset, set predictor(dataset, normalizedata=FALSE). Note that normalization is highly recommended. 

The dataset object contains CpG names in the rows and samples in the columns. The first column should be named "CpGName"

Results output to the working directory. 
