# Aggregating Knockoffs for False Discovery Rate Control With an Application to Gut Microbiome Data

This repository provides the data and implementations of the methods described in the paper
*Aggregating Knockoffs for False Discovery Rate Control With an Application to Gut Microbiome Data* by
Fang Xie and Johannes Lederer.

## Usage

The file `lib/AggregatedKnockoffs.R` contains a function `AKO` to implement the `AKO` and `AKO+` methods described in the  paper.
 
The file `lib/Lasso_max_lambda.R` contains two functions `Lasso_max_lambda` and `stat.Lasso_lambdasmax` that can obtain the statistic W of lasso described in the paper.

The file `lib/SCAD_max_lambda.R` contains two functions `SCAD_max_lambda` and `stat.SCAD_lambdasmax` that can obtain the statistic W of scad described in the paper.

The file `lib/Sigma.R` contains a function `Sig` to generate the covariance matrix of design matrix.

The file `lib/myPlot.R` contains two functions `myPlot.FP` and `myPlot.x.Target` to display the relationships between target FDR and selection accuarcy, and between target FDR and actual FDR.

## Simulation
We provide the file `SCAD_Lasso.R` for comparing the FDR and selection accuracy of `KO` and `KO+` with those of `AKO` and `AKO+` respectively in the cases of lasso and scad.

## Real Data Analysis
The processed data in `RealData/ag-cleaned_L2.txt` is download from the website of American Gut Project (http://americangut.org).

The files `lib/CenteringByGMM.R`,`lib/CLRTransformation.R`,`lib/DataUS012018.R`,`lib/ImportData.R`, `SelectIndLinCols.R` and `lib/TransformData.R` are used to import and transform the bacteria phyla and BMI data of United State from `RealData/ag-cleaned_L2.txt`. The data we used contains 8404 samples and 58 phyla.

We provide the code `RealDataAnalysis_AGP.R` to estimate the bacteria phyla that may have an influence on the BMI.

## Repository Authors

* **[Fang Xie](fang.xie@rub.de)** &mdash; Postdoc in Statistics, Ruhr-University Bochum &mdash; *methodology and `R` implementation*
* **[Johannes Lederer](johannes.lederer@rub.de)** &mdash; Professor in Statistics, Ruhr-University Bochum, &mdash; *methodology and `R` implementation*


## Acknowledgements

The raw AGP data was downloaded from http://americangut.org. 

## Reference

[Aggregating Knockoffs for False Discovery Rate Control With an Application to Gut Microbiome Data](https://arxiv.org/abs/1907.03807)

Cite as "Xie and Lederer, *Aggregating Knockoffs for False Discovery Rate Control With an Application to Gut Microbiome Data,* arXiv:1907.03807, 2019"


