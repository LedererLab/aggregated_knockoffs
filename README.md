# Aggregating Knockoffs for False Discovery Rate Control With an Application to Gut Microbiome Data

This repository provides the data and implementations of the methods described in the paper
*Aggregating Knockoffs for False Discovery Rate Control With an Application to Gut Microbiome Data* by
Fang Xie and Johannes Lederer.

## Usage

The file `lib/AggregatedKnockoffs.R` contains the functions `AKO` and `AKO.m` to implement the `AKO` and `modified AKO` methods described in the  paper.
 
The file `lib/stat.R` contains the functions `stat.ToDoFDRCon`, which is to compute the test statistic W described in the paper. The argument `stat.type` in the function `stat.ToDoFDRCon` gives four choices `max_lambda_lasso`, `max_lambda_scad`, `max_lambda_logistic` and `estimator_logistic`, which stand for four different test statistics.


The file `lib/myPlot.R` contains the function `myPlot` to display the relationships between target FDR and actual FDR and between target FDR and power.

## Simulations
We provide the files `Simulation1.R` and `Simulation2.R` for comparing the FDR and power of `KO` with those of our pipelines `AKO` and `modified AKO`. `Simulation1.R` is designed for lasso and `Simulation2.R` is designed for l1 regularized logistic regression.

## Real Data Analysis
The processed data in `RealData/ag-cleaned_L2.txt` is download from the website of the American Gut Project (http://americangut.org).

The files `lib/DataUS012018.R`,`lib/ImportData.R` and `SelectIndLinCols.R` are used to import and transform the bacteria phyla and BMI data of United State from `RealData/ag-cleaned_L2.txt`. The data we used contains 8404 samples and 56 phyla.

We provide the code `AGPAnalysis.R` to estimate the bacteria phyla that may have an influence on the obesity.

## Repository Authors

* **[Fang Xie](fang.xie@rub.de)** &mdash; Postdoc in Statistics, Ruhr-University Bochum &mdash; *methodology and `R` implementation*
* **[Johannes Lederer](johannes.lederer@rub.de)** &mdash; Professor in Statistics, Ruhr-University Bochum &mdash; *methodology and `R` implementation*


## Acknowledgements
The two files `lib/ImportData.R` and `lib/TransformData.R` are written by Lun Li and Johannes Lederer.

## Reference

[Aggregating Knockoffs for False Discovery Rate Control With an Application to Gut Microbiome Data](https://arxiv.org/abs/1907.03807)

Cite as "Xie and Lederer, *Aggregating Knockoffs for False Discovery Rate Control With an Application to Gut Microbiome Data,* arXiv:1907.03807, 2019"


