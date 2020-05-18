---
output: html_document
---

# EPNClassifier
EPNClassifier is able to assign ependymal tumors profiled by transcriptomic techniques into their molecular subgroup. This method is based on gene set enrichment analysis and utilizes pre-determined gene signatures to calculate a running sum statistic and p-value. Significant p-values can be used to classify EPN tumors into their most-likely molecular subgroup or classify PFA tumors into their subtypes (PFA_1 and PFA_2). 

## Installation

```{r}
devtools::install_github("CamaraLab/EPNClassifier")
library(EPNClassifier)
```

## Tutorial

Example data for this tutorial is from [Pajtler K.W. et al., Cancer Cell (2015)](https://www.sciencedirect.com/science/article/pii/S153561081500135X?via%3Dihub) and can be
downloaded here ...

### Classify EPN into their molecular subgroups

Calculate the running sum statistic and p-vlaue for each EPN tumor

```{r}
classification <- ClassifyEPN(bulk = data, permutations = 10000)
classification$GSM1570843$pvalue
```

![](examples/figure-markdown_github/ClassifyEPN_output.png)

Plot the running sum statistic (ie. enrichment score ES) for a specific EPN sample.

```{r}
ES <- EnrichmentScore_EPN(bulk_sample = data[,1,drop=F])
plot(x = ES$ST_EPN_RELA$x, y = ES$ST_EPN_RELA$y, xlab = "Rank List of Genes", ylab = "Running Sum Statistic", main = "ST RELA Enrichment Score")
```

![](examples/figure-markdown_github/ESPlot_EPN.png)

Assign each sample to a single molecular subgroup

```{r}
molec_gps <- Classify(classification, min_pvalue = 0.35)
molec_gps
```
![](examples/figure-markdown_github/molec_gps.png)

### Classify PF_A tumors into their subgroups (PFA_1 and PFA_2)

```{r}
#Restrict to PFA tumors
pfa_data <- data[,names(molec_gps)[molec_gps %in% "PF_EPN_A"]]

classification_pf <- ClassifyPFA(bulk = pfa_data, permutatins = 10000)
classification_pf$GSM1570848$pvalue
```
![](examples/figure-markdown_github/ClassifyPFA_output.png)

Plot the running sum statistic (ie. enrichment score ES) for a specific EPN sample.

```{r}
ES <- EnrichmentScore_PFA(bulk_sample = pfa_data[,1,drop=F])
plot(x = ES$PFA_1$x, y = ES$PFA_1$y, xlab = "Rank List of Genes", ylab = "Running Sum Statistic", main = "PFA_1 Enrichment Score")
```

![](examples/figure-markdown_github/ESPlot_PFA.png)

Assign each sample to a single PFA subgroup

```{r}
pf_gps <- Classify(classification, min_pvalue = 0.35)
pf_gps
```
![](examples/figure-markdown_github/pfa_gps.png)



