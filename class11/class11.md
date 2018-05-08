---
title: "Bioinformatics Class 11"
author: "Barry Grant"
date: "5/8/2018"
output: 
  html_document: 
    keep_md: yes
---



## PDB statistics

Downlad CSV file from PDB database <http://www.rcsb.org/stats/summary>. Read this into R and determine fraction of X-ray structures.


```r
pdb.stats <- read.csv("Data Export Summary.csv")
```

Lets claculate something


```r
precent <- (pdb.stats$Total / sum(pdb.stats$Total) ) * 100
names(precent) <- pdb.stats$Experimental.Method
precent
```

```
##               X-Ray                 NMR Electron Microscopy 
##         89.51253323          8.72181096          1.50770286 
##               Other        Multi Method 
##          0.17006317          0.08788979
```

