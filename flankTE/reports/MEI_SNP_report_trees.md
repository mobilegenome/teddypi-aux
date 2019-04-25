Report: SNVs in TE loci: Phylogenetic analysis
========================================================

Research report on the analysis of flanking regions of TE insertion loci identified using the TeddyPi pipeline. 


 1. Load the code 


```
## Error in eval(expr, envir, enclos): could not find function "eval"
```


2. Load the data and store some variables

```r
genomefile <- data.frame(read.table("~/teddypi/data/reference_genome/BGI.scaf.filtered_1Mb.genome"))
colnames(genomefile) <- c("chrom", "end")
start_flank = 10000
end_flank = 10000
window_size = 10000
test_types <- vector()

test_types["tree_SINE_asb_sun"] = "/home/newuser/teddypi/data/2015-12-15/MEI_SNPs_analysis/new072016/10kb_trees_2/SINE_asb_sun/"
test_types["tree_SINE_asb_slo"] = "/home/newuser/teddypi/data/2015-12-15/MEI_SNPs_analysis/new072016/10kb_trees_2/SINE_asb_slo/"
test_types["tree_SINE_asian"] = "/home/newuser/teddypi/data/2015-12-15/MEI_SNPs_analysis/new072016/10kb_trees_2/SINE_asian/"
test_types["tree_SINE_asb_amb"] = "/home/newuser/teddypi/data/2015-12-15/MEI_SNPs_analysis/new072016/10kb_trees_2/SINE_asb_amb/"
test_types["tree_SINE_sun_slo"] = "/home/newuser/teddypi/data/2015-12-15/MEI_SNPs_analysis/new072016/10kb_trees_2/SINE_sun_slo/"
test_types["tree_SINE_amb_bro"] = "/home/newuser/teddypi/data/2015-12-15/MEI_SNPs_analysis/new072016/10kb_trees_2/SINE_amb_bro/"
#test_types["tree_SINE_amb_pol"] = "/home/newuser/teddypi/data/2015-12-15/MEI_SNPs_analysis/new072016/10kb_trees_2/SINE_amb_pol/"
test_types["tree_SINE_amb_bro_pol"] = "/home/newuser/teddypi/data/2015-12-15/MEI_SNPs_analysis/new072016/10kb_trees_2/SINE_amb_bro_pol/"
#test_types["tree_SINE_bro_pol"] = "/home/newuser/teddypi/data/2015-12-15/MEI_SNPs_analysis/new072016/10kb_trees_2/SINE_bro_pol/"


df_collection <- list()
for (in_dir in names(test_types)) {
  df_collection[[in_dir]] <- load_data(test_types[[in_dir]], start_flank, end_flank, end_flank, tree_analysis = TRUE)
}
element <- "SINE"
```


3. Plots

```r
for (dset in names(df_collection)) {
  df <- df_collection[[dset]]
    summarySE(df, measurevar="insSNP",groupvars="flank")
  
  p1 <- do.monophyly.histogram(df) + ggtitle(dset)
  p2 <- do.monophyly.barplot(df)
  multiplot(p1,p2)
  
}
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

```
## Warning: Removed 1 rows containing non-finite values (stat_bin).
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-2.png)

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-3.png)

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-4.png)

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-5.png)

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

```
## Warning: Removed 3 rows containing non-finite values (stat_bin).
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-6.png)

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-7.png)


4. MS plot

```r
for (dset in names(df_collection)) {
df <- df_collection[[dset]]
p <- do.monophyly.barplot.percent(df, title = dset)
ggsave(paste("~/teddypi/results//MS/MEI_SNP/plot",dset,"pdf", sep="."),plot=p, device="pdf", width=100, height=50, units="mm")
}
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-2.png)![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-3.png)![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-4.png)![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-5.png)![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-6.png)![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-7.png)

```r
df <- df_collection[["tree_SINE_asb_amb"]]
p2 <- do.monophyly.barplot.percent(df)
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-8.png)

```r
multiplot(p1,p2)
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-9.png)

```r
df <- df_collection[["tree_SINE_sun_asb"]]
print(do.monophyly.barplot.percent(df))
```

```
## Error in factor(ins_monophyletic, levels = c("True", "False"), ordered = TRUE): object 'ins_monophyletic' not found
```

```r
df <- df_collection[["tree_SINE_slo_asb"]]
print(do.monophyly.barplot.percent(df))
```

```
## Error in factor(ins_monophyletic, levels = c("True", "False"), ordered = TRUE): object 'ins_monophyletic' not found
```

```r
df <- df_collection[["tree_SINE_sun_slo"]]
head(df)
```

```
##     seqname infSNP              taxa allSNP insDiv   stop absDiv  start
## 1 scaffold1     29 sun_bear,slo_bear    194    2.5  92109   12.5  92057
## 2 scaffold1     38 sun_bear,slo_bear    165    3.0  92109   17.0  92057
## 3 scaffold1     67 sun_bear,slo_bear    359    5.5  92109   29.5  92057
## 4 scaffold1     19 sun_bear,slo_bear    146    1.0 838175    9.5 838091
## 5 scaffold1     31 sun_bear,slo_bear    154    5.0 838175   15.0 838091
## 6 scaffold1     50 sun_bear,slo_bear    300    6.0 838175   24.5 838091
##      flank othSNP insSNP noSNP ins_monophyletic sptreeSNP baseSNP
## 1   5prime    165      4  8155             True         8       0
## 2   3prime    127      4  8440             True        23       0
## 3 combined    292      8 16595             True        31       0
## 4   5prime    127      0  7813             True         3       0
## 5   3prime    123      1  8082            False         5       0
## 6 combined    250      1 15895             True         8       0
##                     locus   dset clade_support spinsratio
## 1   scaffold1:92057-92109 -10000       sun_slo  0.3333333
## 2   scaffold1:92057-92109  10000       sun_slo  0.7037037
## 3   scaffold1:92057-92109  10000       sun_slo  0.5897436
## 4 scaffold1:838091-838175 -10000       sun_slo  1.0000000
## 5 scaffold1:838091-838175  10000       sun_slo  0.6666667
## 6 scaffold1:838091-838175  10000       sun_slo  0.7777778
```

```r
print(do.monophyly.barplot.percent(df))
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-10.png)![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-11.png)

```r
print(do.monophyly.vertical.barplot.percent(df))
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-12.png)![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-13.png)

```r
multiplot(p1,p2)
```

```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-14.png)
