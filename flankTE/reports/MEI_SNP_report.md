Report: SNVs in TE loci 
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
test_types <- vector()
test_types["SINE"] = "/home/newuser/teddypi/data/2015-12-15/MEI_SNPs_analysis/new072016/window_based2/combined/"

df_collection <- list()
for (in_dir in names(test_types)) {
  df_collection[[in_dir]] <- load_data(test_types[[in_dir]])
}
element <- "SINE"
```

## Show a histogram


```
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![plot of chunk unnamed-chunk-3](figure/unnamed-chunk-3-1.png)

|flank  |clade_support |    N|   allSNP|       sd|        se|        ci|
|:------|:-------------|----:|--------:|--------:|---------:|---------:|
|3prime |american      | 97.7| 16.21699| 5.998634| 0.1919133| 0.3766102|
|3prime |asb_amb       | 99.7| 17.65496| 6.276215| 0.1987697| 0.3900555|
|3prime |asian         | 97.5| 17.15487| 6.121931| 0.1960587| 0.3847461|
|3prime |bro_amb       | 95.7| 16.74922| 6.095973| 0.1970548| 0.3867100|
|3prime |other         | 97.6| 16.65164| 6.208807| 0.1987391| 0.3900056|
|3prime |pol_bro       | 99.4| 16.36922| 6.127998| 0.1943683| 0.3814198|
|3prime |slo_asb       | 99.0| 16.99697| 6.066283| 0.1927991| 0.3783424|
|3prime |sun_asb       | 99.7| 17.55166| 6.448949| 0.2042403| 0.4007906|
|3prime |sun_slo       | 99.1| 16.98890| 6.198801| 0.1969114| 0.3864116|
|5prime |american      | 99.2| 16.03831| 6.206223| 0.1970478| 0.3866788|
|5prime |asb_amb       | 99.3| 16.81873| 6.150611| 0.1951837| 0.3830204|
|5prime |asian         | 98.4| 16.97154| 6.143354| 0.1958430| 0.3843184|
|5prime |bro_amb       | 98.4| 16.86789| 5.814251| 0.1853516| 0.3637302|
|5prime |other         | 97.8| 16.65644| 5.978344| 0.1911664| 0.3751439|
|5prime |pol_bro       | 98.6| 16.82252| 5.846943| 0.1862046| 0.3654033|
|5prime |slo_asb       | 97.7| 16.92426| 6.162120| 0.1971437| 0.3868743|
|5prime |sun_asb       | 99.6| 17.50602| 6.165634| 0.1953656| 0.3833758|
|5prime |sun_slo       | 99.1| 17.10999| 6.525708| 0.2072959| 0.4067899|


## Plot everything
![plot of chunk fig](figure/fig-1.png)

```
## Scale for 'y' is already present. Adding another scale for 'y', which
## will replace the existing scale.
```

```
## Warning: Removed 360 rows containing non-finite values (stat_smooth).
```

```
## Warning: Removed 360 rows containing missing values (geom_point).
```

![plot of chunk fig](figure/fig-2.png)

```
## Error in paste("~/teddypi/results//MS/MEI_SNP/insratio", dset, "pdf", : object 'dset' not found
```


## Regression analysis
![plot of chunk unnamed-chunk-4](figure/unnamed-chunk-4-1.png)

|    |flank    |clade    |    r_sqrd|    Pvalue|
|:---|:--------|:--------|---------:|---------:|
|2   |5prime   |american | 0.0050274| 0.0255354|
|21  |3prime   |american | 0.0155182| 0.0000946|
|3   |combined |american | 0.0202634| 0.0000062|
|4   |5prime   |asb_amb  | 0.0110978| 0.0008851|
|5   |3prime   |asb_amb  | 0.0200030| 0.0000074|
|6   |combined |asb_amb  | 0.0255337| 0.0000004|
|7   |5prime   |asian    | 0.0334343| 0.0000000|
|8   |3prime   |asian    | 0.0279007| 0.0000002|
|9   |combined |asian    | 0.0519884| 0.0000000|
|10  |5prime   |bro_amb  | 0.0308198| 0.0000000|
|11  |3prime   |bro_amb  | 0.0146336| 0.0001760|
|12  |combined |bro_amb  | 0.0329226| 0.0000000|
|13  |5prime   |other    | 0.0089678| 0.0030329|
|14  |3prime   |other    | 0.0028388| 0.0961944|
|15  |combined |other    | 0.0076178| 0.0058456|
|16  |5prime   |pol_bro  | 0.0046819| 0.0316858|
|17  |3prime   |pol_bro  | 0.0029263| 0.0882689|
|18  |combined |pol_bro  | 0.0071172| 0.0076024|
|19  |5prime   |slo_asb  | 0.0065519| 0.0113743|
|20  |3prime   |slo_asb  | 0.0306894| 0.0000000|
|211 |combined |slo_asb  | 0.0369938| 0.0000000|
|22  |5prime   |sun_asb  | 0.0036479| 0.0567190|
|23  |3prime   |sun_asb  | 0.0111715| 0.0008302|
|24  |combined |sun_asb  | 0.0110393| 0.0008760|
|25  |5prime   |sun_slo  | 0.0254934| 0.0000004|
|26  |3prime   |sun_slo  | 0.0242567| 0.0000008|
|27  |combined |sun_slo  | 0.0438093| 0.0000000|

```
## Using r_sqrd as value column: use value.var to override.
```



========  =========  =========  =========
clade        3prime     5prime   combined
========  =========  =========  =========
american  0.0155182  0.0050274  0.0202634
asb_amb   0.0200030  0.0110978  0.0255337
asian     0.0279007  0.0334343  0.0519884
bro_amb   0.0146336  0.0308198  0.0329226
other     0.0028388  0.0089678  0.0076178
pol_bro   0.0029263  0.0046819  0.0071172
slo_asb   0.0306894  0.0065519  0.0369938
sun_asb   0.0111715  0.0036479  0.0110393
sun_slo   0.0242567  0.0254934  0.0438093
========  =========  =========  =========

# Manuscript plot


```r
df <- subset(df_collection[[element]], dset != "" & flank != "combined" , select=c("dset", "insratio", "clade_support", "flank"))
dfm <- melt(df, id.vars=c("dset", "clade_support", "flank"))
insSNP_summary <- summarySE(dfm,groupvars=c("dset","variable", "clade_support", "flank"), measurevar="value")

plot_context_MS(insSNP_summary, "SINE")
```

```
## Scale for 'y' is already present. Adding another scale for 'y', which
## will replace the existing scale.
```

![plot of chunk unnamed-chunk-5](figure/unnamed-chunk-5-1.png)

## Manuscript plot - alternative topologies

```r
pdf(file=paste("~/teddypi//results/MS/MEI_SNP/alternative_topos.pdf", sep=""),
    width=11.69, height=3.5 )
  
df <- subset(df_collection[[element]], dset != "" & flank != "combined" , select=c("dset", "alt1_SNV", "alt2_SNV", "alt3_SNV", "alt4_SNV", "alt5_SNV", "clade_support", "flank"))
for (signal in levels(as.factor(df$clade_support))) {
  
  dfs <- df[df$clade_support == signal, ]
  dfm <- melt(dfs, id.vars=c("dset", "clade_support", "flank"))
  insSNP_summary <- summarySE(dfm,groupvars=c("dset","variable", "clade_support", "flank"), measurevar="value")

  print(plot_context_alt_topologies_MS(insSNP_summary, "SINE"))
}
```

```
## Scale for 'y' is already present. Adding another scale for 'y', which
## will replace the existing scale.
## Scale for 'y' is already present. Adding another scale for 'y', which
## will replace the existing scale.
## Scale for 'y' is already present. Adding another scale for 'y', which
## will replace the existing scale.
## Scale for 'y' is already present. Adding another scale for 'y', which
## will replace the existing scale.
## Scale for 'y' is already present. Adding another scale for 'y', which
## will replace the existing scale.
## Scale for 'y' is already present. Adding another scale for 'y', which
## will replace the existing scale.
## Scale for 'y' is already present. Adding another scale for 'y', which
## will replace the existing scale.
## Scale for 'y' is already present. Adding another scale for 'y', which
## will replace the existing scale.
## Scale for 'y' is already present. Adding another scale for 'y', which
## will replace the existing scale.
```

```r
dev.off()
```

```
## png 
##   2
```
