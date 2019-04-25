#
# (c) Fritjof Lammers
#

# Load librares
library(scales)
library(knitr)
library(reshape2)
library(ggplot2)
library(ggthemes)
library(plyr)



#
# Change sign of flank distance depending on origin (3'prime, 5' prime or combined)
#
change_flanksign <- function (x, flank_pos) {
  if ( x == "5prime")
  {
    return(flank_pos*(-1))
  }
  else if ( x == "3prime")
  {
    return(flank_pos)
  }
  else if ( x == "combined")
  {
    return(flank_pos)
  }
}

#
# identify NaN in dataframe
#

is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))


# Function: load_data
# Loads datafiles from different flank-specific directories and saves as dataframe
#

load_data <- function(in_dir, start_flank = 1000, end_flank = 10000, window_size = 1000, tree_analysis = FALSE) {
  data <- list()
  for (flank_pos in as.integer(seq(from = start_flank, to = end_flank, by = window_size))) {
    infile <- paste(in_dir,as.character(flank_pos),"combined_sites.tsv", sep="/") # non_ref_ins/4_taxa

#    infile <- paste(in_dir,as.character(flank_pos),"sites.csv", sep="/") # non_ref_ins/4_taxa
    
    df.temp <- read.table(infile, header=TRUE, sep="\t", stringsAsFactors=FALSE)
    df.temp$locus <- paste(df.temp$seqname, ":", df.temp$start, "-", df.temp$stop, sep="")
    df.temp$dset <- lapply(df.temp$flank, function (x) change_flanksign(x, flank_pos))
    df.temp$dset <- as.integer(df.temp$dset)
    df.temp$clade_support <- as.character(lapply(df.temp$taxa, function (x) get_clade_support(x, FALSE)))
    df.temp <- df.temp[complete.cases(df.temp), ]
    # df.temp$spTratio <- (df.temp$sptreeSNP -  df.temp$incSNP) / (df.temp$incSNP + df.temp$sptreeSNP)
    df.temp$spinsratio <- as.numeric(df.temp$sptreeSNP -  df.temp$insSNP) / as.numeric(df.temp$insSNP + df.temp$sptreeSNP)
    if (tree_analysis == FALSE) {
        df.temp$insratio <- df.temp$insSNP / df.temp$seg_sites
        df.temp$alt1_SNV <- df.temp$p1SNP / df.temp$seg_sites
        df.temp$alt2_SNV <- df.temp$p2SNP / df.temp$seg_sites
        df.temp$alt3_SNV <- df.temp$p3SNP / df.temp$seg_sites
        df.temp$alt4_SNV <- df.temp$p4SNP / df.temp$seg_sites
        df.temp$alt5_SNV <- df.temp$p5SNP / df.temp$seg_sites
        }
    df.temp <- rbind(df.temp, data.frame())
    df.temp[is.nan(df.temp)] <- 0
    data[[flank_pos]] <- df.temp
  }

  compare_df <- data.frame(do.call("rbind", data))
  compare_df <- compare_df[compare_df$start >= (end_flank+1000), ] # what does this do?
  compare_df <- compare_df[apply(compare_df, 1, function(x) genomefile[genomefile$chrom == x[["seqname"]][[1]], "end"] > as.numeric(x[["stop"]][[1]])+end_flank+1000), ]

  return(compare_df)
}

# Function: kb_formatter
# format axes as kilobasepairs

kb_formatter <- function(x)
  {
  lab <- sprintf('%01dkb',(x/1000))
  }

# Plot 1: "Explore" A histogram SNp counts

plot_explore <- function(data, data_label) {
  p <- ggplot(data, aes(x=nuc_div,  group=flank))
  p <- p +   geom_histogram(alpha=.5,  position= "dodge",  aes(group=flank, fill=flank, color=flank)) +
    ggtitle(label=paste("Histograms of SNPs (",data_label,")")) + theme_bw()
  return(p)
}

# Plot 2:  Plot frequencies
#

plot_freqs <- function(data, data_label) {
  p <- ggplot(dfm, aes(x=dset, y=value))
  p <- p + geom_boxplot(aes(group=dset, fill=variable), notch=FALSE) +
    scale_y_continuous() +
    scale_x_continuous(breaks=seq((-1)*max_length,max_length, by=1000), label=kb_formatter) +
     theme(axis.text.x  = element_text(angle=90, vjust=0.5, size=10),
          axis.text.y  = element_text(angle=0, vjust=0.5, size=14)) +
  ggtitle(label=paste("Mean SNP frequencies per nucleotide", data_label)) +
  theme_bw() + facet_grid(variable ~ clade_support, scales="free")
  return(p)
}


# Plot 3: Scatterplot of means
#

plot_context <- function(data, data_label) {
  p <- ggplot(data, aes(x=dset, y=value, color=variable, fill=variable))
  p <- p + geom_point() +
            geom_smooth(aes(group=flank),method=loess) +
            theme(strip.text = element_text(angle=0, vjust=0.5, size=9)) +
    ggtitle(label=paste("Flanking regions of TE insertions", data_label)) +
    facet_wrap(clade_support ~ variable, scales="free_y", ncol=3)

  return(p)
  }


# Plot 3b: "MANUSCRIPT" Scatterplot of means
#

plot_context_MS <- function(data, data_label) {
p <- ggplot(data, aes(x=dset, y=value, color=variable, fill=variable))
p <- p + geom_point() +
            theme_base() +
          geom_smooth(aes(group=flank), fill="#00821f", alpha=.5, method=loess) +
          theme(strip.text = element_text(angle=0, vjust=0.5, size=9), panel.grid = element_line(color="#CCCCCC", size=0.3)) +
  ggtitle(label=paste("Flanking regions of TE insertions", data_label)) +
  facet_wrap(clade_support ~ variable, scales="free", ncol=3) +
  scale_fill_grey() + scale_color_grey() +
  scale_y_continuous(name="Frequency", breaks = c(0.0,0.02, 0.04, 0.06)) +
  scale_x_continuous(labels= function (x) x/1000, name="Distance from insertion site (kb)") +
  ylim(c(-0.01,0.06))
return(p)
}


plot_context_alt_topologies_MS <- function(data, data_label) {
p <- ggplot(data, aes(x=dset, y=value, fill=variable))
p <- p + geom_point() +
            theme_base() +
          geom_smooth(aes(group=flank), color="#000000", fill="#00821f", alpha=.5, method=loess) +
          theme(strip.text = element_text(angle=0, vjust=0.5, size=9), panel.grid = element_line(color="#CCCCCC", size=0.3),
                legend.position =  "none")+
  ggtitle(label=paste("Flanking regions of TE insertions", data_label)) +
  facet_wrap(clade_support ~ variable, scales="free", ncol=5) +
  scale_fill_grey() + scale_color_grey() +
  ylab("TE-supporting SNV frequency") +
  scale_y_continuous(name="Frequency", breaks = c(0.0,0.02, 0.04, 0.06)) +
  scale_x_continuous(labels= function (x) x/1000, name="Distance from insertion site (kb)") +
  ylim(c(-0.01,0.06))
return(p)
}


# Plot 4:  Tree plots scatterplot
do.monophyly.scatterplot <- function(data) {
  p <- ggplot( data, aes(x=allSNP,  y=insSNP, color=ins_monophyletic))
  p <- p + geom_point()  + facet_wrap(~ flank) + geom_smooth(method=lm, aes(group=ins_monophyletic)) +
    theme_classic()
  return(p)
}

do.monophyly.histogram <- function(data) {
  p <- ggplot( data, aes(x=(insSNP/allSNP), fill=ins_monophyletic))
  p <- p + geom_histogram()  + facet_wrap(~ flank) +
    theme_classic() + scale_fill_tableau()
  return(p)
}

# Plot 5: Barplot
do.monophyly.barplot <- function(data) {
  p <- ggplot(data, aes(x=ins_monophyletic, fill=ins_monophyletic))
  p <- p +  geom_bar(color="black") + facet_wrap(~ flank) +
    theme_classic() + scale_fill_grey() + theme(axis.text.x  = element_text(angle=90))

  return(p)
}


# Barplot with percentages
translate_flanks <- function(x) {
    if (x == "combined") {
        return("5' + 3'")
    }
    else {
        return(sub("prime", "'", x))
    }

}


do.monophyly.barplot.percent <- function(data, title = NULL) {
    data$flank <- as.factor(as.character((lapply(data$flank, translate_flanks))))
    data <- transform(data, flank = factor(flank,
        levels = c("5'", "3'", "5' + 3'"), ordered = TRUE))
    data <- transform(data, ins_monophyletic = factor(ins_monophyletic,
        levels = c("True", "False"), ordered = TRUE))
    head(data)
    df <- count(data, vars=c('ins_monophyletic', 'flank'))
    # next: compute percentages per group
    df <- ddply(df, .(flank), transform, p = freq/sum(freq))
    head(df)
    p <- ggplot(df, aes(x=flank, y=p, fill=ins_monophyletic))
    p <- p +    geom_bar(stat="identity", color="black") +
                theme_classic() +
                scale_fill_manual(values = c("True" = "#00821f", "False" = "#FFFFFF", "#82001c"))+
                scale_y_continuous(labels = percent_format()) +
                theme(axis.text  = element_text(size=8)) +
                ggtitle(title)
    print(p)

  return(p)
}


do.monophyly.vertical.barplot.percent <- function(data, title = NULL) {
    data$flank <- as.factor(as.character((lapply(data$flank, translate_flanks))))
    data <- transform(data, flank = factor(flank,
        levels = c("5' + 3'","3'","5'" ), ordered = TRUE))
    data <- transform(data, ins_monophyletic = factor(ins_monophyletic,
        levels = c("True", "False"), ordered = TRUE))
    head(data)
    df <- count(data, vars=c('ins_monophyletic', 'flank'))
    # next: compute percentages per group
    df <- ddply(df, .(flank), transform, p = freq/sum(freq))
    head(df)
    p <- ggplot(df,aes(x = flank, y=p,fill=ins_monophyletic))
    p <- p +    geom_bar(stat="identity", color="black") +
                theme_classic() +
                scale_fill_manual(values = c("True" = "#00821f", "False" = "#FFFFFF", "#82001c"))+
                scale_y_continuous(labels = percent_format()) +
                coord_flip() +
                theme(axis.text  = element_text(size=16)) +
                ggtitle(title)
    print(p)
    #scale_fill_grey(start=0.3,end=0.8) +

  return(p)
}
