# flankTE

Scripts to extract SNPs in the flanking regions of TE insertions. 
SNPs are extracted from whole-genome consensus sequences of the sample WGS datasets used for the TE calling analyses.
 

There are two types of the scripts. 

1. `flankTE-SNP.py` extract various types of SNP that support specific groupings of the samples, that need to be
configured in the script itself. 
2. `flankTE-tree.py` Extract SNPs that according to various branches on a species tree. 
This also requires configuration in the script file. 
 
## Dependencies

Resolve dependencies i.e. using *pip*


```bash
pip install biopython pyfaidx dendropy ete2
```

flankTE runs on Python 2.7


## Usage

```bash
python flankTE-tree.py 
usage: flankTE-tree.py [-h] -i FASTA -f FLANK [-o OVERLAPPING]
                                [-w WINDOWSIZE] -m MATRIX
```

### Command line arguments

- `-i ` Whole-genome consensus sequence (fasta format), supplied for each taxa that is analysed,
 e.g. `-i A.cons.fa -i B.cons.fa ...`
- `-f ` Size of flanking regions in base pairs (i.e. 10000). Note that the first 500 bp surrounding the breakpoint will
not be analysed in order to account for spuriously mapped reads. 
- `-o` True/[False]. Indicate whether the window size should be extended in each round (True) or
 if windows should be distinct regions of size [windowssize] (This is the default) 
- `-w` Window size in base pairs [1000]
- `-m` Tab-separated files with presence information on TE insertions (see below for details).


#### Overlapping windows

The default behaviour is to count SNPs in distinct windows of size W surrounding the insertion site (breakpoint).


```
|<--W-4-->|<--W-3-->|<--W-2-->|<--W-1-->|<500bp omitted><BREAKPOINT><500bp omitted>|<--W+1-->|<--W+2-->|<--W+3-->|<--W+4-->|

```

With `-o True` the program returns the increasing window sizes. This is most meaningful if 
the window size equals the flanking size, i.e. only a specific window around the breakpoint is investigated.

```
Window 1:         |<--W-->|<500bp omitted><BREAKPOINT><500bp omitted>|<--W-->|

Window 2: |<------W------>|<500bp omitted><BREAKPOINT><500bp omitted>|<------W------>|

etc...
```

#### Format of matrix file

```bash
chr1    10000  10001    sampleA,sampleB
chr1    20000  20001    sampleA,sampleC
...

