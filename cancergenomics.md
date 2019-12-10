Untitled
================

## Epitope Database (IEDB)

``` r
library(bio3d)
library(maftools)
```

``` r
alig <- read.fasta("lecture18_sequences.fa")

cons <- conserv(alig, method="identity")
mut <- which(cons<1)

gaps <- gap.inspect(alig)
mut <- mut[mut %in% gaps$f.inds]

mut
```

    ## [1]  41  65 213 259

``` r
alig$ali[,mut]
```

    ##            [,1] [,2] [,3] [,4]
    ## P53_wt     "D"  "R"  "R"  "D" 
    ## P53_mutant "L"  "W"  "V"  "V"

``` r
alt1<-alig$ali[2,c((mut-8):(mut+8))]
```

    ## Warning in (mut - 8):(mut + 8): numerical expression has 4 elements: only
    ## the first used
    
    ## Warning in (mut - 8):(mut + 8): numerical expression has 4 elements: only
    ## the first used

``` r
write.fasta(seq=alt1, file="alt1.fa")
```
