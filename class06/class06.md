Class 6 R Functions
================
Terry Lin
10/18/2019

``` r
#can do code here
```

text editor

# This is h1 big

## this is H2 small

### this is h3 smaller plain txt **bold**

Test 1

``` r
test1<-read.table("C:/Users/Terry/Desktop/R/class06/test1.txt", sep=",",header=TRUE)
View(test1)
```

Test 2

``` r
test2<-read.table("C:/Users/Terry/Desktop/R/class06/test2.txt", sep="$",header=TRUE)
View(test2)
```

Test 3

``` r
test3<-read.table("C:/Users/Terry/Desktop/R/class06/test3.txt", sep="",row.names="V1")
View(test3)
```

Function

``` r
add<- function(x, y=1) {
  #sum of x and y
  x+y
}
add(5,5)
```

    ## [1] 10

rescale

``` r
rescale <- function(x){
  rng<-range(x, na.rm=TRUE)
  (x-rng[1])/(rng[2]-rng[1])
}
rescale(c(1,10,NA,3))
```

    ## [1] 0.0000000 1.0000000        NA 0.2222222

``` r
rescale3 <- function(x, na.rm=TRUE, plot=FALSE) {
 if(na.rm) {
 rng <-range(x, na.rm=TRUE)
 } else {
 rng <-range(x)
 }
 print("Hello")
 answer <- (x - rng[1]) / (rng[2] - rng[1])
 
 print("is it me you are looking for?")
 if(plot) {
 plot(answer, typ="b", lwd=4)
 }
 print("I can see it in ...")
return (answer)
}

rescale3(1:10)
```

    ## [1] "Hello"
    ## [1] "is it me you are looking for?"
    ## [1] "I can see it in ..."

    ##  [1] 0.0000000 0.1111111 0.2222222 0.3333333 0.4444444 0.5555556 0.6666667
    ##  [8] 0.7777778 0.8888889 1.0000000
