Rdocumentation: model\_diffnet
================

\#Install required R-packages.

``` r
install.packages("glmnet",repos = "http://cran.us.r-project.org");library("glmnet")
install.packages("matrixStats",repos = "http://cran.us.r-project.org");library("matrixStats")
install.packages("qgraph",repos = "http://cran.us.r-project.org");library("qgraph")
```

\#Example of the sign-adjusted dCCN method.

\#Simulate a simple model consisting of two main effects, two type I
interactions, and one type II interaction.

``` r
#At first, simulate 1000 i.i.d samples of ten explanatory variables.
X <- matrix(0,1000,10)
for(i in 1:10){
  X[,i] <- rnorm(1000,0,1)}

#Let us consider a scenario, where the resistance of an individual for a particular disease is controlled, to some extent, by the sum of two interaction terms of type I.
resistance <- X[,1]*X[,5] + X[,2]*X[,8] + rnorm(1000,0,0.5)

#Suppose that there exist a regulatory relationship between genes X[,9] and X[,10] that is associated only with high resistance-levels - let's say all values above the median of "resistance" - and can be modeled as follows. 

X[which(resistance > quantile(resistance, 0.5)),10] <- X[which(resistance > quantile(resistance, 0.5)),9] + rnorm(500,0,0.2)

#Let us also assume that, among individuals with low resistance-levels, this regulatory mechanism is somehow disrupted and the expression levels of these genes become independent.

X[which(resistance < quantile(resistance, 0.5)),10] <-  rnorm(500,0,1)
X[which(resistance < quantile(resistance, 0.5)),9] <-  rnorm(500,0,1)


# Finally, simulate the trait associated with that particular disease as the sum of the previous resistance-levels and the main effects of two genes. NOTE that the "resistance" vector already contains the random error term.
y <- X[,3] + X[,7] + resistance

#Now the proposed model guided differential network estimation method can be used to find these disrupted regulation patterns, the interaction terms behind them, as well as the additional main effects (although not shown in the estimated differential network as they are estimated separately in the first step).
```

\#Apply the proposed model\_diffnet function

``` r
source("./model_diffnet.R")
dCCN_diff_net <- model_diffnet(X,y,cut_off_a = 1/2,alpha = 1/3,r = 0.1, corr_metric = "CCN", main_effect = TRUE, CV = TRUE,symmetric = TRUE, signAdj = TRUE, adj_matrices = FALSE)
qgraph(dCCN_diff_net)
```

![](Toy_example_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
dCCN_diff_net
```

    ##            [,1]      [,2]       [,3] [,4]       [,5]       [,6] [,7]      [,8]
    ##  [1,] 0.0000000 0.0000000  0.0000000    0  0.8767788  0.0000000    0 0.0000000
    ##  [2,] 0.0000000 0.0000000  0.0000000    0  0.0000000  0.0000000    0 0.9927316
    ##  [3,] 0.0000000 0.0000000  0.0000000    0  0.0000000 -0.1369188    0 0.0000000
    ##  [4,] 0.0000000 0.0000000  0.0000000    0  0.0000000  0.0000000    0 0.0000000
    ##  [5,] 0.8767788 0.0000000  0.0000000    0  0.0000000 -0.1110961    0 0.0000000
    ##  [6,] 0.0000000 0.0000000 -0.1369188    0 -0.1110961  0.0000000    0 0.0000000
    ##  [7,] 0.0000000 0.0000000  0.0000000    0  0.0000000  0.0000000    0 0.0000000
    ##  [8,] 0.0000000 0.9927316  0.0000000    0  0.0000000  0.0000000    0 0.0000000
    ##  [9,] 0.0000000 0.0000000  0.0000000    0  0.0000000  0.0000000    0 0.0000000
    ## [10,] 0.0000000 0.0000000  0.0000000    0  0.0000000  0.0000000    0 0.0000000
    ##             [,9]      [,10]
    ##  [1,]  0.0000000  0.0000000
    ##  [2,]  0.0000000  0.0000000
    ##  [3,]  0.0000000  0.0000000
    ##  [4,]  0.0000000  0.0000000
    ##  [5,]  0.0000000  0.0000000
    ##  [6,]  0.0000000  0.0000000
    ##  [7,]  0.0000000  0.0000000
    ##  [8,]  0.0000000  0.0000000
    ##  [9,]  0.0000000 -0.9642571
    ## [10,] -0.9642571  0.0000000

\#Example of the sign-adjusted dPCCN method. (a more complex type II
interaction scenario) \#Simulate a simple model consisting of two main
effects, two type I interactions, and one type II interaction.

``` r
#Simulate 1000 i.i.d samples of ten explanatory variables.
X <- matrix(0,1000,10)
for(i in 1:10){
  X[,i] <- rnorm(1000,0,1)}

#Let us consider the same scenario as above, and simulate the resistance vector as the sum of two interaction terms of type I.
resistance <- X[,1]*X[,5] + X[,2]*X[,8] + rnorm(1000,0,0.5)

#Now let us consider a more complex type II interaction scenario, where genes X[,9] and X[,10] among individuals with high resistance-levels (above the median value) are expressed independently from each other. Let us also assume that low or high expression levels are extremely rare among high-resistant individuals and use relatively small variance to simulate them around zero from the following normal distribution.
X[which(resistance < quantile(resistance, 0.5)),10] <-  rnorm(500,0,0.1)
X[which(resistance < quantile(resistance, 0.5)),9] <-  rnorm(500,0,0.1)

#Now suppose  that this "stable", low-variance behaviour of gene X[,9] is assumed to be somehow disrupted among individuals with low resistance-levels (below the median value), and as a more "active" version,  we simulate X[,9] in the low-resistance group with larger variance than in the high-resistant group. 
X[which(resistance > quantile(resistance, 0.5)),9] <- rnorm(500,0,1)

#Let us also assume that, among individuals with low resistance-levels, gene X[,9] is regulating the expression levels of gene X[,10] which can be seen as a regulation pattern that has an an adverse effect on the resistance-level (backwards simulation).
X[which(resistance > quantile(resistance, 0.5)),10] <- X[which(resistance > quantile(resistance, 0.5)),9] + rnorm(500,0,0.2)


#Simulate the final trait associated with that particular disease as the sum of the previous resistance-levels and the main effects of two genes. NOTE that the "resistance" vector again contains the random error term.
y <- X[,3] + X[,7] + resistance
```

\#Apply the proposed model\_diffnet function

``` r
source("./model_diffnet.R")
diff_net <- model_diffnet(X,y,cut_off_a = 1/2,alpha = 1/3,r = 0.1, corr_metric = "PCCN", main_effect = TRUE, CV = TRUE,symmetric = TRUE, signAdj = TRUE, adj_matrices = TRUE)
qgraph(diff_net)
```

![](Toy_example_files/figure-gfm/unnamed-chunk-5-1.png)<!-- -->

``` r
diff_net
```

    ##       [,1] [,2] [,3] [,4] [,5] [,6] [,7] [,8] [,9] [,10]
    ##  [1,]    0    0    0    0    1    0    0    0    0     0
    ##  [2,]    0    0    0    0    0    0    0    1    0     0
    ##  [3,]    0    0    0    0    0    0    0    0    0     0
    ##  [4,]    0    0    0    0    0    0    0    0    0     0
    ##  [5,]    1    0    0    0    0    0    0    0    0     0
    ##  [6,]    0    0    0    0    0    0    0    0    0     0
    ##  [7,]    0    0    0    0    0    0    0    0    0     0
    ##  [8,]    0    1    0    0    0    0    0    0    0     0
    ##  [9,]    0    0    0    0    0    0    0    0    0    -1
    ## [10,]    0    0    0    0    0    0    0    0   -1     0
