\name{model_diffnet}
\alias{model_diffnet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Model guided trait-specific co-expression network estimation
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{model_diffnet(X,y,cut_off_a = 1/2,alpha = alpha,r = r, 

corr_metric = c("CCN","PCCN"), signAdj = TRUE, main_effect = TRUE, 

lambda = lambda, CV = TRUE,symmetric = TRUE)}
 

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{A n x p-dimensional model matrix, where p is the number of explanatory variables and n is the number of samples. 
}
\item{y}{ Response variable. 
}
\item{cut_off_a}{Defines how samples are divided into the high and low groups in the estimation process through empirical quantiles Q(cut_off_a) and Q(1-cut_off_a) of the estimated residual vector from the main effect model. Default is 1/2 (median) which uses the whole dataset by dividing it into two equally sized parts.}

\item{alpha}{The elastic net penalty-mixing parameter ([0,1]) -- it yields the lasso penalty with alpha = 1, and the ridge penalty with alpha = 0. Default is 1/3.}
\item{r}{Dependencies between random variables are deemed to be zero in the estimation process if their magnitudes are below the value of r.}

\item{corr.metric}{Defines which correlation metric is used: "CCN" applies standard Pearson correlation coefficient and "PCCN" is a part-correlation coefficient presented in Kontio et al. 2021. NOTE that the "PCCN" should not be confused with the partial correlation coefficient.}
\item{SignAdj}{Logical argument: The sign-adjusted version is used when SignAdj = TRUE (default).
}

\item{main_effect}{Logical argument. With main_effect = "TRUE", model_diffnet first estimates the main effects and removes their impact from the estimation of interaction effects.}
\item{lambda}{A user-specified penalty parameter lambda value for estimating the main effects with the elastic-net algorithm using \textit{glmnet} R-package (Friedman et al. 2010). Using
CV = "TRUE" overrides this.}
\item{CV}{If CV = "TRUE", model_diffnet uses 5-fold cross-validation to find an optimal lambda for estimating the main effects with the elastic-net algorithm.}
\item{symmetric}{Logical argument: Defines whether (TRUE) or not (FALSE) the output is symmetric. Default is TRUE.}

}



}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
  \item{dPCCN}{If SignAdj = "FALSE" and corr_metric = "PCCN", the output is a differential part-correlation matrix.}
  \item{dCCN}{If SignAdj = "FALSE" and corr_metric = "CCN", the output is a differential correlation matrix.}
  \item{signdPCCN}{If SignAdj = "TRUE" and corr_metric = "PCCN", the output is a sign-adjusted differential part-correlation matrix.}
  \item{signdCCN}{If SignAdj = "TRUE" and corr_metric = "CCN", the output is a sign-adjusted differential correlation matrix.}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{


### SIMULATE A SET OF EXPLANATORY VARIABLES ### 

     X <- matrix(0,1000,10) 
     for(i in 1:10){ 
     X[,i] <- rnorm(1000,0,1)} 

 
### SIMULATE THE INTERACTION PART (OF TYPE I) OF THE MODEL ### 

     interactions_and_rand.error <- X[,1]*X[,5] + X[,2]*X[,8] + rnorm(1000,0,0.5) 

### TYPE II INTERACTIONS ARE SIMULATED  
     c <- interactions_and_rand.error
     X[which(c < median(c)),9] <- rnorm(500,0,0.1) 
     X[which(c < median(c)),10] <- rnorm(500,0,0.1) 
     X[which(c > median(c)),9] <- X[which(c > median(c)),10] + rnorm(500,0,0.2) 

  

### SIMULATE THE MAIN EFFECTS AND THE RESPONSE VARIABLE

   main_effects <- X[,3] + X[,7] 

   y <- main_effects + interactions_and_rand.error 

  
  
### APPLY MODEL_DIFFNET
  
dPCCN <- model_diffnet(X,y, corr_metric = "PCCN", main_effect = TRUE, signAdj = TRUE) 

qgraph(dPCCN) 


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory (show via RShowDoc("KEYWORDS")):
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }
% Use only one keyword per line.
% For non-standard keywords, use \concept instead of \keyword:
% \concept{ ~cpt1 }
% \concept{ ~cpt2 }
% Use only one concept per line.
