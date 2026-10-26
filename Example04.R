library(cmna)
library(actuar)
library(univariateML)

# Step 1: Preparation of Wind Speed Data

A                    <- read.table('C:/Users/benja/Desktop/Forschung/2026-WKS-OptimizationTwoParameterWeibullInverseWeibull/ReadData-01639.txt', header=TRUE, sep=",")
B                    <- A[,4]
B                    <- B[which(B>=-0.1)]
C                    <- B
C                    <- C[C>=0.05]

# Step 2: Computation of scale and shape

test_wei <- mlweibull(C)

x1    <- seq(0, 35, length = 351)

shape_wei <- 1.51
scale_wei <- 2.93

eq_wei <- function(x){
  (shape_wei/scale_wei)*(x/scale_wei)^(shape_wei-1)*exp(-(x/scale_wei)^(shape_wei))
}

y1     <- eq_wei(x1)
y1[1]  <- 0

# Step 3: Plotting of histogram and Weibull graph

hist(C,prob = TRUE, col = "white",
     breaks=seq(min(C),max(C),l=35),
     main="Example 3: Histogram and Estimated PDF",
     xlab="Wind Speed", ylab="Relative Frequencies", 
     xlim = c(0,20), ylim = c(0,0.45))
     lines(x1,y1,col="green4",lty="solid",lwd="3")