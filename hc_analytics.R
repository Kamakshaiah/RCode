# HEALTHCARE ANALYTICS - NONCLINICAL DATA ANALYSIS
# 
# service variability of patient discharge services
# statistical diagnosis

# STEPS:
# 1. DATA SIMULATIN
# 2. TEESTING SIMILARITY BETWEEN ACTUAL DATA AND SIMULATED DATA
# 3. STATITICAL DIAGNOSIS THROGH MARDIA, BARTLETT'S, LEVENES and BROWN-FORSYTHES tests

# DATA SIMULATION

set.seed(1) 
prod.lines <- matrix(abs(round((rnorm(150)*10*0.2)+0.01, 0)), 30, 5) 
colnames(prod.lines) <- c("Line1", "Line2", "Line3", "Line4", "Line5") 
prod.lines

# testing of simulated data sets with actual data 

# ACTUAL DISCHARGES
# you need to give your actual data with in paranthesis i.e. in c()
Actual <- c(27, 28, 14, 26, 12, 24, 27, 22, 22, 27, 13, 22, 25, 26, 19, 22, 22, 15, 24, 27, 21, 26, 23, 25, 19, 22, 33, 17, 21, 27) 
Actual <- as.numeric(Actual) 

# set.seed(2) 
# PLANNED DISCHARGES
Planed <- matrix(NA, dim(prod.lines)[1]) 
for(i in 1:dim(prod.lines)[1]) {Planed[i] <- sum(prod.lines[i,])} 
cor(Actual, as.vector(Planed)) 

# MARDIA TEST OF MULTIVARIATE NORMALITY

mardia.test <- function(x1) {
  x1 <- matrix(abs(round(rnorm(16)*10, 0)), 4, 4)
  col_sums <- colSums(x1) 
  row_sums <- rowSums(x1) 
  sigma_hat <- (1/length(x1))*sum((col_sums - mean(colMeans(x1)))*t(col_sums - mean(colMeans(x1)))) 
  
  # chi-square statistic "A" 
  A <- (1/6*(length(x1)))*sum(t(row_sums - mean(colMeans(x1)))*(1/sigma_hat)*(row_sums - mean(colMeans(x1))))^3 
  B_nume_1 <- sqrt(length(x1)/(8*(dim(x1)[2])*(dim(x1)[2]+2))) 
  B_nume_2 <- (1/length(x1))*sum(t(row_sums - mean(colMeans(x1)))*(1/sigma_hat)*(row_sums - mean(colMeans(x1))))^2 
  B_nume_3 <- ((dim(x1)[2])*(dim(x1)[2]+2)) 
  B <- B_nume_1*(B_nume_2 - B_nume_3) 
  
  # DoF 
  dof <- (1/6)*(dim(x1)[2])*(dim(x1)[2]+1)*(dim(x1)[2]+2) 
  
  # pvalue 
  pvalue1 <- pchisq(abs(A), dof) 
  pvalue2 <- pnorm(abs(B), mean(colMeans(x1)), sd(x1)) 
  
  mardia.test <- list(A = A, pvalue1 = pvalue1, B=B,  pvalue2 = pvalue2)
  #mardia.test <- list(chisquare.measure = A, DoF = dof, p-value = pvalue) 
  } 

  ## testing mardia's test of multivariate normality

my.mar.test <- mardia.test(prod.lines) 
my.mar.test 

# BARTLETT'S TEST OF SPHERECITY

bartletts.test <- function(x) {
  
  # Demo data set 
  #prod.lines <- matrix(abs(round(rnorm(30)*100, 0)), 6, 5) #colnames(prod.lines) <- c("Line1", "Line2", "Line3", "Line4", "Line5") 
  #prod.lines 
  # variance vector 
  
  if(is.data.frame(x)){x=x}else{x = as.data.frame(x)} 
  var.vec <- diag(var(x)) 
  
  # Alternative 
  # var.vec <- matrix(NA, 1, 5) 
  #for(i in 1:5) {var.vec[i] <- var(prod.lines[,i])} 
  #var.vec 
  # pooled variance 
  
  n <- dim(x)[1] 
  num <- (n-1)*(sum(var.vec)) 
  N <- prod(dim(x)) 
  a <- dim(x)[2] 
  pooled.var <- num/(N-a) 
  
  # Numerator 
  log(var.vec) 
  m <- (n-1)*sum(log(var.vec)) 
  sh <- (N-a)*log(pooled.var) 
  Q <- m - sh 
  
  # Denominator 
  
  C <- 1+(1/3*(a-1))*(sum(for(i in 1:5){1/length(x[,1])})-(1/(N-a))) 
  
  # Chi-squared statistic 
  
  chi2 <- 2.3026*(Q/C)
  
  # p-value 
  
  p <- pchisq(abs(chi2), df = dim(x)[2]-1) 
  
  bartletts.test <- list(Q = Q, C = C, test.static = chi2, p.value = p) 
  
} 

  ## testing bartlett's test of spherecity

my.bart.res <- bartletts.test(prod.lines) 
my.bart.res

# LEVENE'S TEST

levene.test <- function(x) {
  
  x <- matrix(round(rnorm(16)*10, 0), 4, 4) 
  
  # check if "x" is a numeric 
  # x <- as.numeric(x) 
  # column means.... 
  y_bar <- colMeans(x) 
  z_i_j <- x-y_bar 
  z_i_j <- abs(z_i_j) 
  Z_i_dot <- (1/dim(x)[1])*colMeans(z_i_j) 
  z_dot_dot <- (1/length(x))*sum(z_i_j) 
  
  # numerator
  
  numer_part <- dim(x)[1]*(Z_i_dot - z_dot_dot) 
  sum_of_num_part <- sum(numer_part) 
  numer <- (length(x)-dim(x)[2])*sum_of_num_part 
  
  # denominator 
  
  denom_part <- sum(z_i_j - Z_i_dot) 
  denom <- (dim(x)[2]-1)*denom_part 
  W <- numer/denom 
  p.value <- pf(abs(W), (dim(x)[2]-1), length(x)-dim(x)[2]) 
  levene.test <- list(means = y_bar, p.value = p.value) 
  
  # return(p.value) 
  } 


  ## testing Levene's Test of variability 

lev.res <- levene.test(prod.lines) 
lev.res 

# BROWN-FORSYTHE'S TEST

bf.test <- function(x) {
  x <- matrix(abs(round(rnorm(30)*10, 0)), 5, 6) 
  
  # check if "x" is a numeric 
  # x <- as.numeric(x) 
  # column means.... 
  
  med_vec <- matrix(NA, dim(x)[2]) 
  for(i in 1:dim(x)[2]) {med_vec[i, ] <- median(x[,i])} 
  z_i_j <- x-as.vector(med_vec) 
  z_i_j <- abs(z_i_j) 
  z_i_j_med <- matrix(NA, dim(x)[2]) 
  
  for(i in 1:dim(x)[2]) {z_i_j_med[i, ] <- median(z_i_j[,i])} 
  Z_i_dot <- (1/dim(x)[1])*z_i_j_med 
  z_dot_dot <- (1/length(x))*sum(z_i_j_med)
  
  # numerator 
  
  numer_part <- dim(x)[1]*(Z_i_dot - z_dot_dot) 
  sum_of_num_part <- sum(numer_part) 
  numer <- (length(x)-dim(x)[2])*sum_of_num_part 
  
  # denominator 
  
  denom_part <- sum(z_i_j - as.vector(Z_i_dot)) 
    denom <- (dim(x)[2]-1)*denom_part 
  
  W <- numer/denom 
  p.value <- pf(W, (dim(x)[2]-1), (length(x)-dim(x)[2])) 
  
  levene.test <- list(medians = med_vec, p.value = p.value) 
  
  # return(p.value) 
  } 

  ## testing Brown-Forsythe's test 

my.brown.test <- bf.test(prod.lines) 
my.brown.test 



