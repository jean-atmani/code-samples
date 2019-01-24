install.packages('sn')
install.packages("readxl")
install.packages("R.matlab")

library('sn')
library("readxl")
library("R.matlab")


% import data %
  
data1 <- read_excel("R inputs.xlsx")
data2 <- read_excel("xi2, alpha2.xlsx")
data3 <- read_excel("xi3, alpha3.xlsx")
Omega <- as.matrix(read.csv("Omega.csv", header=FALSE, sep=","))

% simulate multivariate random numbers from a skew-t distribution and conversion into Matlab scripts %

% simulate it 1000 times%
  
  set.seed(5)
  RetSimMskt1 = array(0, dim=c(1237,108,1000))

  for (i in 1:1000){
    
    RndMsktData<-rmst(n=1237, data1$xi, Omega,data1$alpha,7)
    RetSimMskt1[ , ,i] <- RndMsktData
  }

writeMat("RetSimMskt1.mat",RetSimMskt1=RetSimMskt1)

set.seed(23)
RetSimMskt2 = array(0, dim=c(1237,108,1000))

for (i in 1:1000){
  
  RndMsktData<-rmst(n=1237, data2$xi2, Omega,data2$alpha2,7)
  RetSimMskt2[ , ,i] <- RndMsktData
}

writeMat("RetSimMskt2.mat",RetSimMskt2=RetSimMskt2)

set.seed(11)
RetSimMskt3 = array(0, dim=c(1237,108,1000))

for (i in 1:1000){
  
  RndMsktData<-rmst(n=1237, data3$xi3, Omega,data3$alpha3,7)
  RetSimMskt3[ , ,i] <- RndMsktData
}

writeMat("RetSimMskt3.mat",RetSimMskt3=RetSimMskt3)

