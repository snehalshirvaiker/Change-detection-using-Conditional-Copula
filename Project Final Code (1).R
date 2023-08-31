#Loading Required Package

library(imager)
library(stats)
library(methods)
library(prevtoinc)   
library(copula)
library(VineCopula)
library(copula)
library(VineCopula)
library(VC2copula)
library(EBImage)
library(seewave)

#Loading Images

t1 = load.image("C:/Users/sneha/OneDrive/Desktop/data/p2/t1.jpg")
t2 = load.image("C:/Users/sneha/OneDrive/Desktop/data/p2/t2.jpg")

#Getting dimension of image

dimt1<-dim(t1) 

dim(t1)
dim(t2)
  
#Converting Image into Gray Scale

t1 = grayscale(t1)
t2 = grayscale(t2)
par(mfrow=c(1,2))
plot(t1)
plot(t2)

#difference
diff=t2-t1
plot(diff)
#resizing the image
t2=resize(t2,dimt1[1],dimt1[2])


#Getting Matrix from imported image

oim1 = as.matrix(t1)
oim2 = as.matrix(t2)

#Converting Image into treshold image and getting its matrix

th1 = threshold(t1,thr="80%" , approx = T , adjust = 1) %>% plot
th2 = threshold(t2,thr="80%" , approx = T , adjust = 1) %>% plot
thm1 = as.matrix(th1)
thm2 = as.matrix(th2)

#Getting estimation of Copula Parameter

#Getting Frank Copula Parameters
cpf = BiCopEst(
    oim1,
    oim2,
    family=5,
    method = "mle",
    se = FALSE
)
summary(cpf)

#Getting Clayton Copula Parameters
cpc = BiCopEst(
  oim1,
  oim2,
  family=3,
  method = "mle",
  se = FALSE
)
summary(cpc)

#Getting Gumbel Copula Parameters
cpg = BiCopEst(
  oim1,
  oim2,
  family=4,
  method = "mle",
  se = FALSE
)

summary(cpg)

#Getting Student t Copula Parameters
cps = BiCopEst(
  oim1,
  oim2,
  family=2,
  method = "mle",
  se = FALSE
)
summary(cps)

#Getting Gaussian Copula Parameters
cps = BiCopEst(
  oim1,
  oim2,
  family=1,
  method = "mle",
  se = FALSE
)

summary(cps)

#Getting ind Copula Parameters
cps = BiCopEst(
  oim1,
  oim2,
  family=0,
  method = "mle",
  se = FALSE
)

summary(cps)

#Transforming data into unifrom distribution using 
#Probability Integral transformation

u1 = pobs(oim1)
u2 = pobs(oim2)

cpc = BiCopEst(
  u1,
  u2,
  family=3,
  method = "mle",
  se = FALSE
)
summary(cpc)
ye = BiCopHfunc1(u2, u1, cpc)
yem1 = matrix(ye , nrow = dimt1[1] , ncol=dimt1[2])
yyem=as.Image(yem)
image(yem)

cpc = BiCopEst(
  t1,
  t2,
  family=3,
  method = "mle",
  se = FALSE
)
summary(cpc)
ye = BiCopHfunc1(t2, t1, cpc)
yem1 = matrix(ye , nrow = dimt1[1] , ncol=dimt1[2])
yyem=as.Image(yem)
image(yem)
#Fitting Copula Function 
#Since clayton Copula is best fit for this data we uses frank parameter
#Estimation and fitting clayton copula using those parameters

claytonc = BiCop(family=3, par=cpc$par ,tau=cpc$tau, check.pars = TRUE)
ye = BiCopHfunc1(u2, u1, cpc)
yem = matrix(ye , nrow = dimt1[1] , ncol=dimt1[2])
yyem=as.Image(yem)
image(yem)
th = threshold(yem,thr="22%" , approx = T , adjust = 1) %>% plot
summary(claytonc)



try = BiCop(family=3, par=cpf$par ,tau=cpf$tau, check.pars = TRUE)
ye = BiCopHfunc1(u2, u1, try)
yem = matrix(ye , nrow = dimt1[1] , ncol=dimt1[2])
image(yem)

#Getting Difference of two images

diff = t2-t1

#Plotting diff image

difi=as.Image(diff)
image(difi)

#Plotting Copula Estimated image

yem1=as.Image(yem)
image(yem1)
