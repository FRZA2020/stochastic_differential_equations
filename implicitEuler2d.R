# Date 25 March, Author: Franziska Nitsch
# This script implements the implicit Euler-method and tests it for a two dimensional differential equation dX_t=A *X_t *dt + B *X_t *dW_t, where A and B are matrices and W is a one dimensional Brownian motion. The example for A and B used here are 

# A = [-1.  1 	and B = [0.1  0
#		1   -1]  		 0    0.1]

# The iteration for the implicit Euler method looks as follows
# Y_n+1 = Y_n + A*Y_n+1 *Delta + B W_n+Y_n
# and can be rearranged to
# Y_n+1 = Ainv *Y_n + Ainv * B * W_n *Y_n, with Ainv = inv(Identity matrix - A*delta). 
# Since A is a 2x2 matrix, the inverse can be calculated via a simple formula. 

# the library sde contains the function BM to simulate paths of a Brownian motion
library(sde)

# set number of intervalls
N<-1000
delta<-1/N


# entries for Ainv
v<-c((1+delta)/((1+delta)^2-delta^2),delta/((1+delta)^2-delta^2),delta/((1+delta)^2-delta^2),(1+delta)/((1+delta)^2-delta^2))
Ainv<-matrix(data=v,2,2,byrow=FALSE)
# entries of B matrix
w<-c(0.1,0,0,0.1)
B<-matrix(data=w,2,2,byrow=FALSE)

# initialise Y for the 2 dimensional solution of the differential equation
y<-vector(mode="double",length=2*(N+1))
Y<-matrix(data=y,2,N+1,byrow=FALSE)

# start value for Y is set to 1
Y0 <- c(1,1)
Y[1,1] = Y0[1]
Y[2,1] = Y0[2]

# We solve the differential equation 3 times, using 3 different samples for the Brownian Motion W
 
for(j in 1:3){
# two dimensional Brownian Motion
W_1 <- BM(x=0, t0=0, T=1, N)
W_2 <- BM(x=0, t0=0, T=1, N)

for(i in 2:(N+1)){
 W <- matrix(c(W_1[i]-W_1[i-1],W_2[i]-W_2[i-1],W_1[i]-W_1[i-1],W_2[i]-W_2[i-1]),2,2,byrow=FALSE)
# iteration for implicit Euler
Y[,i] = Ainv%*%Y[,i-1]+Ainv%*%B%*%W%*%Y[,i-1]
}
# store 3 different solutions in matrices P1, P2, P3
 if (j==1) {
# save BM so that we can compare the solution to the exact solution
W_1s <- W_1
W_2s <- W_2 	
P1 = Y
} 
if (j == 2){
P2 = Y
}
if (j == 3){
P3 = Y
}
}





# two dimensional brownian motion for calculation of the exact solution
W_1 <- W_1s
W_2 <- W_2s

# vector A
A<-matrix(c(-1,1,1,-1),2,2,byrow=FALSE)

# initialise matrix for exact solution
expmat<-vector(mode="double",length=2*(N+1))
EXPMA<-matrix(data=expmat,2,N+1,byrow=FALSE)

# produce the exact solution with a vector "time" for plotting it 
time<-seq(from = 0, to = 1, by =1/N)
# m is used for a matrix multiplication
m<-matrix(c(1,1),2,1,byrow = FALSE)

for (i in 1:(N+1)){
W<-matrix(c(W_1[i],W_2[i],W_1[i],W_2[i]),2,2,byrow=FALSE)
temp<-(expm((A-0.5*B%*%B)*time[i]+B%*%W)%*%m)
EXPMA[1,i]<-temp[1]
EXPMA[2,i]<-temp[2]
}

# component11 is the first, component21 the second component of the solution calculated via the implicit euler (3 samples). EXPMA is the exact solution and should be identical to the first solution, because we have reused its sample of the Brownian Motion in the exact solution. The fifth plot is the difference between the exact solution and the first solution and should be close to zero.

png('component11.png')
par(mfrow=c(2,3))
plot(time,P1[1,],  cex=.2)
#par(new=TRUE)
plot(time,P2[1,],  cex=.2)
#par(new=TRUE)
plot(time,P3[1,],  cex=.2)
#par(new=TRUE)
plot(time,EXPMA[1,] , cex=.2)
#par(new=TRUE)
plot(time,EXPMA[1,]-P1[1,], cex=.2,  col = "red")
dev.off()

png('component21.png')
par(mfrow=c(2,3))
plot(time,P1[2,],  cex=.2)
#par(new=TRUE)
plot(time,P2[2,],  cex=.2)
#par(new=TRUE)
plot(time,P3[2,],  cex=.2)
#par(new=TRUE)
plot(time,EXPMA[2,],  cex=.2)
plot(time,EXPMA[2,]-P1[2,], cex=.2, col = "red")
dev.off()