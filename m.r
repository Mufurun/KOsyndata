# Just for someone who does not like Rmd. 
# I do not know if there is anyone like that.

library("rethinking")
library(dplyr)
library(dagitty)


## Causal effects of the DAG 

kdag<-dagitty("dag {
  S-> K <- T 
  T -> S
  T -> SL -> K
}")
adjustmentSets(kdag, exposure="SL", outcome = "K")# This shows the backdoor path and how to stratify


#####################
## SYnthetic Data  ##
#####################


# Function of SL by Y 
myFunc <- function(x) {
case_when(
    -16.000 <= x && x <= -10   ~ (x + 16) * 30 / 6 - 80,
    -10.000 < x  && x <= -7.300 ~ (x + 10) * 60 / 3 - 50,
    -7.300 < x   && x <= -6    ~ 4,
    -6.000 < x   && x <= -3.000 ~ (x + 6) * (-5) / 3 + 4,
    -3.000 < x   && x <= -2.400 ~ (x + 3) * (-5) / 0.5 - 1,
    TRUE                      ~ 1 
  )
}
curve(myFunc(x), col = 2, xlim=c(-16.000, -2.400), ylim=c(-85, 10))

# Initialize Y 
Y<- runif(200, -16, -2.4)

# T and SL caused by Y
T<- NULL
SL<-NULL
for (i in 1:length(Y)) {
    T[i]<-rpois(1,  (Y[i]+16)/13.6 * 25)
    if(T[i] < 1) {T[i]<- 1}
    else if(T[i] > 25) {T[i]<- 25} # Make sure T is between 1 and 25
    SL[i]<- rnorm(1, myFunc(Y[i]), 5) 
}

SL <- (SL - min(SL)) / (max(SL) - min(SL)) # to estimate easily, make it between 0 and 1 

# Check the data
plot(T, SL)
plot(Y, T)
plot(Y, SL)

# initialize dS with random numbers
dS<-rbeta(length(T), 1.2,1) * 200 
# S is affected by Y and dS (omit T as T and S are correlated by Y)
S<-NULL
for(i in 1:length(T))S[i] <- rnorm(1, 0.3*Y[i] -dS[i]/50, 0.1)/25 + 0.4

# Check the data
range(S) # try to make less than 0.5
plot(density(S))
plot(dS, S)
points(dS[T==2], S[T==2], col = 2, lwd =4)
points(dS[T==13], S[T==13], col = 3, lwd =4)
points(dS[T==7], S[T==7], col = 4, lwd =4)
plot(T, S)
plot(Y, S)
plot(T, dS) # no correlation 
plot(SL, S)

# Initialize dK with random numbers
dK<-rbeta(length(T), 1.4,1)*400
# K is affected by SL, Y, S, dK
K<-NULL
for(i in 1:length(T))K[i] <- rnorm(1, SL[i] + 0.1* Y[i] + S[i]- dK[i]/250, 0.2)/9 + 0.4

# Check the data
range(K)# try to make less than 0.5
plot(density(K))
plot(dK, K)
plot(SL, K)
plot(Y, K)
plot(S, K)
plot(dS, K)

# Standardize data without SL, S, K (those are ratio)
d<-list(T=T, dS=standardize(dS), S=S, K=K, dK=standardize(dK))

# check the data
plot(density(d$K))
plot(d$dK, d$K)
plot(SL, d$K)
plot(d$T, d$K)
plot(d$S, d$K)
plot(d$dS, d$S)
plot(Y, d$S)
plot(d$dS, d$K) # no correlation



## Prior Prediction 


# Model 1 
m1<- ulam(
    alist (
        K ~ normal(mu, U), 
        logit(mu) <- a[T] + b_SL[T] * SL[i] + b_S[T] * S + b_dK * dK + b_Y_K[T] * Y[i],

        # We do not know anything about this. 
        # So use partial pooling to regularize the estimate. 
        # mean should start at 0 and has small variation.
        vector[25]: a <- a_bar + siga * za[T],
        a_bar ~ normal(0,0.3),
        siga ~ exponential(1),
        za[T] ~ normal(0,0.1),


        # In this case, we do not include correlation
        # But we have to regularize. 
        vector [25]: b_S~normal(mub_S, sigmab_S),
        vector [25]: b_SL~normal(mub_SL, sigmab_SL),
        vector [25]: b_Y_K~normal(mub_Y_K, sigmab_Y_K),
        vector[25]: mub_S <- b_S_bar + sigb_S * zb_S[T],
        vector[25]: mub_SL <- b_SL_bar + sigb_SL * zb_SL[T],
        vector[25]: mub_Y_K <- b_Y_K_bar + sigb_Y_K * zb_Y_K[T],
        c(b_S_bar,  b_SL_bar, b_Y_K_bar) ~ normal(0,0.3),
        c(sigb_S, sigb_SL, sigb_Y_K)~exponential(6),
        zb_S[T]~ normal(0,0.01),
        zb_SL[T]~ normal(0,0.01),
        zb_Y_K[T] ~ normal(0,0.01),
        

        c(sigmab_S, sigmab_SL, sigmab_Y_K)~exponential(6),
        

        # don't need regularization
        # We are hoping/believing b_dK is positive 
        b_dK~normal(0, 0.5),
        U~exponential(1),


        # Model S 
        # I am hoping/believing b_dS is negative and b_Y_S is positive 

        S~normal(muS, sigmaS),
        logit(muS) <- a2[T] + b_dS * dS + b_Y_S[T] * Y[i],

        # partial pooling for the same reason for a2 and b_Y_S
        vector[25]: a2<-a2_bar + siga2 * za2[T],
        a2_bar~normal(0,0.1), 
        siga2 ~ exponential(6),
        za2[T]~normal(0,0.05),

        vector[25]: b_Y_S<-b_Y_S_bar + sigb_Y_S * zb_Y_S[T],
        b_Y_S_bar~normal(0,0.1), 
        sigb_Y_S ~ exponential(6),
        zb_Y_S[T]~normal(0,0.05),


        b_dS~normal(0,0.1),

        sigmaS~exponential(1),


        # Model SL by Y 
        # The log function and prior was refered to the previous research of Sea Level Change over time. 
        # This function is abstruct so I will have to actually stan code to provide accurate estimates
        # ulam code in r does not allow us to hartd code the conditional function. 
        # vector[200]: SL ~ normal(muSL + aSL2, sigmaSL),
        vector[200]: SL ~ normal(muSL, sigmaSL),
        vector[200]: logit(muSL) <- aSL + Y[i] * bSL,
        aSL~normal(8, 0.01), # very low variation since this part is hard coded. 
        bSL~normal(0.8, 0.01), # very low variation since this part is hard coded. 

        sigmaSL ~ exponential(1), 

         # Model T 
        T~poisson((Y + 16)/ ( -2.4 + 16) * 25),

        # Model T by Y 
        vector[200]:Y ~ uniform( -16 , -2.4 )


    ), dat = d, constraints = list(
        SL    = "lower=0,upper=1" # SL is between 0 and 1 
        #, T    = "lower=1,upper=25" this is unnecessary since T is provided in the dat 
    ),chains=4, cores=4, log_lik=TRUE, 
#    custom_block = stan_functions
)


dashboard(m1)
precis(m1, depth=1)


# Posterior Check
post <- extract.samples(m1)

# Posterior Y changes 
plot(Y, d$K)
points(post$Y[1,], d$K, col = 2)
for(i in 1 : 200) lines(c(post$Y[1,i], Y[i]), c(d$K[i], d$K[i]), col = 2)

plot(SL, d$K)
points(post$SL[1,], d$K, col = 2)
for(i in 1 : 200) lines(c(post$SL[1,i], SL[i]), c(d$K[i], d$K[i]), col = 2)
points(post$SL[2,], d$K, col = 3)
for(i in 1 : 200) lines(c(post$SL[2,i], SL[i]), c(d$K[i], d$K[i]), col = 3)

plot(Y, SL)
points(post$Y[1,], post$SL[1,], col=2)
for(i in 1 : 200) lines( c(Y[i], post$Y[1,i]), c(SL[i], post$SL[1,i]), col = 2)


# Posterior Simulation 
plot(SL, d$K)
SLseq<-seq(0, 1, len = 20)
KbySL<- sapply(SLseq, # Tricky but t and t to calculate estimates by group, though dK does not have group so do not need 
    function (i) inv_logit(post$a + (post$b_SL) * i + t(t(post$b_S) * d$S) + rep(post$b_dK,25) * d$K + t(t(post$b_Y_K) * colMeans(post$Y)))) 
means<- apply( KbySL , 2 , mean )
PIs<- apply( KbySL , 2 , PI )
lines(SLseq, means, col = 2, lwd = 3)
lines(SLseq, PIs[1, ], col = 2, lty = 2, lwd = 3)
lines(SLseq, PIs[2, ], col = 2, lty = 2, lwd = 3)



plot(d$S, d$K)
Sseq<-seq(0, 0.5, len = 20)
dim(post$SL)
KbyS<- sapply(Sseq, 
    function (i) inv_logit(post$a + t(t(post$b_SL) * colMeans(post$SL)) + post$b_S * i + rep(post$b_dK,25) * d$K + t(t(post$b_Y_K) * colMeans(post$Y)))) 
means<- apply( KbyS , 2 , mean )
PIs<- apply( KbyS , 2 , PI )
lines(Sseq, means, col = 2, lwd = 3)
lines(Sseq, PIs[1, ], col = 2, lty = 2, lwd = 3)
lines(Sseq, PIs[2, ], col = 2, lty = 2, lwd = 3)


plot(d$dK, d$K)
dKseq<-seq(-2.5, 2, len = 20)
KbydK<- sapply(dKseq, 
    function (i) inv_logit(post$a + t(t(post$b_SL) * colMeans(post$SL)) + t(t(post$b_S) * d$S) + rep(post$b_dK,25) * i + t(t(post$b_Y_K) * colMeans(post$Y)))) 
PIs<- apply( KbydK , 2 , PI )
means<- apply( KbydK , 2 , mean )
lines(dKseq, means, col = 2, lwd = 3)
lines(dKseq, PIs[1, ], col = 2, lty = 2, lwd = 3)
lines(dKseq, PIs[2, ], col = 2, lty = 2, lwd = 3)

plot(Y, d$K)
dKseq<-seq(-16, -2.4, len = 20)
KbydK<- sapply(dKseq, 
    function (i) inv_logit(post$a + t(t(post$b_SL) * colMeans(post$SL)) + t(t(post$b_S) * d$S) + rep(post$b_dK,25) * d$K + post$b_Y_K * i)) 
PIs<- apply( KbydK , 2 , PI )
means<- apply( KbydK , 2 , mean )
lines(dKseq, means, col = 2, lwd = 3)
lines(dKseq, PIs[1, ], col = 2, lty = 2, lwd = 3)
lines(dKseq, PIs[2, ], col = 2, lty = 2, lwd = 3)



# Model S by dS and Y 

# By the way, there is no post$S exists as it seems absolutely same as post$SL
# would be because it is one of the predictors that i provided data 
post$SL[1,] == post$S[1,] 



dim(post$Y)
Y
plot(d$dS, d$S)
dSseq<-seq(-2, 2, len = 20)
SbydS<- sapply(dSseq, 
    function (i) inv_logit(post$a2 + rep(post$b_dS, 25) * i + post$b_Y_S * colMeans(post$Y))) # Used mean of each simulated data
PIs<- apply( SbydS , 2 , PI )
means<- apply( SbydS , 2 , mean )
lines(dSseq, means, col = 2, lwd = 3)
lines(dSseq, PIs[1, ], col = 2, lty = 2, lwd = 3)
lines(dSseq, PIs[2, ], col = 2, lty = 2, lwd = 3)

plot(Y, d$S)
Yseq<-seq(-16, -2.4, len = 20)
SbyY<- sapply(Yseq, 
    function (i) inv_logit(post$a2 + rep(post$b_dS, 25) * d$dS + post$b_Y_S * i))
means<- apply( SbyY , 2 , mean )
PIs<- apply( SbyY , 2 , PI )
lines(Yseq, means, col = 2, lwd = 3)
lines(Yseq, PIs[1, ], col = 2, lty = 2, lwd = 3)
lines(Yseq, PIs[2, ], col = 2, lty = 2, lwd = 3)




# Model SL by Y 

plot(Y, SL)
points(post$Y[1,], post$SL[1,], col = 2)
Yseq<-seq(-16, -2.4, len = 20)
SLbyY<- sapply(Yseq, 
    function (i) inv_logit(post$aSL + post$bSL * i))
PIs<- apply( SLbyY , 2 , PI )
means<- apply( SLbyY , 2 , mean )
lines(Yseq, means, col = 2, lwd = 3)
lines(Yseq, PIs[1, ], col = 2, lty = 2, lwd = 3)
lines(Yseq, PIs[2, ], col = 2, lty = 2, lwd = 3)

rm(list=ls())
