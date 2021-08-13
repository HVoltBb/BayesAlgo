## Forward backward algorithm for a hidden markov chain
## version: 0.1
## author: HVoltBb
## email: cannoic at gmail.com
##
## Premises
## 1, state space is finite
## 2, transistion distributions are known
## 3, emission distribution is known
## Objective
## 4, to compute the marginal distribution of the hidden states

T <- 50 # length of the chain
S <- 1:10 # sample space
n <- length(S) # size of the sample space
x0 <- S[1] # initial state
trans <- matrix(0, nrow = n, ncol = n) # transition matrix
trans <- apply(trans, 2, function(x) {
    tmp <- runif(n)
    tmp / sum(tmp)
})

x <- rep(NA, T) # hidden states
y <- rep(NA, T) # observations

x[1] <- x0 # initial hiden state
y[1] <- rpois(1, x[1]) # initial observed state

trans_cdf <- apply(trans, 2, cumsum)

for (i in 2:T) {
    x[i] <- S[which(trans_cdf[, which(S == x[i - 1])] > runif(1))[1]]
    y[i] <- rpois(1, x[i])
}

## Forward-backward algorithm
pf0 <- rep(1, n) / n
pf <- matrix(0, nrow = n, ncol = T)
tmp <- pf0 * dpois(y[1], S)
pf[, 1] <- tmp / sum(tmp)

for (i in 2:T) {
    pf[, i] <- dpois(y[i], S) * t(pf[, i - 1]) %*% trans
}

pb0 <- rep(1, n) / n
pb <- matrix(0, nrow = n, ncol = T)
pb[, T] <- pb0

for (i in (T - 1):1) {
    pb[, i] <- t(pb[, i + 1] * dpois(y[i + 1], S)) %*% trans
}

tmp <- pf * pb
pi <- apply(tmp, 2, function(x) {
    x / sum(x)
})

Ex <- S %*% pi

## Visualization
plot(y, type = "l")
points(1:T, Ex, col = "red")