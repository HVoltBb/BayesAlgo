## Viterbi algorithm
## version: 0.1
## auther: HVoltBb
## email: cannoic at gmail.com
##
## Premises
## 1, finite state space
## 2, state transition matrix is known
## 3, emmision probabilities are known
## Objective
## 4, to compute the maximum a posteriori estiamte of the hidden states

T <- 50 # length of the chain
S <- 1:10 # state space
n <- length(S)
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

## Viterbi algorithm
p0 <- rep(1, n) / n
pi <- matrix(0, nrow = n, ncol = T)
path <- pi

pi[, 1] <- p0 * dpois(y[1], S)

for (i in 2:T) {
    tmp <- t(c(pi[, i - 1]) * t(trans))
    path[, i - 1] <- apply(tmp, 1, function(x) {
        S[which.max(x)]
    })
    pi[, i] <- apply(dpois(y[i], S) * tmp, 1, max)
}

path_map <- rep(S[which.max(pi[, T])], T)

for (i in (T - 1):1) {
    path_map[i] <- path[path_map[i + 1], i]
}

plot(y, type = "l")
lines(path_map, col = "red")