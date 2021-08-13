## Baum-Welch algorithm
## version: 0.1
## author: HVoltBb
## email: cannoic at gmail.com
##
## Premises
## 1, state space is finite
## 2, emission distribution is known
## Objective
## 4, to compute the max a posterior estimates of the
##    transition distributions and initial states

T <- 1000 # length of the chain
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

## Baum-Welch algorithm
## E step using the forward-backward algorithm

iter_max <- 100
rel_conv <- 1e-20

## random initialization
tmp <- runif(n)
q0 <- tmp / sum(tmp)
tmp <- matrix(runif(n * n), n, n)
p0 <- apply(tmp, 2, function(x) {
    x / sum(x)
})
trans_old <- p0
pf0_old <- q0
iter <- 1
check <- rep(NA, iter_max)
trace <- array(0, dim = c(n, n, iter_max))
trace[, , 1] <- trans_old
pf <- matrix(0, n, T)
pb <- matrix(1, nrow = n, ncol = T)

while (iter < iter_max) {
    # E step
    tmp <- pf0_old * dpois(y[1], S)
    pf[, 1] <- tmp / sum(tmp)

    for (i in 2:T) {
        tmp <- dpois(y[i], S) * t(pf[, i - 1]) %*% trans_old
        pf[, i] <- tmp / sum(tmp)
    }


    for (i in (T - 1):1) {
        tmp <- t(pb[, i + 1] * dpois(y[i + 1], S)) %*% trans_old
        pb[, i] <- tmp / sum(tmp)
    }

    # M step
    tmp <- pf[, 1] * pb[, 1]
    pf0_new <- tmp / sum(tmp)

    trans_new <- matrix(0, n, n)
    for (i in 1:(T - 1)) {
        tmp <- trans_old * matrix(pf[, i], n, n) *
            t(matrix(dpois(y[i + 1], S), n, n)) *
            t(matrix(pb[, i], n, n))
        tmp <- tmp / sum(tmp)

        trans_new <- trans_new + tmp
    }

    trans_new <- apply(trans_new, 1, function(x) {
        x / sum(x)
    })

    # Relative convergence check
    check[iter] <- max(abs(trans_new - trans_old))

    pf0_old <- pf0_new
    trans_old <- (trans_new)
    iter <- iter + 1
    trace[, , iter] <- (trans_new)
}

tmp <- pf * pb
pi <- apply(tmp, 2, function(x) {
    x / sum(x)
})
Ey <- S %*% pi

par(mfrow = c(2, 2))
plot(y[1:50], type = "l")
points(1:50, Ey[1:50], col = "red")