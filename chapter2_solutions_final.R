## Statistical Computing 2026 - Chapter 2 R Exercises (FINAL version)
## HOW TO RUN: open in RStudio, click "Source" (top right of editor) or press Ctrl+Shift+S.
## Results print in the Console; plots appear in the Plots pane (use the arrows to flip through them).
## To run one exercise only: select from its "#### Exercise" header down to just before the next one, press Ctrl+Enter.

rm(list = ls())                        # start clean
graphics.off()                         # reset plot devices
say <- function(...) cat("\n[Explanation]", ..., "\n")
out <- function(label, x) { cat("\n--", label, "--\n"); print(x) }

###############################################################
#### Exercise 2.1
###############################################################
set.seed(2026)
values <- c(0, 1, 2, 4)
p      <- c(0.15, 0.25, 0.35, 0.25)
out("sum of probabilities (must be 1)", sum(p))
cp <- cumsum(p)                                   # cumulative probabilities
out("cumulative probabilities", cp)

# (a) pmf and right-continuous cdf
par(mfrow = c(1, 2))
plot(values, p, type = "h", lwd = 3, main = "pmf", ylab = "P(X = x)")
Fx <- stepfun(values, c(0, cp))                   # right-continuous step function
plot(Fx, main = "cdf", verticals = TRUE)
par(mfrow = c(1, 1))

out("P(X<=2)",     cp[values == 2])
out("P(1<X<=4)",   cp[values == 4] - cp[values == 1])
out("P(X=2)",      p[values == 2])
# The jump of F at 2 is F(2) - F(2-) = cp[3] - cp[2] = p[3] = P(X = 2)
out("jump of cdf at 2", cp[3] - cp[2])
say("F is a right-continuous step function. At x = 2 it jumps from F(2-) = 0.40 to F(2) = 0.75,",
    "so the jump size F(2) - F(2-) = 0.35 equals the point probability P(X = 2).")

# (b) moments
EX   <- sum(values * p)
EX2  <- sum(values^2 * p)
EX3  <- sum(values^3 * p)
VarX <- EX2 - EX^2
out("E[X], E[X^2], E[X^3], Var(X)", c(EX, EX2, EX3, VarX))
yv   <- 3 * values - 2                            # Y = 3X - 2
EY   <- sum(yv * p)
VarY <- sum(yv^2 * p) - EY^2
out("E[Y] vs 3E[X]-2",    c(EY, 3 * EX - 2))
out("Var(Y) vs 9Var(X)",  c(VarY, 9 * VarX))

# (c) generalized inverse (quantile function)
qfun <- function(q) values[which(cp >= q - 1e-12)[1]]
qs   <- c(0.15, 0.40, 0.50, 0.75, 0.90)
xq   <- sapply(qs, qfun)
Fxq  <- sapply(xq, function(z) Fx(z))                                  # F(x_q)
Fxq_left <- sapply(xq, function(z) { i <- which(values == z); if (i == 1) 0 else cp[i - 1] })  # F(x_q-)
out("quantile table", data.frame(q = qs, xq = xq, F_xq = Fxq, F_left = Fxq_left,
                                 F_xq_ge_q = Fxq >= qs, F_left_lt_q = Fxq_left < qs))

# simulation check
sim <- sample(values, size = 1e5, replace = TRUE, prob = p)
out("simulated relative frequencies", table(factor(sim, levels = values)) / length(sim))

###############################################################
#### Exercise 2.2
###############################################################
set.seed(2026)
J <- matrix(c(0.10, 0.15,
              0.20, 0.10,
              0.15, 0.30), nrow = 3, byrow = TRUE)
x <- 0:2; y <- 0:1
out("total of J", sum(J))
pX <- rowSums(J); pY <- colSums(J)
out("marginal pmf of X", pX)
out("marginal pmf of Y", pY)
out("P(X<=1, Y=1)", sum(J[x <= 1, y == 1]))

condX <- sweep(J, 2, pY, "/")                     # divide each column by its column sum
out("conditional pmfs (col 1: X|Y=0, col 2: X|Y=1)", condX)
cm <- colSums(x * condX)                          # E[X | Y = y]
cv <- colSums(x^2 * condX) - cm^2                 # Var(X | Y = y)
out("conditional means", cm)
out("conditional variances", cv)
out("P(X>=1 | Y=1)", sum(condX[x >= 1, 2]))

# (b)
EX  <- sum(x * pX); EY <- sum(y * pY)
EXY <- sum(outer(x, y, "*") * J)
covXY <- EXY - EX * EY
VarX  <- sum(x^2 * pX) - EX^2
VarY  <- sum(y^2 * pY) - EY^2
corXY <- covXY / sqrt(VarX * VarY)
out("E[X], E[Y], E[XY], Cov, Cor", c(EX, EY, EXY, covXY, corXY))
out("product of marginals", outer(pX, pY))
out("J minus product of marginals (nonzero => NOT independent)", J - outer(pX, pY))

out("total expectation: E[X] vs sum E[X|Y]P(Y)", c(EX, sum(cm * pY)))
VarEcond <- sum(cm^2 * pY) - sum(cm * pY)^2       # Var(E[X|Y])
out("total variance: Var(X) vs E[Var(X|Y)] + Var(E[X|Y])", c(VarX, sum(cv * pY) + VarEcond))

out("Var(2X-Y) via covariance formula", 4 * VarX + VarY - 2 * 2 * covXY)
W <- outer(2 * x, y, "-")                         # value of 2x - y in each cell
out("Var(2X-Y) directly", sum(W^2 * J) - sum(W * J)^2)

# (c) f(x,y) = 2 on 0 < y < x < 1:  f_X(x) = 2x,  f(y | x) = 1/x on (0, x)  (uniform on (0,x))
fjoint <- function(y, x0) ifelse(y > 0 & y < x0, 2, 0)
fX <- function(x0) integrate(function(t) fjoint(t, x0), 0, x0)$value   # marginal: integrate joint over 0<y<x
out("marginal f_X(x) via integrate at x = 0.3, 0.8 (theory 2x = 0.6, 1.6)", c(fX(0.3), fX(0.8)))
out("f_X integrates to 1 over (0,1)", integrate(function(t) sapply(t, fX), 0, 1)$value)
out("conditional density f(y|x=0.8) = joint / marginal, at y = 0.5 (theory 1/0.8 = 1.25)", fjoint(0.5, 0.8) / fX(0.8))
fYx <- function(y, x0) ifelse(y > 0 & y < x0, 1 / x0, 0)
out("conditional density integrates to 1", integrate(function(t) fYx(t, 0.8), 0, 0.8)$value)
out("E[Y | X=0.8] (should be 0.4)", integrate(function(t) t * fYx(t, 0.8), 0, 0.8)$value)
out("P(Y<0.3 | X=0.8) (should be 0.375)", integrate(function(t) fYx(t, 0.8), 0, 0.3)$value)

# (d) uncorrelated but dependent
u <- runif(10000, -1, 1); v <- u^2
out("cor(u, v)", cor(u, v))
plot(u, v, pch = 16, cex = 0.4, main = "v = u^2: cor ~ 0 but strongly dependent")
# Correlation measures only LINEAR association; v is a deterministic function of u.
say("cor(u, v) is about 0 because the relation v = u^2 is symmetric (U-shaped), not linear.",
    "But v is a deterministic function of u, so they are strongly dependent. Zero correlation does NOT imply independence.")

###############################################################
#### Exercise 2.3   (n = 5000)
###############################################################
set.seed(2026)
n <- 5000

# (a) Bernoulli
b <- rbinom(n, size = 1, prob = 0.3)
out("Bernoulli mean: sample vs p",       c(mean(b), 0.3))
out("Bernoulli variance: sample vs p(1-p)", c(var(b), 0.3 * 0.7))
# Binomial(20, 0.15)
k <- 0:20
plot(k, dbinom(k, 20, 0.15), type = "h", lwd = 3, main = "Binomial(20, 0.15) pmf", ylab = "P(B = k)")
out("P(2<=B<=5) exact", pbinom(5, 20, 0.15) - pbinom(1, 20, 0.15))
Bs <- rbinom(n, 20, 0.15)
out("P(2<=B<=5) simulated", mean(Bs >= 2 & Bs <= 5))

# (b) Multinomial
pr <- c(0.2, 0.5, 0.3)
out("multinomial P(6,15,9)", dmultinom(c(6, 15, 9), prob = pr))
counts <- rmultinom(n, size = 30, prob = pr)      # each COLUMN is one count vector
out("mean counts: simulated", rowMeans(counts))
out("mean counts: theory 30*p", 30 * pr)
out("covariance: simulated", cov(t(counts)))
out("covariance: theory (diag 30p(1-p), off-diag -30 pi pj)", 30 * (diag(pr) - outer(pr, pr)))

# (c) Geometric, p = 0.25 (G = failures before first success)
pg <- 0.25
out("P(G=3)",  dgeom(3, pg))
out("P(G>=5)", pgeom(4, pg, lower.tail = FALSE))  # P(G>=5) = P(G>4)
g <- rgeom(n, pg)
out("P(G=3): simulated vs exact",  c(mean(g == 3), dgeom(3, pg)))
out("P(G>=5): simulated vs exact", c(mean(g >= 5), pgeom(4, pg, lower.tail = FALSE)))
out("mean: sample vs (1-p)/p",     c(mean(g), (1 - pg) / pg))
out("var: sample vs (1-p)/p^2",    c(var(g), (1 - pg) / pg^2))
# T = G + 1 counts TRIALS up to and including the first success: E[T] = 1/p, Var unchanged.
say("T = G + 1 counts trials (failures plus the final success): E[T] = 1/p = 4, Var(T) = (1-p)/p^2 (unchanged, a shift only).")
out("memoryless: P(G>=8 | G>=3)", pgeom(7, pg, lower.tail = FALSE) / pgeom(2, pg, lower.tail = FALSE))
out("memoryless: P(G>=5)",        pgeom(4, pg, lower.tail = FALSE))

# Negative binomial: K = failures before the 4th success, p = 0.4
r <- 4; pn <- 0.4
out("P(K=6)",  dnbinom(6, size = r, prob = pn))
out("P(K<=8)", pnbinom(8, size = r, prob = pn))
kk <- rnbinom(n, size = r, prob = pn)
out("P(K=6): simulated vs exact",  c(mean(kk == 6), dnbinom(6, r, pn)))
out("P(K<=8): simulated vs exact", c(mean(kk <= 8), pnbinom(8, r, pn)))
out("mean: sample vs r(1-p)/p",    c(mean(kk), r * (1 - pn) / pn))
out("var: sample vs r(1-p)/p^2",   c(var(kk), r * (1 - pn) / pn^2))
# K + 4 is the variable that counts total trials.
say("Total number of trials until the 4th success is K + 4 (K failures + 4 successes); its mean is r/p, variance unchanged.")

# (d) Poisson(4)
out("P(N>=6) exact", ppois(5, 4, lower.tail = FALSE))
N <- rpois(n, 4)
out("P(N>=6): simulated vs exact", c(mean(N >= 6), ppois(5, 4, lower.tail = FALSE)))
out("mean: sample vs 4", c(mean(N), 4))
out("var: sample vs 4",  c(var(N), 4))
S <- rpois(n, 2) + rpois(n, 3)                    # independent sum ~ Poisson(5)
rel <- table(factor(S, levels = 0:15)) / n
bp <- barplot(rel, main = "Poisson(2)+Poisson(3) relative frequencies vs Poisson(5) pmf")
points(bp, dpois(0:15, 5), pch = 19, col = "red")
legend("topright", c("simulated sum", "Poisson(5) pmf"), fill = c("grey", NA), border = c("black", NA),
       pch = c(NA, 19), col = c(NA, "red"))

###############################################################
#### Exercise 2.4
###############################################################
set.seed(2026)

# (a) Normal(65, 10^2)
xn <- rnorm(1000, 65, 10)
hist(xn, probability = TRUE, main = "N(65, 10^2)", xlab = "x")
curve(dnorm(x, 65, 10), add = TRUE, col = "red", lwd = 2)
lines(density(xn), col = "blue", lwd = 2)
legend("topright", c("theoretical density", "kernel density"), col = c("red", "blue"), lwd = 2)
out("P(60<X<80): exact vs sample", c(pnorm(80, 65, 10) - pnorm(60, 65, 10), mean(xn > 60 & xn < 80)))
out("90th percentile: exact vs sample", c(qnorm(0.90, 65, 10), quantile(xn, 0.90)))

# (b) Gamma(shape 3, rate 0.5)
out("P(W>8)", pgamma(8, shape = 3, rate = 0.5, lower.tail = FALSE))
out("95th percentile", qgamma(0.95, shape = 3, rate = 0.5))
w <- rgamma(5000, shape = 3, rate = 0.5)
out("mean: sample vs shape/rate",     c(mean(w), 3 / 0.5))
out("var: sample vs shape/rate^2",    c(var(w), 3 / 0.5^2))
ssum <- rexp(5000, rate = 0.5) + rexp(5000, rate = 0.5) + rexp(5000, rate = 0.5)
par(mfrow = c(1, 2))
hist(ssum, probability = TRUE, main = "Sum of 3 Exp(0.5)", xlab = "value")
curve(dgamma(x, shape = 3, rate = 0.5), add = TRUE, col = "red", lwd = 2)
plot(ecdf(ssum), main = "ecdf of sum vs gamma cdf")
curve(pgamma(x, shape = 3, rate = 0.5), add = TRUE, col = "red", lwd = 2)
par(mfrow = c(1, 1))

# (c) Exponential(0.2): memorylessness
out("P(T>15 | T>5) exact", pexp(15, 0.2, lower.tail = FALSE) / pexp(5, 0.2, lower.tail = FALSE))
tt <- rexp(5000, 0.2)
out("P(T>15 | T>5) simulated", mean(tt[tt > 5] > 15))
out("P(T>10)", pexp(10, 0.2, lower.tail = FALSE))   # same value: memoryless
say("Memorylessness: P(T > s+t | T > s) = P(T > t). Having already waited 5 units does not change the chance of waiting 10 more.")

# (d) Lognormal
mu_l <- 0.5; sd_l <- 0.7
Lg <- exp(rnorm(5000, mu_l, sd_l))
out("mean: sample vs exp(mu+s^2/2)",  c(mean(Lg), exp(mu_l + sd_l^2 / 2)))
out("var: sample vs (e^{s^2}-1)e^{2mu+s^2}", c(var(Lg), (exp(sd_l^2) - 1) * exp(2 * mu_l + sd_l^2)))
out("P(L>3): exact vs sample", c(plnorm(3, meanlog = mu_l, sdlog = sd_l, lower.tail = FALSE), mean(Lg > 3)))
out("median: exact vs sample", c(qlnorm(0.5, mu_l, sd_l), median(Lg)))
# meanlog / sdlog = mean and sd of log(L) (the underlying normal), NOT of L itself.
say("meanlog = mean of log(L) (= mu = 0.5), sdlog = standard deviation of log(L) (= sigma = 0.7); they are NOT the mean and sd of L.")

###############################################################
#### Exercise 2.5
###############################################################
set.seed(2026)

# (a) chi-square(6) from squared standard normals
M  <- matrix(rnorm(5000 * 6), nrow = 5000, ncol = 6)
Cs <- rowSums(M^2)
hist(Cs, probability = TRUE, main = "Sum of 6 squared normals", xlab = "value")
curve(dchisq(x, df = 6), add = TRUE, col = "red", lwd = 2)
out("mean: sample vs 6",  c(mean(Cs), 6))
out("var: sample vs 12",  c(var(Cs), 12))
out("P(C>12): exact vs sample", c(pchisq(12, 6, lower.tail = FALSE), mean(Cs > 12)))
out("95th percentile: exact vs sample", c(qchisq(0.95, 6), quantile(Cs, 0.95)))

# (b) t = Z / sqrt(C/6)
Zn <- rnorm(5000); C2 <- rchisq(5000, df = 6)
T6 <- Zn / sqrt(C2 / 6)
hist(T6, probability = TRUE, breaks = 50, main = "T = Z/sqrt(C/6)", xlab = "t")
curve(dt(x, df = 6), add = TRUE, col = "red", lwd = 2)
out("P(|T|>2): exact vs sample", c(2 * pt(2, 6, lower.tail = FALSE), mean(abs(T6) > 2)))
out("97.5th percentile (df=6)", qt(0.975, 6))
dfs <- c(1, 4, 30)
curve(dnorm(x), -5, 5, lwd = 2, ylab = "density", main = "Normal vs t densities")
for (i in seq_along(dfs)) curve(dt(x, df = dfs[i]), add = TRUE, col = i + 1, lwd = 2)
legend("topright", c("normal", paste("t, df =", dfs)), col = 1:4, lwd = 2)
out("tail P(|T|>2) for df = 1, 4, 30", sapply(dfs, function(d) 2 * pt(2, d, lower.tail = FALSE)))
out("tail P(|Z|>2) normal", 2 * pnorm(2, lower.tail = FALSE))
# Mean exists for df > 1 (df = 4, 30; NOT df = 1). Variance exists for df > 2 (df = 4, 30; NOT df = 1).
say("Mean exists for df = 4 and 30 (needs df > 1), not df = 1. Variance exists for df = 4 and 30 (needs df > 2), not df = 1.",
    "Smaller df gives heavier tails: P(|T|>2) is largest for df = 1 and approaches the normal value as df grows.")

# (c) Beta(2, 5)
curve(dbeta(x, 2, 5), 0, 1, lwd = 2, main = "Beta(2, 5) density", ylab = "density")
Bt <- rbeta(5000, 2, 5)
out("P(B<0.3): exact vs sample", c(pbeta(0.3, 2, 5), mean(Bt < 0.3)))
out("90th percentile: exact vs sample", c(qbeta(0.9, 2, 5), quantile(Bt, 0.9)))
out("mean: sample vs a/(a+b)", c(mean(Bt), 2 / 7))
out("var: sample vs ab/((a+b)^2(a+b+1))", c(var(Bt), 2 * 5 / (7^2 * 8)))

# (d) Uniform(2, 8)
out("P(3<U<6) exact", punif(6, 2, 8) - punif(3, 2, 8))
out("75th percentile exact", qunif(0.75, 2, 8))
Un <- runif(5000, 2, 8)
out("mean: sample vs (a+b)/2",    c(mean(Un), 5))
out("var: sample vs (b-a)^2/12",  c(var(Un), 36 / 12))
U01 <- (Un - 2) / 6                               # transform to Uniform(0,1)
hist(U01, probability = TRUE, main = "Uniform(0,1) (blue) vs Beta(2,5) (red)", xlab = "value")
curve(dbeta(x, 2, 5), add = TRUE, col = "red", lwd = 2)
abline(h = 1, col = "blue", lwd = 2)

###############################################################
#### Exercise 2.6
###############################################################
set.seed(2026)

# (a) bivariate normal built from independent U, V
U <- rnorm(5000); V <- rnorm(5000)
X <- 2 + 3 * U
Y <- 5 + 2 * (0.6 * U + sqrt(1 - 0.6^2) * V)
out("means of X, Y (theory 2, 5)", c(mean(X), mean(Y)))
out("sds of X, Y (theory 3, 2)",   c(sd(X), sd(Y)))
out("cov(X,Y) (theory 3.6)",       cov(X, Y))
out("cor(X,Y) (theory 0.6)",       cor(X, Y))
plot(X, Y, pch = 16, cex = 0.4, main = "Bivariate normal, rho = 0.6")
muX <- 2; sX <- 3; muY <- 5; sY <- 2; rho <- 0.6
mu_c  <- muY + rho * sY / sX * (4 - muX)          # E[Y | X=4]
var_c <- sY^2 * (1 - rho^2)                       # Var(Y | X=4)
out("E[Y|X=4], Var(Y|X=4)", c(mu_c, var_c))
out("P(Y>6 | X=4)", pnorm(6, mean = mu_c, sd = sqrt(var_c), lower.tail = FALSE))
out("lm prediction at X=4 (compare with E[Y|X=4])", predict(lm(Y ~ X), data.frame(X = 4)))

# (b) rho = 0  (repeat of part (a) with Y0 = 5 + 2*(0*U + 1*V))
rho0 <- 0
Y0 <- 5 + 2 * (rho0 * U + sqrt(1 - rho0^2) * V)
out("means of X, Y0 (theory 2, 5)", c(mean(X), mean(Y0)))
out("sds of X, Y0 (theory 3, 2)",   c(sd(X), sd(Y0)))
out("cov(X,Y0) (theory 0)",         cov(X, Y0))
out("cor(X,Y0) (theory 0)",         cor(X, Y0))
plot(X, Y0, pch = 16, cex = 0.4, main = "Bivariate normal, rho = 0")
mu_c0  <- muY + rho0 * sY / sX * (4 - muX)
var_c0 <- sY^2 * (1 - rho0^2)
out("E[Y0|X=4], Var(Y0|X=4) (= unconditional 5, 4)", c(mu_c0, var_c0))
out("P(Y0>6 | X=4)", pnorm(6, mean = mu_c0, sd = sqrt(var_c0), lower.tail = FALSE))
out("lm prediction at X=4 (compare with E[Y0|X=4] = 5)", predict(lm(Y0 ~ X), data.frame(X = 4)))
# For a JOINTLY normal pair the joint density factorizes iff the covariance is 0,
# so uncorrelated => independent. In 2.2(d), (u, u^2) is not jointly normal.
say("For a jointly normal pair the density factorizes into f_X(x) f_Y(y) exactly when rho = 0 (the conditional law of Y|X=x no longer depends on x),",
    "so zero covariance implies independence. In 2.2(d), (u, u^2) is not jointly normal, so zero correlation proves nothing.")

# (c) three-dimensional normal via Cholesky
mu <- c(1, 2, 3)
Sigma <- matrix(c(1,   0.4, 0.2,
                  0.4, 2,   0.3,
                  0.2, 0.3, 1.5), nrow = 3, byrow = TRUE)
out("eigenvalues of Sigma (all > 0 => positive definite)", eigen(Sigma)$values)
Zs <- matrix(rnorm(5000 * 3), nrow = 5000, ncol = 3)
R  <- chol(Sigma)                                 # upper triangular, t(R) %*% R = Sigma
Xs <- sweep(Zs %*% R, 2, mu, "+")
out("column means (theory 1, 2, 3)", colMeans(Xs))
out("sample covariance (theory Sigma)", cov(Xs))
pairs(Xs, pch = 16, cex = 0.3, main = "Simulated 3-D normal")

# (d) linear combination A = X1 + 2 X2 - X3
a <- c(1, 2, -1)
mean_th <- sum(a * mu)
var_th  <- drop(t(a) %*% Sigma %*% a)
out("theoretical mean and variance of A", c(mean_th, var_th))
A <- Xs %*% a
out("simulated mean and variance of A", c(mean(A), var(A)))
out("P(A>4): exact vs simulated", c(pnorm(4, mean_th, sqrt(var_th), lower.tail = FALSE), mean(A > 4)))

dmvn <- function(z, mu, Sigma) {
  d  <- length(mu)
  dz <- z - mu
  as.numeric((2 * pi)^(-d / 2) * det(Sigma)^(-1 / 2) * exp(-0.5 * t(dz) %*% solve(Sigma) %*% dz))
}
out("density at mu (the peak)", dmvn(mu, mu, Sigma))
out("density at (2,3,4)", dmvn(c(2, 3, 4), mu, Sigma))
# The second value is smaller: (2,3,4) is farther from the mean in Mahalanobis distance.
say("The density is maximal at mu. (2,3,4) lies farther from mu in Mahalanobis distance (z-mu)' Sigma^-1 (z-mu), so its density is much smaller.")

cat("\nDone. Flip through the Plots pane with the arrows to see all figures.\n")
