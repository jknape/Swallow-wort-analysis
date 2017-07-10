library(rjags)
load.module("glm")
load.module("mix")

# Load data
dat = read.csv("seedData.csv", sep=",",  nrow=40 , row.names=1, flush=TRUE)

V = dat[["Pods"]]
A = dat[["Area"]]
E = dat[["No.E.attacked.pods"]]


# Function returning design matrix for a thin plate spline from Crainiceanu et al.
getZ = function(y, num.knots) {
  knots = quantile(unique(as.numeric(as.matrix(y))), na.rm=TRUE, seq(0,1,length = num.knots + 2))[-c(1,num.knots)]
  Z_K = (abs(outer(X[,2], knots, "-")))^3
  OMEGA_all = (abs(outer(knots, knots,"-")))^3
  svd.OMEGA_all = svd(OMEGA_all)
  sqrt.OMEGA_all = t(svd.OMEGA_all$v %*% (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d)))
  Z = t(solve(sqrt.OMEGA_all, t(Z_K))) 
  Z
}

# Grid over which to compute spline for resource feedback
grid = seq(-7,7, length.out = 40)
num.knots = 7
X = cbind(1, grid)
Z = getZ(log(V/A), num.knots) 

# Grid over which to compute spline for predation rate
gridP = seq(-12,12, length.out = 40)
XP = cbind(1, gridP)
ZP = getZ(log(V[2:length(V)]/E[1:(length(V)-1)]), num.knots) 

# Standardize rainfall
rain = as.numeric(scale(log(dat[["Rain.JJ"]])))

model = jags.model("swallowWortModel.bugs", 
                   data = list(V = as.vector(V),
                               A = as.vector(A),
                               E = as.vector(E),         
                               Z = Z,
                               X = X,
                               ZP = ZP,
                               XP = XP,
                               rain = rain,
                               grid = grid,
                                gridP = gridP,
                               nG = length(grid),
                               nK = ncol(Z),
                               nYear = length(V)
                   ),  
                   inits = list(a0 = 3, 
                                b = rep(0, ncol(Z)),
                                beta = c(0,-0.5),
                                bP = rep(0, ncol(ZP)),
                                betaP = c(0,-0.5),
                                la = approx(1:length(A),log(A),1:length(A))$y,
                                s.epsilon = 0.50,
                                s.eta = 0.13,
                                s.gamma = 0.04
                   ),
                   n.adapt = 4000,
                  n.chains = 4)

# This will take quite some time to run
t1 = proc.time()
samp = coda.samples(model, variable.names = c("tp","tpP", "ls.pred","ls.eps", "s.epsilon", 
                                              "c1",  "ls.mu", "la", "s.gamma", "s.eta", "beta", 
                                              "b", "a0","s.b", "s.bP", "pe", "th.bb", "d0", "d1"), 
                    n.iter=1000000, thin = 500)
proc.time()-t1

source("mcmc2array.r")

# Convert mcmc chains from coda into arrays for posterior computations
tp = mcmc2array(samp, "tp")
a0 = mcmc2array(samp, "a0")
s.epsilon = mcmc2array(samp, "s.epsilon")
c1 = mcmc2array(samp, "c1")
d0 = mcmc2array(samp, "d0")
d1 = mcmc2array(samp, "d1")
th = mcmc2array(samp, "th.bb")
ls.eps = mcmc2array(samp, "ls.eps")
ls.pred = mcmc2array(samp, "ls.pred")
la = mcmc2array(samp, "la")
tpP = mcmc2array(samp, "tpP")



#################
# Figs
#################

# Fig 2.
par(mfcol = c(3,1))

addCI = function(x, hpd, ...) {
  polygon(x = c(rev(x), x), y = c(rev(hpd[,1]), hpd[,2]), ...)
}

par(mfcol = c(4,1), ps = 10, oma = c(0,0,0,0.5), mai = c(.5,.5,.1,0),mgp=c(2.5,1,0), cex = .8)
plot(1977:2016, A, t = "n", xlab = "Year", ylab = "Plant cover area (m^2)", ylim = c(0,2000), xlim = c(1975,2020), bty = "n", xaxt = "n")
axis(1, at = seq(1975, 2020, by = 5))
hpdA = HPDinterval(as.mcmc(exp(la)))
addCI(1977:2016, hpdA, col = "gray", border = NA)
points(1977:2016, A, pch = 20, cex = 1)
points(1977:2016, exp(apply(la,2,mean)), t = "l")

veg = c(192,192,186,198,200,192,198,193,188,190,194,195,194,191,190,190,192,191,192,192,196,196,196,197,200,205,205,207,210,210,210,210,211,212,213,213)
plot(1976:2011, veg, xlab = "Year", ylab = "Vegetation period (days)", ylim = c(180,215),xlim = c(1975,2020), t = "l", bty = "n", pch = 20, cex = 1, xaxt = "n")
axis(1, at = seq(1975, 2020, by = 5))

plot(1977:2016, V, xlab = "Year", ylab = "Seed pods", xlim = c(1972,2020), t = "b", bty = "n", pch = 20, cex = 1, xaxt = "n")
axis(1, at = seq(1975, 2020, by = 5))

plot(1977:2016, E/V, xlab = "Year", ylab = "Predation rate", xlim = c(1973,2020), t = "b", bty = "n", pch = 20, cex = 1, xaxt = "n")
axis(1, at = seq(1975, 2020, by = 5))

# Fig 3.
par(mfcol = c(1,1), ps = 10, oma = c(0,0,0,0.5), mai = c(.5,.5,.1,0),mgp=c(2.5,.8,0), cex = .8)
Ai = apply(exp(la),2,mean)
plot(grid, apply(tp,2,mean)+ mean(a0), type = "l", xlim = c(-5,6), ylim = c(-5,6), xlab = expression(exp(s[t-1])), ylab = expression(exp(f(s[t]))), bty = "n", axes = FALSE)
hpdTP = HPDinterval(as.mcmc(tp + a0 %*% rep(1, 40)))
addCI(grid, hpdTP, col = "gray", border = NA)
points(grid, apply(tp,2,mean)+mean(a0), type = "l")

points(log(V[1:19]/Ai[1:19]), log(V[2:20]/Ai[2:20]), pch = 1)
points(log(V[20:39]/Ai[20:39]), log(V[21:40]/Ai[21:40]), pch = 20)
axis(1, at = log(10^(-2:2)), labels = 10^(-2:2))
axis(2, at = log(10^(-2:2)), labels = 10^(-2:2))



# Fig 4.
par(mfcol = c(1,1), oma = c(0,0,0,0.5), mai = c(.5,.5,.1,0),mgp=c(2.5,.8,0), cex = .8)

plot(log(V[2:40]/E[1:39]), (E[2:40]/V[2:40]),  bty = "n", type = "n", xlab = expression(S[t]/E[t-1]), ylab = "E. connexa predation rate in year t", axes= FALSE)
hpdE = HPDinterval(as.mcmc(plogis(d0 %*% rep(1, 40) + tpP)))

addCI(gridP, hpdE, col = "gray", border = NA)
points(gridP, plogis(mean(d0) + apply(tpP,2,mean)), type = "l")
points(log(V[2:40]/E[1:39]), (E[2:40]/V[2:40]),  pch = 20)
axis(1, at = log(10^(seq(-4,4,by=2))), labels = paste(formatC(10^(seq(-4,4,by=2)), digits = 4, format = "fg", width = -1)))
axis(2)