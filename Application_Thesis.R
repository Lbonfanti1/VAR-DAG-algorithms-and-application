getwd()

##### libraries and datasets #####

library(igraph)
library(network)
library(vars)
library(mvtnorm)
library(BCDAG)
library(lubridate)
library(urca)
library(readr)
library(corrplot)
library(ggplot2)
library(ggcorrplot) 
library(gridExtra)

# Load handmade codes

irate <- read.csv("FEDFUNDS.csv")
xrate <- read.csv("DEXCHUS.csv")
oil <- read.csv("POILBREUSDM.csv")
unemp <- read.csv("UNRATE.csv")
pcet <- read.csv("PCE.csv")
GDP <- read.csv("GDP.csv")
mon <- read.csv("M3.csv")

# Convert variables into time series

par(mfrow = c(3,3), mar = c(2,2,2,2), oma = c(1,1,1,1))

PCE <- ts(pcet$PCE_CHG, start = c(2001,7), end = c(2022,7), frequency = 3)
fedf <- ts(irate$FEDFUNDS_CHG, start = c(2001,7), end = c(2022,7), frequency = 3)
xr <- ts(xrate$DEXCHUS, start = c(2001,7), end = c(2022,7), frequency = 3)
brent <- ts(oil$POILBREUSDM_CHG, start = c(2001,7), end = c(2022,7), frequency = 3)
gdp <- ts(GDP$GDP_PCH, start = c(2001,7), end = c(2022,7), frequency = 3)
un <- ts(unemp$UNRATE_CHG, start = c(2001,7), end = c(2022,7), frequency = 3)
money <- ts(mon$MABMM301USM189S_CHG/10^11, start = c(2001,7), end = c(2022,7), frequency = 3)


## GRAFICI

plot(xr, type = "l", main = "Exchange rate")
plot(gdp, type = "l", main = "GDP")
plot(un, type = "l", main = "Unemployment")
plot(money, type = "l", main = "M3")
plot(fedf, type = "l", main = "interest rate")
plot(brent, type = "l", main = "Oil Global price")
plot(PCE, type = "l", main = "PCE")



#### Exploratory analysis ####


acf(ir)
pacf(ir)
head(xrate)
scale_x_date()

#ACF AND PACF

plot(xchgr, type = "l")

par(mfrow = c(3,3), mar = c(3,3,3,3), oma = c(1,1,1,1))
acf(pce, main = "PCE")
acf(ir, main = "interest rate")
acf(brent, main = "brent")
acf(xchgr, main = "xchgr")
acf(gdp, main = "gdp")
acf(un, main = "un")
acf(money, main = "money")

# Scelgo lag 2

par(mfrow = c(1,2))
acf(brent)
pacf(brent)

acf(un)
pacf(un)

acf(pce)
pacf(pce)

acf(xchg)
pacf(xchg)

# UNIT ROOT TEST

unitir <- ur.df(ir, lags = 10, selectlags = "AIC")
summary(unitir)

unitbr <- ur.df(pce, lags = 10, selectlags = "AIC")
summary(unitbr)

unitun <- ur.df(un, lags = 10, selectlags = "AIC")
summary(unitun)

unitx <- ur.df(xchgr, lags = 10, selectlags = "AIC")
summary(unitx)

# example of non stationary process

rd <- arima.sim(model = list(order = c(0,1,0)), n = 200)
plot(rd)

rdts <- ts(rd)
test <- ur.df(rdts, lags = 10, selectlags = "AIC")
summary(test)

plot(moneym, type = "l")

#differencing 

moneym <- diff(money)
test1 <- ur.df(moneym, lags = 10, selectlags = "AIC")
summary(test1)

#### graphical representation ####


Df <- data.frame(brent = oil$POILBREUSDM_CHG,
                 xchgr = xrate$DEXCHUS_CHG,
                 pce = pcet$PCE_CHG, 
                 un = unemp$UNRATE_CHG, 
                 ir = irate$FEDFUNDS_CHG,
                 gdp =  GDP$GDP_PCH,
                 money = (mon$MABMM301USM189S_CHG)/(10^11))

head(Df)

corrplot(cor(Df),method = "color",tl.col = "black", order = "alphabet")

my_cols <- c("#00AFBB")  
pairs(Df[1:4], pch = 19,  cex = 0.5,
      col = my_cols,
      lower.panel=NULL)
pairs(Df, pch = 19,  cex = 0.5,
      col = my_cols,
      lower.panel=NULL)

pairs(Df[4:7], pch = 19,  cex = 0.5,
      col = my_cols,
      lower.panel=NULL)


### DAG CREATION ####

X <- as.matrix(Df)

n = nrow(X)
q = ncol(X)


# MCMC INPUT

  # Lag of choice : 2

Y = X[-c(1:2),]
Z = cbind(1,X[-c(1,n),], X[-c(n,n-1),])
p = ncol(Z)
dim(Y); dim(Z)


source("DAG_functions_library.R")

S = 25000
burn = 5000

out_mcmc = mcmc_dag_revised(Y = Y, Z = Z, w = 0.5, S = S)

str(out_mcmc)

plot(out_mcmc$num_edges, type = "l", ylab = "Edges",xlab = "Iterations", col = "#00008B", main = "Number of Edges")

post_probs = apply(out_mcmc$Graph_post[,,burn:S], MARGIN = c(1,2), FUN = mean)

post_probs

med_dag = round(post_probs)

med_dag


labelc = c("brent","xchgr","pce","un","ir","gdp","money")
G = network(med_dag, label = labelc)
G$val
vertex_col = "#00008B"
edge_col = "#CDC8B1"

##################
## Plot 1 (DAG) ##
##################

par(mfrow = c(1,1))

pdf("dag_est.pdf", height = 6, width = 6)

plot.network(G, displaylabels = TRUE, vertex.col = vertex_col,
                   label = labelc,
                   mode = "circle",
                   label.pos = 3,
                   vertex.cex = 3.5,
                   label.cex = 0.9,
                   edge.lwd = 0.1, arrowhead.cex = 1.5)




plot.network(G, displaylabels = TRUE, vertex.col = vertex_col,
             label = labelc,
             mode = "circle",
             label.pos = 1,
             vertex.cex = 3.5,
             label.cex = 0.8,
             edge.lwd = 0.1, arrowhead.cex = 1.5)
plot.network(G, displaylabels = TRUE)
dev.off()


####################
## Plot 2 (probs) ##
####################

library(fields)

par(mar = c(4,4,1,1), oma = c(1,1,1,1), mfrow = c(1,1))

seq_pos = seq(0, 1, length = 7)
post_probs
colori = colorRampPalette(c('white','black'))

pdf("probsright.pdf", height = 6, width = 6)

image.plot(t(post_probs), col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "v", ylab = "u", axes = F,
           horizontal = F, legend.shrink = 1, axis.args = list(cex.axis = 1), cex.lab = 1)
axis(1, at = seq_pos, lab = labelc, las = 1)
axis(2, at = seq_pos, lab = labelc, las = 1)

dev.off()

colnames(post_probs) = rownames(post_probs) = labelc

post_probs

image(post_probs)

## Diagnostics of convergence

plot(out_mcmc$num_edges[burn:S], type = "l")


out_mcmc4 = mcmc_dag_revised(Y = Y, Z = Z, w = 0.5, S = S)


out_mcmc4$Graph_post

str(out_mcmc)

post_probs4 = apply(out_mcmc4$Graph_post[,,burn:S], MARGIN = c(1,2), FUN = mean)



#### plots heatmap #####

par(mar = c(4,4,1,1), oma = c(1,1,1,1), mfrow = c(1,1))

seq_pos = seq(0, 1, length = 7)
post_probs
colori = colorRampPalette(c('white','#8B0000'))

pdf("probs11.pdf", height = 6, width = 6)

image.plot(t(post_probs), col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "v", ylab = "u", axes = F,
           horizontal = F, legend.shrink = 1, axis.args = list(cex.axis = 1), cex.lab = 1)
axis(1, at = seq_pos, lab = labelc, las = 1)
axis(2, at = seq_pos, lab = labelc, las = 1)

dev.off()

pdf("probs22.pdf", height = 6, width = 6)

image.plot(t(post_probs4), col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "v", ylab = "u", axes = F,
           horizontal = F, legend.shrink = 1, axis.args = list(cex.axis = 1), cex.lab = 1)
axis(1, at = seq_pos, lab = labelc, las = 1)
axis(2, at = seq_pos, lab = labelc, las = 1)



dev.off()



post_probs11 = apply(out_mcmc$Graph_post[,,5000:25000], MARGIN = c(1,2), FUN = mean)
post_probs22 = apply(out_mcmc4$Graph_post[,,5000:25000], MARGIN = c(1,2), FUN = mean)

par(mfrow = c(1,2))

image(post_probs1, main = "out_mcmc1")

image(post_probs2, main = "out_mcmc2")


par(mfrow = c(1,1))
pdf("probsconfronto.pdf", height = 6, width = 6)
plot(post_probs11,post_probs22,
     pch = 25,
     bg = "#00008B",
     col = "#00008B",
     ylab = "posterior probabilties m2",
     xlab = "posterior probabilities m1",
     main = "posterior probabilities comparison") + abline(a = 0, b = 1, col = "red", lty = "dashed")
dev.off()



# S = 25000
# burn = 5000

# MCMC convergence of diagnostics


S = 25000
burn = 5000
out_mcmc = mcmc_dag_revised(Y = Y, Z = Z, w = 0.5, S = S)

str(out_mcmc)
plot(out_mcmc$num_edges, type = "l", ylab = "Edges",xlab = "Iterations", col = "#00008B", main = "Number of Edges")
post_probs = apply(out_mcmc$Graph_post[,,burn:S], MARGIN = c(1,2), FUN = mean)

par(mfrow = c(1,1))

pdf("probs25.pdf", height = 6, width = 6)

image.plot(t(post_probs), col = colori(100), zlim = c(0,1), cex.sub = 1, xlab = "v", ylab = "u", axes = F,
           horizontal = F, legend.shrink = 1, axis.args = list(cex.axis = 1), cex.lab = 1)
axis(1, at = seq_pos, lab = labelc, las = 1)
axis(2, at = seq_pos, lab = labelc, las = 1)

dev.off()


med_dag25 = round(post_probs)
G25 = network(med_dag25, label = labelc)
vertex_col = "#00008B"
edge_col = "#CDC8B1"

##################
## Plot 2 (DAG) ##
##################

par(mfrow = c(1,1))

pdf("dag_now.pdf", height = 6, width = 6)

plot.network(G25, displaylabels = TRUE, vertex.col = vertex_col,
             label = labelc,
             mode = "circle",
             label.pos = 6,
             usecurve = T, edge.curve = 0, vertex.cex = 1, edge.col = "gray40",
             label.cex = 0.75, edge.lwd = 0.01, arrowhead.cex = 0.8)
dev.off()

post_probs

med_dag = round(post_probs)

med_dag

