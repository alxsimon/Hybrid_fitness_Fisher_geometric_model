# Analysis of simulations
# to run on computing station
library(caTools)

set_opt <- c(1,11,12,13,14)
res <- list()
for (i in 1:length(set_opt)) {
    res[[i]] <- readRDS(paste0("rand_samples_", set_opt[i], ".rds"))
    # produced by script simulations_save_randomizations.R
}

params <- read.table("source_code/Fraisse2016/parameters", row.names = NULL)[,-1]
params$d <- NA
disp_params <- c(1,2,4,5,6,9)
text_params <- c("N",
                 "n",
                 "Moving optimum P1",
                 "Moving optimum P2",
                 expression(alpha),
                 "MVN model",
                 "d")

p <- which(params[set_opt[1],] != params[set_opt[2],])
dat <- list()
plot_max <- 0
plot_max_V <- 0
for (i in 1:length(set_opt)) {
    xi <- set_opt[i]
    N <- params$N[xi]
    n <- params$n[xi]
    s <- params$s[xi]
    diff1 <- params$diff1[xi]
    diff2 <- params$diff2[xi]
    k <- params$k[xi]
    nmuts <- params$nmuts[xi]
    ncycles <- params$ncycles[xi]
    MVN <- params$MVN[xi]
    sigma <- params$sigma[xi]
    
    file <- paste0("simulation_results/SelCoef_N", N,
                   "_n", n,
                   "_sigma", sigma,
                   "_s", s,
                   "_k", k,
                   "_Ncycles", ncycles,
                   "_Nmuts", nmuts,
                   "_opt1_", diff1,
                   "_opt2_", diff2,
                   "_MVN", MVN,
                   "_cycle", ncycles,
                   "_1.txt")
    sim <- read.table(file, header = T)

    coltraits <- 10:(10+n-1)
    sim_P0 <- sim[sim$pop == 0,]
    sim_P1 <- sim[sim$pop == 1,]
    sim_P2 <- sim[sim$pop == 2,]
    
    if (n > 1){
        common_anc <- apply(sim_P0[,coltraits], 2, sumexact) 
        zp1 <- apply(rbind(sim_P0[,coltraits], sim_P1[,coltraits]), 2, sumexact)
        zp2 <- apply(rbind(sim_P0[,coltraits], sim_P2[,coltraits]), 2, sumexact)
        m <- rbind(-sim_P1[nrow(sim_P1):1,coltraits],
                   sim_P2[,coltraits])
        walk_opt <- c(rev(sim_P1$o), 0, sim_P2$o)
        d <- nrow(m)
        xx <- (0:d)/d

        sim_walk <- t(apply(rbind(rep(0,n), apply(m, 2, cumsumexact)), 1, FUN = function(x, zp1){x + zp1}, zp1 = zp1))
        walk_opt <- c(rev(sim_P1$o), 0, sim_P2$o)
        sim_walk[,1] <- sim_walk[,1] - walk_opt
        z_sim <- apply(sim_walk, 1, FUN = function(x) sum(x^2))
        #var_step <- var(sqrt(rowSums(m^2)))
        var_m <- mean(apply(m, 2, var))
        v_sum <- sum(apply(m, 2, var))
        
    } else {
        common_anc <- sum(sim_P0[,coltraits])
        zp1 <- sum(c(sim_P0[,coltraits], sim_P1[,coltraits]))
        zp2 <- sum(c(sim_P0[,coltraits], sim_P2[,coltraits]))
        m <- c(-sim_P1[nrow(sim_P1):1,coltraits],
               sim_P2[,coltraits])
        d <- length(m)
        sim_walk <- zp1 + c(0,cumsum(m)) - c(rev(sim_P1$o), 0, sim_P2$o)
        z_sim <- sim_walk^2
        var_m <- var(m)
        v_sum <- var_m
        #var_step <- var(m)
    }
    
    z_bridge <- res[[i]]
    z_bridge_mean <- rowMeans(z_bridge)
    z_bridge_var <- apply(z_bridge, 1, var)
    
    if (max(z_bridge_mean)/d/v_sum > plot_max) {
        plot_max <- max(z_bridge_mean)/d/v_sum
    }
    
    if (max(z_bridge_var)/d^2/v_sum^2 * n > plot_max_V) {
        plot_max_V <- max(z_bridge_var)/d^2/v_sum^2 * n
    }

    dat[[i]] <- list(xi, m = m,
                     zs = z_sim,
                     var_m = var_m,
                     v_sum = v_sum,
                     zb = z_bridge_mean,
                     zbv = z_bridge_var,
                     d = d)
    params$d[xi] <- d
}
#===========================
# plot
fig_file <- "figures/simulations_opt_10k_sapply.pdf"

pdf(file = fig_file, width = 9, height = 5)
layout(mat=matrix(c(1,1,1,2,3,3,4,4),2,4,byrow=T),heights=c(1,0.6))
mmar <- c(4,5,1,0.5)
par(mar=mmar,
    cex.axis = 1.5,
    cex.lab = 1.5,
    lwd = 1.5)

mylims <- c(0, plot_max)*1.1
cols <- rainbow(length(set_opt))
t <- seq(0, 1, length.out = 1000)
f0 <- t*(1-t)

par(mar = c(4,5,3,0.5))
plot(c(0,1), mylims, col='white',
     xlab=expression(paste('hybrid index, ', italic(h))), 
     ylab=expression(paste('relative breakdown, ', italic(f[h]))))
for (j in 1:length(set_opt)){
    xx <- (0:dat[[j]]$d)/dat[[j]]$d
    lines(xx, dat[[j]]$zs/dat[[j]]$d/dat[[j]]$v_sum, col = cols[j])
    lines(xx, dat[[j]]$zb/dat[[j]]$d/dat[[j]]$v_sum, col = cols[j])
}
lines(t,f0,col='black',lty='dotted')
#legend(-0.03, mylims[2], legend = params[set_opt, p],
#       col = cols, lty = rep(1, length(set_opt)), title = expression(alpha))

par(mar = c(0,0,3,0.5))
plot.new()
legend(0, 0.8,
       legend = paste0(params[set_opt, 4], " / ", params[set_opt, 5], " (d = ", params[set_opt,12], ")"),
       col = cols, lty = rep(1, length(set_opt)), title = "opt 1 / opt 2",
       cex = 1.5)
#legend(0.5, 0.8, legend = params[set_opt, 12],
#       col = cols, lty = rep(1, length(set_opt)), title = "d",
#       cex = 1.5)
#text2add <- matrix(c(text_params, c(as.character(params[i, disp_params])), d), 7, 2)
#sapply(1:2, function(i) text((i/2)-0.2, (1:7)/7, rev(text2add[c(1,2,5,7,3,4,6),i])))

par(mar=mmar)
for (j in 1:length(set_opt)){
    if(j == 1) append = F else append = T
    if (n > 1){
        hist(c(as.matrix(dat[[j]]$m))/sqrt(dat[[j]]$var_m),
             border=F, main='',
             breaks = 40,
             xlab=expression(paste('mutational effect, ',italic(m[i]/root(v)))),
             ylab='Pr. density', freq=F, xlim=4*c(-1,1), ylim=c(0,0.6),
             col = sub("FF$","40",cols[j]), add = append)
    } else {
        hist(dat[[j]]$m/sqrt(dat[[j]]$var_m),
             border=F, main='',
             breaks = 40,
             xlab=expression(paste('mutational effect, ',italic(m[i]/root(v)))),
             ylab='Pr. density', freq=F, xlim=4*c(-1,1),
             ylim=c(0,0.6),
             col = sub("FF$","40",cols[j]), add = append)
    }
}
lines(seq(-4,4,0.01),dnorm(seq(-4,4,0.01),mean=0,sd=1),col='black',lty='dotted')
box()

plot(c(0,1), c(0,plot_max_V)*1.1, col='white',
     xlab=expression(italic(h)), ylab=expression(paste('Var(',italic(f[h]),')')))
for (j in 1:length(set_opt)){
    xx <- (0:dat[[j]]$d)/dat[[j]]$d
    V <- dat[[j]]$zbv/dat[[j]]$d^2/dat[[j]]$v_sum^2 * params$n[set_opt[j]]
    lines(xx, V, col = cols[j])
}
lines(t,2*f0*f0, col='black', lty='dotted')

mtext("Environmental change", side = 3, line = -2, outer = T, cex = 1.5)
mtext("f)", side = 3, line = -2.2, adj = 0.05, outer = T, cex = 2)

dev.off()

# > devtools::session_info()
# Session info ------------------------------------------------------------------
# setting  value
# version  R version 3.4.1 (2017-06-30)
# system   x86_64, linux-gnu
# ui       X11
# language (EN)
# collate  en_US.UTF-8
# tz       Europe/Paris
# date     2018-06-28
# 
# Packages ----------------------------------------------------------------------
# package   * version date       source
# base      * 3.4.1   2017-11-15 local
# bitops      1.0-6   2013-08-17 CRAN (R 3.4.1)
# caTools   * 1.17.1  2014-09-10 CRAN (R 3.4.1)
# compiler    3.4.1   2017-11-15 local
# datasets  * 3.4.1   2017-11-15 local
# devtools    1.13.6  2018-06-27 CRAN (R 3.4.1)
# digest      0.6.15  2018-01-28 CRAN (R 3.4.3)
# graphics  * 3.4.1   2017-11-15 local
# grDevices * 3.4.1   2017-11-15 local
# memoise     1.1.0   2017-04-21 CRAN (R 3.4.3)
# methods   * 3.4.1   2017-11-15 local
# stats     * 3.4.1   2017-11-15 local
# tools       3.4.1   2017-11-15 local
# utils     * 3.4.1   2017-11-15 local
# withr       2.1.1   2017-12-19 CRAN (R 3.4.3)