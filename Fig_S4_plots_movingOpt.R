# Script to produce Fig. S4

setwd("./simulation_results_jumping_opt/")
library(caTools)
library(shape)


#==================
# Functions

draw_z2 <- function(i, m, zp1){
    if (!is.null(ncol(m))){
        p_rand <- apply(m[sample(1:nrow(m)),], 2, caTools::cumsumexact)
        z_rand <- apply(sweep(rbind(rep(0,n), p_rand), 2, zp1, "+"),
                        1, FUN = function(x) sum(x^2))
        return(z_rand)
    } else {
        m_rand <- sample(m)
        z_rand <- (zp1 + c(0, caTools::cumsumexact(m_rand)))^2
        return(z_rand)
    }
}

draw_z2_opt <- function(i, m, zp1, opt){
    if (!is.null(ncol(m))){
        p_rand <- apply(m[sample(1:nrow(m)),], 2, caTools::cumsumexact)
        p_rand <- sweep(rbind(rep(0,n), p_rand), 2, opt, "-")
        z_rand <- apply(sweep(p_rand, 2, zp1, "+"),
                        1, FUN = function(x) sum(x^2))
        return(z_rand)
    } else {
        m_rand <- sample(m)
        z_rand <- (zp1 + c(0, caTools::cumsumexact(m_rand)))^2
        return(z_rand)
    }
}

plot_onetrait <- function(m, zp1, opt1){
    sim_walk <- c(0, cumsumexact(m[,2])) + zp1[2]
    z_sim <- sim_walk^2
    v_sum <- var(m[,2])
    
    rep <- 10000
    z_bridge <- sapply(1:rep, draw_z2, m = m[,2], zp1 = zp1[2])
    z_bridge_mean <- rowMeans(z_bridge)
    plot_max <- max(c(max(z_bridge_mean)/d/v_sum, max(z_sim)/d/v_sum))
    
    #============
    # Brownian bridge approx
    t <- seq(0, 1, length.out = 1000)
    
    fp1 <- sum((zp1[2] - opt1[2])^2)/d/v_sum
    fp2 <- sum((zp2[2] - opt1[2])^2)/d/v_sum
    fmal <- (1 - t)^2*fp1 + t^2*fp2
    f <- t*(1 - t) + fmal
    #============
    
    par(mar = c(4,7,1,0.5))
    mylims <- c(0, plot_max)
    plot(c(0,1), mylims, col='white',
         xlab=expression(paste('hybrid index, ', italic(h))), 
         ylab=expression(atop('relative breakdown, '*italic(f[h])*',', 'on the intrinsic trait')),
         cex.axis = 1.5,
         cex.lab = 1.5)
    lines(xx, z_sim/d/v_sum)
    lines(xx, z_bridge_mean/d/v_sum, col = "blue")
    lines(t,f,col='red',lty='dotted', lwd = 2)
}

plot_twotraits <- function(m, zp1, opt1){
    sim_walk <- t(apply(rbind(rep(0,n), apply(m, 2, cumsumexact)), 1, FUN = function(x, zp1){x + zp1}, zp1 = zp1))
    sim_walk[,1] <- sim_walk[,1] - opt1[1]
    z_sim <- apply(sim_walk, 1, FUN = function(x) sum(x^2))
    v_sum <- sum(apply(m, 2, var))
    
    rep <- 10000
    z_bridge <- sapply(1:rep, draw_z2_opt, m = m, zp1 = zp1, opt = opt1)
    z_bridge_mean <- rowMeans(z_bridge)
    plot_max <- max(c(max(z_bridge_mean)/d/v_sum, max(z_sim)/d/v_sum))
    
    #============
    # Brownian bridge approx
    t <- seq(0, 1, length.out = 1000)
    
    fp1 <- sum((zp1 - opt1)^2)/d/v_sum
    fp2 <- sum((zp2 - opt1)^2)/d/v_sum
    fmal <- (1 - t)^2*fp1 + t^2*fp2
    f <- t*(1 - t) + fmal
    #============
    
    par(mar = c(4,7,1,0.5))
    mylims <- c(0, plot_max)
    plot(c(0,1), mylims, col='white',
         xlab=expression(paste('hybrid index, ', italic(h))), 
         ylab=expression(atop('relative breakdown, '*italic(f[h])*',', 'on both traits')),
         cex.axis = 1.5,
         cex.lab = 1.5)
    lines(xx, z_sim/d/v_sum)
    lines(xx, z_bridge_mean/d/v_sum, col = "blue")
    lines(t,f,col='red',lty='dotted', lwd = 2)
}

get_arrows <- function(sim, start, coltraits){
    tmp <- sweep(rbind(rep(0,length(coltraits)), apply(sim[,coltraits], 2, caTools::cumsumexact)), 2, start, "+")
    A <- cbind(tmp[-nrow(tmp),], tmp[-1,])
    return(A)
}

#=========================
#circler <- c(
#    0.5,
#    0.3,
#    0.18,
#    0.10)
circler <- c(
    0.7,
    0.35,
    0.22,
    0.15)
circlecol <- c(
    "#40404030",
    "#40404030",
    "#40404040",
    "#40404040"
)


sims <- list()
sims[[1]] <- list(read.table("SelCoef_N5000_n2_sigma0.5_s0.01_k2_Ncycles5_Nmuts100000_opt1_0_opt2_0_MVN0_cycle5_1.txt",
                             header = T),
                  opt1 = c(0,0),
                  opt2 = c(0,0),
                  xlim = c(-0.5,0.5))
sims[[2]] <- list(read.table("SelCoef_N5000_n2_sigma0.5_s0.01_k2_Ncycles5_Nmuts100000_opt1_0.5_opt2_0.5_MVN0_cycle5_1.txt",
                             header = T),
                  opt1 = c(0.5,0),
                  opt2 = c(0.5,0),
                  xlim = c(-0.1,1))
sims[[3]] <- list(read.table("SelCoef_N5000_n2_sigma0.5_s0.01_k2_Ncycles5_Nmuts100000_opt1_0.5_opt2_-0.5_MVN0_cycle5_1.txt",
                             header = T),
                  opt1 = c(0.5,0),
                  opt2 = c(-0.5,0),
                  xlim = c(-1,1))


titles <- c(
    "a)    Drift around a fixed optimum",
    "b)    Adaptation to the same shifted optimum",
    "c)    Adaptation to different shifted optima"
)


#=========================
pdf("Fig_jumping_opt.pdf", width = 15, height = 8)
layout(mat = matrix(1:9, ncol = 3, byrow=T),
       widths = c(0.7,1,1))
for (i in 1:3){
    sim <- sims[[i]][[1]]
    n <- 2
    coltraits <- 10:(10+n-1)
    trait <- 11
    sim_P0 <- sim[sim$pop == 0,]
    sim_P1 <- sim[sim$pop == 1,]
    sim_P2 <- sim[sim$pop == 2,]
    zmrca <- apply(sim_P0[,coltraits], 2, sumexact)
    zp1 <- apply(rbind(zmrca, sim_P1[,coltraits]), 2, sumexact)
    zp2 <- apply(rbind(zmrca, sim_P2[,coltraits]), 2, sumexact)
    m <- rbind(-sim_P1[nrow(sim_P1):1,coltraits],
               sim_P2[,coltraits])
    d <- nrow(m)
    xx <- (0:d)/d
    
    #============
    par(mar = c(0.5,0.5,2,0.5))
    arrows_p1 <- get_arrows(sim_P1, zmrca, coltraits)
    arrows_p2 <- get_arrows(sim_P2, zmrca, coltraits)
    
    emptyplot(xlim = sims[[i]]$xlim, ylim = c(-0.5,0.65))
    mapply(plotcircle, r = circler, col = circlecol, lcol =list(NULL), mid = list(sims[[i]][[2]]))
    if (i == 3) {
        mapply(plotcircle, r = circler, col = c("#cc00ff30", "#cc00ff30", "#cc00ff40", "#cc00ff40"),
               lcol =list(NULL), mid = list(sims[[i]][[3]]))
    }
    arrows(sims[[i]]$xlim[1],0,sims[[i]]$xlim[2],0, length = 0.1)
    arrows(0,-0.5,0,0.5, length = 0.1)
    points(zmrca[1], zmrca[2], col = "blue", pch = 16, cex = 1.5)
    mapply(arrows,
           arrows_p1[,1], arrows_p1[,2],
           arrows_p1[,3], arrows_p1[,4],
           length = 0.05)
    mapply(arrows,
           arrows_p2[,1], arrows_p2[,2],
           arrows_p2[,3], arrows_p2[,4],
           length = 0.05, col = "red")
    text(0, 0.55, "Intrinsic", cex = 1.5)
    if(i == 1) {
        text(0.6, -0.05, "Extrinsic", cex = 1.5)
    } else if (i == 2) {
        text(1.1, -0.05, "Extrinsic", cex = 1.5)
    } else {
        text(0.8, -0.05, "Extrinsic", cex = 1.5)
    }
    mtext(titles[i], side = 3)
    
    #============
    plot_twotraits(m, zp1, sims[[i]]$opt1)
    plot_onetrait(m, zp1, sims[[i]]$opt1)
    
}
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