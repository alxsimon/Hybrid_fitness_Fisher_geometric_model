# Randomization for 10,000 Brownian bridges,
# saved as .rds for easier later use by plotting functions.
library(caTools)

params <- read.table("source_code/Fraisse2016/parameters", row.names = NULL)[,-1]
disp_params <- c(1,2,4,5,6,9)
text_params <- c("N",
                 "n",
                 "Moving optimum P1",
                 "Moving optimum P2",
                 expression(alpha),
                 "MVN model",
                 "d")

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

#for (i in 1:nrow(params)){
for (i in 1:3){
    N <- params$N[i]
    n <- params$n[i]
    s <- params$s[i]
    diff1 <- params$diff1[i]
    diff2 <- params$diff2[i]
    k <- params$k[i]
    nmuts <- params$nmuts[i]
    ncycles <- params$ncycles[i]
    MVN <- params$MVN[i]
    sigma <- params$sigma[i]
    
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
        d <- nrow(m)
    } else {
        common_anc <- sum(sim_P0[,coltraits])
        zp1 <- sum(c(sim_P0[,coltraits], sim_P1[,coltraits]))
        zp2 <- sum(c(sim_P0[,coltraits], sim_P2[,coltraits]))
        m <- c(-sim_P1[nrow(sim_P1):1,coltraits],
               sim_P2[,coltraits])
        d <- length(m)
    }
    
    rep <- 1000
    z_bridge <- sapply(1:rep, draw_z2, m = m, zp1 = zp1)
    saveRDS(z_bridge, file = paste0("rand_samples_", i, ".rds"))
}

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