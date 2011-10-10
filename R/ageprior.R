find.prior <- function(id,ty,to,r,p_q,lag=0,corr='none',plot=TRUE)
{
  if (to<ty) {print("ERROR: to < ty")}
  else {
    # calculate the mean of exponential distribution P't (tm) and from that, calculate lambda
    if (lag>0) {
      Pdasht <- function(t) {exp(((-2*r)/(p_q)) * (exp((p_q)*t) -1)) * (1-exp(-(1/lag)*t))}
      area <- integrate(Pdasht,0,Inf)$value
      Ptt <- function(t) {t * (1/area) * exp(((-2*r)/(p_q)) * (exp((p_q)*t) -1)) * (1-exp(-(1/lag)*t))} #only needed to calculate the mean of Pt (tm)
    }
    else {
      Pdasht <- function(t) {exp(((-2*r)/(p_q)) * (exp((p_q)*t) -1))}      
      area <- integrate(Pdasht,0,Inf)$value
      Ptt <- function(t) {t * (1/area) * exp(((-2*r)/(p_q)) * (exp((p_q)*t) -1))} #only needed to calculate the mean of Pt (tm)
    }
    tm <- integrate(Ptt,0,Inf)$value
    tmm <- (to-ty)/2 + tm
    lambda <- 1/tm

    # specify the precision (increases computation time)
    stepsPerUnit <- 10.0

    # calculate, up to which t calculations should be performed
    endY <- 0.001
    if (to==0) {endX <- (-1*log(endY/lambda))/lambda}
    else {
      if (to==ty) {endX <- (-1*log(endY/lambda))/lambda + ty}
      else {
        logNom <- endY*(to-ty)
        logDenom <- exp(lambda*(to-ty))-1
        nom <- log(logNom/logDenom)
        denom <- -lambda
        endX <- round(nom/denom + ty)        
      }
    }
    if (corr != 'none') {endX <- 5*endX}


    # produce a correction curve (Ct)
    if (corr != 'none') {
      data(list=c(corr))
      if (corr=="global_marine") {
        pys <- global_marine$py
        relrs <- global_marine$relr
      }
      if (corr=="global_terrestrial") {
        pys <- global_terrestrial$py
        relrs <- global_terrestrial$relr
      }
      if (corr=="global") {
        pys <- global$py
        relrs <- global$relr
      }
      if (corr=="NorthAmerica_marine") {
        pys <- NorthAmerica_marine$py
        relrs <- NorthAmerica_marine$relr
      }
      if (corr=="NorthAmerica_terrestrial") {
        pys <- NorthAmerica_terrestrial$py
        relrs <- NorthAmerica_terrestrial$relr
      }
      if (corr=="NorthAmerica") {
        pys <- NorthAmerica$py
        relrs <- NorthAmerica$relr
      }
      if (corr=="WesternEurope_marine") {
        pys <- WesternEurope_marine$py
        relrs <- WesternEurope_marine$relr
      }
      if (corr=="WesternEurope_terrestrial") {
        pys <- WesternEurope_terrestrial$py
        relrs <- WesternEurope_terrestrial$relr
      }
      if (corr=="WesternEurope") {
        pys <- WesternEurope$py
        relrs <- WesternEurope$relr
      }
      if (corr=="Australia_marine") {
        pys <- Australia_marine$py
        relrs <- Australia_marine$relr
      }
      if (corr=="Australia_terrestrial") {
        pys <- Australia_terrestrial$py
        relrs <- Australia_terrestrial$relr
      }
      if (corr=="Australia") {
        pys <- Australia$py
        relrs <- Australia$relr
      }
      cts <- c()
      for (i in 1:(endX*stepsPerUnit)) {
        t <- i/stepsPerUnit
        for (z in 1:(length(relrs)-1)) {
          if (t>=pys[z] & t<pys[z+1]) {ct <- relrs[z]}
        }
        if (t>=pys[length(relrs)]) {ct <- relrs[length(relrs)]}
        cts <- c(cts,ct)
      }

      # calculate the scale factor for all cts depending on ty
      for (z in 1:(length(relrs)-1)) {
        if (ty>=pys[z] & ty<pys[z+1]) {scale <- 1/relrs[z]}
      }
      if (ty>=pys[length(relrs)]) {scale <- 1/relrs[length(relrs)]}
    }


    if (corr=="none") {
      # calculate all values of Rt
      rts <- c()
      for (i in 1:(endX*stepsPerUnit)) {
        t <- i/stepsPerUnit
        if (to==0) {
          if (lag>0) {rt <- ( exp(((-2*r)/(p_q)) * (exp((p_q)*t) -1)) * (1-exp(-(1/lag)*t)) )/area}
          else {rt <- ( exp(((-2*r)/(p_q)) * (exp((p_q)*t) -1)) )/area}
        }
        else {
          if (to==ty) {
            if (t<ty) {rt <- 0} else {
              if (lag>0) {rt <- ( exp(((-2*r)/(p_q)) * (exp((p_q)*(t-ty)) -1)) * (1-exp(-(1/lag)*(t-ty))) )/area}
              else {rt <- ( exp(((-2*r)/(p_q)) * (exp((p_q)*(t-ty)) -1)) )/area}
            }
          }
          else {
            if (t<ty) {rt <- 0} else {
              if (t<to) {
                if (lag>0) {tmp <- function(dt) {( exp(((-2*r)/(p_q)) * (exp((p_q)*(t-ty-dt)) -1)) * (1-exp(-(1/lag)*(t-ty-dt))) )/area}}
                else {tmp <- function(dt) {( exp(((-2*r)/(p_q)) * (exp((p_q)*(t-ty-dt)) -1)) )/area}}
                rt <- integrate(tmp,0,t-ty)$value / (to-ty)
              }
              else {
                if (lag>0) {tmp <- function(dt) {( exp(((-2*r)/(p_q)) * (exp((p_q)*(t-ty-dt)) -1)) * (1-exp(-(1/lag)*(t-ty-dt))) )/area}}
                else {tmp <- function(dt) {( exp(((-2*r)/(p_q)) * (exp((p_q)*(t-ty-dt)) -1)) )/area}}
                rt <- integrate(tmp,0,to-ty)$value / (to-ty)
              }
            }
          }
        }
        rts <- c(rts,rt)
      }      
    } else {
      rt_corrs <- c()
      ctOld <- 1
      ctNew <- 1
      noRtsCalculatedYet <- TRUE
      for (i in 1:(endX*stepsPerUnit)) {
        t <- i/stepsPerUnit
        if (t < ty) {rt_corr <- 0} else {
          ctNew <- cts[i]*scale
          if (noRtsCalculatedYet | ctNew != ctOld) {

            # unless it's the very first Rt curve, calculate rtHereOld for later adjustment to the previous one
            if (noRtsCalculatedYet == FALSE) {rtHereOld <- rts[i]}

            # calculate all values of Rt with scaled preservation rates (r_scaled instead of r)
            r_scaled <- r*ctNew
            rts <- c()
            for (j in 1:(endX*stepsPerUnit)) {
              t_ <- j/stepsPerUnit
              if (to==0) {
                if (lag>0) {rt <- ( exp(((-2*r_scaled)/(p_q)) * (exp((p_q)*t_) -1)) * (1-exp(-(1/lag)*t_)) )/area}
                else {rt <- ( exp(((-2*r_scaled)/(p_q)) * (exp((p_q)*t_) -1)) )/area}
              }
              else {
                if (to==ty) {
                  if (t_<ty) {rt <- 0} else {
                    if (lag>0) {rt <- ( exp(((-2*r_scaled)/(p_q)) * (exp((p_q)*(t_-ty)) -1)) * (1-exp(-(1/lag)*(t_-ty))) )/area}
                    else {rt <- ( exp(((-2*r_scaled)/(p_q)) * (exp((p_q)*(t_-ty)) -1)) )/area}
                  }
                }
                else {
                  if (t_<ty) {rt <- 0} else {
                    if (t_<to) {
                      if (lag>0) {tmp <- function(dt) {( exp(((-2*r_scaled)/(p_q)) * (exp((p_q)*(t_-ty-dt)) -1)) * (1-exp(-(1/lag)*(t_-ty-dt))) )/area}}
                      else {tmp <- function(dt) {( exp(((-2*r_scaled)/(p_q)) * (exp((p_q)*(t_-ty-dt)) -1)) )/area}}
                      rt <- integrate(tmp,0,t_-ty)$value / (to-ty)
                    }
                    else {
                      if (lag>0) {tmp <- function(dt) {( exp(((-2*r_scaled)/(p_q)) * (exp((p_q)*(t_-ty-dt)) -1)) * (1-exp(-(1/lag)*(t_-ty-dt))) )/area}}
                      else {tmp <- function(dt) {( exp(((-2*r_scaled)/(p_q)) * (exp((p_q)*(t_-ty-dt)) -1)) )/area}}
                      rt <- integrate(tmp,0,to-ty)$value / (to-ty)
                    }
                  }
                }
              }
              rts <- c(rts,rt)
            }

            # unless it's the very first Rt curve, adjust it to the previous one
            if (noRtsCalculatedYet == FALSE) {
              rtHereNew <- rts[i]
              if (rtHereNew > 0) {rts <- rts*(rtHereOld/rtHereNew)}
            }

            ctOld <- ctNew
            noRtsCalculatedYet <- FALSE
          }
          rt_corr <- rts[i]
        }
        rt_corrs <- c(rt_corrs,rt_corr)
      }
      rts <- rt_corrs


      # calculate new mean of rts and adjust tmm
      rts_total <- 0
      for (i in 1:(endX*stepsPerUnit)) {
        t <- i/stepsPerUnit
        rts_total <- rts_total + rts[i]*t # multiplication with t is important here
      }
      rts_mean <- rts_total/sum(rts)
      tmm <- rts_mean-ty
      
      # scale rts so that its volume is 1.0
      rts <- rts/(sum(rts)/stepsPerUnit)

      # correct endX if unnecessarily large
      for (i in 2:(endX*stepsPerUnit)) {
        if (rts[i] < 0.001 & rts[i-1] >= 0.001) {endX <- i/stepsPerUnit}
      }
      
      # truncate rts and cts to the shorter endX
      rts <- rts[1:(endX*stepsPerUnit)]
      cts <- cts[1:(endX*stepsPerUnit)]
      
      if (plot==TRUE) {
        # plot cts
        x <- c(1:(endX*stepsPerUnit))
        plot(x/stepsPerUnit,cts,type="l",xlab="Time (my)",ylab="Relative preservation rate",ylim=c(0,max(cts)),main="Preservation rate correction")        
      }
    }    


    #calculate all values of Lt (a lognormal distribution with offset ty, and mean tmm), and find the best sigma value
    ltss <- matrix(nrow=(endX*stepsPerUnit),ncol=(130))
    rmsd_ls <- c()
    sigmas <- c()
    for (a in 1:130) {
      sigma <- a/100
      lts <- c()
      for (i in 1:(endX*stepsPerUnit)) {
        t <- i/stepsPerUnit
        if (t<=ty) {lt <- 0} else {lt <- dlnorm((t-ty), meanlog = log(tmm)-(sigma^2)/2, sdlog = sigma, log = FALSE)}
        lts <- c(lts,lt)
      }
      ltss[,a] <- lts
      sd_l <- 0
      for (i in 1:(endX*stepsPerUnit)) {sd_l <- sd_l+(rts[i]-lts[i])^2}
      msd_l <- sd_l/(endX*stepsPerUnit)
      rmsd_l <- sqrt(msd_l)
      rmsd_ls <- c(rmsd_ls,rmsd_l)
      sigmas <- c(sigmas,sigma)
    }
    bestPos_l <- order(rmsd_ls)[1]
    lts <- ltss[,bestPos_l]
    sigma <- sigmas[bestPos_l]
    mu <- log(0.5*(to-ty)+tm) - ((sigma**2)/2.0)


    # calculate all values of Gt (a gamma distribution with offset ty, and mean tmm), and find the best shape parameter
    gtss <- matrix(nrow=(endX*stepsPerUnit),ncol=(500))
    rmsd_gs <- c()
    shapes <- c()
    for (a in 1:500) {
      shape <- a/100
      scale <- tmm/shape #to make sure the gamma distribution Gt has the same mean as Rt
      gts <- c()
      for (i in 1:(endX*stepsPerUnit)) {
        t <- i/stepsPerUnit
        if (t<ty) {gt <- 0} else {gt <- dgamma((t-ty), shape=shape, scale=scale)}
        gts <- c(gts,gt)
      }
      gtss[,a] <- gts
      sd_g <- 0
      for (i in 1:(endX*stepsPerUnit)) {sd_g <- sd_g+(rts[i]-gts[i])^2}
      msd_g <- sd_g/(endX*stepsPerUnit)
      rmsd_g <- sqrt(msd_g)
      rmsd_gs <- c(rmsd_gs,rmsd_g)
      shapes <- c(shapes,shape)
    }
    bestPos_g <- order(rmsd_gs)[1]
    gts <- gtss[,bestPos_g]
    shape <- shapes[bestPos_g]
    scale <- tmm/shape


    # calculate all values of Et (an exponential distribution with offset ty, and mean tmm)
    ets <- c()
    lambda_e <- 1/tmm
    for (i in 1:(endX*stepsPerUnit)) {
      t <- i/stepsPerUnit
      if (t<ty) {et <- 0} else {et <- lambda_e*exp(-lambda_e*(t-ty))}
      ets <- c(ets,et)
    }
    sd_e <- 0
    for (i in 1:(endX*stepsPerUnit)) {sd_e <- sd_e+(rts[i]-ets[i])^2}
    msd_e <- sd_e/(endX*stepsPerUnit)
    rmsd_e <- sqrt(msd_e)


    # check whether lts, gts, or ets is the best approximation of rts
    rmsds <- c(rmsd_ls[bestPos_l],rmsd_gs[bestPos_g],rmsd_e)
    if (order(rmsds)[1]==1) {
      # Lt is chosen with corresponding sigma value
      if (plot==TRUE) {
        # plot both Rt and Lt
        x <- c(1:(endX*stepsPerUnit))
        if (max(rts)>max(lts)) {
          plot(x/stepsPerUnit,rts,type="l",xlab="Time (my)",ylab="Probability",main="Rt (solid) and Lt (dashed)")
          points(x/stepsPerUnit,lts,type="l",lty=2)
          lines(x=c(ty,ty),y=c(0,max(rts)))
        } else {
          plot(x/stepsPerUnit,lts,type="l",lty=2,xlab="Time (my)",ylab="Probability",main="Rt (solid) and Lt (dashed)")
          points(x/stepsPerUnit,rts,type="l")
          lines(x=c(ty,ty),y=c(0,max(lts)))
        }
        text(x=(ty-1),y=max(rts)/40,label=paste("ty"))
        lines(x=c(to,to),y=c(0,max(rts)))
        text(x=(to-1),y=max(rts)/40,label=paste("to"))
      }
      # This text should go into the BEAST XML input file
      cat(sprintf("				<!-- Age Prior for node \"%s\"                                              -->\n",id))
      cat(sprintf("				<logNormalPrior mean=\"%s\" stdev=\"%s\" offset=\"%s\" meanInRealSpace=\"false\">\n",mu,sigma,ty))
      cat(sprintf("				    <statistic idref=\"tmrca(%s)\"/>\n",id))
      cat(sprintf("				</logNormalPrior>\n\n"))
    } else {
      if (order(rmsds)[1]==2) {
        # Gt is chosen with corresponding shape parameter
        if (plot==TRUE) {
          # plot both Rt and Gt
          x <- c(1:(endX*stepsPerUnit))
          if (max(rts)>max(gts)) {
            plot(x/stepsPerUnit,rts,type="l",xlab="Time (my)",ylab="Probability",main="Rt (solid) and Gt (dashed)")
            points(x/stepsPerUnit,gts,type="l",lty=2)
            lines(x=c(ty,ty),y=c(0,max(rts)))  
          } else {
            plot(x/stepsPerUnit,gts,type="l",lty=2,xlab="Time (my)",ylab="Probability",main="Rt (solid) and Gt (dashed)")
            points(x/stepsPerUnit,rts,type="l")
            lines(x=c(ty,ty),y=c(0,max(gts)))  
          }
          text(x=(ty-1),y=max(rts)/40,label=paste("ty"))
          lines(x=c(to,to),y=c(0,max(rts)))
          text(x=(to-1),y=max(rts)/40,label=paste("to"))
        }
        # This text should go into the BEAST XML input file
        cat(sprintf("				<!-- Age Prior for node \"%s\"                                              -->\n",id))
        cat(sprintf("				<gammaPrior shape=\"%s\" scale=\"%s\" offset=\"%s\">\n",shape,scale,ty))
        cat(sprintf("				    <statistic idref=\"tmrca(%s)\"/>\n",id))
        cat(sprintf("				</gammaPrior>\n\n"))
      }
      else {
        # Et is chosen
        if (plot==TRUE) {
          x <- c(1:(endX*stepsPerUnit))
          plot(x/stepsPerUnit,ets,type="l",lty=2,xlab="Time (my)",ylab="Probability",main="Rt (solid) and Et (dashed)")
          points(x/stepsPerUnit,rts,type="l")
          lines(x=c(ty,ty),y=c(0,max(ets)))
          text(x=(ty-1),y=max(rts)/40,label=paste("ty"))
          lines(x=c(to,to),y=c(0,max(ets)))
          text(x=(to-1),y=max(rts)/40,label=paste("to"))
        }
        # This text should go into the BEAST XML input file
        cat(sprintf("				<!-- Age Prior for node \"%s\"                                              -->\n",id))
        cat(sprintf("				<exponentialPrior mean=\"%s\" offset=\"%s\">\n",tmm,ty))
        cat(sprintf("				    <statistic idref=\"tmrca(%s)\"/>\n",id))
        cat(sprintf("				</exponentialPrior>\n\n"))
      }
    }
  }    
}


adjust.pres.rate <- function(r,ty,t1,t2,corr,plot=TRUE) {
  if (t2<t1) {print("ERROR: t2 < t1")}
  else {
    data(list=c(corr))
    if (corr=="global_marine") {
      pys <- global_marine$py
      relrs <- global_marine$relr
    }
    if (corr=="global_terrestrial") {
      pys <- global_terrestrial$py
      relrs <- global_terrestrial$relr
    }
    if (corr=="global") {
      pys <- global$py
      relrs <- global$relr
    }
    if (corr=="NorthAmerica_marine") {
      pys <- NorthAmerica_marine$py
      relrs <- NorthAmerica_marine$relr
    }
    if (corr=="NorthAmerica_terrestrial") {
      pys <- NorthAmerica_terrestrial$py
      relrs <- NorthAmerica_terrestrial$relr
    }
    if (corr=="NorthAmerica") {
      pys <- NorthAmerica$py
      relrs <- NorthAmerica$relr
    }
    if (corr=="WesternEurope_marine") {
      pys <- WesternEurope_marine$py
      relrs <- WesternEurope_marine$relr
    }
    if (corr=="WesternEurope_terrestrial") {
      pys <- WesternEurope_terrestrial$py
      relrs <- WesternEurope_terrestrial$relr
    }
    if (corr=="WesternEurope") {
      pys <- WesternEurope$py
      relrs <- WesternEurope$relr
    }
    if (corr=="Australia_marine") {
      pys <- Australia_marine$py
      relrs <- Australia_marine$relr
    }
    if (corr=="Australia_terrestrial") {
      pys <- Australia_terrestrial$py
      relrs <- Australia_terrestrial$relr
    }
    if (corr=="Australia") {
      pys <- Australia$py
      relrs <- Australia$relr
    }
    if (t1>=pys[length(relrs)]) {r_corr <- r} else {
      sum <- 0
      for (z in 1:(length(relrs)-1)) {
        # t1 smaller, t2 inside
        if (t1<=pys[z] & t2>pys[z] & t2<=pys[z+1]) {sum <- sum + (t2-pys[z])*relrs[z]}
        # t1 smaller, t2 larger
        if (t1<=pys[z] & t2>pys[z+1]) {sum <- sum + (pys[z+1]-pys[z])*relrs[z]}
        # t1 inside, t2 inside
        if (t1>pys[z] & t2<=pys[z+1]) {sum <- sum + (t2-t1)*relrs[z]}
        # t1 inside, t2 larger
        if (t1>pys[z] & t1<=pys[z+1] & t2>pys[z+1]) {sum <- sum + (pys[z+1]-t1)*relrs[z]}
      }
      if (t1<pys[length(relrs)] & t2>pys[length(relrs)]) {sum <- sum + (t2-pys[length(relrs)])*relrs[length(relrs)]}
      avg_relr <- sum/(t2-t1)
      for (z in 1:length(relrs)) {
        if (ty>=pys[z]) {ty_relr <- relrs[z]}
      }
      r_corr <- r*(ty_relr/avg_relr)
    }
    cts <- c()
    for (i in 1:pys[length(relrs)]*100) {
      t <- i/100
      for (z in 1:length(relrs)) {
        if (t>=pys[z]) {ct <- relrs[z]}
      }
      cts <- c(cts,ct)
    }
    if (plot==TRUE) {
      rs <- cts*r*(1/avg_relr)
      x <- c(1:pys[length(relrs)]*100)
      plot(x/100,rs,type="l",xlab="Time (my)",ylab="r",ylim=c(0,max(rs)),main="Preservation rate")
      lines(x=c(ty,ty),y=c(0,max(rs)),lty=2)
      text(x=(ty-7),y=max(rs)/40,label=paste("ty"))
      lines(x=c(t1,t1),y=c(0,max(rs)),lty=2)
      text(x=(t1-7),y=max(rs)/40,label=paste("t1"))
      lines(x=c(t2,t2),y=c(0,max(rs)),lty=2)
      text(x=(t2-7),y=max(rs)/40,label=paste("t2"))
      lines(x=c(t1,t2),y=c(r,r),lty=2)
      text(x=(t2-7),y=max(rs)/40,label=paste("r"))      
    }
    cat(sprintf("The corrected preservation rate r, at time ty is: %s",r_corr))
    if (ty>t2 | ty<t1) {cat(sprintf("\nWARNING: The fossil age may be outside the period for which the given preservation rate had been calculated!"))}
  }
}


net.div.rate <- function(n,dt,sat=FALSE) {
  if (sat) {
    p_q <- log(n/2)/(dt/2)
    cat(sprintf("Assuming exponential net diversification with saturation, the initial rate (p-q) is: %s",p_q))
  } else {
    p_q <- log(n)/dt
    cat(sprintf("Assuming exponential net diversification, the rate (p-q) is: %s",p_q))
  }
}