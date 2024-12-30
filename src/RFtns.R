simulateNHPP = function(k, intensity_ftn, maxintensity, maxT){
  points = c()
  newTime = 0
  
  while(TRUE){
    newTime = newTime + rexp(1, maxintensity)
    if(newTime > maxT){
      break
    }
    if(runif(1) < intensity_ftn(newTime, k) / maxintensity){
      points = rbind(points, c(newTime, k))
    }
  }
  
  return( points )
}


simulateTrig = function(l, intensity_ftn, maxintensity, maxT, distmat, soundspeed){
  points = c()
  
  for(k in 1:K){
    
    newTime = 0
    while(TRUE){
      newTime = newTime + rexp(1, maxintensity)
      if(newTime > maxT){
        break
      }

      if( (newTime > distmat[l,k] / soundspeed) & (runif(1) < intensity_ftn(l, newTime, distmat[l,k]) / maxintensity) ) { # sound speed
        points = rbind(points, c(newTime, k))
      }
    }
    
  }
  
  return( points )
}


simulateMultiselfex = function(
    Xm, distmat, maxT, lam0m, knts, 
    beta0, beta, Wm = NULL, delta = NULL, zeta, eta, phi, soundspeed,
    displayOutput = TRUE){
  
  K = nrow(distmat)
  
  maxh = 1
  maxtrig = zeta * maxh
  maxtrig[maxtrig == 0] = min(maxtrig[maxtrig != 0])
  maxlam0 = sapply(1:K, function(k) max(lam0m[,k]))
  
  if( is.null(Wm) ){
    fn_lam0 = function(t, k){
      lam0 = compLam0k(k, t, k, maxT, beta0[k], Xm[[k]], beta[,k], 0, rep(0, nrow(lam0m)), which(knts >= t)[1] - 1 - 1, knts)
      return( as.vector(lam0) ) 
    } 
    
  } else {
    fn_lam0 = function(t, k){
      lam0 = compLam0k(k, t, k, maxT, beta0[k], Xm[[k]], beta[,k], delta[k], Wm, which(knts >= t)[1] - 1 - 1, knts)
      return( as.vector(lam0) )
    }
  }
  
  
  x = foreach(k = 1:K, .combine = 'rbind') %do% {
    simulateNHPP(k, fn_lam0, maxlam0[k], maxT)
  }
  
  if (displayOutput == TRUE) { 
    print( sprintf(paste0("%s events generated from the background process"), nrow(x)) )
    print( sprintf(paste0("(", 1:10, ") %s"), table(x[,2])) )
  }
  
  if( all(zeta == 0) ) {
    ord = order(x[,1])
    x = x[ord,]
    caused = rep(0, nrow(x))
    
  } else {
    count = 1
    caused = rep(0, nrow(x))
    while (TRUE) {
      if ( (count > nrow(x)) ) { break }
      pt = x[count,]
      
      fn_trig = function(l, tdiff, sdiff) {
        return( zeta[l] * compH(tdiff, sdiff, soundspeed, eta, phi)  )
      }
      
      pts = simulateTrig(pt[2], fn_trig, maxtrig[pt[2]], maxT - pt[1], distmat, soundspeed)
      if( !is.null(pts) ){
        pts = matrix(pts, ncol = 2)
        x = rbind(x, cbind(pts[,1] + pt[1], pts[,2]))
        if (displayOutput == TRUE) { print( sprintf("%s events generated so far", nrow(x)) ) }
        caused = c(caused, rep(pt[1], nrow(pts)))
      }
      count = count + 1
    }
    ord = order(x[,1])
    x = x[ord,]
    caused = caused[ord]
    for (i in 1:length(caused)) {
      if (caused[i] == 0) {
        next
      }
      caused[i] = which(x[,1] == caused[i])
    } 
  }
  
  return(list(ts = x, branching = caused))
}

