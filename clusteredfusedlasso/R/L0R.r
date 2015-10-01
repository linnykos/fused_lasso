L0Msg = function() {
  # a message is piecewise quadratic with knot locations stored in 'x'
  return(list(x = -Inf, const = 0, lin = 0, quad = 0))
}

L0MsgNumKnots = function(msg) {
  return(length(msg$x))
}

L0MsgFindKnot = function(msg, x) {
  # Finds the index of the segment of 'msg' which contains 'x'
  
  nKnot = L0MsgNumKnots(msg)
  for (j in 1:nKnot) {
    if (msg$x[j] <= x && (nKnot == j || msg$x[j+1] >= x)) {
      return(j)
    }
  }
  stop("knot not found")
}

L0MsgEval = function(msg, x, knotIdx = NULL) {
  # 'msg' is a piecewise quadratic function and this evaluates it
  # at the point 'x'  i.e.  msg(x)
  
  if (is.null(knotIdx)) knotIdx = L0MsgFindKnot(msg, x)
  return(msg$const[knotIdx] + msg$lin[knotIdx] * x + msg$quad[knotIdx] * x * x)
}

L0DPPlotMsg = function(msgObj, xlim = NULL, nPts = 100, add = FALSE,
                       plotKnots = TRUE, ...) {
  nKnot = L0MsgNumKnots(msgObj)

  if (is.null(xlim)) {
    xLwr = min(msgObj$x[msgObj$x > -Inf])
    xUpr = max(msgObj$x[msgObj$x < Inf])

    if (sum(is.finite(msgObj$x)) == 0 || nKnot <= 2) xLwr = -5
    if (sum(is.finite(msgObj$x)) == 0 || nKnot <= 2) xUpr = 5
    xlim = c(xLwr - 0.2 * (xUpr - xLwr),
             xUpr + 0.2 * (xUpr - xLwr))
  }

  evalPts = seq(from = xlim[1], to = xlim[2], length = nPts)
  y = sapply(evalPts, L0MsgEval, msg = msgObj)

  if (add) {
    lines(evalPts, y, ...)
  } else {
    plot(evalPts, y, type = "l", ...)
  }
  
  if (plotKnots) abline(v = msgObj$x, col = "gray")
  return(rbind(evalPts, y))
}

L0MsgMax = function(msg) {
  # Performs  \argmax_x msg(x) and returns the value and segment index as well
  nKnot = L0MsgNumKnots(msg)
  if (nKnot == 1 && msg$quad[1] == 0) {
    if (msg$lin[1] != 0) stop("invalid message (1)")
    return(msg$y[1])
  }
  
  yMax = xMax = knotIdx = -Inf    
  for (j in 1:nKnot) {
    xNext = ifelse(j == nKnot, Inf, msg$x[j+1])
    lin = msg$lin[j]  
    quad = msg$quad[j]
    x = msg$x[j]
    
    if (quad == 0) stop("LS segmentation only")
    
    yLeft = msg$const[j] + lin * x + quad * x * x
    yRight = msg$const[j] + lin * xNext + quad * xNext * xNext
                                                             
    if (!is.finite(x)) yLeft = -Inf
    if (!is.finite(xNext)) yRight = -Inf
    
    xBest = -lin / (2 * quad)
    if (xBest < x) xBest = x 
    if (xBest > xNext) xBest = xNext
    yBest = msg$const[j] + lin * xBest + quad * xBest * xBest
    
    checkX = c(x, xNext, xBest)
    checkY = c(yLeft, yRight, yBest)
    
    for (k in 1:3) {
      if (checkY[k] > yMax) {
        yMax = checkY[k]
        xMax = checkX[k]
        knotIdx = j
      }
    }
  }
  return(list(knotIdx = knotIdx, x = xMax, y = yMax))
}

L0MsgSegCrosses = function(msg, knotIdx, val) {
  # Returns where this segment of 'msg' cross 'val'
  # Also indicates if this is an up or down-crossing
  # i.e.  finds A := \{x : msg(x) = val\}
  #   and evalutes the derivative msg'(x) for each x \in A
  # We assume that the the segments of msg are concave, so if a segment
  # crosses 'val' twice the first is an upcrossing and the second is
  # a downcrossing.
  #
  # Some extra care is needed if the quadratic part of the message 'quad' is zero
  # (e.g. LAD segmentation).  To keep the code below simpler, we do not check
  # for these cases
  
  nKnot = L0MsgNumKnots(msg)
  xNext = ifelse(knotIdx == nKnot, Inf, msg$x[knotIdx+1])
   
  c0 = msg$const[knotIdx] 
  lin = msg$lin[knotIdx]  
  quad = msg$quad[knotIdx]
  x = msg$x[knotIdx]
  
  if (quad == 0) stop("LS segmentation only")
  
  up = xr = c()
  
  z = (lin * lin - 4 * quad * (c0 - val))
  if (z < 0) {
    return(list(up = c(), x = c()))
  } else if (quad > 0) {
    solns = c((-lin - sqrt(z)) / (2 * quad),
              (-lin + sqrt(z)) / (2 * quad))
  } else {
    solns = c((-lin + sqrt(z)) / (2 * quad),
              (-lin - sqrt(z)) / (2 * quad))
  }
  
  for (j in 1:2) {
    if (solns[j] >= x && solns[j] <= xNext) {
      up = c(up, (j == 1))
      xr = c(xr, solns[j])
    }
  }
  return(list(up = up, x = xr))
}

L0MsgArgMax = function(msg, lambda2) {
  # Creates the function  f(x) = \max_y \{msg(y) + \lambda_2 * 1\{y \neq x\}\}
  # and in 'backSegs' stores the intervals such that the argmax was
  #  \max_y msg(y)
  AddKnot = function(m, x0, const0, lin0, quad0) {
    m$x = c(m$x, x0)
    m$const = c(m$const, const0)
    m$lin = c(m$lin, lin0)
    m$quad = c(m$quad, quad0)
    return(m)
  }
   
  maxPt = L0MsgMax(msg)
  waterLev = maxPt$y - lambda2
  nKnot = L0MsgNumKnots(msg)

  retMsg = list(x = -Inf, const = waterLev, lin = 0, quad = 0)

  belowWater = TRUE
  backSegs = -Inf
  
  for (knotIdx in 1:nKnot) {
    crs = L0MsgSegCrosses(msg, knotIdx, waterLev)
    nCross = length(crs$up)
    
    if (belowWater) { 
      if (nCross >= 1) {
        retMsg = AddKnot(retMsg, crs$x[1], msg$const[knotIdx], msg$lin[knotIdx],
                         msg$quad[knotIdx])
        backSegs = c(backSegs, crs$x[1])
        belowWater = FALSE
      }
      
      if (nCross == 2) {
        retMsg = AddKnot(retMsg, crs$x[2], waterLev, 0, 0)
        backSegs = c(backSegs, crs$x[2])
        belowWater = TRUE
      }
    } else {  # if (belowWater)
      retMsg = AddKnot(retMsg, msg$x[knotIdx], msg$const[knotIdx],
                       msg$lin[knotIdx], msg$quad[knotIdx])
        
      if (nCross == 1) {
        retMsg = AddKnot(retMsg, crs$x[1], waterLev, 0, 0)
        backSegs = c(backSegs, crs$x[1])
        belowWater = TRUE
      }
    }   
  }
  
  if (length(backSegs) %% 2 != 1) stop("message appears invalid (5)")
  backSegs = matrix(c(backSegs, Inf), nrow = 2)
  
  return(list(backSegs = backSegs, msg = retMsg,
              maxPt = maxPt))
}

LSL0MsgAddObs = function(msg, x, wt = -0.5) {
  # Adds the term   e(h) = wt * (h - x)^2
  # to the message.  The constant shift is  (x * x * wt),
  # the linear shift is (-2 * wt * x)
  # and the quadratic shift is (wt)
  msg$const = msg$const + (x * x * wt)
  msg$lin = msg$lin + (-2.0 * x * wt)
  msg$quad = msg$quad + wt
  return(msg)
}

L0Seg = function(x, lambda2, wts = NULL, addObsFunc = LSL0MsgAddObs) {
  # Performs the dynamic programming algorithm using penalty lambda2
  n = length(x)
  preMsgs = vector(mode = "list", length = n)
  vitMsgs = vector(mode = "list", length = n)
  
  if (is.null(wts)) wts = rep(-0.5, n)
  if (length(wts) != n) wts = rep(wts, n)
  
  curMsg = L0Msg()
  
  for (j in 1:n) {
    preMsgs[[j]] = addObsFunc(curMsg, x[j], wts[j])
    vitMsgs[[j]] = L0MsgArgMax(preMsgs[[j]], lambda2)
    curMsg = vitMsgs[[j]]$msg
  }
  
  fit = rep(NA, n)
  fit[n] = vitMsgs[[n]]$maxPt$x
  for (j in (n-1):1) {
    backSegs = vitMsgs[[j]]$backSegs
    nSegs = ncol(backSegs)
    
    inSeg = FALSE
    for (segIdx in 1:nSegs) {
      if (fit[j+1] >= backSegs[1, segIdx] && fit[j+1] <= backSegs[2, segIdx]) {
        inSeg = TRUE
        fit[j] = vitMsgs[[j]]$maxPt$x
      }
    }
    if (!inSeg) {
      fit[j] = fit[j + 1]
    }
  }
  return(list(fit = fit, lambda2 = lambda2,
              x = x, wts = wts, addObsFunc = addObsFunc,
              vitMsgs = vitMsgs, preMsgs = preMsgs))
}