# Takes a coda object and a variable name and returns an array with the mcmc output along the first dimension
mcmc2array = function(codaObject, varName, stack = TRUE) {
  stopifnot(is.mcmc(codaObject) | is.mcmc.list(codaObject))
  stopifnot(is.character(varName) & length(varName) == 1)
  if (is.mcmc(codaObject)) {
    invisible(.mcArr(codaObject, varName))
  }
  out = lapply(codaObject, .mcArr, varName = varName)
  if (!stack)
    invisible(out)
  invisible(do.call(rbind, lapply(out, as.matrix)))
}


.mcArr = function(codaObject, varName) {
  nIter = niter(codaObject)
  arrayNames = varnames(codaObject)[grep(paste0("^", varName, "(\\[|$)"), varnames(codaObject))]  
  if (length(arrayNames) == 0) {
    stop(sprintf("Variable %s not found.", varName))
  }
  if(length(arrayNames) == 1) 
    return(as.array(codaObject[, arrayNames]))  
  coords = parseDim(gsub(varName, "", arrayNames))
  if(!is.matrix(coords)) {
    coords = t(coords)
  }
  dims = apply(coords, 1, max)
  stopifnot(prod(dims) == length(arrayNames))
  varArray = array(dim = c(nIter, dims))
  for (i in 1:length(arrayNames)) {
    if (length(dims) == 1) {
      varArray[1:nIter, i] = codaObject[, arrayNames[i]]
    } else {
      varArray[1:nIter + cumprod(c(nIter, dims[seq.pos(1,length(dims) - 1)])) %*% (coords[, i, drop = FALSE] - 1) ] = codaObject[, arrayNames[i]]
    }
  }
  invisible(varArray)
}

seq.pos = function(start, end) {
  if(start <= end)
    return(start:end)
  else 
    return(integer(0))
}

parseDim = function(strDim) {
  sapply(strDim, function(str) as.integer(unlist(strsplit(substr(str, 2, nchar(str)-1), ","))))
}

