synth_fun = function (time, values, nens = 1, sameTrend = TRUE, index.to.model = NA) 
{
  if (is.list(time)) {
    time = time$values
  }
  if (is.list(values)) {
    values = values$values
  }
  if (any(is.na(index.to.model))) {
    index.to.model = seq_along(time)
  }
  orig.time = time
  time = as.matrix(time[index.to.model])
  values = as.matrix(values[index.to.model])
  tnai = which(is.na(time))
  vnai = which(is.na(values))
  trend = predict(lm(values ~ time))
  trendmodel = lm(values ~ time)
  longtrend = orig.time * trendmodel$coefficients[2] + trendmodel$coefficients[1]
  notrend = values - trend
  m = mean(notrend, na.rm = TRUE)
  s = sd(notrend, na.rm = TRUE)
  a = acf(notrend, na.action = na.pass, plot = FALSE)
  ar = max(0, as.numeric(unlist(a[1])[1]))
  fit = arima(x = notrend, order = c(1, 0, 0), method='ML')
  synValues = matrix(NA, nrow = length(orig.time), ncol = nens)
  for (jj in 1:nens) {
    rdata = arima.sim(model = list(ar = ar), n = length(orig.time))
    if (sameTrend) {
      rtrend = predict(lm(rdata ~ orig.time))
      rdata = rdata - rtrend
      rdata = scale(rdata) * s + m
      withTrend = rdata + longtrend
      withTrend[vnai] = NA
      synValues[, jj] = withTrend
    }
    else {
      synValues[, jj] = rdata
    }
  }
  return(synValues)
}