cci.reproducibility = function(dta_cci){
  
  output = data.frame(matrix(nrow = 0, ncol = 5))
  
  for (i in unique(c(dta_cci$spectrum.in, dta_cci$spectrum.out))) {
    
    temp = dta_cci[dta_cci$spectrum.in == i | dta_cci$spectrum.out == i, ]
    n = nrow(temp)
    cci_vals = log10(quantile(temp$cci, probs = c(0,0.5,1)))
    
    newline = c(i, n, cci_vals)
    output = rbind(output, newline)
    
  }
  
  names(output) = c("spectrum", "n", "min.cci", "median.cci", "max.cci")
  output$n = as.numeric(output$n)
  output$min.cci = as.numeric(output$min.cci)
  output$median.cci = as.numeric(output$median.cci)
  output$max.cci = as.numeric(output$max.cci)
  return(output)

}