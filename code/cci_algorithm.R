cci.algorithm = function(spectrum1, spectrum2, outfile, min.mass, max.mass, interval){
  
  # Declare mass.class and ccf.val
  
  mass.class = sort(unique(cut(min.mass:max.mass, breaks = seq(min.mass,max.mass, by=interval), include.lowest = T, right = F)))
  ccf.val = vector()
  
  # Process spectra
  
  spectra = MALDIquant::alignSpectra(c(spectrum1, spectrum2))
  spectra = trim(spectra,  range = c(min.mass, max.mass))
  metaData(spectra[[1]])$mass.class = cut(spectra[[1]]@mass, breaks=seq(min.mass,max.mass, by=interval), include.lowest = T, right = F)
  metaData(spectra[[2]])$mass.class = cut(spectra[[2]]@mass, breaks=seq(min.mass,max.mass, by=interval), include.lowest = T, right = F)
  
  # Calculate ccf values
  
  for (j in mass.class) {
    
    ss.temp1 =  spectra[[1]][metaData(spectra[[1]])$mass.class == j]
    ss.temp2 =  spectra[[2]][metaData(spectra[[2]])$mass.class == j]
    
    temp = ccf(ss.temp1@intensity,ss.temp2@intensity, plot = FALSE)
    temp = temp$acf[which(diff(sign(diff(temp$acf)))==-2)+1]
    temp = max(temp[temp>0])
    
    ccf.val = c(ccf.val, temp)
    
  }
  
  outvec = as.list(c(names(spectrum1), names(spectrum2), ccf.val))
  
  # Write result to outputfile
  
  fwrite(outvec, outfile, append = TRUE, row.names =F, col.names = FALSE, sep = ',', qmethod = 'escape')
  
}