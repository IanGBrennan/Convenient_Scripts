# Write shell scripts for IQtree, RAxML, whatever

shell.IQ <- function(path.to.files, bootstraps, cores, num.sets, file.ext){
  infiles <- dir(path.to.files, pattern=file.ext)
  
  sets <- split(infiles, ceiling(seq_along(infiles)/num.sets))
  for (j in 1:sets){
    shscript <- NULL
    for (k in 1:length(sets[j]))
    
  }
}

doo <- seq(1, 100, 0.5)
testo <- split(doo, ceiling(seq_along(doo)/20))

seq_along(doo)/20
split(doo, 10)
