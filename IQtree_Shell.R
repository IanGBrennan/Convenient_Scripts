setwd("/Users/Ian/Desktop/UCE_Alignments/Trimmed_Alignments/Top500_Loci")
#alignment.files <- align.files
alignment.files <- dir(getwd(), pattern =".phy")
iq.path <- "/Applications/iqtree-1.7.11/bin/iqtree"

iqtree.shell <- function(in.files, batch.size=10){
  curr.dir <- getwd()
  dir.create(paste0(curr.dir,"/IQTREE_Shells")) # make a directory for the new files
  
  splits <- split(in.files, ceiling(seq_along(in.files)/batch.size))
  for (k in 1:length(splits)){
    cd.call <- paste("cd", curr.dir) 
    write.table(cd.call, file=paste0(curr.dir,"/IQTREE_Shells/", "IQTREE_Shell",k,".sh"), append=T, row.names=F, col.names=F, quote=F)
    for (g in 1:length(splits[[k]])){
      curr.call <- paste(iq.path, "-s", paste0(curr.dir,"/",splits[[k]][g]), "-bb 1000 -nt AUTO")
      write.table(curr.call, file=paste0(curr.dir,"/IQTREE_Shells/", "IQTREE_Shell",k,".sh"), append=T, row.names=F, col.names=F, quote=F)
      permission.call <- paste("chmod =rwx,g+s", paste0(curr.dir,"/IQTREE_Shells/", "IQTREE_Shell",k,".sh"))
      system(permission.call)
    }
  }
  shell.files <- dir(paste0(curr.dir,"/IQTREE_Shells"), pattern=".sh")
  cd.call <- paste("cd", curr.dir) 
  shell.calls <- NULL
  for (q in 1:length(shell.files)){shell.calls <- append(shell.calls, paste("source", paste0(curr.dir, "/IQTREE_Shells/", shell.files[q]), "&"))}
  combo <- c(cd.call, shell.calls)
  write.table(combo, file=paste0(curr.dir,"/IQTREE_Shells/", "MASTER_IQTREE_Shell.sh"), append=F, row.names=F, col.names=F, quote=F)
  permission.call <- paste("chmod =rwx,g+s", paste0(curr.dir,"/IQTREE_Shells/", "MASTER_IQTREE_Shell.sh"))
  system(permission.call)
}
iqtree.shell(alignment.files, batch.size=65)

