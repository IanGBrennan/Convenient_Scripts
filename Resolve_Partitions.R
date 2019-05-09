# Short script to resolve the partitions of a concatenated file

positions <- read.csv("/Users/Ian/Desktop/Geckomics_Partition_Positions.csv", header=F)

pos.bin <- NULL
for (k in 1:nrow(positions)){
  if(k==1){pos.bin <- data.frame(start=1, to="-", end=positions[1,1])}
  else if (k>1){
    start.pos <- pos.bin[k-1,"end"] + 1
    end.pos <- pos.bin[k-1, "end"] + positions[k,1]
    pos.bin <- rbind(pos.bin, data.frame(start=start.pos, to="-", end=end.pos))
  }
}
write.csv(pos.bin, file="/Users/Ian/Desktop/Geckomics_Partition.csv")
