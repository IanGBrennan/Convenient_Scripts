#install.packages("purrr"); install.packages("magick")
library(purrr); library(magick);

setwd("PATH_TO_YOUR_DIRECTORY")

# Step 2: List those Plots, Read them in, and then make animation
list.files(path = paste0(getwd(),"/"), pattern = "*.jpeg", full.names = T) %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=0.25) %>% # animates, can set the Frames Per Second or number of loops
  image_write(paste0(getwd(),"/","FILE_NAME",".gif")) # write to current dir