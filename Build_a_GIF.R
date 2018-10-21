library(purrr); library(magick);

# Step 2: List those Plots, Read them in, and then make animation
list.files(path = paste0(getwd(),"/","Desktop/PRSB_gif","/"), pattern = "*.jpeg", full.names = T) %>% 
  map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=0.25) %>% # animates, can opt for number of loops
  image_write(paste0(getwd(),"/Desktop/PRSB_gif/","PRSB",".gif")) # write to current dir

getwd()
