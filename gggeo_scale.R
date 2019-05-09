gggeo_scale <- function(gg, fill = NULL, color = "black", alpha = 1, height = .05, size = 5, quat = FALSE, pos = "bottom", abbrv = TRUE, periods = NULL, neg = FALSE) {
  #This is a function to add a geologic time scale to a ggplot object.
  #gg: the ggplot object
  #fill: the fill color of the boxes; default are the colors from the Commission for the Geological Map of the World (CGMW);
  #   custom fill colors can be provided and will be recycled if necessary
  #   if a custom dataset is provided with periods without color and without fill, a greyscale will be used
  #color: the outline color of the boxes
  #alpha: transparency of the fill colors
  #height: the proportional height of the plot to use for the scale 
  #size: the size of the text in the scale
  #quat: specifies whether the Quaternary should be labelled
  #pos: which side to add the scale to (left, right, top, or bottom)
  #abbrv: whether to use abbreviations instead of full period names
  #periods: a custom data set of time interval boundaries, with the following columns:
  #   period: name of each period (will be used as labels if no abbreviations are provided)
  #   max_age: the oldest boundary of each time interval
  #   min_age: the youngest boundary of each time interval
  #   abbr: (optional) abbreviations that will be used as labels
  #   color: (optional) a hex color code (which can be obtained with rgb()) for each time interval
  #neg: set this to true if your x-axis is actually negative values
  
  require(ggplot2)
  if(is.null(periods)){
    periods <- data.frame(period = c("Quaternary", "Neogene", "Paleogene", "Cretaceous", "Jurassic", "Triassic", "Permian", "Carboniferous", "Devonian", "Silurian", "Ordovician", "Cambrian", "Ediacaran", "Cryogenian", "Tonian"),
                          max_age = c(2.588, 23.03, 66, 145, 201.3, 252.2, 298.9, 358.9, 419.2, 443.4, 485.4, 541, 635, 720, 1000),
                          min_age = c(0, 2.588, 23.03, 66, 145, 201.3, 252.2, 298.9, 358.9, 419.2, 443.4, 485.4, 541, 635, 720),
                          abbr = c("Q", "N", "Pg", "K", "J", "Tr", "P", "C", "D", "S", "O", "Cm","E","Cr","To"),
                          color = c(rgb(249, 249, 127, maxColorValue = 255),rgb(255, 230, 25, maxColorValue = 255),rgb(253, 154, 82, maxColorValue = 255),rgb(127, 198, 78, maxColorValue = 255),rgb(52, 178, 201, maxColorValue = 255),rgb(129, 43, 146, maxColorValue = 255),rgb(240, 64, 40, maxColorValue = 255),rgb(103, 165, 153, maxColorValue = 255),rgb(203, 140, 55, maxColorValue = 255),rgb(179, 225, 182, maxColorValue = 255),rgb(0, 146, 112, maxColorValue = 255),rgb(127, 160, 86, maxColorValue = 255),rgb(254, 217, 106, maxColorValue = 255),rgb(254, 204, 92, maxColorValue = 255),rgb(254, 191, 78, maxColorValue = 255)),
                          stringsAsFactors = FALSE)
  }
  if(neg){
    periods$max_age <- -1 * (periods$max_age)
    periods$min_age <- -1 * (periods$min_age)
  }
  periods$mid_age <- (periods$max_age + periods$min_age)/2
  if(!is.null(fill)){
    periods$color <- rep(fill, length.out = nrow(periods))
  }else if(!("color" %in% colnames(periods))){
    periods$color <- rep(c("grey60","grey80"), length.out = nrow(periods))
  }
  lims <- ggplot_build(gg)$layout$panel_params[[1]]
  if(abbrv & "abbr" %in% colnames(periods)){
    periods$names <- periods$abbr
  }else{
    periods$names <- periods$period
  }
  if(!quat){
    periods$names[periods$abbr=="Q"] <- ""
  }
  if(pos %in% c("bottom", "top", "b", "t")){
    if(pos %in% c("top","t")){
      ymax <- max(lims$y.range)
      ymin <- max(lims$y.range) - height * (max(lims$y.range) - min(lims$y.range))
    }else{
      ymin <- min(lims$y.range)
      ymax <- min(lims$y.range) + height * (max(lims$y.range) - min(lims$y.range))
    }
    gg <- gg +
      annotate("rect", xmin = periods$min_age, xmax = periods$max_age, ymin = ymin, ymax = ymax,
               fill = periods$color, color = color, alpha = alpha) +
      annotate("text", x = periods$mid_age, label = periods$names, y = (ymin+ymax)/2,
               vjust = "middle", hjust = "middle", size = size)
  }else if(pos %in% c("left", "right","l","r")){
    if(pos %in% c("right","r")){
      xmax <- max(lims$x.range)
      xmin <- max(lims$x.range) - height * (max(lims$x.range) - min(lims$x.range))
    }else{
      xmin <- min(lims$x.range)
      xmax <- min(lims$x.range) + height * (max(lims$x.range) - min(lims$x.range))
    }
    gg <- gg +
      annotate("rect", ymin = periods$min_age, ymax = periods$max_age, xmin = xmin, xmax = xmax,
               fill = periods$color, color = color, alpha = alpha) +
      annotate("text", y = periods$mid_age, label = periods$names, x = (xmin+xmax)/2,
               vjust = "middle", hjust = "middle", size = size, angle = 90)
  }
  gg
}
