#############################################################################################
# Now let's plot the evolutionary rates estimated from models preferred models (>0.5 aiccweight)
#############################################################################################
# Adjust the RPANDA plotting function so we can fix the axes, and do a bunch of plots
plot.fixed_t_env <- function (x, steps = 100, xlim=c(-10,0), ylim=c(0,1), linecol="red", ...) 
{
  if (is.function(x$model)) {
    fun_temp <- function(x, temp, model, param) {
      rate_fun <- function(x) {
        model(x, temp, param)
      }
      rate <- rate_fun(x)
      return(rate)
    }
  }
  else if (x$model == "EnvExp") {
    fun_temp <- function(x, temp, model, param) {
      sig <- param[1]
      beta <- param[2]
      rate <- (sig * exp(beta * temp(x)))
      return(rate)
    }
  }
  else if (x$model == "EnvLin") {
    fun_temp <- function(x, temp, model, param) {
      sig <- param[1]
      beta <- param[2]
      rate <- sig + (beta - sig) * temp(x)
      return(rate)
    }
  }
  t <- seq(0, x$tot_time, length.out = steps)
  rates <- fun_temp(x = t, temp = x$env_func, model = x$model, 
                    param = x$param)
  plot(-t, rates, type = "l", xlab = "Times", ylab = bquote(paste("Evolutionary rates ", 
                                                                  sigma)), xlim=xlim, ylim=ylim, col=linecol, ...)
  results <- list(time_steps = t, rates = rates)
  invisible(results)
}