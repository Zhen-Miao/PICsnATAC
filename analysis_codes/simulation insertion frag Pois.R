library('ggplot2')
library(fitdistrplus)

## data simulation procedure

## parameters
peak_length = 1000 # peak length
# insertion_rate = 0.0005 ## the rate of insertion at any location
# n_cells = 5000 ## number of cells
max_frag_len = 1001
min_frag_len = 0
pad_length = 0 ## pad region length, the insertion can happen between these regions
count_type = 'fragment'
dirs = c(1,-1)
# min_frag_length = 20

dist = rep(list(), length = 10)
# mean_obs = vector(mode='numeric', length = 10)
# var_obs = vector(mode='numeric', length = 10)
sum_sim = matrix(ncol = 8, nrow = 32)
colnames(sum_sim) <- c('n_obs', 'rate', 'mean', 'mean_sd','variance', 'variance_sd','theoretical_mean', 'theoretical_variance')
n_cells = c(5000, 2000, 1000, 500)
insertion_rate = c(0.0005, 0.001, 0.002,0.004,0.005, 0.008, 0.01, 0.02) /2

for(n_cell_iter in 1:4){
  for(rate_iter in 1:8){
    row_id = (n_cell_iter-1)*8+rate_iter
    mean_obs = vector(mode='numeric', length = 10)
    var_obs = vector(mode='numeric', length = 10)
    for(iter in 1:10){
      dist[[iter]] = get_n_frag(n_cells = n_cells[n_cell_iter],
                             peak_length = peak_length,
                             insertion_rate = insertion_rate[rate_iter],
                             pad_length =  pad_length,
                             count_type = 'fragment',
                             min_frag_length = 35,
                             max_frag_length = 600
                              )
      mean_obs[iter] <- mean(dist[[iter]])
      var_obs[iter] <- var(dist[[iter]])
    }
    sum_sim[row_id,'n_obs'] <- n_cells[n_cell_iter]
    sum_sim[row_id,'rate'] <- insertion_rate[rate_iter] * peak_length
    sum_sim[row_id,'mean'] <- mean(mean_obs)
    sum_sim[row_id,'mean_sd'] <- sd(mean_obs)
    sum_sim[row_id,'variance'] <- mean(var_obs)
    sum_sim[row_id,'variance_sd'] <- sd(var_obs)
    sum_sim[row_id, 'theoretical_mean'] <- mean_fun(insertion_rate[rate_iter] * peak_length)
    sum_sim[row_id, 'theoretical_variance'] <- var_fun(insertion_rate[rate_iter] * peak_length)

  }


}


## for ploting the mean variance relationship
var_fun <- function(lam){
  varw = (2 * lam + 2 * exp(-1*lam) - 2 * lam * exp(-1*lam) - exp(-2 * lam) -1) / 4
  return(varw)
}
mean_fun <- function(lam){
  meanw = 0.5 * (lam - 1 + exp(lam * -1))
  return(meanw)
}


get_n_frag <- function(
  n_cells,
  peak_length,
  pad_length,
  insertion_rate,
  dirs = c(-1,1),
  count_type = 'fragment',
  min_frag_length = 0,
  max_frag_length = 1000
){

  frags = vector(mode = 'numeric', length = n_cells)

  for(j in 1:n_cells){
    sim_insertion <- rbinom(n = peak_length + pad_length * 2, size = 1, prob = insertion_rate)
    sim_direction <- sample(dirs,size = peak_length + pad_length * 2,replace = T,prob = c(0.5,0.5)) *
      sim_insertion

    names(sim_insertion) <- names(sim_direction) <- c(1:(peak_length + pad_length * 2))

    ## counting
    if(sum(sim_insertion != 0) <= 1){
      frags[j] <- 0
      next
    }
    sim_event = sim_direction[which(sim_insertion != 0)]
    sim_pos = as.numeric(names(sim_event))
    frag = 0
    if(count_type == 'fragment'){
      for(i in 1:(length(sim_event)-1)){
        if(sim_pos[i+1] < pad_length |
           sim_pos[i] > peak_length + pad_length ){
          next
        }
        if(sim_event[i+1] != sim_event[i] &
           sim_pos[i+1] - sim_pos[i] > min_frag_length &
           sim_pos[i+1] - sim_pos[i] < max_frag_length){
          frag = frag + 1
        }
      }
    }else if(count_type == 'PIC'){
      for(i in 1:(length(sim_event)-1)){
        if(sim_pos[i+1] < pad_length |
           sim_pos[i] > peak_length + pad_length |
            (sim_pos[i] < pad_length & sim_pos[i+1] > peak_length + pad_length )
           ){
          next
        }
        if(sim_event[i+1] != sim_event[i] &
           sim_pos[i+1] - sim_pos[i] > min_frag_length &
           sim_pos[i+1] - sim_pos[i] < max_frag_length){
          frag = frag + 1
        }
      }
    }


    frags[j] <- frag
  }
  return(frags)

}





