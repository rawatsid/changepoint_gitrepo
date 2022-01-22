library(dplyr)
library(utils)
library(readr)
library(invgamma)
library(goftest)
library(matlib)
library(geosphere)
library(ape)
library(zoo)
library(matrixcalc)
library(MASS)
library(mnormt)
library(matrixStats)
library(truncnorm)
library(mcmc)
#library(tmvtnorm)
library(sandwich) 
#library(gmm)
library(Matrix)
library(evd)
library(truncdist)
library(stats4)
library(parallel)
library(caret)
library(Rcpp)
library(RcppArmadillo)

df = read.csv("D:/Year_3/Term2/Research/changepoint_analysis/US_COVID_weekly_NY_data.csv", row.names = NULL)
df$week = as.Date(df$week)
set.seed(10)
test_loc = sample(unique(df$id),30)
training_df = df %>% filter(id %in% test_loc)
#test_loc = sample(unique(df$key),2)
# start_training_wk = max(df$week) - 52*7
# test_wk_min = max(df$week) - 0
# validation_wk_min = max(df$week) - 0
# if(is_all_loc){
#   training_df = df %>% filter(id %in% test_loc)
#   validation_df = df %>% filter(week > validation_wk_min & week <= test_wk_min)
# }else{
#   training_df = df %>% filter(! state %in% test_loc & week <= validation_wk_min)
#   validation_df = df %>% filter(! state %in% test_loc & week > validation_wk_min & week <= test_wk_min)
# }


unique_county_lat_long = unique(training_df %>% dplyr::select(key,lat,long))
distance_mat_training = as.matrix(distm(cbind(unique_county_lat_long$long,unique_county_lat_long$lat), 
                                        fun = distHaversine))
#training values

n_total_training = nrow(training_df)
n_sp_training = nrow(unique_county_lat_long)
n_tmp_training = n_total_training/n_sp_training
# #total dataset values
# n_total = nrow(df)
# n_sp_total = length(unique(df$key))
# n_tmp_total = n_total/ n_sp_total
# #validation set values
# y_validation = validation_df$category
# n_total_validation = length(y_validation)
# n_sp_validation = length(unique(validation_df$key))
# n_tmp_validation = n_total_validation/n_sp_validation

#No of categories
no_of_categ = length(unique(training_df$category))
#X matrix for the model
# X_mat_pre_changepoint_training = cbind(rep(1,n_total_training), log(training_df$population), training_df$time, training_df$time_sq, 
#                                        training_df$prev_log_new_death)
#Number of covariates

#Getting week difference mat for Sigma_T
# week_vec = 1:n_tmp_training
# week_diff_mat = as.matrix(dist(week_vec,diag  =  TRUE,upper = TRUE))

#Error Vectors
# mape_vec_training_training = vector()
# mape_vec_validation_training = vector() 
accuracy_pred_vec = vector()
pred_category_mat = vector()

burn_period_training = 5000
sample_size_training = 1000
diff_in_random_draws_training = 10

#Getting Artificial Values############################################################################
training_df_clone = training_df

t0_art = 35
beta_art =  c(1.5, 4, 7, 3, 9) #rep(0,p+1)
beta_star_art = c(3.5, 6, 8, 5, 11) #rep(0,p+1)
sigma1_sq_art = 1
sigma1_sq_star_art = 1
sigma2_sq_art = 1
sigma2_sq_star_art = 1


rho_sp = 0.01

rho_tmp = 0.1

week_diff_mat_pre_t0_art = as.matrix(dist(1:t0_art, 
                                          diag =  TRUE,upper = TRUE))
week_diff_mat_post_t0_art = as.matrix(dist((t0_art+1):n_tmp_training, 
                                           diag =  TRUE,upper = TRUE))
SIGMA_sp_art = exp(-rho_sp * distance_mat_training/(1000)) # dividing by 80000 to get matrix values similar to temporal values
# SIGMA_tmp_pre_changepoint_training = exp(-rho_tmp * week_diff_mat_pre_changepoint)
# SIGMA_tmp_post_changepoint_training = exp(-rho_tmp * week_diff_mat_post_changepoint)
#inv_SIGMA_sp_training = chol2inv(chol(SIGMA_sp_training))
SIGMA_tmp_pre_t0_art = exp(-rho_tmp * week_diff_mat_pre_t0_art)
SIGMA_tmp_post_t0_art = exp(-rho_tmp * week_diff_mat_post_t0_art)

SIGMA_v_pre_t0_art = SIGMA_sp_art %x% SIGMA_tmp_pre_t0_art

SIGMA_v_post_t0_art = SIGMA_sp_art %x% SIGMA_tmp_post_t0_art
set.seed(1)
v_vec_art = sqrt(sigma1_sq_art) * t(chol(SIGMA_v_pre_t0_art)) %*% rnorm(t0_art * n_sp_training)
set.seed(2)
v_vec_star_art =sqrt(sigma1_sq_star_art) * t(chol(SIGMA_v_post_t0_art)) %*% rnorm((n_tmp_training - t0_art) * n_sp_training)

cp_t0_week_art = unique(training_df_clone$week)[t0_art]

training_pre_t0_art_df = training_df_clone %>% filter(week <= cp_t0_week_art)
training_post_t0_art_df = training_df_clone %>% filter(week > cp_t0_week_art)

X_mat_pre_t0_art = cbind(rep(1,n_sp_training * t0_art), log(training_pre_t0_art_df$population), 
                         training_pre_t0_art_df$time, training_pre_t0_art_df$time_sq, 
                         training_pre_t0_art_df$prev_log_new_death)

X_mat_post_t0_art = cbind(rep(1,n_sp_training * (n_tmp_training - t0_art)), 
                          log(training_post_t0_art_df$population), training_post_t0_art_df$time, 
                          training_post_t0_art_df$time_sq, training_post_t0_art_df$prev_log_new_death)

p = ncol(X_mat_pre_t0_art) - 1
set.seed(3)
pi_vec_art = (X_mat_pre_t0_art %*% beta_art + v_vec_art + 
                sqrt(sigma2_sq_art)* rnorm(t0_art * n_sp_training))
set.seed(4)
pi_vec_star_art = (X_mat_post_t0_art %*% beta_star_art + v_vec_star_art + 
                     sqrt(sigma2_sq_star_art)* rnorm((n_tmp_training - t0_art) * n_sp_training))

training_pre_t0_art_df$pi_vec_art = as.numeric(pi_vec_art)

training_post_t0_art_df$pi_vec_art = as.numeric(pi_vec_star_art)

training_df_clone = (rbind(training_pre_t0_art_df, training_post_t0_art_df) %>% 
                       arrange(country, state,key, lat, long, id,  week))

delta_art = c(-Inf, 40, 60, 95, Inf)

training_df_clone = training_df_clone %>% 
  mutate(category_art = ifelse(pi_vec_art > delta_art[4],4,
                               ifelse(pi_vec_art > delta_art[3], 3,
                                      ifelse(pi_vec_art > delta_art[2],2,1))))
y_training = training_df$category_art
#constant parameters values###########################################################################
a_training = 2
lambda_training = 0.1
#Saving sample of posterior values for each rho values
pi_mat_sample_list_training = list()
delta_mat_sample_list_training = list()
v_vec_mat_sample_list_training = list()

beta_sample_list_training = list()
beta_star_sample_list_training = list()
t_0_sample_vec_training = vector()
sigma1_sq_sample_vec_training = vector()
sigma1_sq_star_sample_vec_training = vector()
sigma2_sq_sample_vec_training = vector()
sigma2_sq_star_sample_vec_training = vector()
#y_pred_val_mat_list_training = list()

# rho_vec_tmp_training = c(0.10, 0.50, 1, 3)#c(0.10, 0.50, 0.75, 1, 1.5, 3) 
# rho_vec_sp_training = c(0.01, 0.05, 0.25)#c(0.006, 0.01, 0.05, 0.1, 0.25, 0.5) c(0.006, 0.01, 0.05, 0.25)
rho_idx = 1

gap_bw_changepoints = 5
changepoint_vec = seq(2*gap_bw_changepoints, n_tmp_training - gap_bw_changepoints, 
                      by = gap_bw_changepoints)# - gap_bw_changepoints
#Getting y in matrix form
#y_mat_training = matrix(y_training,nrow = n_tmp_training,ncol = n_sp_training)

#Giving values to different parameters of our model
n_tmp_pre_changepoint_training = changepoint_vec[floor(length(changepoint_vec)/2)]
t_0_training = n_tmp_pre_changepoint_training
n_tmp_post_changepoint_training = n_tmp_training - n_tmp_pre_changepoint_training
changepoint_t_0_week = unique(training_df$week)[n_tmp_pre_changepoint_training]

set.seed(5)
delta_training = runif(no_of_categ -1, -100, 100)
delta_training = sort(delta_training)
delta_training = c(-Inf , delta_training ,Inf)
alpha_vec = rep(1,no_of_categ) # shape and scale parameters of trunc Beta .
stdnorm = 10
#training_df_clone = training_df
# training_df_clone$pi_vec = rep(0,n_total_training)
# pi_mat_pre_changepoint_training = matrix(((training_df_clone %>% 
#                                              filter(week <= changepoint_t_0_week))$pi_vec),
#                                          nrow = n_tmp_pre_changepoint_training, ncol = n_sp_training)
# 
# pi_mat_post_changepoint_training = matrix(((training_df_clone %>% 
#                                               filter(week > changepoint_t_0_week))$pi_vec),
#                                           nrow = (n_tmp_post_changepoint_training), ncol = n_sp_training)

#selecting midpoint as the initial changepoint
training_df_clone$v_vec1 = rep(0, n_total_training)
training_df_clone$v_vec2 = rep(0, n_total_training)
training_df_clone$v_vec3 = rep(0, n_total_training)



beta_training = matrix(rep(0,p+1), ncol = 1)
beta_star_training = matrix(rep(0,p+1), ncol = 1)
sigma1_sq_training = 1
sigma1_sq_star_training = 1
sigma2_sq_training = 1
sigma2_sq_star_training = 1


#Setting matrices
# week_diff_mat_pre_changepoint = as.matrix(dist(1:n_tmp_pre_changepoint_training, 
#                                                diag =  TRUE,upper = TRUE))
# week_diff_mat_post_changepoint = as.matrix(dist((n_tmp_pre_changepoint_training+1):n_tmp_training, 
#                                                 diag =  TRUE,upper = TRUE))
SIGMA_sp_training = exp(-rho_sp * distance_mat_training/(1000)) # dividing by 80000 to get matrix values similar to temporal values

inv_SIGMA_sp_training = chol2inv(chol(SIGMA_sp_training))

#time taken analysis for parameters
total_tt_v = 0
total_tt_beta = 0
total_tt_sigma1_sq = 0
total_tt_sigma2_sq = 0
total_tt_t_0 = 0
total_tt_pi = 0
total_tt_delta = 0

#C++ code##########################################################################################
arma_code <- 
  "arma::vec arma_mm(const arma::mat& m, const arma::vec& v, const arma::rowvec& tv) {
       return tv * m * v;
   };"

arma_mm = cppFunction(code = arma_code, depends = "RcppArmadillo")
#Setting values for markov chains before convergence########################################################

beta_mat_preconvergence_training = matrix(c(seq(1,5,1),
                                            seq(1,5,1) + rnorm(p+1, sd =0.4),
                                            seq(1,5,1) + rnorm(p+1, sd =0.4)), 
                                          ncol = 3, nrow = p+1)

beta_star_mat_preconvergence_training = matrix(c(seq(1,5,1),
                                                 seq(1,5,1)+ rnorm(p+1, sd =0.6),
                                                 seq(1,5,1) + rnorm(p+1, sd =0.6)), 
                                               ncol = 3, nrow = p+1)
sigma1_sq_preconvergence_training = c(1, 1, 1)

sigma1_sq_star_preconvergence_training = c(1, 1, 1)

sigma2_sq_preconvergence_training = c(1, 1, 1)

sigma2_sq_star_preconvergence_training = c(1, 1, 1)

t0_preconvergence_training = c(changepoint_vec[floor(length(changepoint_vec)/2)], 
                               changepoint_vec[floor(length(changepoint_vec)/2) + 1],
                               changepoint_vec[floor(length(changepoint_vec)/2) - 1])

changepoint_t0_week_preconvergence_training = unique(training_df$week)[t0_preconvergence_training]

delta_mat_preconvergence_training = cbind(sort(runif(no_of_categ -1, -100, 100)), sort(runif(no_of_categ -1, -100, 100)),
                                          sort(runif(no_of_categ -1, -100, 100)))

delta_mat_preconvergence_training = rbind(c(-Inf, -Inf, -Inf),delta_mat_preconvergence_training,
                                          c(Inf, Inf, Inf))

############################################################################################################

#Getting inv Sigma_sp_22 to use for every iteration
inv_SIGMA_sp_trans_22_matlist_training = list()
# inv_SIGMA_cond_pre_changepoint_matlist_training = list()
# inv_SIGMA_cond_post_changepoint_matlist_training = list()
SIGMA_sp_trans_21_mat_training = vector()
for (sp in (1:n_sp_training)) {
  idx = 1:n_sp_training
  
  if(sp != 1){
    idx[1]=sp
    idx[sp]=1
  }
  SIGMA_sp_trans_training = (diag(n_sp_training)[idx,] %*% SIGMA_sp_training %*% 
                               diag(n_sp_training)[idx,])
  inv_SIGMA_sp_trans_training = (diag(n_sp_training)[idx,] %*% inv_SIGMA_sp_training %*% 
                                   diag(n_sp_training)[idx,])
  
  inv_SIGMA_sp_trans_a_training = inv_SIGMA_sp_trans_training[1,1]
  inv_SIGMA_sp_trans_b_training = inv_SIGMA_sp_trans_training[1,-1]
  inv_SIGMA_sp_trans_c_training = inv_SIGMA_sp_trans_training[-1,1]
  inv_SIGMA_sp_trans_d_training = inv_SIGMA_sp_trans_training[-1,-1]
  
  inv_SIGMA_sp_trans_22_training = (inv_SIGMA_sp_trans_d_training - (inv_SIGMA_sp_trans_c_training %*%
                                                                       t(inv_SIGMA_sp_trans_b_training))/
                                      inv_SIGMA_sp_trans_a_training)
  
  SIGMA_sp_trans_21_training = SIGMA_sp_trans_training[-1,1]
  
  # SIGMA_cond_pre_changepoint_training = (1 - t(SIGMA_sp_trans_21_training) %*% inv_SIGMA_sp_trans_22_training %*% 
  #                                          SIGMA_sp_trans_21_training) %x% SIGMA_tmp_pre_changepoint_training
  # 
  # SIGMA_cond_post_changepoint_training = (1 - t(SIGMA_sp_trans_21_training) %*% inv_SIGMA_sp_trans_22_training %*%
  #                                           SIGMA_sp_trans_21_training) %x% SIGMA_tmp_post_changepoint_training
  # 
  # inv_SIGMA_cond_pre_changepoint_matlist_training[[sp]] = chol2inv(chol(SIGMA_cond_pre_changepoint_training))
  # 
  # inv_SIGMA_cond_post_changepoint_matlist_training[[sp]] = chol2inv(chol(SIGMA_cond_post_changepoint_training))
  
  inv_SIGMA_sp_trans_22_matlist_training[[sp]] = inv_SIGMA_sp_trans_22_training
  SIGMA_sp_trans_21_mat_training = cbind(SIGMA_sp_trans_21_mat_training, SIGMA_sp_trans_21_training)
  
}

#Saving  SIGMA_v matrix for different t0
inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list = list()
inv_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_mat_list = list()
log_det_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_vec = vector()
log_det_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_vec = vector()
inv_SIGMA_tmp_pre_t_mat_list = list()
inv_SIGMA_tmp_post_t_mat_list = list()

for (tp in (changepoint_vec)){
  week_diff_mat_pre_t = as.matrix(dist(1:tp, diag =  TRUE,upper = TRUE))
  week_diff_mat_post_t = as.matrix(dist((tp+1):n_tmp_training, diag =  TRUE,upper = TRUE))
  SIGMA_sp_training = exp(-rho_sp * distance_mat_training/(1000)) # dividing by 80000 to get matrix values similar to temporal values
  SIGMA_tmp_pre_t = exp(-rho_tmp * week_diff_mat_pre_t)
  SIGMA_tmp_post_t = exp(-rho_tmp * week_diff_mat_post_t)
  
  det_log_obj_SIGMA_sp_training = determinant(SIGMA_sp_training, logarithm = TRUE)
  det_log_SIGMA_sp_training = det_log_obj_SIGMA_sp_training$sign * det_log_obj_SIGMA_sp_training$modulus[1]
  
  det_log_obj_SIGMA_tmp_pre_t = determinant(SIGMA_tmp_pre_t, logarithm = TRUE)
  det_log_SIGMA_tmp_pre_t = det_log_obj_SIGMA_tmp_pre_t$sign * det_log_obj_SIGMA_tmp_pre_t$modulus[1]
  
  det_log_obj_SIGMA_tmp_post_t = determinant(SIGMA_tmp_post_t, logarithm = TRUE)
  det_log_SIGMA_tmp_post_t = det_log_obj_SIGMA_tmp_post_t$sign * det_log_obj_SIGMA_tmp_post_t$modulus[1]
  
  inv_SIGMA_tmp_pre_t = chol2inv(chol(SIGMA_tmp_pre_t))
  inv_SIGMA_tmp_post_t = chol2inv(chol(SIGMA_tmp_post_t))
  
  inv_SIGMA_tmp_pre_t_mat_list[[tp/gap_bw_changepoints - 1]] = inv_SIGMA_tmp_pre_t
  inv_SIGMA_tmp_post_t_mat_list[[tp/gap_bw_changepoints - 1]] = inv_SIGMA_tmp_post_t
  
  inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t = (inv_SIGMA_sp_training %x% inv_SIGMA_tmp_pre_t)
  inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list[[tp/gap_bw_changepoints - 1]] = inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t
  log_det_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t = (tp * det_log_SIGMA_sp_training + 
                                                n_sp_training * det_log_SIGMA_tmp_pre_t)
  log_det_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_vec = c(log_det_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_vec,
                                                   log_det_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t)
  
  inv_SIGMA_sp_kro_inv_SIGMA_tmp_post_t = (inv_SIGMA_sp_training %x% inv_SIGMA_tmp_post_t)
  inv_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_mat_list[[tp/gap_bw_changepoints - 1]] = inv_SIGMA_sp_kro_inv_SIGMA_tmp_post_t
  log_det_SIGMA_sp_kro_inv_SIGMA_tmp_post_t = ((n_tmp_training - tp) * det_log_SIGMA_sp_training + 
                                                 n_sp_training * det_log_SIGMA_tmp_post_t)
  log_det_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_vec = c(log_det_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_vec,
                                                    log_det_SIGMA_sp_kro_inv_SIGMA_tmp_post_t)
}


pi_vec_sample_training = vector()
delta_sample_training = vector()
beta_sample_training = vector()
beta_star_sample_training = vector()
sigma1_sq_sample_training = vector()
sigma1_sq_star_sample_training = vector()
sigma2_sq_sample_training = vector()
sigma2_sq_star_sample_training = vector()
v_vec_sample_training = vector()
t_0_sample_training = vector()
#Initialtzing for different markov chain samples#########################################
beta_sample_chains_list_training = list()
beta_star_sample_chains_list_training = list()
sigma1_sq_sample_chains_list_training = list()
sigma1_sq_star_sample_chains_list_training = list()
sigma2_sq_sample_chains_list_training = list()
sigma2_sq_star_sample_chains_list_training = list()
v_vec_sample_chains_list_training = list()
pi_vec_sample_chains_list_training = list()
t_0_sample_chains_list_training = list()
for (init in 1:3) {
  beta_sample_chains_list_training[[init]] = vector()
  beta_star_sample_chains_list_training[[init]] = vector()
  sigma1_sq_sample_chains_list_training[[init]] = vector()
  sigma1_sq_star_sample_chains_list_training[[init]] = vector()
  sigma2_sq_sample_chains_list_training[[init]] = vector()
  sigma2_sq_star_sample_chains_list_training[[init]] = vector()
  v_vec_sample_chains_list_training[[init]] = vector()
  pi_vec_sample_chains_list_training[[init]] = vector()
  t_0_sample_chains_list_training[[init]] = vector()
  # v_mat_preconvergence_training[[init]] = matrix(rep(0,n_training),nrow = n_tmp_training,
  #                                                ncol = n_sp_training)
}

#Iteration Starts here
converged = FALSE
converged_i = -1
while (i > 0){
  
  for (mc in (1:3)) {
    v_vec_mc_col = paste("v_vec",as.character(mc), sep = "")
    pi_vec_mc_col = paste("pi_vec",as.character(mc), sep = "")
    #training Df pre and post change point
    # training_pre_changepoint_df = training_df_clone %>% filter(week <= changepoint_t_0_week)
    # training_post_changepoint_df = training_df_clone %>% filter(week > changepoint_t_0_week)
    
    #pi posterior
    st_pi = Sys.time()
    pi_vec = vector()
    for (st in 1:n_total_training) {
      st_week = training_df_clone$week[st]
      y_val_st = training_df_clone$category_art[st]
      v_val_st = training_df_clone[[v_vec_mc_col]][st]
      x_vec_st = c(1, log(training_df_clone$population[st]), training_df_clone$time[st], 
                   training_df_clone$time_sq[st], training_df_clone$prev_log_new_death[st])
      if(st_week <= changepoint_t_0_week){
        pi_st = truncnorm::rtruncnorm(1, a = delta_mat_preconvergence_training[,mc][y_val_st], 
                                      b = delta_mat_preconvergence_training[,mc][y_val_st + 1],
                                      mean = t(x_vec_st) %*% beta_mat_preconvergence_training[,mc] + v_val_st, 
                                      sd = sqrt(sigma2_sq_preconvergence_training[mc]))
      }else{
        pi_st = truncnorm::rtruncnorm(1, a = delta_mat_preconvergence_training[,mc][y_val_st], 
                                      b = delta_mat_preconvergence_training[,mc][y_val_st + 1],
                                      mean = t(x_vec_st) %*% beta_star_mat_preconvergence_training[,mc] + v_val_st, 
                                      sd = sqrt(sigma2_sq_preconvergence_training[mc]))
      }
      
      pi_vec = c(pi_vec, pi_st)
    }
    training_df_clone$pi_vec = NULL
    training_df_clone[[pi_vec_mc_col]] = pi_vec
    
    pi_vec_training = pi_vec
    
    training_pre_changepoint_df = training_df_clone %>% filter(week <= changepoint_t0_week_preconvergence_training[mc])
    training_post_changepoint_df = training_df_clone %>% filter(week > changepoint_t0_week_preconvergence_training[mc])
    
    #Getting X matrix for pre and post change point
    X_mat_pre_changepoint_training = cbind(rep(1,n_sp_training * t0_preconvergence_training[mc]), 
                                           log(training_pre_changepoint_df$population), 
                                           training_pre_changepoint_df$time, 
                                           training_pre_changepoint_df$time_sq, 
                                           training_pre_changepoint_df$prev_log_new_death)
    
    X_mat_post_changepoint_training = cbind(rep(1,n_sp_training * (n_tmp_training - t0_preconvergence_training[mc])), 
                                            log(training_post_changepoint_df$population), 
                                            training_post_changepoint_df$time, 
                                            training_post_changepoint_df$time_sq, 
                                            training_post_changepoint_df$prev_log_new_death)
    
    p = ncol(X_mat_pre_changepoint_training) - 1
    
    total_tt_pi = total_tt_pi + as.numeric(difftime(Sys.time(), st_pi), units="secs")
    ###############################################################################################################
    #v posterior#############################################################
    st_v = Sys.time()
    #Xbeta matrix
    v_mat_training = matrix((training_pre_changepoint_df[[v_vec_mc_col]]),
                             ncol = n_sp_training)
    
    v_star_mat_training = matrix((training_post_changepoint_df[[v_vec_mc_col]]),
                                  ncol = n_sp_training)
    
    pi_mat_pre_changepoint_training = matrix((training_pre_changepoint_df[[pi_vec_mc_col]]),
                                             ncol = n_sp_training)
    
    pi_mat_post_changepoint_training = matrix((training_post_changepoint_df[[pi_vec_mc_col]]),
                                              ncol = n_sp_training)
    
    
    xbeta_pre_changepoint_mat_training = matrix((X_mat_pre_changepoint_training %*% 
                                                   beta_mat_preconvergence_training[,mc]),
                                                ncol = n_sp_training)
    
    xbeta_post_changepoint_mat_training = matrix((X_mat_post_changepoint_training %*% 
                                                    beta_star_mat_preconvergence_training[,mc]),
                                                 ncol = n_sp_training)
    
    week_diff_mat_pre_changepoint = as.matrix(dist(1:t0_preconvergence_training[mc], 
                                                   diag =  TRUE,upper = TRUE))
    week_diff_mat_post_changepoint = as.matrix(dist((t0_preconvergence_training[mc]+1):n_tmp_training, 
                                                    diag =  TRUE,upper = TRUE))
    
    # SIGMA_tmp_pre_changepoint_training = exp(-rho_tmp * week_diff_mat_pre_changepoint)
    # SIGMA_tmp_post_changepoint_training = exp(-rho_tmp * week_diff_mat_post_changepoint)
    # inv_SIGMA_tmp_pre_changepoint_training = chol2inv(chol(SIGMA_tmp_pre_changepoint_training))
    # inv_SIGMA_tmp_post_changepoint_training = chol2inv(chol(SIGMA_tmp_post_changepoint_training))
    
    # inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_changepoint_training = (inv_SIGMA_sp_training %x%
    #                                                              inv_SIGMA_tmp_pre_changepoint_training)
    # inv_SIGMA_sp_kro_inv_SIGMA_tmp_post_changepoint_training = (inv_SIGMA_sp_training %x%
    #                                                               inv_SIGMA_tmp_post_changepoint_training)
    
    for (sp in (1:n_sp_training)){
      idx = 1:n_sp_training
      
      if(sp!=1){
        idx[1] = sp
        idx[sp] = 1
      }
      #Switching first and jth index column for V matrix and then taking values from the 2nd column
      v_cond_training = (v_mat_training %*% diag(n_sp_training)[idx,])[(1+t0_preconvergence_training[mc]):
                                                                         (t0_preconvergence_training[mc] * 
                                                                            n_sp_training)]
      
      v_cond_star_training = (v_star_mat_training %*% diag(n_sp_training)[idx,])[(1+(n_tmp_training - 
                                                                                       t0_preconvergence_training[mc])):
                                                                                   ((n_tmp_training - 
                                                                                       t0_preconvergence_training[mc]) * 
                                                                                      n_sp_training)]
      #SIGMA' conditional without sigma1_sq as it's included in the posterior 
      # SIGMA_cond_training = (1 -t(SIGMA_sp_trans_21_mat_training[,j])%*%inv_SIGMA_sp_trans_22_matlist_training[[j]] %*%
      #                          SIGMA_sp_trans_21_mat_training[,j])%x% SIGMA_tmp_training
      mu_cond_training = ((t(SIGMA_sp_trans_21_mat_training[,sp]) %*% 
                             inv_SIGMA_sp_trans_22_matlist_training[[sp]]) %x% 
                            diag(t0_preconvergence_training[mc])) %*% v_cond_training
      
      mu_cond_star_training = ((t(SIGMA_sp_trans_21_mat_training[,sp]) %*% 
                                  inv_SIGMA_sp_trans_22_matlist_training[[sp]]) %x% 
                                 diag((n_tmp_training - 
                                        t0_preconvergence_training[mc]))) %*% v_cond_star_training
      
      
      #Can be made faster here by saving these values
      mat_mult_temp = as.numeric(1 - t(SIGMA_sp_trans_21_mat_training[,sp]) %*% 
                                   inv_SIGMA_sp_trans_22_matlist_training[[sp]] %*% 
                                   SIGMA_sp_trans_21_mat_training[,sp])
      
      inv_SIGMA_cond_pre_changepoint_training = (inv_SIGMA_tmp_pre_t_mat_list[[t0_preconvergence_training[mc]/gap_bw_changepoints - 1]]/
                                                   mat_mult_temp)
      
      inv_SIGMA_cond_post_changepoint_training = (inv_SIGMA_tmp_post_t_mat_list[[t0_preconvergence_training[mc]/gap_bw_changepoints - 1]]/
                                                    mat_mult_temp)
      
      v_covar_training = chol2inv(chol(inv_SIGMA_cond_pre_changepoint_training/sigma1_sq_preconvergence_training[mc] + 
                                         diag(n_tmp_pre_changepoint_training)/sigma2_sq_preconvergence_training[mc]))
      
      v_star_covar_training = chol2inv(chol(inv_SIGMA_cond_post_changepoint_training/sigma1_sq_star_preconvergence_training[mc] + 
                                              diag(n_tmp_post_changepoint_training)/sigma2_sq_star_preconvergence_training[mc]))
      
      v_mean_training  = v_covar_training %*% (inv_SIGMA_cond_pre_changepoint_training %*% 
                                                 mu_cond_training/sigma1_sq_preconvergence_training[mc] + 
                                                 (pi_mat_pre_changepoint_training[,sp] - 
                                                    xbeta_pre_changepoint_mat_training[,sp])/
                                                 sigma2_sq_preconvergence_training[mc])
      
      v_star_mean_training  = v_star_covar_training %*% (inv_SIGMA_cond_post_changepoint_training %*% 
                                                           mu_cond_star_training/sigma1_sq_star_preconvergence_training[mc] + 
                                                           (pi_mat_post_changepoint_training[,sp] - 
                                                              xbeta_post_changepoint_mat_training[,sp])/
                                                           sigma2_sq_star_preconvergence_training[mc] )
      
      v_mat_training[,sp] = v_mean_training + t(chol(v_covar_training))  %*% rnorm(t0_preconvergence_training[mc])
      v_star_mat_training[,sp] = (v_star_mean_training + t(chol(v_star_covar_training))  %*% 
                                    rnorm((n_tmp_training - t0_preconvergence_training[mc])))
      
    }
    training_pre_changepoint_df[[v_vec_mc_col]] = v_mat_training[1:(n_sp_training * t0_preconvergence_training[mc])]
    training_post_changepoint_df[[v_vec_mc_col]] = v_star_mat_training[1:(n_sp_training * (n_tmp_training - 
                                                                                   t0_preconvergence_training[mc]))]
    
    training_df_clone = rbind(training_pre_changepoint_df, training_post_changepoint_df) %>% 
      arrange(country, state,key, lat, long, id,  week)
    
    v_vec_training = training_df_clone[[v_vec_mc_col]]
    
    total_tt_v = total_tt_v + as.numeric(difftime(Sys.time(), st_v), units="secs")
    ########################################################################################################################
    #beta posterior###################################
    st_beta = Sys.time()
    pi_pre_changepoint_training = training_pre_changepoint_df[[pi_vec_mc_col]]#pi_mat_pre_changepoint_training[1:(n_sp_training * n_tmp_pre_changepoint_training)]
    pi_post_changepoint_training = training_post_changepoint_df[[pi_vec_mc_col]]#pi_mat_pre_changepoint_training[1:(n_sp_training * n_tmp_post_changepoint_training)]
    
    
    v_training = v_mat_training[1:(n_sp_training * t0_preconvergence_training[mc])]
    v_star_training = v_star_mat_training[1:(n_sp_training * (n_tmp_training - t0_preconvergence_training[mc]))]
    #dividing by 10^4 for large variance for beta prior
    beta_covar_training = solve( (t(X_mat_pre_changepoint_training) %*% 
                                    X_mat_pre_changepoint_training)/sigma2_sq_preconvergence_training[mc]) #+ diag(p+1)/10^4
    
    beta_star_covar_training = solve( (t(X_mat_post_changepoint_training) %*% 
                                         X_mat_post_changepoint_training)/sigma2_sq_star_preconvergence_training[mc])
    
    beta_mean_training = beta_covar_training %*% ( t(X_mat_pre_changepoint_training)%*% 
                                                     (pi_pre_changepoint_training - v_training))/sigma2_sq_preconvergence_training[mc]
    
    beta_star_mean_training = (beta_star_covar_training %*% ( t(X_mat_post_changepoint_training)%*% 
                                                                (pi_post_changepoint_training - v_star_training))/
                                 sigma2_sq_star_preconvergence_training[mc])
    
    beta_mat_preconvergence_training[,mc] = beta_mean_training + t(chol(beta_covar_training)) %*% rnorm(p+1)
    
    beta_star_mat_preconvergence_training[,mc] = beta_star_mean_training + t(chol(beta_star_covar_training)) %*% rnorm(p+1)
    
    total_tt_beta = total_tt_beta + as.numeric(difftime( Sys.time(), st_beta), units="secs")
    #######################################################################################################################
    #Sigma1^2 posterior######################################################
    st_sigma1_sq = Sys.time()
    
    # sigma1_sq_a_training = a_training + (n_sp_training * n_tmp_pre_changepoint_training)/2
    # 
    # sigma1_sq_star_a_training = a_training + (n_sp_training * n_tmp_post_changepoint_training)/2
    # 
    # sigma1_sq_lambda_training = (t(v_training) %*% 
    #                                inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list[[t_0_training/gap_bw_changepoints - 1]] 
    #                              %*% v_training)/2 + lambda_training
    # 
    # sigma1_sq_star_lambda_training = (t(v_star_training) %*% 
    #                                     inv_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_mat_list[[t_0_training/gap_bw_changepoints - 1]] 
    #                                   %*% v_star_training)/2 + lambda_training
    
    # sigma1_sq_training = 1#invgamma::rinvgamma(1,sigma1_sq_a_training,rate = sigma1_sq_lambda_training)
    # 
    # sigma1_sq_star_training = 1 #invgamma::rinvgamma(1, sigma1_sq_star_a_training,
    #rate = sigma1_sq_star_lambda_training)
    
    total_tt_sigma1_sq = total_tt_sigma1_sq + as.numeric(difftime(Sys.time(), st_sigma1_sq), units="secs")
    ##########################################################################################################################
    #Sigma2^2 posterior############################################
    st_sigma2_sq = Sys.time()
    
    # sigma2_sq_a_training = a_training + (n_sp_training * n_tmp_pre_changepoint_training)/2
    # 
    # sigma2_sq_star_a_training = a_training + (n_sp_training * n_tmp_post_changepoint_training)/2
    # 
    # temp_vec = pi_pre_changepoint_training - X_mat_pre_changepoint_training%*%beta_training - v_training
    # 
    # temp_vec2 = (pi_post_changepoint_training - X_mat_post_changepoint_training %*% beta_star_training - 
    #                v_star_training)
    # 
    # sigma2_sq_lambda_training = sum(temp_vec^2)/2 + lambda_training
    # 
    # sigma2_sq_star_lambda_training = sum(temp_vec2^2)/2 + lambda_training
    
    # sigma2_sq_training = 1#invgamma::rinvgamma(1,sigma2_sq_a_training,rate = sigma2_sq_lambda_training)
    # 
    # sigma2_sq_star_training = 1#invgamma::rinvgamma(1, sigma2_sq_star_a_training,
    #rate = sigma2_sq_star_lambda_training)
    
    total_tt_sigma2_sq = total_tt_sigma2_sq + as.numeric(difftime(Sys.time(), st_sigma2_sq), units="secs")
    ###################################################################################################################
    #delta posterior###############################
    st_delta = Sys.time()
    for (j in 1: (no_of_categ -1)){
      pi_cat1 = pi_vec[y_training == j] #because j starts from 1 hence category values are okay but delta vec values are shifted
      c1 = max(pi_cat1) 
      pi_cat2 = pi_vec[y_training == (j + 1)] 
      c2 = min(pi_cat2)
      a = pnorm( delta_mat_preconvergence_training[,mc][j], mean =0, stdnorm )
      b = pnorm( delta_mat_preconvergence_training[,mc][j+2] , mean =0, stdnorm )
      w1 = (pnorm(c1 ,0, stdnorm)-a)/(b-a)
      w2 = (pnorm(c2 ,0, stdnorm)-a)/(b-a)
      w = rtrunc (1, spec = "beta",shape1 = alpha_vec[j],
                  shape2 = alpha_vec[j+1] ,a=w1 ,b=w2)
      
      delta_mat_preconvergence_training[,mc][j +1] = qnorm ((b-a)*w+a, mean =0, stdnorm )
    }
    total_tt_delta = total_tt_delta + as.numeric(difftime(Sys.time(), st_delta), units="secs")
    ################################################################################################################
    #t_0 posterior####################################
    st_t_0 = Sys.time()
    log_prop_pmf_t_0_vec_training = vector()
    for (tp in (changepoint_vec)) {
      
      week_t = unique(training_df_clone$week)[tp]
      
      X_mat_pre_t = cbind(rep(1,n_sp_training * tp), 
                          log((training_df_clone %>% filter(week <= week_t))$population), 
                          (training_df_clone %>% filter(week <= week_t))$time, 
                          (training_df_clone %>% filter(week <= week_t))$time_sq, 
                          (training_df_clone %>% filter(week <= week_t))$prev_log_new_death)
      
      X_mat_post_t = cbind(rep(1,n_sp_training * (n_tmp_training - tp)), 
                           log((training_df_clone %>% filter(week > week_t))$population),
                           (training_df_clone %>% filter(week > week_t))$time, 
                           (training_df_clone %>% filter(week > week_t))$time_sq, 
                           (training_df_clone %>% filter(week > week_t))$prev_log_new_death)
      
      training_df_pre_t = training_df_clone %>% filter(week <= week_t)
      training_df_post_t = training_df_clone %>% filter(week > week_t)
      
      
      v_vec_pre_t = training_df_pre_t[[v_vec_mc_col]]
      v_vec_post_t = training_df_post_t[[v_vec_mc_col]]
      pi_vec_pre_t = training_df_pre_t[[pi_vec_mc_col]]
      pi_vec_post_t = training_df_post_t[[pi_vec_mc_col]]
      
      temp_vec = pi_vec_pre_t - X_mat_pre_t %*% beta_mat_preconvergence_training[,mc] - v_vec_pre_t
      temp_vec2 = pi_vec_post_t - X_mat_post_t %*% beta_star_mat_preconvergence_training[,mc] - v_vec_post_t
      
      log_prop_pmf_t_0_training = ((-n_sp_training*tp/2)*log(2*pi*sigma2_sq_training) + 
                                     ((-1/(2*sigma2_sq_training)) * sum(temp_vec^2)) + 
                                     ((-n_sp_training*(n_tmp_training - tp)/2) * 
                                        log(2*pi*sigma2_sq_star_training)) +
                                     ((-1/(2*sigma2_sq_star_training)) * sum(temp_vec2^2)) + 
                                     ((-n_sp_training*tp/2) * log(2*pi*sigma1_sq_training)) + 
                                     ((-1/2) * 
                                        log_det_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_vec[tp/gap_bw_changepoints - 1]) + 
                                     ((-1/(2 * sigma1_sq_training)) * 
                                        (arma_mm(inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list[[tp/gap_bw_changepoints - 1]],
                                                 v_vec_pre_t, v_vec_pre_t))) + 
                                     ((-n_sp_training*(n_tmp_training - tp)/2) * 
                                        log(2*pi*sigma1_sq_star_training)) + 
                                     (-1/2) * log_det_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_vec[tp/gap_bw_changepoints - 1] + 
                                     ((-1/(2*sigma1_sq_star_training)) * 
                                        (arma_mm(inv_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_mat_list[[tp/gap_bw_changepoints - 1]],
                                                 v_vec_post_t, v_vec_post_t))))
      
      log_prop_pmf_t_0_vec_training = c(log_prop_pmf_t_0_vec_training, log_prop_pmf_t_0_training)
    }
    #min_prop_pmf_t_0_vec_training = min(prop_pmf_t_0_vec_training)
    max_log_prop_pmf_t_0_vec_training = max(log_prop_pmf_t_0_vec_training)
    # if(min_prop_pmf_t_0_vec_training < 0){
    #   prop_pmf_t_0_vec_training = prop_pmf_t_0_vec_training + (-min_prop_pmf_t_0_vec_training)
    # }
    log_prop_pmf_t_0_vec_training = log_prop_pmf_t_0_vec_training - max_log_prop_pmf_t_0_vec_training
    prop_pmf_t_0_vec_training = exp(log_prop_pmf_t_0_vec_training)
    
    const_prop = 1/sum(prop_pmf_t_0_vec_training)
    
    t0_preconvergence_training[mc] = sample(x = (changepoint_vec),1, replace = FALSE, 
                                            prob = round(prop_pmf_t_0_vec_training * const_prop, digits = 3))
    
    # n_tmp_pre_changepoint_training = t0_preconvergence_training
    # n_tmp_post_changepoint_training = n_tmp_training - n_tmp_pre_changepoint_training
    changepoint_t0_week_preconvergence_training = unique(training_df$week)[t0_preconvergence_training]
    
    total_tt_t_0 = total_tt_t_0 + as.numeric(difftime(Sys.time(), st_t_0), units="secs")
    
    #Collecting samples for prediction
    #Collecting samples for Gelman Rubin stat#####################################################################
    
    v_vec_sample_chains_list_training[[mc]] = cbind(v_vec_sample_chains_list_training[[mc]], 
                                                    v_vec_training)
    
    pi_vec_sample_chains_list_training[[mc]] = cbind(pi_vec_sample_chains_list_training[[mc]],
                                                     pi_vec_training)
    
    beta_sample_chains_list_training[[mc]] = cbind(beta_sample_chains_list_training[[mc]], 
                                                   beta_mat_preconvergence_training[,mc])
    
    beta_star_sample_chains_list_training[[mc]] = cbind(beta_star_sample_chains_list_training[[mc]],
                                                        beta_star_mat_preconvergence_training[,mc])
    
    sigma1_sq_sample_chains_list_training[[mc]] = c(sigma1_sq_sample_chains_list_training[[mc]],
                                                    sigma1_sq_preconvergence_training[mc])
    
    sigma1_sq_star_sample_chains_list_training[[mc]] = c(sigma1_sq_star_sample_chains_list_training[[mc]], 
                                                         sigma1_sq_star_preconvergence_training[mc])
    
    sigma2_sq_sample_chains_list_training[[mc]] = c(sigma2_sq_sample_chains_list_training[[mc]],
                                                    sigma2_sq_preconvergence_training[mc])
    
    sigma2_sq_star_sample_chains_list_training[[mc]] = c(sigma2_sq_star_sample_chains_list_training[[mc]],
                                                         sigma2_sq_star_preconvergence_training[mc])
    
    t0_sample_chains_list_training[[mc]] = c(t0_sample_chains_list_training[[mc]], 
                                                 beta_L_preconvergence_training[mc])
    
    
  }
  
  #Post convergence################################################################################
  #training Df pre and post change point
  # training_pre_changepoint_df = training_df_clone %>% filter(week <= changepoint_t_0_week)
  # training_post_changepoint_df = training_df_clone %>% filter(week > changepoint_t_0_week)
  
  #pi posterior
  st_pi = Sys.time()
  pi_vec = vector()
  for (st in 1:n_total_training) {
    st_week = training_df_clone$week[st]
    y_val_st = training_df_clone$category_art[st]
    v_val_st = training_df_clone$v_vec[st]
    x_vec_st = c(1, log(training_df_clone$population[st]), training_df_clone$time[st], 
                 training_df_clone$time_sq[st], training_df_clone$prev_log_new_death[st])
    if(st_week <= changepoint_t_0_week){
      pi_st = truncnorm::rtruncnorm(1, a = delta_training[y_val_st], b = delta_training[y_val_st + 1],
                                    mean = t(x_vec_st) %*% beta_training + v_val_st, 
                                    sd = sqrt(sigma2_sq_training))
    }else{
      pi_st = truncnorm::rtruncnorm(1, a = delta_training[y_val_st], b = delta_training[y_val_st + 1],
                                    mean = t(x_vec_st) %*% beta_star_training + v_val_st, 
                                    sd = sqrt(sigma2_sq_star_training))
    }
    
    pi_vec = c(pi_vec, pi_st)
  }
  training_df_clone$pi_vec = NULL
  training_df_clone$pi_vec = pi_vec
  
  pi_vec_training = pi_vec
  
  training_pre_changepoint_df = training_df_clone %>% filter(week <= changepoint_t_0_week)
  training_post_changepoint_df = training_df_clone %>% filter(week > changepoint_t_0_week)
  
  #Getting X matrix for pre and post change point
  X_mat_pre_changepoint_training = cbind(rep(1,n_sp_training * n_tmp_pre_changepoint_training), 
                                         log(training_pre_changepoint_df$population), 
                                         training_pre_changepoint_df$time, 
                                         training_pre_changepoint_df$time_sq, 
                                         training_pre_changepoint_df$prev_log_new_death)
  
  X_mat_post_changepoint_training = cbind(rep(1,n_sp_training * n_tmp_post_changepoint_training), 
                                          log(training_post_changepoint_df$population), 
                                          training_post_changepoint_df$time, 
                                          training_post_changepoint_df$time_sq, 
                                          training_post_changepoint_df$prev_log_new_death)
  
  p = ncol(X_mat_pre_changepoint_training) - 1
  
  total_tt_pi = total_tt_pi + as.numeric(difftime(Sys.time(), st_pi), units="secs")
  ###############################################################################################################
  #v posterior#############################################################
  st_v = Sys.time()
  #Xbeta matrix
  v_mat_training = matrix((training_pre_changepoint_df$v_vec),
                          nrow = n_tmp_pre_changepoint_training, ncol = n_sp_training)
  
  v_star_mat_training = matrix((training_post_changepoint_df$v_vec),
                               nrow = (n_tmp_post_changepoint_training), ncol = n_sp_training)
  
  pi_mat_pre_changepoint_training = matrix((training_pre_changepoint_df$pi_vec),
                                           nrow = n_tmp_pre_changepoint_training, ncol = n_sp_training)
  
  pi_mat_post_changepoint_training = matrix((training_post_changepoint_df$pi_vec),
                                            nrow = (n_tmp_post_changepoint_training), ncol = n_sp_training)
  
  
  
  xbeta_pre_changepoint_mat_training = matrix((X_mat_pre_changepoint_training %*% beta_training),
                                              nrow = n_tmp_pre_changepoint_training, ncol = n_sp_training)
  
  xbeta_post_changepoint_mat_training = matrix((X_mat_post_changepoint_training %*% beta_star_training),
                                               nrow = (n_tmp_training - n_tmp_pre_changepoint_training), 
                                               ncol = n_sp_training)
  week_diff_mat_pre_changepoint = as.matrix(dist(1:n_tmp_pre_changepoint_training, 
                                                 diag =  TRUE,upper = TRUE))
  week_diff_mat_post_changepoint = as.matrix(dist((n_tmp_pre_changepoint_training+1):n_tmp_training, 
                                                  diag =  TRUE,upper = TRUE))
  
  # SIGMA_tmp_pre_changepoint_training = exp(-rho_tmp * week_diff_mat_pre_changepoint)
  # SIGMA_tmp_post_changepoint_training = exp(-rho_tmp * week_diff_mat_post_changepoint)
  # inv_SIGMA_tmp_pre_changepoint_training = chol2inv(chol(SIGMA_tmp_pre_changepoint_training))
  # inv_SIGMA_tmp_post_changepoint_training = chol2inv(chol(SIGMA_tmp_post_changepoint_training))
  
  # inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_changepoint_training = (inv_SIGMA_sp_training %x%
  #                                                              inv_SIGMA_tmp_pre_changepoint_training)
  # inv_SIGMA_sp_kro_inv_SIGMA_tmp_post_changepoint_training = (inv_SIGMA_sp_training %x%
  #                                                               inv_SIGMA_tmp_post_changepoint_training)
  
  for (sp in (1:n_sp_training)){
    idx=1:n_sp_training
    
    if(sp!=1){
      idx[1] = sp
      idx[sp] = 1
    }
    #Switching first and jth index column for V matrix and then taking values from the 2nd column
    v_cond_training = (v_mat_training %*% diag(n_sp_training)[idx,])[(1+n_tmp_pre_changepoint_training):
                                                                       (n_tmp_pre_changepoint_training * 
                                                                          n_sp_training)]
    
    v_cond_star_training = (v_star_mat_training %*% diag(n_sp_training)[idx,])[(1+n_tmp_post_changepoint_training):
                                                                                 (n_tmp_post_changepoint_training * 
                                                                                    n_sp_training)]
    #SIGMA' conditional without sigma1_sq as it's included in the posterior 
    # SIGMA_cond_training = (1 -t(SIGMA_sp_trans_21_mat_training[,j])%*%inv_SIGMA_sp_trans_22_matlist_training[[j]] %*%
    #                          SIGMA_sp_trans_21_mat_training[,j])%x% SIGMA_tmp_training
    mu_cond_training = ((t(SIGMA_sp_trans_21_mat_training[,sp]) %*% 
                           inv_SIGMA_sp_trans_22_matlist_training[[sp]]) %x% 
                          diag(n_tmp_pre_changepoint_training)) %*% v_cond_training
    
    mu_cond_star_training = ((t(SIGMA_sp_trans_21_mat_training[,sp]) %*% 
                                inv_SIGMA_sp_trans_22_matlist_training[[sp]]) %x% 
                               diag(n_tmp_post_changepoint_training)) %*% v_cond_star_training
    # inv_SIGMA_cond_training = chol2inv(chol(SIGMA_cond_training))
    # SIGMA_cond_pre_changepoint_training = (1 - t(SIGMA_sp_trans_21_mat_training[,sp]) %*% inv_SIGMA_sp_trans_22_matlist_training[[sp]] %*% 
    #                                          SIGMA_sp_trans_21_mat_training[,sp]) %x% SIGMA_tmp_pre_changepoint_training
    # 
    # SIGMA_cond_post_changepoint_training = (1 - t(SIGMA_sp_trans_21_mat_training[,sp]) %*% inv_SIGMA_sp_trans_22_matlist_training[[sp]] %*%
    #                                           SIGMA_sp_trans_21_mat_training[,sp]) %x% SIGMA_tmp_post_changepoint_training
    
    mat_mult_temp = as.numeric(1 - t(SIGMA_sp_trans_21_mat_training[,sp]) %*% 
                                 inv_SIGMA_sp_trans_22_matlist_training[[sp]] %*% 
                                 SIGMA_sp_trans_21_mat_training[,sp])
    
    inv_SIGMA_cond_pre_changepoint_training = (inv_SIGMA_tmp_pre_t_mat_list[[t_0_training/gap_bw_changepoints - 1]]/
                                                 mat_mult_temp)
    
    inv_SIGMA_cond_post_changepoint_training = (inv_SIGMA_tmp_post_t_mat_list[[t_0_training/gap_bw_changepoints - 1]]/
                                                  mat_mult_temp)
    
    v_covar_training = chol2inv(chol(inv_SIGMA_cond_pre_changepoint_training/sigma1_sq_training + 
                                       diag(n_tmp_pre_changepoint_training)/sigma2_sq_training))
    
    v_star_covar_training = chol2inv(chol(inv_SIGMA_cond_post_changepoint_training/sigma1_sq_star_training + 
                                            diag(n_tmp_post_changepoint_training)/sigma2_sq_star_training))
    
    v_mean_training  = v_covar_training %*% (inv_SIGMA_cond_pre_changepoint_training %*% 
                                               mu_cond_training/sigma1_sq_training + 
                                               (pi_mat_pre_changepoint_training[,sp] - 
                                                  xbeta_pre_changepoint_mat_training[,sp])/
                                               sigma2_sq_training)
    
    v_star_mean_training  = v_star_covar_training %*% (inv_SIGMA_cond_post_changepoint_training %*% 
                                                         mu_cond_star_training/sigma1_sq_star_training + 
                                                         (pi_mat_post_changepoint_training[,sp] - 
                                                            xbeta_post_changepoint_mat_training[,sp])/
                                                         sigma2_sq_star_training )
    
    v_mat_training[,sp] = v_mean_training + t(chol(v_covar_training))  %*% rnorm(n_tmp_pre_changepoint_training)
    v_star_mat_training[,sp] = (v_star_mean_training + t(chol(v_star_covar_training))  %*% 
                                  rnorm(n_tmp_post_changepoint_training))
    
  }
  training_pre_changepoint_df$v_vec = v_mat_training[1:(n_sp_training * n_tmp_pre_changepoint_training)]
  training_post_changepoint_df$v_vec = v_star_mat_training[1:(n_sp_training * n_tmp_post_changepoint_training)]
  
  training_df_clone = rbind(training_pre_changepoint_df, training_post_changepoint_df) %>% 
    arrange(country, state,key, lat, long, id,  week)
  
  v_vec_training = training_df_clone$v_vec
  
  total_tt_v = total_tt_v + as.numeric(difftime(Sys.time(), st_v), units="secs")
  ########################################################################################################################
  #beta posterior###################################
  st_beta = Sys.time()
  pi_pre_changepoint_training = (training_df_clone %>% filter(week <= changepoint_t_0_week))$pi_vec#pi_mat_pre_changepoint_training[1:(n_sp_training * n_tmp_pre_changepoint_training)]
  pi_post_changepoint_training = (training_df_clone %>% filter(week > changepoint_t_0_week))$pi_vec#pi_mat_pre_changepoint_training[1:(n_sp_training * n_tmp_post_changepoint_training)]
  v_training = v_mat_training[1:(n_sp_training * n_tmp_pre_changepoint_training)]
  v_star_training = v_star_mat_training[1:(n_sp_training * n_tmp_post_changepoint_training)]
  #dividing by 10^4 for large variance for beta prior
  beta_covar_training = solve( (t(X_mat_pre_changepoint_training) %*% 
                                  X_mat_pre_changepoint_training)/sigma2_sq_training) #+ diag(p+1)/10^4
  
  beta_star_covar_training = solve( (t(X_mat_post_changepoint_training) %*% 
                                       X_mat_post_changepoint_training)/sigma2_sq_star_training)
  
  beta_mean_training = beta_covar_training %*% ( t(X_mat_pre_changepoint_training)%*% 
                                                   (pi_pre_changepoint_training - v_training))/sigma2_sq_training
  
  beta_star_mean_training = (beta_star_covar_training %*% ( t(X_mat_post_changepoint_training)%*% 
                                                              (pi_post_changepoint_training - v_star_training))/
                               sigma2_sq_star_training)
  
  beta_training = beta_mean_training + t(chol(beta_covar_training)) %*% rnorm(p+1)
  
  beta_star_training = beta_star_mean_training + t(chol(beta_star_covar_training)) %*% rnorm(p+1)
  
  total_tt_beta = total_tt_beta + as.numeric(difftime( Sys.time(), st_beta), units="secs")
  #######################################################################################################################
  #Sigma1^2 posterior######################################################
  st_sigma1_sq = Sys.time()
  
  # sigma1_sq_a_training = a_training + (n_sp_training * n_tmp_pre_changepoint_training)/2
  # 
  # sigma1_sq_star_a_training = a_training + (n_sp_training * n_tmp_post_changepoint_training)/2
  # 
  # sigma1_sq_lambda_training = (t(v_training) %*% 
  #                                inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list[[t_0_training/gap_bw_changepoints - 1]] 
  #                              %*% v_training)/2 + lambda_training
  # 
  # sigma1_sq_star_lambda_training = (t(v_star_training) %*% 
  #                                     inv_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_mat_list[[t_0_training/gap_bw_changepoints - 1]] 
  #                                   %*% v_star_training)/2 + lambda_training
  
  sigma1_sq_training = 1#invgamma::rinvgamma(1,sigma1_sq_a_training,rate = sigma1_sq_lambda_training)
  
  sigma1_sq_star_training = 1 #invgamma::rinvgamma(1, sigma1_sq_star_a_training,
  #rate = sigma1_sq_star_lambda_training)
  
  total_tt_sigma1_sq = total_tt_sigma1_sq + as.numeric(difftime(Sys.time(), st_sigma1_sq), units="secs")
  ##########################################################################################################################
  #Sigma2^2 posterior############################################
  st_sigma2_sq = Sys.time()
  
  # sigma2_sq_a_training = a_training + (n_sp_training * n_tmp_pre_changepoint_training)/2
  # 
  # sigma2_sq_star_a_training = a_training + (n_sp_training * n_tmp_post_changepoint_training)/2
  # 
  # temp_vec = pi_pre_changepoint_training - X_mat_pre_changepoint_training%*%beta_training - v_training
  # 
  # temp_vec2 = (pi_post_changepoint_training - X_mat_post_changepoint_training %*% beta_star_training - 
  #                v_star_training)
  # 
  # sigma2_sq_lambda_training = sum(temp_vec^2)/2 + lambda_training
  # 
  # sigma2_sq_star_lambda_training = sum(temp_vec2^2)/2 + lambda_training
  
  sigma2_sq_training = 1#invgamma::rinvgamma(1,sigma2_sq_a_training,rate = sigma2_sq_lambda_training)
  
  sigma2_sq_star_training = 1#invgamma::rinvgamma(1, sigma2_sq_star_a_training,
  #rate = sigma2_sq_star_lambda_training)
  
  total_tt_sigma2_sq = total_tt_sigma2_sq + as.numeric(difftime(Sys.time(), st_sigma2_sq), units="secs")
  ###################################################################################################################
  #delta posterior###############################
  st_delta = Sys.time()
  for (j in 1: (no_of_categ -1)){
    pi_cat1 = pi_vec[y_training == j] #because j starts from 1 hence category values are okay but delta vec values are shifted
    c1 = max(pi_cat1) 
    pi_cat2 = pi_vec[y_training == (j + 1)] 
    c2 = min(pi_cat2)
    a = pnorm( delta_training [j], mean =0, stdnorm )
    b = pnorm( delta_training [j+2] , mean =0, stdnorm )
    w1 = (pnorm(c1 ,0, stdnorm)-a)/(b-a)
    w2 = (pnorm(c2 ,0, stdnorm)-a)/(b-a)
    w = rtrunc (1, spec = "beta",shape1 = alpha_vec[j],
                shape2 = alpha_vec [j+1] ,a=w1 ,b=w2)
    
    delta_training[j +1] = qnorm ((b-a)*w+a, mean =0, stdnorm )
  }
  total_tt_delta = total_tt_delta + as.numeric(difftime(Sys.time(), st_delta), units="secs")
  ################################################################################################################
  #t_0 posterior####################################
  st_t_0 = Sys.time()
  log_prop_pmf_t_0_vec_training = vector()
  for (tp in (changepoint_vec)) {
    
    week_t = unique(training_df_clone$week)[tp]
    
    X_mat_pre_t = cbind(rep(1,n_sp_training * tp), 
                        log((training_df_clone %>% filter(week <= week_t))$population), 
                        (training_df_clone %>% filter(week <= week_t))$time, 
                        (training_df_clone %>% filter(week <= week_t))$time_sq, 
                        (training_df_clone %>% filter(week <= week_t))$prev_log_new_death)
    
    X_mat_post_t = cbind(rep(1,n_sp_training * (n_tmp_training - tp)), 
                         log((training_df_clone %>% filter(week > week_t))$population),
                         (training_df_clone %>% filter(week > week_t))$time, 
                         (training_df_clone %>% filter(week > week_t))$time_sq, 
                         (training_df_clone %>% filter(week > week_t))$prev_log_new_death)
    
    v_vec_pre_t = (training_df_clone %>% filter(week <= week_t))$v_vec
    v_vec_post_t = (training_df_clone %>% filter(week > week_t))$v_vec
    pi_vec_pre_t = (training_df_clone %>% filter(week <= week_t))$pi_vec
    pi_vec_post_t = (training_df_clone %>% filter(week > week_t))$pi_vec
    
    temp_vec = pi_vec_pre_t - X_mat_pre_t %*% beta_training - v_vec_pre_t
    temp_vec2 = pi_vec_post_t - X_mat_post_t %*% beta_star_training - v_vec_post_t
    
    log_prop_pmf_t_0_training = ((-n_sp_training*tp/2)*log(2*pi*sigma2_sq_training) + 
                                   ((-1/(2*sigma2_sq_training)) * sum(temp_vec^2)) + 
                                   ((-n_sp_training*(n_tmp_training - tp)/2) * 
                                      log(2*pi*sigma2_sq_star_training)) +
                                   ((-1/(2*sigma2_sq_star_training)) * sum(temp_vec2^2)) + 
                                   ((-n_sp_training*tp/2) * log(2*pi*sigma1_sq_training)) + 
                                   ((-1/2) * 
                                      log_det_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_vec[tp/gap_bw_changepoints - 1]) + 
                                   ((-1/(2 * sigma1_sq_training)) * 
                                      (arma_mm(inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list[[tp/gap_bw_changepoints - 1]],
                                               v_vec_pre_t, v_vec_pre_t))) + 
                                   ((-n_sp_training*(n_tmp_training - tp)/2) * 
                                      log(2*pi*sigma1_sq_star_training)) + 
                                   (-1/2) * log_det_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_vec[tp/gap_bw_changepoints - 1] + 
                                   ((-1/(2*sigma1_sq_star_training)) * 
                                      (arma_mm(inv_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_mat_list[[tp/gap_bw_changepoints - 1]],
                                               v_vec_post_t, v_vec_post_t))))
    
    log_prop_pmf_t_0_vec_training = c(log_prop_pmf_t_0_vec_training, log_prop_pmf_t_0_training)
  }
  #min_prop_pmf_t_0_vec_training = min(prop_pmf_t_0_vec_training)
  max_log_prop_pmf_t_0_vec_training = max(log_prop_pmf_t_0_vec_training)
  # if(min_prop_pmf_t_0_vec_training < 0){
  #   prop_pmf_t_0_vec_training = prop_pmf_t_0_vec_training + (-min_prop_pmf_t_0_vec_training)
  # }
  log_prop_pmf_t_0_vec_training = log_prop_pmf_t_0_vec_training - max_log_prop_pmf_t_0_vec_training
  prop_pmf_t_0_vec_training = exp(log_prop_pmf_t_0_vec_training)
  
  const_prop = 1/sum(prop_pmf_t_0_vec_training)
  
  t_0_training = sample(x = (changepoint_vec),1, replace = FALSE, 
                        prob = round(prop_pmf_t_0_vec_training * const_prop, digits = 3))
  
  n_tmp_pre_changepoint_training = t_0_training
  n_tmp_post_changepoint_training = n_tmp_training - n_tmp_pre_changepoint_training
  changepoint_t_0_week = unique(training_df$week)[n_tmp_pre_changepoint_training]
  
  total_tt_t_0 = total_tt_t_0 + as.numeric(difftime(Sys.time(), st_t_0), units="secs")
  
  #Collecting samples for prediction
  if(i > burn_period_training & i %% diff_in_random_draws_training == 0){
    
    v_vec_sample_training = cbind(v_vec_sample_training, v_vec_training)
    pi_vec_sample_training = cbind(pi_vec_sample_training, pi_vec_training)
    delta_sample_training = cbind(delta_sample_training, delta_training)
    beta_sample_training = cbind(beta_sample_training, beta_training)
    beta_star_sample_training = cbind(beta_star_sample_training, beta_star_training)
    sigma1_sq_sample_training = c(sigma1_sq_sample_training, sigma1_sq_training)
    sigma1_sq_star_sample_training = c(sigma1_sq_star_sample_training, sigma1_sq_star_training)
    sigma2_sq_sample_training = c(sigma2_sq_sample_training, sigma2_sq_training)
    sigma2_sq_star_sample_training = c(sigma2_sq_star_sample_training, sigma2_sq_star_training)
    t_0_sample_training = c(t_0_sample_training, t_0_training)
    
  }
  
  
  
}
#time analysis
avg_tt_v = total_tt_v/(burn_period_training+sample_size_training * diff_in_random_draws_training)
avg_tt_beta = total_tt_beta/(burn_period_training+sample_size_training*diff_in_random_draws_training)
avg_tt_sigma1_sq = total_tt_sigma1_sq/(burn_period_training+sample_size_training*diff_in_random_draws_training)
avg_tt_sigma2_sq = total_tt_sigma2_sq/(burn_period_training+sample_size_training*diff_in_random_draws_training)
avg_tt_delta = total_tt_delta/(burn_period_training+sample_size_training*diff_in_random_draws_training)
avg_tt_pi = total_tt_pi/(burn_period_training+sample_size_training*diff_in_random_draws_training)
avg_tt_t_0 = total_tt_t_0/(burn_period_training+sample_size_training*diff_in_random_draws_training)

sum(avg_tt_v,avg_tt_beta, avg_tt_sigma1_sq, avg_tt_sigma2_sq, avg_tt_delta, avg_tt_pi, avg_tt_t_0)
#estimation of y with the posterior parameters
# estimated_pi_training = X_mat_pre_changepoint_training %*% rowMeans(beta_sample_training) + 
#   rowMeans(v_sample_training) + rnorm(n_total_training,mean = 0, sd= sqrt(mean(sigma2_sq_sample_training)))
# 
# estimated_pi_star_training = X_mat_post_changepoint_training %*% rowMeans(beta_star_sample_training) + 
#   rowMeans(v_star_sample_training) + rnorm(n_total_training,mean = 0, sd= sqrt(mean(sigma2_sq_star_sample_training)))

#####################################################################################################

training_df_clone$pi_vec_est = rowMeans(pi_vec_sample_training)

mean_delta_training = rowMeans(delta_sample_training)

training_df_clone = training_df_clone %>% 
  mutate(pred_category = ifelse(pi_vec_est > mean_delta_training[4],4,
                                ifelse(pi_vec_est > mean_delta_training[3], 3,
                                       ifelse(pi_vec_est > mean_delta_training[2],2,1))))
pred_category = training_df_clone$pred_category

pred_category_mat = cbind(pred_category_mat, pred_category)

accuracy_pred = sum(training_df_clone$category_art == training_df_clone$pred_category)/n_total_training

accuracy_pred_vec = c(accuracy_pred_vec, accuracy_pred)
# mape_training = (sum(abs( (y_training - estimated_pi_training)/y_training )))/n_total_training
# mape_vec_training_training = c(mape_vec_training_training, mape_training)



#Testing for significance
beta_significance_vec = vector()
for (ts in (1:(p+1))) {
  beta_significance_vec = c(beta_significance_vec, 
                            unname(!(0 >= quantile(beta_sample_training[ts,], 0.025) & 
                                       0 <= quantile(beta_sample_training[ts,], 0.975))))
}
beta_significance_mat_training = cbind(beta_significance_mat_training, beta_significance_vec)

beta_star_significance_vec = vector()
for (ts in (1:(p+1))) {
  beta_star_significance_vec = c(beta_star_significance_vec, 
                                 unname(!(0 >= quantile(beta_star_sample_training[ts,], 0.025) & 
                                            0 <= quantile(beta_star_sample_training[ts,], 0.975))))
}
beta_star_significance_mat_training = cbind(beta_star_significance_mat_training, 
                                            beta_star_significance_vec)

#Saving samples for each rho


pi_mat_sample_list_training[[rho_idx]] = pi_vec_sample_training
delta_mat_sample_list_training[[rho_idx]] = delta_sample_training
v_vec_mat_sample_list_training[[rho_idx]] = v_vec_sample_training

beta_sample_list_training[[rho_idx]] = beta_sample_training
beta_star_sample_list_training[[rho_idx]] = beta_star_sample_training
t_0_sample_vec_training = cbind(t_0_sample_vec_training, t_0_sample_training)
sigma1_sq_sample_vec_training = cbind(sigma1_sq_sample_vec_training, sigma1_sq_sample_training)
sigma1_sq_star_sample_vec_training = cbind(sigma1_sq_star_sample_vec_training, sigma1_sq_star_sample_training)
sigma2_sq_sample_vec_training = cbind(sigma2_sq_sample_vec_training, sigma2_sq_sample_training)
sigma2_sq_star_sample_vec_training = cbind(sigma2_sq_star_sample_vec_training, sigma2_sq_star_sample_training)


rho_idx = rho_idx + 1




min_idx = which(mape_vec_validation_training == min(mape_vec_validation_training))
if(min_idx %% length(rho_vec_sp_training) == 0){
  phi_s = length(rho_vec_sp_training)
  phi_t = min_idx/length(rho_vec_sp_training)
}else{
  phi_s = rho_vec_sp_training[min_idx %% length(rho_vec_sp_training)]
  phi_t = rho_vec_tmp_training[floor(min_idx/length(rho_vec_sp_training)) + 1]
}
return(c(phi_s, phi_t))


#scoring validation################################################################################


len_phi_grid = length(pi_mat_sample_list_training)
pred_df_list = list()
scoring_vec = vector()

training_pred_mat = vector()
training_pred_vec = rep(0, 4)
for (dt_pt in 1:n_total_training) {
  training_pred_vec[training_df_clone$category_art[dt_pt]] = 1
  training_pred_mat = rbind(training_pred_mat, unname(training_pred_vec))
  training_pred_vec = rep(0, 4)
}



for (phi_i in 1:len_phi_grid) {
  
  pi_sample_curr = pi_mat_sample_list_training[[phi_i]]
  delta_sample_curr = delta_mat_sample_list_training[[phi_i]]
  sample_size = ncol(delta_sample_curr)
  
  pred_df = training_df_clone %>% dplyr::select(category_art)
  for (sample_j in 1:sample_size) {
    pred_df_sample_j = as.data.frame(pi_sample_curr[,sample_j])
    names(pred_df_sample_j) = "pi_Values"
    
    pred_df_sample_j = pred_df_sample_j %>% 
      dplyr::mutate(pred_category = ifelse(pi_Values > delta_sample_curr[,sample_j][4],4,
                                           ifelse(pi_Values > delta_sample_curr[,sample_j][3], 3,
                                                  ifelse(pi_Values > delta_sample_curr[,sample_j][2],2,1))))
    
    colname = paste0("pred_category_sample_",sample_j)
    pred_df[colname] = pred_df_sample_j$pred_category
  }
  pred_df_list[[phi_i]] = pred_df
  
  pred_cat_1_vec = rowSums(pred_df[,2:1001] == 1)
  
  pred_cat_2_vec = rowSums(pred_df[,2:1001] == 2)
  
  pred_cat_3_vec = rowSums(pred_df[,2:1001] == 3)
  
  pred_cat_4_vec = rowSums(pred_df[,2:1001] == 4)
  
  pred_cat_prob_mat = cbind(unname(pred_cat_1_vec), unname(pred_cat_2_vec),
                            unname(pred_cat_3_vec), unname(pred_cat_4_vec))/(sample_size)
  
  error_score = sum((pred_cat_prob_mat - training_pred_mat)^2)
  
  scoring_vec = c(scoring_vec, error_score)
}











