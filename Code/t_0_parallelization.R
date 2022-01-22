library(parallel)
library(doParallel)
library(foreach)

st =Sys.time()

for (tp in (changepoint_vec)) {
  st =Sys.time()
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
  Sys.time() -st
  log_prop_pmf_t_0_training = ((-n_sp_training*tp/2)*log(2*pi*sigma2_sq_training) + 
                                 ((-1/(2*sigma2_sq_training)) * sum(temp_vec^2)) + 
                                 ((-n_sp_training*(n_tmp_training - tp)/2) * 
                                    log(2*pi*sigma2_sq_star_training)) +
                                 ((-1/(2*sigma2_sq_star_training)) * sum(temp_vec2^2)) + 
                                 ((-n_sp_training*tp/2) * log(2*pi*sigma1_sq_training)) + 
                                 ((-1/2) * 
                                    log_det_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_vec[tp/gap_bw_changepoints]) + 
                                 ((-1/(2 * sigma1_sq_training)) * 
                                    (crossprod(v_vec_pre_t,inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list[[tp/gap_bw_changepoints]]) %*%  
                                      v_vec_pre_t)) + 
                                 ((-n_sp_training*(n_tmp_training - tp)/2) * 
                                    log(2*pi*sigma1_sq_star_training)) + 
                                 (-1/2) * log_det_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_vec[tp/gap_bw_changepoints] + 
                                 ((-1/(2*sigma1_sq_star_training)) * 
                                    (crossprod(v_vec_post_t, inv_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_mat_list[[tp/gap_bw_changepoints]]) 
                                        %*% v_vec_post_t)))
  
  #log_prop_pmf_t_0_vec_training = c(log_prop_pmf_t_0_vec_training, log_prop_pmf_t_0_training)
}


Sys.time() -st

#######################################################################################################################
st =Sys.time()
# numCores = detectCores() - 2
# cl = makePSOCKcluster(numCores)
# registerDoParallel(cl)
for(tp in (changepoint_vec)){
  
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
  
  v_vec_pre_t = matrix((training_df_clone %>% filter(week <= week_t))$v_vec, ncol = 1)
  v_vec_post_t = matrix((training_df_clone %>% filter(week > week_t))$v_vec, ncol = 1)
  pi_vec_pre_t = matrix((training_df_clone %>% filter(week <= week_t))$pi_vec, ncol = 1)
  pi_vec_post_t = matrix((training_df_clone %>% filter(week > week_t))$pi_vec, ncol = 1)
  
  temp_vec = pi_vec_pre_t - X_mat_pre_t %*% beta_training - v_vec_pre_t
  temp_vec2 = pi_vec_post_t - X_mat_post_t %*% beta_star_training - v_vec_post_t
  
  log_prop_pmf_t_0_training = ((-n_sp_training*tp/2)*log(2*pi*sigma2_sq_training) + 
                                 ((-1/(2*sigma2_sq_training)) * sum(temp_vec^2)) + 
                                 ((-n_sp_training*(n_tmp_training - tp)/2) * 
                                    log(2*pi*sigma2_sq_star_training)) +
                                 ((-1/(2*sigma2_sq_star_training)) * sum(temp_vec2^2)) + 
                                 ((-n_sp_training*tp/2) * log(2*pi*sigma1_sq_training)) + 
                                 ((-1/2) * 
                                    log_det_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_vec[tp/gap_bw_changepoints]) + 
                                 ((-1/(2 * sigma1_sq_training)) * 
                                    (crossprod(v_vec_pre_t,inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list[[tp/gap_bw_changepoints]])  
                                     %*% v_vec_pre_t)) + 
                                 ((-n_sp_training*(n_tmp_training - tp)/2) * 
                                    log(2*pi*sigma1_sq_star_training)) + 
                                 (-1/2) * log_det_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_vec[tp/gap_bw_changepoints] + 
                                 ((-1/(2*sigma1_sq_star_training)) * 
                                    (crossprod(v_vec_post_t, inv_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_mat_list[[tp/gap_bw_changepoints]]) 
                                     %*% v_vec_post_t)))
  
  #log_prop_pmf_t_0_vec_training = c(log_prop_pmf_t_0_vec_training, log_prop_pmf_t_0_training)
}

# stopCluster(cl)
Sys.time() -st 

