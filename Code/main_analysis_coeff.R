phi_pos_1 = 1

phi_pos_10 = 10

beta_est = rowMeans(beta_sample_list_training[[phi_pos_1]])

beta_est_LB = apply(beta_sample_list_training[[phi_pos_1]],1, FUN = function(x) quantile(x,0.025))

beta_est_UB = apply(beta_sample_list_training[[phi_pos_1]],1, FUN = function(x) quantile(x,0.975))

sigma1_sq_est = mean(sigma1_sq_sample_vec_training[,phi_pos_1])

sigma1_sq_est_LB = unname(quantile(sigma1_sq_sample_vec_training[,phi_pos_1], 0.025))

sigma1_sq_est_UB = unname(quantile(sigma1_sq_sample_vec_training[,phi_pos_1], 0.975))

sigma2_sq_est = mean(sigma2_sq_sample_vec_training[,phi_pos_1])

sigma2_sq_est_LB = unname(quantile(sigma2_sq_sample_vec_training[,phi_pos_1], 0.025))

sigma2_sq_est_UB = unname(quantile(sigma2_sq_sample_vec_training[,phi_pos_1], 0.975))

beta_star_est = rowMeans(beta_star_sample_list_training[[phi_pos_1]])

beta_star_est_LB = apply(beta_star_sample_list_training[[phi_pos_1]],1, FUN = function(x) quantile(x,0.025))

beta_star_est_UB = apply(beta_star_sample_list_training[[phi_pos_1]],1, FUN = function(x) quantile(x,0.975))

sigma1_sq_star_est = mean(sigma1_sq_star_sample_vec_training[,phi_pos_1])

sigma1_sq_star_est_LB = unname(quantile(sigma1_sq_star_sample_vec_training[,phi_pos_1], 0.025))

sigma1_sq_star_est_UB = unname(quantile(sigma1_sq_star_sample_vec_training[,phi_pos_1], 0.975))

sigma2_sq_star_est = mean(sigma2_sq_star_sample_vec_training[,phi_pos_1])

sigma2_sq_star_est_LB = unname(quantile(sigma2_sq_star_sample_vec_training[,phi_pos_1], 0.025))

sigma2_sq_star_est_UB = unname(quantile(sigma2_sq_star_sample_vec_training[,phi_pos_1], 0.975))

t_0_est = mean(t_0_sample_vec_training[,phi_pos_1])

t_0_est_LB = unname(quantile(t_0_sample_vec_training[,phi_pos_1], 0.025))

t_0_est_UB = unname(quantile(t_0_sample_vec_training[,phi_pos_1], 0.975))

est_vec = c(beta_est, sigma1_sq_est, sigma2_sq_est, beta_star_est, sigma1_sq_star_est, sigma2_sq_star_est, t_0_est)

est_LB_vec = c(beta_est_LB, sigma1_sq_est_LB, sigma2_sq_est_LB, beta_star_est_LB, sigma1_sq_star_est_LB, 
                   sigma2_sq_star_est_LB, t_0_est_LB)

est_UB_vec = c(beta_est_UB, sigma1_sq_est_UB, sigma2_sq_est_UB, beta_star_est_UB, sigma1_sq_star_est_UB, 
                   sigma2_sq_star_est_UB, t_0_est_UB)

Var_names = c("Intercept", "Population", "time", "time_sq", "prev_log_new_death","sigma_sq_v","sigma_sq_epsilon",
              "Intercept*", "Population*", "time*", "time_sq*", "prev_log_new_death*","sigma_sq_v*","sigma_sq_epsilon*",
              "change point")

temp_df = as.data.frame(cbind(Var_names,est_vec,est_LB_vec, est_UB_vec))

temp_df$est_vec = as.numeric(temp_df$est_vec)

temp_df$est_LB_vec = as.numeric(temp_df$est_LB_vec)

temp_df$est_UB_vec = as.numeric(temp_df$est_UB_vec)

temp_df %>% mutate_if(is.numeric,round,digits = 4) %>% View()
########################################################################################################################

beta_est2 = rowMeans(beta_sample_list_training[[phi_pos_10]])

beta_est2_LB = apply(beta_sample_list_training[[phi_pos_10]],1, FUN = function(x) quantile(x,0.025))

beta_est2_UB = apply(beta_sample_list_training[[phi_pos_10]],1, FUN = function(x) quantile(x,0.975))

sigma1_sq_est2 = mean(sigma1_sq_sample_vec_training[,phi_pos_10])

sigma1_sq_est2_LB = unname(quantile(sigma1_sq_sample_vec_training[,phi_pos_10], 0.025))

sigma1_sq_est2_UB = unname(quantile(sigma1_sq_sample_vec_training[,phi_pos_10], 0.975))

sigma2_sq_est2 = mean(sigma2_sq_sample_vec_training[,phi_pos_10])

sigma2_sq_est2_LB = unname(quantile(sigma2_sq_sample_vec_training[,phi_pos_10], 0.025))

sigma2_sq_est2_UB = unname(quantile(sigma2_sq_sample_vec_training[,phi_pos_10], 0.975))

beta_star_est2 = rowMeans(beta_star_sample_list_training[[phi_pos_10]])

beta_star_est2_LB = apply(beta_star_sample_list_training[[phi_pos_10]],1, FUN = function(x) quantile(x,0.025))

beta_star_est2_UB = apply(beta_star_sample_list_training[[phi_pos_10]],1, FUN = function(x) quantile(x,0.975))

sigma1_sq_star_est2 = mean(sigma1_sq_star_sample_vec_training[,phi_pos_10])

sigma1_sq_star_est2_LB = unname(quantile(sigma1_sq_star_sample_vec_training[,phi_pos_10], 0.025))

sigma1_sq_star_est2_UB = unname(quantile(sigma1_sq_star_sample_vec_training[,phi_pos_10], 0.975))

sigma2_sq_star_est2 = mean(sigma2_sq_star_sample_vec_training[,phi_pos_10])

sigma2_sq_star_est2_LB = unname(quantile(sigma2_sq_star_sample_vec_training[,phi_pos_10], 0.025))

sigma2_sq_star_est2_UB = unname(quantile(sigma2_sq_star_sample_vec_training[,phi_pos_10], 0.975))

t_0_est2 = mean(t_0_sample_vec_training[,phi_pos_10])

t_0_est2_LB = unname(quantile(t_0_sample_vec_training[,phi_pos_10], 0.025))

t_0_est2_UB = unname(quantile(t_0_sample_vec_training[,phi_pos_10], 0.975))

est2_vec = c(beta_est2, sigma1_sq_est2, sigma2_sq_est2, beta_star_est2, sigma1_sq_star_est2, sigma2_sq_star_est2, t_0_est2)

est2_LB_vec = c(beta_est2_LB, sigma1_sq_est2_LB, sigma2_sq_est2_LB, beta_star_est2_LB, sigma1_sq_star_est2_LB, 
                sigma2_sq_star_est2_LB, t_0_est2_LB)

est2_UB_vec = c(beta_est2_UB, sigma1_sq_est2_UB, sigma2_sq_est2_UB, beta_star_est2_UB, sigma1_sq_star_est2_UB, 
                sigma2_sq_star_est2_UB, t_0_est2_UB)

Var_names = c("Intercept", "Population", "time", "time_sq", "prev_log_new_death","sigma_sq_v","sigma_sq_epsilon",
              "Intercept*", "Population*", "time*", "time_sq*", "prev_log_new_death*","sigma_sq_v*","sigma_sq_epsilon*",
              "change point")

temp_df2 = as.data.frame(cbind(Var_names,est2_vec,est2_LB_vec, est2_UB_vec))

temp_df2$est2_vec = as.numeric(temp_df2$est2_vec)

temp_df2$est2_LB_vec = as.numeric(temp_df2$est2_LB_vec)

temp_df2$est2_UB_vec = as.numeric(temp_df2$est2_UB_vec)

temp_df2 %>% mutate_if(is.numeric,round,digits = 4) %>% View()
#pi vec analysis##############################################################################################################


pi_vec_df = as.data.frame(cbind(apply(pi_vec_sample_training,1, FUN = function(x) quantile(x, 0.025)),
                          apply(pi_vec_sample_training,1, FUN = function(x) quantile(x, 0.975)),
                          training_df_clone$pi_vec_art,
                          training_df_clone$category_art))



