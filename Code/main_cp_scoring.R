len_phi_grid = length(pi_mat_sample_list_training)
pred_df_list = list()
scoring_vec = vector()

training_pred_mat = vector()
training_pred_vec = rep(0, 4)
for (dt_pt in 1:n_total_training) {
  training_pred_vec[training_df$category[dt_pt]] = 1
  training_pred_mat = rbind(training_pred_mat, unname(training_pred_vec))
  training_pred_vec = rep(0, 4)
}



for (phi_i in 1:len_phi_grid) {
  
  pi_sample_curr = pi_mat_sample_list_training[[phi_i]]
  delta_sample_curr = delta_mat_sample_list_training[[phi_i]]
  sample_size = ncol(delta_sample_curr)
  
  pred_df = training_df %>% dplyr::select(category)
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
