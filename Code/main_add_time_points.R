new_df = training_df[,1:8]

id_vec = unique(new_df$id)

sp_i = 1

for (id_num in id_vec) {
  
  id_df = new_df %>%  filter(id == id_num)
  
  min_date = min(id_df$week)
  
  max_date = max(id_df$week)
  
  uniq_df = unique(id_df[,1:7])
  
  for (tp in 1:30) {
    tmp_df = uniq_df
    tmp_df$week = max_date + 7 * tp
    
    id_df = rbind(id_df, tmp_df)
  }
  
  if(sp_i == 1){
    
    new_training_df = id_df
    
  } else{
    
    new_training_df = rbind(new_training_df, id_df)
  }
  
  
  sp_i = sp_i + 1
}

training_df = new_training_df %>% arrange(country, state,key, lat, long, id,  week)

tot_dp = nrow(training_df)

set.seed(100)
training_df$X1 = rnorm(tot_dp)

set.seed(101)
training_df$X2 = rnorm(tot_dp)

set.seed(102)
training_df$X3 = rnorm(tot_dp)

set.seed(103)
training_df$X4 = rnorm(tot_dp)
