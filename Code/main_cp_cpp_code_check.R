library(Rcpp)
library(RcppArmadillo)

st = Sys.time() 
( crossprod(v_vec_pre_t,inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list[[tp/gap_bw_changepoints - 1]]) 
  %*% v_vec_pre_t)
Sys.time() - st

st = Sys.time() 
(t(v_vec_pre_t) %*% 
    inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list[[tp/gap_bw_changepoints - 1]] 
  %*% v_vec_pre_t)
Sys.time() - st

st = Sys.time() 
(v_vec_pre_t %*% 
    inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list[[tp/gap_bw_changepoints - 1]] 
  %*% v_vec_pre_t)
Sys.time() - st

st = Sys.time() 
(as.vector(v_vec_pre_t) %*% 
    inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list[[tp/gap_bw_changepoints - 1]] 
  %*% as.vector(v_vec_pre_t))
Sys.time() - st

#cpp function################################################################################################
#cppFunction('double matrixmult(NumericMatrix m, )
#             ')

arma_code <- 
  "arma::vec arma_mm(const arma::mat& m, const arma::mat& v, const arma::mat& tv) {
       return tv * m * v;
   };"

arma_code2 <- 
  "arma::vec arma_mm(const arma::mat& m, const arma::vec& v, const arma::rowvec& tv) {
       return tv * m * v;
   };"

arma_mm = cppFunction(code = arma_code, depends = "RcppArmadillo")

arma_mm2 = cppFunction(code = arma_code2, depends = "RcppArmadillo")
#arma_mm##############################################################
st = Sys.time()
tmp = matrix(v_vec_pre_t, ncol = 1)
arma_mm(inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list[[tp/gap_bw_changepoints - 1]], 
        tmp, t(tmp))
rm(tmp)
Sys.time() - st
#arma_mm2##############################################################
st = Sys.time() 
arma_mm2(inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list[[tp/gap_bw_changepoints - 1]], 
        v_vec_pre_t, v_vec_pre_t)
Sys.time() - st

#post mult######################################
st = Sys.time() 
arma_mm(inv_SIGMA_sp_kro_inv_SIGMA_tmp_post_t_mat_list[[tp/gap_bw_changepoints - 1]], 
        matrix(v_vec_post_t, ncol = 1), matrix(v_vec_post_t, nrow = 1))
Sys.time() - st


###############################################################################
st = Sys.time()
as.numeric(1 - t(SIGMA_sp_trans_21_mat_training[,sp]) %*% 
             inv_SIGMA_sp_trans_22_matlist_training[[sp]] %*% 
             SIGMA_sp_trans_21_mat_training[,sp])
Sys.time() - st


st = Sys.time()
# as.numeric(1 - t(SIGMA_sp_trans_21_mat_training[,sp]) %*% 
#              inv_SIGMA_sp_trans_22_matlist_training[[sp]] %*% 
#              SIGMA_sp_trans_21_mat_training[,sp])
1 - as.numeric(arma_mm(inv_SIGMA_sp_trans_22_matlist_training[[sp]], 
                       matrix(SIGMA_sp_trans_21_mat_training[,sp], ncol = 1), matrix(SIGMA_sp_trans_21_mat_training[,sp], nrow = 1)))
Sys.time() - st
####################################################################################
st = Sys.time()
(t(v_training) %*% 
    inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list[[t_0_training/gap_bw_changepoints - 1]] 
  %*% v_training)
Sys.time() - st


st = Sys.time() 
arma_mm2(inv_SIGMA_sp_kro_inv_SIGMA_tmp_pre_t_mat_list[[t_0_training/gap_bw_changepoints - 1]], 
         v_training, v_training)
Sys.time() - st
