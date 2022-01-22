# Initial values of delta #
set.seed(1211)
delta_training = runif (no_of_categ -1 , -20 ,20)
delta_training = sort (delta_training)
delta_training = c(-Inf , delta_training ,Inf)
alpha_vec =rep(1,no_of_categ) # shape and scale parameters of trunc Beta .
stdnorm = 10 # SD value for finding CDF of the delta posterior .
deltamat = vector ( mode =" numeric ",length =0)

for (j in 1: (no_of_categ -1)){
  pi_cat1 = pi_vec[y_training==j]
  c1=max( pi_cat1 ); 
  pi_cat2 = pi_vec[y_training== (j + 1)] 
  c2=min( pi_cat2 )
  a= pnorm ( delta_training [j], mean =0, stdnorm )
  b= pnorm ( delta_training [j+2] , mean =0, stdnorm )
  w1 =( pnorm (c1 ,0, stdnorm )-a)/(b-a)
  w2 =( pnorm (c2 ,0, stdnorm )-a)/(b-a)
  w= rtrunc (1, spec =" beta ",shape1 = alpha_vec [j],
             shape2 = alpha_vec [j+1] ,a=w1 ,b=w2)
  delta_training [j +1]= qnorm ((b-a)*w+a, mean =0, stdnorm )
}
deltamat = cbind ( deltamat , delta_training )
