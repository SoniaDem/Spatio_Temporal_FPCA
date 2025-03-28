library('fda')

nn = 100
M = 50
N = 50

#sd of all image = 1.83 or 1.65 with omega 1
noise_level = 0.32
noise_level_name = 1
iterations = 2
#this can only be 0 (no omega 1) or 1 (with)
omega_check = 1

VE_proc_iter = array(NA, dim=c(iterations,5))
VE_cov_iter = array(NA, dim=c(iterations,5))
VE_pca_iter = array(NA, dim=c(iterations,5))
MSE_proc_iter = array(NA, dim=c(iterations,5))
MSE_cov_iter = array(NA, dim=c(iterations,5))
MSE_pca_iter = array(NA, dim=c(iterations,5))
rproc_iter = rep(NA, iterations)
rcov_iter = rep(NA, iterations)
rpca_iter= rep(NA, iterations)
score_diff = array(NA, dim=c(iterations,3))
PC_diff = array(NA, dim=c(iterations,3))


# set seed for consistency 

set.seed(1)

seed_vector = sample.int(1000, iterations, replace = TRUE)
# functions

omega1 <- function(ss_1, ss_2, M, N)
{
  return(sin(4*pi*(ss_1*ss_2/M/N)))
}

omega2 <- function(ss_1, M)
{
  return(cos(3*pi*(ss_1/M)))
}

omega3 <- function(ss_1, M)
{
  return(sin(pi*(ss_1/M)))
}

xi2 <- function(ss_2, d, N)
{
  return(d*cos(5.5*ss_2/N))
}

xi3 <- function(ss_2, d, N)
{
  return(d*sin(5*ss_2/N))
}

generate_equal_points = function(number){
  sequence = c((1:number)/number - 1/2/number)
  return(sequence)
}

cube_mean22 = function(Y, m, n) {
  M = dim(Y)[1]
  N = dim(Y)[2]
  Y.new = array(0, dim = c(m, n))
  Y.new2 = array(0, dim = c(m, n))
  b1 = floor(M/m)
  b2 = floor(N/n)
  for (mm in 1:m) {
    for (nn in 1:n) {
      Y.new[mm, nn] = mean(Y[((mm - 1) * b1 + 1):(mm * b1), ((nn - 1) * b2 + 1):(nn * b2)])
      Y.new2[mm, nn] = sd(Y[((mm - 1) * b1 + 1):(mm * b1), ((nn - 1) * b2 + 1):(nn * b2)])^2
    }
  }
  result = list(Y.new = Y.new, Y.var = Y.new2, va = mean(Y.new2))
}

for (iter in 1:iterations){
  
  set.seed(seed_vector[iter])
  omega1a = array(dim = c(M, N))
  omega2a = array(dim = c(M))
  omega3a = array(dim = c(M))
  xi2a = array(dim = c(nn, N))
  xi3a = array(dim = c(nn, N))
  aa <- rnorm(n = nn, mean = 1, sd = 1.3)
  bb <- rnorm(n = nn, mean = 1.5, sd = 1.25)
  cc <- rnorm(n = nn, mean = 2, sd = 1.5)
  
  for (m1 in 1:M)
    for (n1 in 1:N){
      
      if(omega_check == 0){omega1a[m1, n1] = 0}
      if(omega_check == 1){omega1a[m1, n1] = omega1(m1, n1, M, N)}
      omega2a[m1] = omega2(m1, M)
      omega3a[m1] = omega3(m1, M)
      for (i in 1:nn){
        xi2a[i, n1] =  xi3(n1, bb[i], N)
        xi3a[i, n1] =  xi3(n1, cc[i], N)
      }
      
    }
  
  
  sim_images = array(0, dim = c(nn, M, N))
  for (i in 1:nn){
    sim_images[i, , ] = aa[i]*omega1a
    for (m1 in 1:M)
      for (n1 in 1:N){
        sim_images[i,m1 ,n1 ] = sim_images[i,m1 ,n1 ] + omega2a[m1]*xi2a[i, n1] + omega3a[m1]*xi3a[i, n1]
      }
  }
  
  noise <- array(rnorm(nn*M*N, mean = 0, sd = noise_level), dim=c(nn,M,N))
  sim_images = sim_images + noise
  
  par(mar=c(0.5,0.5,0.5,0.5))
  par(mfrow = c(1, 10))
  for(rr in 1:(10)){
    image(as.matrix(sim_images[rr, , ]),axes = FALSE)
  }
  #remove mean
  
  mu_st = array(0, c(M, N))
  
  for (d1 in 1:M)
    for (d2 in 1:N){
      #take the mean of these voxels to be mu_st
      mu_st[d1, d2] = mean(sim_images[, d1, d2])      
    }
  
  for (i in 1:nn){
    sim_images[i, , ] = sim_images[i,,] - mu_st
  }
  
  
  #marginal process
  mar_proc = array(dim = c(nn, M))
  
  for (i in 1:nn)
    for(m1 in 1:M){
      mar_proc[i,m1] = sum(sim_images[i, m1,])/M
    }
  
  cov_mar_proc = array(dim = c(M, M))
  for(m1 in 1:M)
    for (m2 in 1:M){
      cov_mar_proc[m1,m2] = cov(mar_proc[,m1],mar_proc[,m2])
    }
  
  # normal covariance
  cov_func = array(dim = c(M, N, M, N))
  for(m1 in 1:M){
    print(m1)
    for(m2 in 1:M)
      for(n1 in 1:N)
        for(n2 in 1:N){
          cov_func[m1, n1, m2, n2] = cov(sim_images[ , m1, n1],sim_images[ , m2, n2])
        }
  }
  
  
  #integrate the covariance
  mar_cov = array(dim = c(M, M))
  for(m1 in 1:M)
    for (m2 in 1:M){
      temp_mar = 0
      for (n1 in 1:N){
        temp_mar = temp_mar + cov_func[m1, n1, m2, n1]
      }
      mar_cov[m1, m2] = temp_mar/N
    }
  
  
  #marginal process covariance
  mar_proc_eigen = eigen(cov_mar_proc)
  
  #Plot the scores for everyone by component
  # par(mar=c(2,2,2,2))
  # par(mfrow = c(2,5))
  # for(j in 1:5){
  #   name <- paste('Marginal Process PC ', j ,  sep = "", collapse=NULL)
  #   plot(c(1:M),(mar_proc_eigen$vectors[,j]), type='l', pch = 19,)
  #   title(name)
  # }
  

  #marginal covariance
  mar_cov_eigen = eigen(mar_cov)
  
  #Plot the scores for everyone by component
  
  # for(j in 1:5){
  #   name <- paste('Marginal Covariance PC ', j ,  sep = "", collapse=NULL)
  #   plot(c(1:M),(mar_cov_eigen$vectors[,j]), type='l', pch = 19,)
  #   title(name)
  # }
  # 

  sum_val_proc = rep(NA, length(mar_proc_eigen$values))
  sum_val_cov = rep(NA, length(mar_cov_eigen$values))
  
  for (idx in 1:(length(mar_proc_eigen$values))){
    sum_val_proc[idx] = sum(mar_proc_eigen$values[1:idx])/sum(mar_proc_eigen$values)
    sum_val_cov[idx] = sum(mar_cov_eigen$values[1:idx])/sum(mar_cov_eigen$values)
  }
  
  diff_val_proc = rep(NA, length(mar_proc_eigen$values)-1)
  diff_val_cov = rep(NA, length(mar_cov_eigen$values)-1)
  
  for (idx in 1:(length(diff_val_cov))){
    diff_val_proc[idx] = sum_val_proc[idx+1]-sum_val_proc[idx]
    diff_val_cov[idx] = sum_val_cov[idx+1]-sum_val_cov[idx]
  }
  
  r_proc = length(which(diff_val_proc>0.01))+1
  r_cov = length(which(diff_val_cov>0.01))+1
  rproc_iter[iter] = r_proc
  rcov_iter[iter] = r_cov
  
  eval_points = generate_equal_points(N)
  bspline_basis = fda::create.bspline.basis(0:1, nbasis = 12, norder=3)
  
  choice_bs = bspline_basis
  
  discrete_basis = fda::eval.basis(eval_points, choice_bs)
  
  #rows go over time, columns go over the basis vectors
  basis_all = discrete_basis
  
  #set up an empty matrix for all the scores for all the people
  scores_process = array(0, c(nn, r_proc, nrow(basis_all)))
  scores_covar = array(0, c(nn, r_cov, nrow(basis_all)))
  
  b_coefs_process = array(0, c(r_proc, ncol(basis_all)))
  b_coefs_covar = array(0, c(r_cov, ncol(basis_all)))
  
  #### marginal process
  
  choice_eF = mar_proc_eigen$vectors
  
  XX = array(0, c(M, N))
  
  for (subject_nr in 1:nn){
    
    for (iscore in 1:r_proc){
      
      for (idxx in 1:ncol(basis_all)){
        #project the principal components onto a time basis - add a 4th dimension
        for (idx in 1:nrow(basis_all)) {
          XX[ , idx] = basis_all[idx, idxx]*choice_eF[, iscore]
        }
        XX = XX/nrow(basis_all)
        #calculate the b coefficient dependant on the person
        FF = sum(XX*XX)^(-1)
        b_coefs_process[iscore, idxx] = FF*(sum(sim_images[subject_nr, , ] * XX)/(nrow(basis_all)))
      }
    }
    
    score_func_ij = b_coefs_process %*% t(basis_all)
    scores_process[subject_nr, ,] = score_func_ij
  }
  
  ##### marginal covariances scores
  
  choice_eF = mar_cov_eigen$vector
  
  XX = array(0, c(M, N))
  
  for (subject_nr in 1:nn){
    
    for (iscore in 1:r_cov){
      
      for (idxx in 1:ncol(basis_all)){
        #project the principal components onto a time basis - add a 4th dimension
        for (idx in 1:nrow(basis_all)) {
          XX[ , idx] = basis_all[idx, idxx]*choice_eF[, iscore]
        }
        XX = XX/nrow(basis_all)
        #calculate the b coefficient dependant on the person
        FF = sum(XX*XX)^(-1)
        b_coefs_covar[iscore, idxx] = FF*(sum(sim_images[subject_nr, , ] * XX)/(nrow(basis_all)))
      }
    }
    
    score_func_ij = b_coefs_covar %*% t(basis_all)
    scores_covar[subject_nr, ,] = score_func_ij
  }
  
  # #### JUST A TEST
  # 
  # scores_process = array(0, c(nn, r_proc, N))
  # choice_eF = mar_proc_eigen$vectors
  # 
  # for (subject_nr in 1:nn){
  #   print(subject_nr)
  #   
  #   for (iscore in 1:r_proc){
  #     
  #     XX = as.vector(choice_eF[, iscore])
  #     FF = solve(t(XX)%*%(XX))
  #     for (n1 in 1:N){
  #       YY = as.vector(sim_images[subject_nr, ,n1])
  #       scores_process[subject_nr, iscore, n1] = (FF)*(XX%*%YY)
  #       
  #     }
  #   }
  # }
  # 
  # 
  # scores_covar = array(0, c(nn, r_cov, N))
  # choice_eF = mar_cov_eigen$vector
  # 
  # for (subject_nr in 1:nn){
  #   print(subject_nr)
  #   
  #   for (iscore in 1:r_cov){
  #     
  #     XX = as.vector(choice_eF[, iscore])
  #     FF = solve(t(XX)%*%(XX))
  #     for (n1 in 1:N){
  #       YY = as.vector(sim_images[subject_nr, ,n1])
  #       scores_covar[subject_nr, iscore, n1] = (FF)*(XX%*%YY)
  #       
  #     }
  #   }
  # }
  # #### END OF TEST
  
  # plot_r = max(r_proc, r_cov)
  # par(mar=c(2,2,2,2))
  # par(mfrow = c(2, plot_r))
  # for(j in 1:r_cov){
  #   name <- paste('Marginal Covariance Score ', j ,  sep = "", collapse=NULL)
  #   matplot(c(1:M),t(scores_covar[1:50,j,]), type='l', pch = 19,)
  #   title(name)
  # }
  # 
  # for(j in 1:r_proc){
  #   name <- paste('Marginal Process Score ', j ,  sep = "", collapse=NULL)
  #   matplot(c(1:M),t(scores_process[1:50,j,]), type='l', pch = 19,)
  #   title(name)
  # }
  
  ### METHOD 3
  
  #full PCA with two dims PCs
  m = floor(M/2)
  n = floor(N/2)
  maraYall = array(NA, c(nn, m, n))
  maraYva = array(NA, c(nn))
  
  #calculate the cubic mean 
  for (i in 1:nn){
    sresult = cube_mean22(sim_images[i, , ], m, n)
    maraYall[i, , ] = sresult$Y.new
    maraYva[i] = sresult$va
  }
  
  
  #V_X over time using the full image
  LAY = matrix(0, nn, nn)
  for (i in 1:nn){
    for (j in i:nn){
      LAY[i, j] = sum(maraYall[i, , ]*maraYall[j, ,])/m/n
    }
  }
  
  LAY = LAY + t(LAY)
  diag(LAY) = diag(LAY)/2
  
  eLAY = eigen(LAY - diag(maraYva))
  r_PCA = length(which(eLAY$values >=0)) 
  
  rpca_iter[iter]=r_PCA
  #eA = t(eLAY$vectors[, 1:r])  #get loadings approach 1
  eA=t(eLAY$vectors[,1:r_PCA] %*% diag(sqrt(eLAY$values[1:r_PCA]))) #get loadings
  
  phi_j_PCA = array(0, c(r_PCA+1,M,N))
  
  vec_score = array(1, c(nn, (r_PCA+1)))
  for(j in 1:r_PCA){
    vec_score[,(j+1)] = eA[j, ]
  }
  XtXinverse = t(vec_score)%*%vec_score
  
  for (i in 1:M)
    for (j in 1:N){
      vec_Y = as.vector(sim_images[,i,j])
      
      XtY = t(vec_score)%*%(vec_Y)
      beta = solve(XtXinverse)%*%XtY
      for (rr in 1:(r_PCA+1)){
        phi_j_PCA[rr, i, j] = beta[rr,]
      }
    }
  
  par(mar=c(0.5,0.5,0.5,0.5))
  par(mfrow = c(1, r_PCA+1))
  for(rr in 1:(r_PCA+1)){
    image(as.matrix(phi_j_PCA[rr, ,  ]),axes = FALSE)
  }
  
  eF_PCA = phi_j_PCA[2:(r_PCA+1), ,  ]
  
  scores_PCA = array(0, dim=c(nn, r_PCA))
  
  vec_PC = array(1, c(M*N, (r_PCA+1)))
  for(j in 1:r_PCA){
    vec_PC[,(j+1)] = as.vector(eF_PCA[j, ,])
  }
  
  XtXinverse = t(vec_PC)%*%vec_PC
  
  for (i in 1:nn){
    vec_Y = as.vector(sim_images[i, , ])
    
    XtY = t(vec_PC)%*%(vec_Y)
    beta = solve(XtXinverse)%*%XtY
    scores_PCA[i, ] = beta[2:length(beta),1]
  }
  
  
  ### RECONSTRUCTION 
  
  ## PROCESS
  recon_mar_proc = array(0, dim = c(nn, M, N))
  explained_variance = rep(NA, r_proc)
  mean_squared_error = rep(NA, r_proc)
  variance_full =  rep(NA, nn)
  
  for (i in 1:nn){
    variance_full[i] =  sum(sim_images[i,,]^2)
  }
  
  for (rr in 1:r_proc){
    for (i in 1:nn){
      for (n1 in 1:N){
        #complete reconstruction
        recon_mar_proc[i, , n1] = recon_mar_proc[i, , n1] + mar_proc_eigen$vector[, rr]*scores_process[i, rr, n1]
        
      }
    }
  }
  
  optimal_c = rep(0, nn)
  for (i in 1:nn){
    optimised_recon <- function(varc){-sum((sim_images[i,,] - varc*recon_mar_proc[i, , ])^2)}
    optimal_c[i] = optimise(optimised_recon, c(0,1), maximum=TRUE)$maximum
  }
  
  recon_mar_proc = array(0, dim = c(nn, M, N))
  for (rr in 1:r_proc){
    difference_i = rep(0, nn)
    sum_of_var = rep(0, nn)
    for (i in 1:nn){
      for (n1 in 1:N){
        #complete reconstruction
        recon_mar_proc[i, , n1] = recon_mar_proc[i, , n1] + mar_proc_eigen$vector[, rr]*scores_process[i, rr, n1]
        
      }
      difference_i[i] = sum((sim_images[i,,] - optimal_c[i]*recon_mar_proc[i, , ])^2)
      sum_of_var[i] = 1 - sum(difference_i[i]/variance_full[i])
    }
    
    explained_variance[rr] =  mean(sum_of_var)
    mean_squared_error[rr] = mean(difference_i)
  }
  
  if(r_proc<=5){
    MSE_proc_iter[iter, 1:r_proc] = mean_squared_error
    VE_proc_iter[iter, 1:r_proc] = explained_variance
  }
  
  if(r_proc>5){
    MSE_proc_iter[iter, 1:4] = mean_squared_error[1:4]
    VE_proc_iter[iter, 1:4] = explained_variance[1:4]
    MSE_proc_iter[iter, 5] = mean_squared_error[r_proc]
    VE_proc_iter[iter, 5] = explained_variance[r_proc]
  }
  
  ## COVARIANCE
  recon_mar_cov = array(0, dim = c(nn, M, N))
  explained_variance = rep(NA, r_cov)
  mean_squared_error = rep(NA, r_cov)

  for (rr in 1:r_cov){
    difference_i = rep(0, nn)
    sum_of_var = rep(0, nn)
    for (i in 1:nn){
      for (n1 in 1:N){
        #complete reconstruction
        recon_mar_cov[i, , n1] = recon_mar_cov[i, , n1] + mar_cov_eigen$vector[, rr]*scores_covar[i, rr, n1]
        
      }
    }
  }
  
  optimal_c = rep(0, nn)
  for (i in 1:nn){
    optimised_recon <- function(varc){-sum((sim_images[i,,] - varc*recon_mar_cov[i, , ])^2)}
    optimal_c[i] = optimise(optimised_recon, c(0,1), maximum=TRUE)$maximum
  }
  
  recon_mar_cov = array(0, dim = c(nn, M, N))
  for (rr in 1:r_cov){
    difference_i = rep(0, nn)
    sum_of_var = rep(0, nn)
    for (i in 1:nn){
      for (n1 in 1:N){
        #complete reconstruction
        recon_mar_cov[i, , n1] = recon_mar_cov[i, , n1] + mar_cov_eigen$vector[, rr]*scores_covar[i, rr, n1]
        
      }
      difference_i[i] = sum((sim_images[i,,] - optimal_c[i]*recon_mar_cov[i, , ])^2)
      sum_of_var[i] = 1 - sum(difference_i[i]/variance_full[i])
    }
    explained_variance[rr] =  mean(sum_of_var)
    mean_squared_error[rr] = mean(difference_i)
    
  }
  
  if(r_cov<=5){
    MSE_cov_iter[iter, 1:r_cov] = mean_squared_error
    VE_cov_iter[iter, 1:r_cov] = explained_variance
  }
  
  if(r_cov>5){
    MSE_cov_iter[iter, 1:4] = mean_squared_error[1:4]
    VE_cov_iter[iter, 1:4] = explained_variance[1:4]
    MSE_cov_iter[iter, 5] = mean_squared_error[r_cov]
    VE_cov_iter[iter, 5] = explained_variance[r_cov]
  }

  
  
  ## PCA OLD SCHOOL
  
  recon_PCA = array(0, dim = c(nn, M, N))
  explained_variance = rep(NA, r_PCA)
  mean_squared_error = rep(NA, r_PCA)
  
  for (rr in 1:r_PCA){
    difference_i = rep(0, nn)
    sum_of_var = rep(0, nn)
    for (i in 1:nn){
      #complete reconstruction
      recon_PCA[i, , ] = recon_PCA[i, , ] + eF_PCA[rr, , ]*scores_PCA[i,rr]
      
      
      difference_i[i] = sum((sim_images[i,,] - recon_PCA[i, , ])^2)
      sum_of_var[i] = 1 - sum(difference_i[i]/variance_full[i])
    }
    
    if(rr<5){
      MSE_pca_iter[iter, rr] = mean(difference_i)
      VE_pca_iter[iter, rr] = mean(sum_of_var)
    }

    if(rr==r_PCA){
      MSE_pca_iter[iter, rr] = mean(difference_i)
      VE_pca_iter[iter, rr] = mean(sum_of_var)
    }
  }
  
  
}


result_df = data.frame(r_proc = rproc_iter, 
                       r_cov = rcov_iter, 
                       r_PCA = rpca_iter,
                       VE1_proc = VE_proc_iter[,1],
                       VE2_proc = VE_proc_iter[,2],
                       VE3_proc = VE_proc_iter[,3],
                       VE4_proc = VE_proc_iter[,4],
                       VE5_proc = VE_proc_iter[,5],
                       VE1_cov = VE_cov_iter[,1],
                       VE2_cov = VE_cov_iter[,2],
                       VE3_cov = VE_cov_iter[,3],
                       VE4_cov = VE_cov_iter[,4],
                       VE5_cov = VE_cov_iter[,5],
                       VE1_pca = VE_pca_iter[,1],
                       VE2_pca = VE_pca_iter[,2],
                       VE3_pca = VE_pca_iter[,3],
                       VE4_pca = VE_pca_iter[,4],
                       VE5_pca = VE_pca_iter[,5],
                       MSE1_proc = MSE_proc_iter[,1],
                       MSE2_proc = MSE_proc_iter[,2],
                       MSE3_proc = MSE_proc_iter[,3],
                       MSE4_proc = MSE_proc_iter[,4],
                       MSE5_proc = MSE_proc_iter[,5],
                       MSE1_cov = MSE_cov_iter[,1],
                       MSE2_cov = MSE_cov_iter[,2],
                       MSE3_cov = MSE_cov_iter[,3],
                       MSE4_cov = MSE_cov_iter[,4],
                       MSE5_cov = MSE_cov_iter[,5],
                       MSE1_pca = MSE_pca_iter[,1],
                       MSE2_pca = MSE_pca_iter[,2],
                       MSE3_pca = MSE_pca_iter[,3],
                       MSE4_pca = MSE_pca_iter[,4],
                       MSE5_pca = MSE_pca_iter[,5])

filename = paste('omega', omega_check, '_sim_nn', nn ,'_noise' , noise_level_name, '.csv',   sep = "", collapse=NULL)
write.csv(result_df, filename)
