# ==== SET UP ====


#variables consistent over the simulation

#dimensions of images and PCs
M = 31
N = 31
L = 31

nn <- 100 # no. subjects
TT <- 20 # no. time points for each subject

#decide the number of smaller cubes to divide the image into
m = 15
n = 15
l = 15

#number of iterations
iterations = 32

#noise level in the images
img_sigma = 0.11


#set up variables to remember

true_outcomes = array(0, c(iterations, nn))
estimated_outcomes = array(0, c(iterations, nn))

#RMSE for outcomes
RMSE_outcomes = array(0, c(iterations))

#MISE (mean integrated squared error) for PCs, scores and beta
MISE_pc = array(0, c(iterations,3))
MISE_score = array(0, c(iterations,3))
MISE_beta = array(0, c(iterations,2))

#every beta estimated during the process
estimated_beta = array(0, c(iterations,2, TT))

#reconstruction error as mean integrated integrated squared error
recon_error = array(0, c(iterations, 3))
recon_error_clean = array(0, c(iterations, 3))

#variance explained true and estimated
VE_true = array(0, c(iterations, 3))
VE_estimated = array(0, c(iterations, 3))
VE_est_clean = array(0, c(iterations, 3))


# ==== FUNCTION  TO GENERATE DATA ====

## PCs
pphi_1_f <- function(ss_1, ss_2, ss_3)
{
  return(sqrt(2)*cos(2*pi*ss_1/M)*sqrt(2)*cos(2*pi*ss_2/M)*sqrt(2)*cos(2*pi*ss_3/M))
}

pphi_2_f <- function(ss_1, ss_2, ss_3)
{
  return(sqrt(2)*sin(2*pi*ss_1/M)*sqrt(2)*sin(2*pi*ss_2/M)*sqrt(2)*sin(2*pi*ss_3/M))
}

#Scores
ppsi_1_f <- function(aa, bb, tt)
{
  return(aa*cos(bb*pi*tt/TT))
}

ppsi_2_f <- function(cc, dd, tt)
{
  return(cc*sin(dd*pi*tt/TT))
}

#Beta
beta1 <- function(tt)
{
  return(1.25*cos(tt*pi/20))
}

beta2 <- function(tt){
  return(2*sin(1.25*tt*pi/20))
}

# ==== ITERATIONS ====

for (iter in 1:iterations){
  
  print(paste('iteration ', iter, sep=''))
  #GENERATE DATA
  #generate PCs
  pphi_1_array <- array(dim = c(M,N,L))
  for (i in 1:M)
    for (j in 1:N)
      for (k in 1:L)
      {pphi_1_array[i,j,k] <- sqrt(2)*pphi_1_f(i,j,k)}
  
  pphi_2_array <- array(dim = c(M,N,L))
  for (i in 1:M)
    for (j in 1:N)
      for (k in 1:L)
      {pphi_2_array[i,j,k] <- sqrt(2)*pphi_2_f(i,j,k)}
  
  #check it meets assumptions
  if(sum(pphi_1_array*pphi_2_array)/31^3 > 0.01) stop('PCs not orthogonal')
  if(1-sum(pphi_1_array^2/(31^3)) > 0.01) stop('PC 1 not normal')
  if(1-sum(pphi_1_array^2/(31^3)) > 0.01) stop('PC 2 not normal')
  
  
  #generate scores
  aa <- rnorm(n = nn, sd = sqrt(2))
  bb <- rnorm(n = nn, mean = 0.85, sd = sqrt(0.25))
  cc <- rnorm(n = nn, sd = sqrt(0.55))
  dd <- rnorm(n = nn, mean = 1, sd = 0.5)
  
  ppsi_1_matrix <- matrix(nrow = nn, ncol = TT)
  for (i in 1:nn)
    for (j in 1:TT)
    {
      ppsi_1_matrix[i, j] <- ppsi_1_f(aa[i], bb[i], j)
    }
  
  ppsi_2_matrix <- matrix(nrow = nn, ncol = TT)
  for (i in 1:nn)
    for (j in 1:TT)
    {
      ppsi_2_matrix[i, j] <- ppsi_2_f(cc[i], dd[i], j)
    }
  
  #check it meets assumptions
  if(sum(ppsi_1_matrix*ppsi_2_matrix)/31^3 > 0.01) stop('scores not orthogonal')
  
  (VE_1 <- sum(ppsi_1_matrix^2)/(sum(ppsi_1_matrix^2)+sum(ppsi_2_matrix^2)))
  (VE_2 <- sum(ppsi_2_matrix^2)/(sum(ppsi_1_matrix^2)+sum(ppsi_2_matrix^2)))
  
  VE_true[iter, 1] = VE_1
  VE_true[iter, 2] = VE_2
  
  beta_outcome = array(0, c(r, TT))
  for (j in 1:TT){
    beta_outcome[1,j] = beta1(j)
    beta_outcome[2,j] = beta2(j)
  }
  
  outcome = array(0, nn)
  for (i in 1:nn){
    outcome[i] = sum(beta_outcome[1, ]*ppsi_1_matrix[i, ])/TT + sum(beta_outcome[2, ]*ppsi_2_matrix[i, ])/TT
  }
  
  true_outcomes[iter, ] = outcome
  
  outcome_noise = rnorm(nn, sd = 0.01)
  outcome = outcome + outcome_noise
  
  
  #put everything in arrays
  pphi_j = array(0, c(r, M, N, L))
  pphi_j[1, , ,] = pphi_1_array
  pphi_j[2, , ,] = pphi_2_array
  
  ppsi_ij = array(0, c(r, nn, TT))
  ppsi_ij[1, ,] = ppsi_1_matrix
  ppsi_ij[2, ,] = ppsi_2_matrix
  
  #Plot the scores for everyone by component
  par(mar=c(2,2,2,2))
  par(mfrow = c(1,2))
  for(j in 1:2){
    name <- paste('True scores for component ', j ,  sep = "", collapse=NULL)
    matplot(c(1:TT),t(ppsi_ij[j, , ]), type='l', ylim = c(-5, 5))
    title(name)
  }
  
  # reconstructing the images
  sim_images = array(0, c(nn, M, N, L, TT)) 
  sim_images_clean = array(0, c(nn, M, N, L, TT)) 
  
  for(i in 1:nn){
    noise = rnorm(M*N*L*TT, mean = 0, sd = img_sigma)
    noise_array = array(noise, c(M,N,L,TT))
    
    for(t in 1:TT)
      for(rr in 1:r){
        sim_images[i, , , , t] = sim_images[i, , , , t] + ppsi_ij[rr, i, t]*pphi_j[rr, , , ]
      }
    
    sim_images_clean[i, , , , ] = sim_images[i, , , , ]
    sim_images[i, , , , ] = sim_images[i, , , , ] + noise_array
  }
  
  #REMOVE GLOBAL MEAN
  
  mean_st = estimate_mu_st(sim_images)
  
  # normalizing the images
  sYall = array(0, c(nn, M, N, L, TT))
  
  for (i in 1:nn){
    sYall[i, , , ,] = sim_images[i, , , ,] - mean_st
  }
  
  mean_st_clean = estimate_mu_st(sim_images_clean)
  sYall_clean = array(0, c(nn, M, N, L, TT))
  for (i in 1:nn){
    sYall_clean[i, , , ,] = sim_images_clean[i, , , ,] - mean_st_clean
  }
  
  #ESTIMATE SCORES
  #set how many components you want to estimate
  r = 3
  
  saYall = array(0, c(nn, 15, 15, 15, TT))
  saYva = array(0, c(nn, TT))
  
  #calculate the cubic mean 
  for (i in 1:nn)
    for (j in 1:TT){
      # divide the image into 2x2x2 cubes of which there are 15*15*15
      sresult = cube_mean22(sYall[i, , , , j], m, n, l)
      saYall[i, , , , j] = sresult$Y.new
      saYva[i,j] = sresult$va
    }
  
  #V_X over time USING FULL IMG
  sLAY = array(0, c(nn, nn, TT))
  for (t in 1:TT)
    for (i in 1:nn) 
      for (j in i:nn) {
        #create upper triangular matrix
        sLAY[i, j, t] = threeDsum(sYall[i, , , , t] * sYall[j, , , , t])
      }
  
  for (t in 1:TT){
    #flip the upper triangular matrix to make it symmetrical
    sLAY[,,t] = sLAY[,,t] + aperm(sLAY[,,t])
    #divide diagonal in half since it was added twice
    diag(sLAY[,,t]) = diag(sLAY[,,t])/2
    #divide by the total number of pixels in the image as 
    #this is an estimation of integration
    sLAY[,,t] = sLAY[,,t]/M/N/L
  }
  
  sim.load.old = array(0, c(r, nn, TT))
  for (t in 1:TT){
    # eigedecompostioin of V_y - V_sigma for every time point
    seLAY = eigen(sLAY[,,t] - diag(saYva[,t]))
    # take the first 3 vectors to be score functions
    seA = t(seLAY$vectors[, 1:r])
    seP = seLAY$values[1:r]
    for(rr in 1:r){
      if(seP[rr]<0){seP[rr] = 0}
    }
    sim.load.old[,,t] = seA*sqrt(seP)
  }
  
  #if(seP[3]==0){r=2}
  for (rr in 1:r){if(seP[rr]==0)warning('Error with eigenvalues when computing scores')}
  
  for (j in 1:r)
    for (tt in 1:(TT-1)) {
      #find which subject has the maximum valued score at this time
      max_idx = which.max(abs(sim.load.old[j, , tt]))
      max_sign = sign(sim.load.old[j, max_idx, tt])
      
      #make sure the scoers don't flip arbitrarily across time points
      if (sign(sim.load.old[j, max_idx, tt+1]) != max_sign) {
        sim.load.old[j, ,tt+1] = (-1)*sim.load.old[j, ,tt+1]
      }
    }
  
  #Plot the scores for everyone by component
  par(mar=c(2,2,2,2))
  par(mfrow = c(1,r))
  for(j in 1:r){
    name <- paste('Raw scores for component ', j ,  sep = "", collapse=NULL)
    matplot(c(1:TT),t(sim.load.old[j, , ]), type='l', ylim = c(-4.5, 4.5))
    title(name)
  }
  
  #ESTIMATING PRINCIAL COMPONENTS
  #for every pixel create basic regression
  
  phi_j = array(0, c(r+1,M,N,L))
  
  vec_score = array(1, c(nn*TT, (r+1)))
  for(j in 1:r){
    vec_score[,(j+1)] = as.vector(sim.load.old[j, , ])
  }
  
  XtXinverse = t(vec_score)%*%vec_score
  
  for (i in 1:M)
    for (j in 1:N)
      for (k in 1:L){
        vec_Y = as.vector(sYall[ ,i, j, k, ])
        XtY = t(vec_score)%*%(vec_Y)
        beta = solve(XtXinverse)%*%XtY
        for (rr in 1:(r+1)){
          phi_j[rr, i, j, k] = beta[rr,]
        }
      }
  
  par(mar=c(0.5,0.5,0.5,0.5))
  par(mfrow = c(1,r+1)) 
  for(rr in 1:(r+1)){
    image(as.matrix(phi_j[rr, , , 15]),axes = FALSE)
    title(paste('Scale ', 
                signif(sqrt(sum(phi_j[rr, , , ]^2/(31^3))), digits = 2), sep=''))
  }
  
  
  seF = array(0, c(r, M, N, L))
  scaling_scores = array(0, r)
  for (rr in 1:r){
    scaling_scores[rr] = sqrt(sum(phi_j[rr+1, , , ]^2/(31^3))) 
    
    if(scaling_scores[rr]>0.97){
      seF[(rr), , , ] = phi_j[rr+1, , , ]
      scaling = FALSE
    }
    else{
      seF[(rr), , , ] = phi_j[rr+1, , , ]/sqrt(sum(phi_j[rr+1, , , ]^2/(31^3)))
      warning('Estimated PCs are not directly normal')
      scaling = TRUE
    }
  }
  
  if(abs(sum(seF[1,,,]*seF[2,,,])/31^3) > 0.15 ) warning('Esrtimated PCs not orthogonal')
  
  #UPDATING SCORES
  
  eval_points = generate_equal_points(TT)
  bspile_basis = create.bspline.basis(0:1, nbasis = 4, norder=3)
  discrete_basis_bs = eval.basis(eval_points, bspile_basis)
  #rows go over time, columns go over the basis vectors
  basis_all = discrete_basis_bs 
  
  #set up an empty matrix for all the scores for all the people
  sim_scores_updated = array(0, c(nn, r, nrow(basis_all)))
  b_coefs = array(0, c(r, ncol(basis_all)))
  
  XX = array(0, c(M, N, L, TT))
  
  for (subject_nr in 1:nn){
    print(subject_nr)
    
    for (iscore in 1:r){
      
      for (idxx in 1:ncol(basis_all)){
        #project the principal components onto a time basis - add a 4th dimension
        for (idx in 1:nrow(basis_all)) {
          XX[ , , , idx] = basis_all[idx, idxx]*seF[iscore, , , ]
        }
        XX = XX/nrow(basis_all)
        #calculate the b coefficient dependant on the person
        FF = sum(XX*XX)^(-1)
        b_coefs[iscore, idxx] = FF*(sum(sYall[subject_nr, , , ,] * XX[ , , , ])/(nrow(basis_all)))
      }
    }
    
    score_func_ij = b_coefs %*% t(basis_all)
    sim_scores_updated[subject_nr, ,] = score_func_ij
  }

  par(mar=c(2,2,2,2))
  par(mfrow = c(1,r))
  for(j in 1:r){
    name <- paste('Updated Scores ', j ,  sep = "", collapse=NULL)
    matplot(c(1:TT),t(sim_scores_updated[,j,]*scaling_scores[j]), type='l')
    title(name)
  }
  
  #ORDERING EVERYTHING FOR THE SAKE OF PLOTTING
  
  ordered_seF = array(0, c(r, M,N,L))
  ordered_scores = array(0, c(nn, r, TT)) 
  
  difference_pc = array(0, c(4,2))
  
  for (rr in 1:2){
    difference_pc[1,rr] = sum(abs(pphi_j[rr, , , 15] - seF[1, , , 15]))
    difference_pc[2,rr] = sum(abs(pphi_j[rr, , , 15] - seF[2, , , 15]))
    difference_pc[3,rr] = sum(abs(pphi_j[rr, , , 15] + seF[1, , , 15])) 
    difference_pc[4,rr] = sum(abs(pphi_j[rr, , , 15] + seF[2, , , 15])) 
  }
  
  for (rr in 1:2){
    if(which(difference_pc[,rr] == min(difference_pc[,rr]))==1){
      ordered_scores[ , rr, ] = sim_scores_updated[ , 1, ]
      ordered_seF[ rr, , , ] = seF[1, , , ]
    }
    
    if(which(difference_pc[,rr] == min(difference_pc[,rr]))==2){
      ordered_scores[ , rr, ] = sim_scores_updated[ , 2, ]
      ordered_seF[ rr, , , ] = seF[2, , , ]
    }
    
    if(which(difference_pc[,rr] == min(difference_pc[,rr]))==3){
      ordered_scores[ , rr, ] = - sim_scores_updated[ , 1, ]
      ordered_seF[ rr, , , ] = - seF[1, , , ]
    } 
    
    if(which(difference_pc[,rr] == min(difference_pc[,rr]))==4){
      ordered_scores[ , rr, ] = - sim_scores_updated[ , 2, ]
      ordered_seF[ rr, , , ] = - seF[2, , , ]
    } 
  }
  
  if(scaling==TRUE){ for(rr in 1:r){
    ordered_scores[ , rr, ] =  ordered_scores[ , rr, ]*scaling_scores[rr]}  }
  
  #RECONSTRUCTION 
  
  recon_Y = array(0, c(nn, M, N, L, TT))
  recon_L_Y = array(0, c(r, nn, M, N, L, TT))
  explained_variance = rep(NA, r)
  mean_squared_error = rep(NA, r)
  
  for (i in 1:nn){
    #recon_Y[i, , , , ] = syall_mean[ , , , ]
    for( t in 1:TT){
      for (j in 1:r){
        #complete reconstruction
        recon_Y[i, , , , t] = recon_Y[i, , , , t] + ordered_seF[j, , , ]*ordered_scores[i, j, t]
        
        #reconstruction given L components
        recon_L_Y[j, i, , , , t] = recon_Y[i, , , , t]
      }
    }
  }
  
  optimal_c = rep(0, nn)
  for (i in 1:nn){
    optimised_recon <- function(varc){-sum((sYall[i, , , ,] - varc*recon_L_Y[j, i, , , ,])^2)}
    optimal_c[i] = optimise(optimised_recon, c(0,1), maximum=TRUE)$maximum
  }
  
  
  for (j in 1:r){
    #print(paste('Component ',j, task_idx, sep='' ))
    variance_full = rep(0, nn)
    difference_i = rep(0, nn)
    sum_of_var = rep(0, nn)
    for (i in 1:nn){
      print(i)
      #variance_full[i] =  sum(sim_images[i, , , , ]^2)
      variance_full[i] =  sum(sYall[i, , , , ]^2)
      difference_i[i] = sum((sYall[i, , , ,] - optimal_c[i]*recon_L_Y[j, i, , , ,])^2)
      sum_of_var[i] = sum(1 - difference_i[i]/variance_full[i])
    }
    
    recon_error[iter, j] = mean(difference_i/31^3)
    VE_estimated[iter, j] = mean(sum_of_var)
  }
  
  
  for (j in 1:2){
    #print(paste('Component ',j, task_idx, sep='' ))
    variance_full = rep(0, nn)
    difference_i = rep(0, nn)
    sum_of_var = rep(0, nn)
    for (i in 1:nn){
      print(i)
      #variance_full[i] =  sum(sim_images[i, , , , ]^2)
      variance_full[i] =  sum(sYall_clean[i, , , , ]^2)
      difference_i[i] = sum((sYall_clean[i, , , ,] - optimal_c[i]*recon_L_Y[j, i, , , ,])^2)
      sum_of_var[i] = sum(1 - difference_i[i]/variance_full[i])
    }
    
    recon_error_clean[iter, j] = mean(difference_i/31^3)
    VE_est_clean[iter, j] = mean(sum_of_var)
  }
  
  
  #CALCULATE MISE FOR PCS AND SCORES
  for(rr in 1:2){
    difference_pc = sqrt((ordered_seF[rr, , , ] - pphi_j[rr, , , ])^2)
    MISE_pc[iter, rr] = sum(difference_pc)/M/N/L
    
    difference_score = sqrt((optimal_c*ordered_scores[, rr, ] - ppsi_ij[rr, , ])^2)
    MISE_score[iter, rr] = mean(apply(difference_score, 1, sum)/TT)
  }
  
  #OUTCOME PREDICTION
  
  score1 <- as.matrix(optimal_c*ordered_scores[,1,])
  score2 <- as.matrix(optimal_c*ordered_scores[,2,])
  
  col_names <- c('col1','col2','col3','col4','col5','col6', 'col7', 'col8','col9',
                 'col10','col11','col12','col13', 'col14','col15', 'col16','col17',
                 'col18','col19', 'col20')
  
  names_list <- vector(mode='list', length = nn)
  for (i in 1:nn){
    subject_numbers <- paste('subj', i ,  sep = "", collapse=NULL)
    names_list[[i]] <- subject_numbers
  }
  
  dimnames(score1)[[1]] <- names_list
  dimnames(score1)[[2]] <- col_names
  dimnames(score2)[[1]] <- names_list
  dimnames(score2)[[2]] <- col_names
  
  #fit the model beta(t)psi(t)
  par(mar=c(2,2,2,2))
  par(mfrow = c(1,2))
  fit.lf <- pfr(outcome ~ lf(score1, k=3)+lf(score2, k=3))
  plot(fit.lf)
  
  pred_fits <- as.data.frame(fit.lf$fitted.values)
  
  RMSE_outcomes[iter] = mean(sqrt((outcome-pred_fits)^2)[,1])
  
  fit_sc1 <- pfr(outcome ~ lf(score1, k=4))
  fit_sc2 <- pfr(outcome ~ lf(score2, k=4))
  
  estimated_beta[iter, 1, ] = coef(fit_sc1)$value
  estimated_beta[iter, 2, ] = coef(fit_sc2)$value
  
  true_outcomes[iter, ] = outcome
  estimated_outcomes[iter, ] = pred_fits[,1]
  
  
  
  #plotting principal components 
  difference = pphi_j[1:2, , , ] - ordered_seF[1:2, , , ]
  mmin <- min(c(pphi_j[, , , 15], seF[, , , 15], difference[, , , 15]))
  mmax <- max(c(as.matrix(pphi_j[, , , 15]), as.matrix(ordered_seF[, , , 15]), as.matrix(difference[ , , , 15])))
  
  plot_r = 2
  
  par(mar=c(2,2,2,2))
  par(mfrow = c(plot_r,3))
  for(rr in 1:plot_r){
    name <- paste('True PC ', rr ,  sep = "", collapse=NULL)
    image(as.matrix(pphi_j[rr, , , 15]),axes = FALSE)
    title(name)
    name <- paste('Est. PC ', rr ,  sep = "", collapse=NULL)
    image(as.matrix(ordered_seF[rr, , , 15]),axes = FALSE)
    title(name)
    name <- paste('Error ', rr ,  sep = "", collapse=NULL)
    image(as.matrix(difference[rr , , , 15]),axes = FALSE,  zlim=c(mmin, mmax))
    title(name)
  }
  
  par(mar=c(2,2,2,2))
  par(mfrow = c(1,3))
  for(rr in 1:2){
    name <- paste('Est. PC ', rr ,  sep = "", collapse=NULL)
    image(as.matrix(ordered_seF[rr, , , 15]),axes = FALSE)
    title(name)
  }
  name <- paste('Est. PC ', 3 ,  sep = "", collapse=NULL)
  image(as.matrix(seF[3, , , 15]),axes = FALSE)
  title(name)
  
  par(mar=c(2,2,2,2))
  par(mfrow = c(2,2))
  
  for(j in 1:2){
    mmin <- min(c(ppsi_ij[j,,], t(ordered_scores[,j,] )))
    mmax <- max(c(ppsi_ij[j,,], t(ordered_scores[,j,] )))
    
    name <- paste('True Score ', j ,  sep = "", collapse=NULL)
    matplot(c(1:TT),t(ppsi_ij[j, ,] ), type='l', ylim = c(mmin, mmax))
    title(name)
    
    name <- paste('Estimated Score ', j ,  sep = "", collapse=NULL)
    matplot(c(1:TT), t(optimal_c*ordered_scores[,j,]), type='l', ylim = c(mmin, mmax))
    title(name)
  }
  
  par(mar=c(2,2,2,2))
  par(mfrow = c(1,3))
  mmin <- min(ordered_scores[,,] )
  mmax <- max(ordered_scores[,,] )
  for(j in 1:2){
    name <- paste('Estimated Score ', j ,  sep = "", collapse=NULL)
    matplot(c(1:TT), t(optimal_c*ordered_scores[,j,]), type='l', ylim = c(mmin, mmax))
    title(name)
  }
  name <- paste('Estimated Score ', 3 ,  sep = "", collapse=NULL)
  matplot(c(1:TT), t(optimal_c*sim_scores_updated[,3,]), type='l', ylim = c(mmin, mmax))
  title(name)
  
  #end of iteration
}

df_simulation_high_noise <- data.frame(recon_error_PC1 = recon_error[,1],
                             recon_error_PC2 = recon_error[,2],
                             True_VE_PC1 = VE_true[,1],
                             True_VE_PC2 = VE_true[,2],
                             VE_PC1 = VE_estimated[,1],
                             VE_PC2 = VE_estimated[,2],
                             VE_PC1_clean = VE_est_clean[,1],
                             VE_PC2_clean = VE_est_clean[,2],
                             MISE_score1 = MISE_score[,1],
                             MISE_score2 = MISE_score[,2],
                             MISE_PC1 = MISE_pc[,1],
                             MISE_PC2 = MISE_pc[,2])

write.csv(ComplexSimulation_LowNoise,'D:\\STFPCA_ComplexSimulation_LowNoiseFull', row.names=FALSE)

df_simulation_high_noise <- data.frame(recon_error_PC1 = recon_error[,1],
                             recon_error_PC2 = recon_error[,2],
                             True_VE_PC1 = VE_true[,1],
                             True_VE_PC2 = VE_true[,2],
                             VE_PC1 = VE_estimated[,1],
                             VE_PC2 = VE_estimated[,2],
                             MISE_score1 = MISE_score[,1],
                             MISE_score2 = MISE_score[,2],
                             MISE_PC1 = MISE_pc[,1],
                             MISE_PC2 = MISE_pc[,2],
                             RMSE_outcome = RMSE_outcomes)


df_outcomes_high_noise <- data.frame(Beta1 = estimated_beta[ ,1,],
                          Beta2 = estimated_beta[ ,2,],
                          True_outcome = true_outcomes,
                          Est_outcome = estimated_outcomes)

write.csv(df_simulation_high_noise, 'D:\\STFPCA_ComplexSimulation_HighNoise_AvgResults', row.names=FALSE)
write.csv(df_outcomes_high_noise, 'D:\\STFPCA_ComplexSimulation_HighNoise_Outcomes', row.names=FALSE)
