
# =========== IMPORT LIBRARIES ===========

libraries = c("R.matlab", "NMOF", "splines", "tictoc", "matrixStats", "fda")
lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x)
})
lapply(libraries, library, quietly = TRUE, character.only = TRUE)

install.packages("BiocManager")
library("BiocManager")
BiocManager::install("rhdf5")
# =========== SET UP FUNCTIONS ===========

cube_mean2 = function(Y, m, n, l) {
  M = dim(Y)[1]
  N = dim(Y)[2]
  L = dim(Y)[3]
  Y.new = array(0, dim = c(m, n, l))
  Y.new2 = array(0, dim = c(m, n, l))
  b1 = floor(M/m)
  b2 = floor(N/n)
  b3 = floor(L/l)
  for (mm in 1:m) {
    for (nn in 1:n) {
      for (ll in 1:l) {
        Y.new[mm, nn, ll] = mean(Y[((mm - 1) * b1 + 1):(mm * b1), ((nn - 
                                                                      1) * b2 + 1):(nn * b2), ((ll - 1) * b3 + 1):(ll * b3)])
        Y.new2[mm, nn, ll] = var(Y[((mm - 1) * b1 + 1):(mm * b1), ((nn - 
                                                                      1) * b2 + 1):(nn * b2), ((ll - 1) * b3 + 1):(ll * b3)])
      }
    }
  }
  result = list(Y.new = Y.new, va = mean(Y.new2))
}

genddd = function(N, m, n, l) {
  mn = m * n
  d3 = floor(N/mn) * (N%%mn == 0) + (floor(N/mn) + 1) * (N%%mn != 0)
  temp = N - (d3 - 1) * mn
  d2 = floor(temp/m) * (temp%%mn == 0) + floor(temp/m + 1) * (temp%%mn != 0)
  d1 = temp - (d2 - 1) * m
  d = c(d1, d2, d3)
  return(d)
}

genB = function(n, K, p) {
  B = bs(c((1:n)/n - 1/2/n), df = K + p, degree = p)
  return(B)
}

genP = function(m, K, p) {
  ms = choose(m, 0:(m)) * (-1)^c(0:(m))
  D_m = matrix(0, K + p - m, K + p)
  for (i in 1:(K + p - m)) {
    D_m[i, i:(i + m)] = ms
  }
  P = t(D_m) %*% D_m
  return(P)
}

gentD = function(B, P) {
  eB = eigen(t(B) %*% B)
  Bhalf = eB$vectors %*% diag(1/sqrt(eB$values)) %*% t(eB$vectors)
  eBP = eigen(Bhalf %*% P %*% t(Bhalf))
  D = eBP$values
  D = D * (D > 1e-05)
  V = eBP$vectors
  temp = t(V) %*% Bhalf %*% t(B)
  result = list(D = D, temp = temp)
  return(result)
}

genlamD = function(lambda, D) {
  lamD = 1/(1 + lambda * D)
  return(lamD)
}

GCV = function(lambda, sY2, dim, D1, D2, D3, tildeY2) {
  lamD1 = genlamD(lambda[1], D1)
  lamD2 = genlamD(lambda[2], D2)
  lamD3 = genlamD(lambda[3], D3)
  
  W = outer(outer(lamD1, lamD2), lamD3)
  WtildeY2 = W * tildeY2
  W2tildeY2 = W^2 * tildeY2
  
  GCVnum = sY2 - 2 * sum(WtildeY2) + sum(W2tildeY2)
  GCVnum/(dim[1] * dim[2] * dim[3] - sum(lamD1) * sum(lamD2) * sum(lamD3))
}

threeDsum = function(array) {
  return(sum(rowSums(colSums(array))))
}




# =========== RETREIVE THE DATA ===========

#subject ID numebrs
subjects = c(1,3,4,5,6,8,9,10,11,12,15,16,17,19,21) #2 & 18 omitted as they couldn't load

# number of instances
I = 15
# number of time points for one task
tasks = 81
ttime = floor(1360/81) #16

# size of small cube
m = 30
n = 35
l = 30
# size of big cube
M = 91
N = 109
L = 91

#set up empty arrays, one for fresh data, one for cubic mean and one for variance
Yall = array(0, c(I, M, N, L, ttime))
aYall = array(0, c(I, m, n, l, ttime))
aYva = array(0, c(I, ttime))

#image directroy, retrieve all images
#image_dir = "D:/R/Spatial_FPCA/FMRI"
image_dir= "/Volumes/T7/R/Spatial_FPCA/FMRI"

fmri_data_names <- list.files(path = image_dir, pattern = "*.mat" )
fmri_data_names <- fmri_data_names[-10][-8]

k=1
task_idx=1
for (i in 1:length(subjects)){
  print(i)
  path <-  paste(image_dir, fmri_data_names[i], sep = "/", collapse=NULL)
  
  fmri_img <- rhdf5::h5read(path, "/temp")
  img_idx = (task_idx-1)*ttime + 1
  fmri_img <- fmri_img[,,, img_idx:(img_idx+(ttime-1))]
  
  
  Yall[k, , , ,  ] = fmri_img[,,,]
  
  k=k+1
}

# normalizing the images
nYall = array(0, c(I, M, N, L, ttime))
yall_mean = array(0, c(M, N, L, ttime))

for (d1 in 1:M) {
  for (d2 in 1:N) {
    for (d3 in 1:L) {
      for (d4 in 1:ttime) {
        
        yall_mean[d1, d2, d3, d4] = mean(Yall[,d1,d2,d3,d4])
        
      }
    }
  }
}

for (i in 1:I){
  nYall[i, , , ,] = Yall[i, , , ,] - yall_mean[ , , , ]
}


#recalculate the cubic mean etc.
for (i in 1:I){
  for (j in 1:ttime){
    result = cube_mean2(nYall[i, , , , j], m, n, l)
    aYall[i, , , , j] = result$Y.new
    aYva[i,j] = result$va
  }
}


# =========== V_X OVEER TIME ===========

#inner product of the images - takes image by image
#had to adjust the way the sum is calculated as it didn't seem to work otherwise
#this is V_y, over time (3rd dimesnion)
LAY = array(0, c(I, I, ttime))
for (t in 1:ttime){
  for (i in 1:I) {
    for (j in i:I) {
      LAY[i, j, t] = threeDsum(aYall[i, , , , t] * aYall[j, , , , t])
      # print(paste(i,',',j)) 
    }
  }
}

#computing the other side of V_Y to make it symmetrical
#the function 'aperm' is the transpose of 3D array
for (t in 1:ttime){
  LAY[,,t] = LAY[,,t] + aperm(LAY[,,t])
  diag(LAY[,,t]) = diag(LAY[,,t])/2
}

#compute the eigendecomposition of V_y - V_\sigma
#aYva is the sum of all the variation in the 30 x 35 x 30 image
#we need to do this separately by time points
r = 8

load.old = array(0, c(ttime, r, 15))
for (t in 1:ttime){
  eLAY = eigen(LAY[,,t]/m/n/l - diag(aYva[,t]))
  eA = t(eLAY$vectors[, 1:r])  #get loadings approach 1
  load.old[t,,] = eA
}

#Plot the scores for everyone by component
par(mar=c(2,2,2,2))
par(mfrow = c(3,2))
for(j in 1:8){
  name <- paste('Loadings for component ', j ,  sep = "", collapse=NULL)
  matplot(c(1:ttime),load.old[,j,], type='l')
  title(name)
}

# get all the scores,
# find the score that is at its highest or furthest away from zero
# track the score over time, when the sign switches, switch all the signs for that time point

for (j in 1:r) {
  
  #find which subject has the maximum valued score at the beginning
  max_idx = which.max(abs(load.old[1,j,]))
  max_sign = sign(load.old[1,j,max_idx])
  
  #if the score for one subject changes signs, change all the signs
  for (tt in 1:ttime) {
    if (sign(load.old[tt,j,max_idx]) != max_sign) {
      
      load.old[tt,j,] = (-1)*load.old[tt,j,]
      
    }
  }
}

par(mar=c(2,2,2,2))
par(mfrow = c(3,2))
for(j in 1:8){
  name <- paste('Loadings for component ', j ,  sep = "", collapse=NULL)
  matplot(c(1:ttime),load.old[,j,], type='l')
  title(name)
}


# =========== SCORES ===========

# METHOD 1 - COVARIANCE FUNCTION
# load.old are the scores over time, we want to now smooth them into a function
# we have 15 subject, 5 eigenfucntions
# this means 15 x 5 score functions over time

#for each component score & subject, generate covariance function G_ij
Gij = array(0, c(I,r,ttime,ttime))
for (i in 1:I){
  for (j in 1:r){
    for (t in 1:ttime){
      for (k in t:ttime){
        #This is calculating 1x17 covariance for ij denoted Gij[i,j,,]
        Gij[i,j,t,k] = load.old[t,j,i]*load.old[k,j,i] 
        #make it symmetrical
        Gij[i,j,,] = Gij[i,j,,] + aperm(Gij[i,j,,])
        #make diagonal 0
        diag(Gij[i,j,,]) = diag(Gij[i,j,,])*0
        
      }
    }
  }
}


# METHOD 2 - DIRECT B-SPLINE BASIS FITTING

# project each score function into a b-spline basis

library('fda')

e.loadings = array(0, c(I, r, ttime))
#create an initial lambda as a guess
lambda_init = 10^-3
# a vector of values along the x axis on which you have data and between which you wish 
# to create a function
argvalues <- c((1:ttime)/ttime - 1/2/ttime)

for (i in 1:I){
  for (j in 1:r){
    # create a bspline basis on a [0,1] interval with 17 equidistant points
    basisobj <- create.bspline.basis(c(0,1),17)
    
    #create an empty array of coefficients and the penalty matrix based on a pre-determined lambda
    fdParobj <- fdPar(basisobj, lambda = lambda_init)
    
    data = load.old[,j,i]
    
    # find optimal smoothing parameter using the GCV criterion
    # the input is the log od initial lambda, argument values along the x axis,
    # the actual data to be fitted and the empty array of coefficients + smoothing parameter
    lambda = lambda2gcv(log(lambda_init, base=10), argvalues, data, fdParobj, wtvec=rep(1,length(argvalues)))
    
    # based on the GCV criterion, create a new lambda and use this in a new fdPar function
    fdParobj <- fdPar(basisobj, lambda = lambda)
    
    # fit the spline to the data on a range of argument vaues using the fdParobj we defined above
    smoothed.loadings <- smooth.basis(argvalues, data, fdParobj )
    #if you inspect smoothed.loadings you will be able to see the coefficients 
    s.loadings <- smoothed.loadings$fd
    
    #evaluate the spline at the 17 points on [0,1] and store it as an array
    e.loadings[i, j, ] <- eval.fd(argvalues, s.loadings)
  }
}

par(mar=c(2,2,2,2))
par(mfrow = c(3,2))
for(j in 1:5){
  name <- paste('Loading functions for component ', j ,  sep = "", collapse=NULL)
  matplot(c(1:17), t(e.loadings[, j, ]), type='l')
  title(name)
}

# =========== COMPONENTS ===========

# METHOD 1
# calculate principal components at each time point separately, normalise the,
# then take the mean across time to be left with a component over space

#for each time point calculate separate raw pricipal components
tepo_eF.raw = array(0, c(ttime, r, M, N, L))
tempo_eF = array(0, c(ttime, r, M, N, L))

for (t in 1:ttime){
  temporary_scores=e.loadings[,,t]
  
  for (k in 1:r) {
    for (i in 1:I) {
      tepo_eF.raw[t, k, , , ] = tepo_eF.raw[t, k, , , ] + temporary_scores[i, k] * nYall[i, , , ,t]
    }
    tempo_eF[t, k, , , ] = tepo_eF.raw[t, k, , , ]/sqrt(threeDsum(tepo_eF.raw[t, k, , , ]^2)) * sqrt(M * N * L)  #normalize eigenfunctions
  }
  
}

par(mar=c(0.5,0.5,2,0.5))
par(mfrow = c(5,5))
for(time in 11:15){
  for(i in 1:5){
    title = paste('PC(',i,'), time=', time, sep = "", collapse=NULL)
    image(as.matrix(tempo_eF[time, i, , ,55]),axes = FALSE, main=title) #, col = grey(seq(0, 1, length = 256)))
  }
}

eF_2 = array(0, c(r, M, N, L))
for (k in 1:r) {
  for (d1 in 1:M) {
    for (d2 in 1:N) {
      for (d3 in 1:L) {
        
        eF_2[k, d1, d2, d3] = mean(tempo_eF[, k, d1, d2, d3])
        
      }
    }
  }
}

# ===============


# METHOD 2

eF.raw = array(0, c(r, M, N, L))
eF = array(0, c(r, M, N, L))
for (k in 1:r) {
  for (t in ttime){
    par(mfrow = c(3,5))
    for (i in 1:I) {
      # eF.raw[k, , , ] = eF.raw[k, , , ] + e.loadings[i, k, t] * nYall[i, , , ,t]
      eF.raw[k, , , ] = eF.raw[k, , , ] + load.old[t, k, i] * nYall[i, , , ,t]
      
      #plot the current component
      title = paste('PC(', k ,'), time=', t, sep = "", collapse=NULL)
      image(as.matrix(eF.raw[k, , , 55]),axes = FALSE, main=title) #, col = grey(seq(0, 1, length = 256)))
    }
    #normalize eigenfunctions
    eF[k, , , ] = eF.raw[k, , , ]/ttime *sqrt(threeDsum(eF.raw[k, , , ]^2)) * sqrt(M * N * L)  
  }
}


#print the PC components
par(mar=c(0.5,0.5,0.5,0.5))
par(mfrow = c(2,5))
for(i in 1:r){
  image(as.matrix(eF[i, , ,55]),axes = FALSE) #, col = grey(seq(0, 1, length = 256)))
}
for(i in 1:5){
  image(as.matrix(eF_2[i, , ,75]),axes = FALSE) #, col = grey(seq(0, 1, length = 256)))
}


## OLD OLD METHOD

eF.raw = array(0, c(r, M, N, L))
eF = array(0, c(r, M, N, L))
for (k in 1:r) {
  par(mfrow = c(3,5))
  for (i in 1:I) {
    for (t in 1:ttime){
      eF.raw[k, , , ] = eF.raw[k, , , ] + e.loadings[i, k, t] * nYall[i, , , ,t]
      eF.raw[k, , , ] = eF.raw[k, , , ] + e.loadings[i, k, t] * nYall[i, , , ,t]
      d
      #plot the current component
      title = paste('PC(', k ,'), time=', t, sep = "", collapse=NULL)
      image(as.matrix(eF.raw[k, , , 55]),axes = FALSE, main=title) #, col = grey(seq(0, 1, length = 256)))
    }
    eF.raw[k, , , ] = eF.raw[k, , , ]/ttime
  }
  #normalize eigenfunctions
  eF_w[k, , , ] = eF.raw[k, , , ]/sqrt(threeDsum(eF.raw[k, , , ]^2)) * sqrt(M * N * L)  
}

# =========== PC COMPONENTS ===========

#### smooth eF

K = c(30, 35, 30)
dim = c(M, N, L)
p = 3
mm = 2

B1 = genB(M, K[1], p)
B2 = genB(N, K[2], p)
B3 = genB(L, K[3], p)
P1 = genP(mm, K[1], p)
P2 = genP(mm, K[2], p)
P3 = genP(mm, K[3], p)

res1 = gentD(B1, P1)
D1 = res1$D
temp1 = res1$temp

res2 = gentD(B2, P2)
D2 = res2$D
temp2 = res2$temp

res3 = gentD(B3, P3)
D3 = res3$D
temp3 = res3$temp


eF.smooth = array(0, dim = c(r, M, N, L))
eF.trim = array(0, dim = c(r, M, N, L))

for (k in 1:r) {
  Y2 = eF[k, , , ]^2
  sY2 = threeDsum(Y2)
  
  oldY2 = array(0, c(K[1] + p, K[2] + p, L))
  for (d3 in 1:L) oldY2[, , d3] = temp1 %*% eF[k, , , d3] %*% t(temp2)
  
  tildeY2 = array(0, c(K[1] + p, K[2] + p, K[3] + p))
  for (d1 in 1:(K[1] + p)) for (d2 in 1:(K[2] + p)) tildeY2[d1, d2, ] = temp3 %*% 
    oldY2[d1, d2, ]
  
  tildeY2 = tildeY2^2
  
  lambda0 = gridSearch(GCV, levels = list(exp(seq(-15, 15, length = 15)), 
                                          exp(seq(-15, 15, length = 15)),
                                          exp(seq(-15, 15, length = 15))), 
                       sY2 = sY2, dim = dim, 
                       D1 = D1, D2 = D2, D3 = D3, tildeY2 = tildeY2, printDetail = FALSE)$minlevels
  lambda = nlm(GCV, lambda0, sY2 = sY2, dim = dim, D1 = D1, D2 = D2, D3 = D3, tildeY2 = tildeY2)$estimate
  
  H1 = B1 %*% solve(t(B1) %*% B1 + lambda[1] * P1) %*% t(B1) #smoothing matrix 1 s_lambda
  H2 = B2 %*% solve(t(B2) %*% B2 + lambda[2] * P2) %*% t(B2) #smoothing matrix 2
  H3 = B3 %*% solve(t(B3) %*% B3 + lambda[3] * P3) %*% t(B3) #smoothing matrix 3
  
  oldfac = array(0, dim = dim)
  for (d3 in 1:L) oldfac[, , d3] = H1 %*% eF[k, , , d3] %*% H2
  
  hatfac = array(0, dim = dim)
  for (d1 in 1:M) for (d2 in 1:N) hatfac[d1, d2, ] = H3 %*% oldfac[d1, d2, ]
  
  eF.smooth[k, , , ] = hatfac
  
  #### trim the factors with the quantile bands for(k in 1:r){
  quantiles = quantile(eF.smooth[k, , , ], probs = c(0.001, 0.999))
  for (d1 in 1:M) {
    for (d2 in 1:N) {
      for (d3 in 1:L) {
        if (quantiles[1] < eF.smooth[k, d1, d2, d3] && eF.smooth[k, d1, d2, 
                                                                 d3] < quantiles[2]) {
          eF.trim[k, d1, d2, d3] = 0
        } else {
          eF.trim[k, d1, d2, d3] = 1
        }
      }
    }
  }
  
}


# =========== UPDATE SCORES ===========


F = matrix(0, r, r)
for (i in 1:r) {
  for (j in 1:r) {
    F[i, j] = threeDsum(eF.smooth[i, , , ] * eF.smooth[j, , , ])
  }
}

Y.tilde = array(0, c(r, I, ttime))
for (i in 1:r) {
  for (j in 1:I) {
    for (t in 1:ttime) {
      #do not separate because then you loose the temporal aspect
      Y.tilde[i, j, t] = threeDsum(nYall[j, , , ,t] * eF.smooth[i, , , ]) 
    }
  }
}

load.new = array(0, c(r,I,ttime))
for (t in 1:ttime){
  load.new[,,t] = solve(F) %*% Y.tilde[,,t]
}

par(mar=c(2,2,2,2))
par(mfrow = c(3,2))
for(j in 1:r){
  name <- paste('Updated loadings for component ', j ,  sep = "", collapse=NULL)
  matplot(c(1:16),t(load.new[j,,]), type='l')
  title(name)
}

#load new are the new loadings/scores (there is r of them)
#eF.smooth are the smoothed PC components in 3D - dimension r, 208, 208, 10 

# =========== UPDATE SCORES ===========


F = matrix(0, r, r)
for (i in 1:r) {
  for (j in 1:r) {
    F[i, j] = threeDsum(eF.raw[i, , , ] * eF.raw[j, , , ])
  }
}

Y.tilde = array(0, c(r, I, ttime))
for (i in 1:r) {
  for (j in 1:I) {
    for (t in 1:ttime) {
      #do not separate because then you loose the temporal aspect
      Y.tilde[i, j, t] = threeDsum(nYall[j, , , ,t] * eF.raw[i, , , ]) 
    }
  }
}

load.new = array(0, c(r,I,ttime))
for (t in 1:ttime){
  load.new[,,t] = solve(F) %*% Y.tilde[,,t]
}

par(mar=c(2,2,2,2))
par(mfrow = c(3,2))
for(j in 1:r){
  name <- paste('Updated loadings for component ', j ,  sep = "", collapse=NULL)
  matplot(c(1:ttime),t(load.new[j,,]), type='l')
  title(name)
}

#===========================================

#print the PC components
par(mar=c(0.5,0.5,0.5,0.5))
par(mfrow = c(2,5))
for(i in 1:5){
  image(as.matrix(eF[i, , ,55]),axes = FALSE, col = grey(seq(0, 1, length = 256)))
}
for(i in 1:5){
  image(as.matrix(eF[i, , ,55]),axes = FALSE)
}
