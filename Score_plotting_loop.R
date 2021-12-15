# =========== DATA ===========

#image directroy, retrieve all images
image_dir = "D:/R/Spatial_FPCA/FMRI"

fmri_data_names <- list.files(path = image_dir, pattern = "*.mat" )
fmri_data_names <- fmri_data_names[-10]
fmri_data_names <- fmri_data_names[-8]

task_total = 81
task_idx=20
k=1


for (task_idx in 2:15){
  ttime = 50
  iter_Yall = array(0, c(I, M, N, L, ttime))
  iter_aYall = array(0, c(I, m, n, l, ttime))
  iter_aYva = array(0, c(I, ttime))
  k=1
  
  for (i in 1:length(subjects)){
    print(i)
    path <-  paste(image_dir, fmri_data_names[i], sep = "/", collapse=NULL)
    
    fmri_img <- rhdf5::h5read(path, "/temp")
    print(dim(fmri_img))
    img_idx = (task_idx-1)*50 + 1
    it_fmri_img <- fmri_img[,,,img_idx:(img_idx+49)]
    
    iter_Yall[k, , , ,  ] = it_fmri_img[,,,]
    
    k=k+1
  }
  
  
  # normalizing the images
  iter_nYall = array(0, c(I, M, N, L, ttime))
  yall_mean = array(0, c(M, N, L, ttime))
  
  for (d1 in 1:M) {
    for (d2 in 1:N) {
      for (d3 in 1:L) {
        for (d4 in 1:ttime) {
          
          yall_mean[d1, d2, d3, d4] = mean(iter_Yall[,d1,d2,d3,d4])
          
        }
      }
    }
  }
  
  for (d1 in 1:15){
    iter_nYall[d1, , , ,] = iter_Yall[d1, , , ,] - yall_mean[ , , , ]
  }
  
  
  # recalculate the cubic mean etc.
  for (i in 1:I){
    for (j in 1:ttime){
      result = cube_mean2(iter_nYall[i, , , , j], m, n, l)
      iter_aYall[i, , , , j] = result$Y.new
      iter_aYva[i,j] = result$va
    }
  }
  
  
  # V_X OVEER TIME 
  
  #inner product of the images - takes image by image
  #had to adjust the way the sum is calculated as it didn't seem to work otherwise
  #this is V_y, over time (3rd dimesnion)
  iter_LAY = array(0, c(I, I, ttime))
  for (t in 1:ttime){
    for (i in 1:I) {
      for (j in i:I) {
        iter_LAY[i, j, t] = threeDsum(iter_aYall[i, , , , t] * iter_aYall[j, , , , t])
        # print(paste(i,',',j)) 
      }
    }
  }
  
  #computing the other side of V_Y to make it symmetrical
  #the function 'aperm' is the transpose of 3D array
  for (t in 1:ttime){
    iter_LAY[,,t] = iter_LAY[,,t] + aperm(iter_LAY[,,t])
    diag(iter_LAY[,,t]) = diag(iter_LAY[,,t])/2
  }
  
  #compute the eigendecomposition of V_y - V_\sigma
  #aYva is the sum of all the variation in the 30 x 35 x 30 image
  #we need to do this separately by time points
  r = 6
  iter_load.old = array(0, c(ttime, r, 15))
  for (t in 1:ttime){
    iter_eLAY = eigen(iter_LAY[,,t]/m/n/l - diag(iter_aYva[,t]))
    iter_eA = t(iter_eLAY$vectors[, 1:r])  #get loadings approach 1
    iter_load.old[t,,] = iter_eA
  }
  
  #Plot the scores for everyone by component
  par(mar=c(2,2,2,2))
  par(mfrow = c(3,2))
  for(j in 1:r){
    name <- paste('T',task_idx,': Loadings for component ', j ,  sep = "", collapse=NULL)
    matplot(c(1:50),iter_load.old[,j,], type='l')
    title(name)
  }
}


#need to make a function to fix the issue that we see with the scores jumping from negative to positive
for (i in 1:15){
  print(min(iter_load.old[,6,i]))
}
# 8 is the person with red graph
print(iter_load.old[,6,8])

for (j in 1:r) {
  score_signs = rep(NA, I)
  for (i in 1:I) {
    
    ij_mean = mean(iter_load.old[,j,i])
    ij_signs = sign(ij_mean)
    score_signs[i] = ij_signs
  }
  
  
  for (tt in 1:ttime) {
    time_sign = rep(1, I)
    for (i in 1:I){
      if (sign(iter_load.old[tt,j,i]) == score_signs[i]) {
        time_sign[i] = 0 
      }
    }
    
    if (sum(time)<(I-7)){ iter_load.old[tt,j,] = (-1)*iter_load.old[tt,j,] }
    
  }
  #sign(iter_load[]) <- gives you whether this is a positive or negative value in a vector
}
#if the mean of a graph is positive ethen we assume this is the right sign for the score at this time and 
# therefore we need to switch the signs for the whole subject group
#this occurs for everyone at the ame time so this way you can se this fact 
# to make sure that you arent taking the wrong scores out f context (some scores can naturally venture between pos and neg)

# get all the scores,
# find the score that is at its highest or furthest away from zero
# track the score over time, when the sign switches, switch all the signs for that time point

for (j in 1:r) {
  
  #find which subject has the maximum valued score at the beginning
  max_idx = which.max(abs(load.old[1,j,]))
  max_sign = sign(load.old[1,j,max_idx])
  
  for (tt in 1:ttime) {
    if (sign(load.old[tt,j,max_idx]) != max_sign) {
      
      load.old[tt,j,] = (-1)*load.old[tt,j,]
      
    }
  }
}