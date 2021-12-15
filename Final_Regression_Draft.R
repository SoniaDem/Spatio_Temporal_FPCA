
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

subj = 1

# ======== CONVERT IMAGE FOR 1 SUBJECT TO VEC =========

vecYall = array(0, c(M*N*L*ttime,1))

idx = 1
for (mm in 1:M) {
  for (nn in 1:N) {
    for (ll in 1:L) {
      for (tt in 1:ttime) {
        
        vecYall[idx] = nYall[subj, mm, nn, ll, tt]
        idx = idx + 1
        
      }
    }
  }
}


# ======== PREDETERMINED BAIS FUNCTION + KRONECKER PRODUCT ========= 
p = 3 

basis = bs(c((1:ttime)/ttime - 1/2/ttime), df = ttime + p, degree = p)

phi_basis = kronecker(eF.smooth[1, , , ], basis)

# ======== 1 BASIS VECTOR CASE ========= 

XX = array(0, c(M, N, L, ttime))

basis_1 = basis[,1]

for (idx in 1:length(basis_1)) {
  XX[ , , , idx] = basis_1[idx]*eF.smooth[1, , , ]
}


b_1_estimate = threeDsum(XX*XX)^(-1)*threeDsum(nYall[1, , , ,] * XX)*basis_1[2]



# ======== ALL BASIS VECTORS CASE - ONE SUBJECT & ONE COMPONENT ========= 

XX = array(0, c(M, N, L, ttime))

#rows go over time, columns go over the basis vectors
basis_all = basis

b_coefs = rep(NA, ncol(basis_all))

for (idxx in 1:ncol(basis_all)){
  #we want to calculate one coefficient at a time so we take each basis vector
  #separately and compute regression that way
  for (idx in 1:nrow(basis_all)) {
    XX[ , , , idx] = basis_all[idx, idxx]*eF.smooth[1, , , ]
  }
  #balculate the b coefficient
  b_coefs[idxx] = threeDsum(XX*XX)^(-1)*threeDsum(Yall[1, , , ,] * XX)
}

#recreate the score function over time
score_func_ij = b_coefs %*% t(basis_all)
matplot(t(score_func_ij), type='l')





# ========= ALL BASIS VECTORS - ONE SUBJECT, ALL COMPONENTS ========= 

XX = array(0, c(M, N, L, ttime))

#rows go over time, columns go over the basis vectors
basis_all = basis

b_coefs = array(0, c(r, ncol(basis_all)))

for (iscore in 1:r){
  
  for (idxx in 1:ncol(basis_all)){
    #we want to calculate one coefficient at a time so we take each basis vector
    #separately and compute regression that way
    for (idx in 1:nrow(basis_all)) {
      XX[ , , , idx] = basis_all[idx, idxx]*eF.smooth[iscore, , , ]
    }
    #balculate the b coefficient
    FF = threeDsum(XX[ , , ,]*XX[ , , , ])^(-1)
    b_coefs[iscore, idxx] = FF*threeDsum(nYall[2, , , ,] * XX[ , , , ])
    
  }
}

#recreate the score function over time
score_func_ij = b_coefs %*% t(basis_all)
matplot(t(score_func_ij), type='l')


#this is a note


##
## ========= ALL BASIS VECTORS - ALL SUBJECTS, ALL COMPONENTS ========= 
##
F = matrix(0, 1, 1)
for (i in 1:1) {
  for (j in 1:1) {
    F[i, j] = threeDsum(eF.smooth[i, , , ] * eF.smooth[j, , , ])
  }
}
Y.tilde = matrix(0, r, I)
for (i in 1:r) {
  for (j in 1:I) {
    Y.tilde[i, j] = threeDsum(nYall[j, , , ] * eF.smooth[i, , , ])
  }
}
load.new = solve(F) %*% Y.tilde
