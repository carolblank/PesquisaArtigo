using Pkg
using Random
using Clustering

PLD_excel = Matrix(DataFrame(XLSX.readtable("Spot_price_scenarios.xlsx", "S_insample (2)")))
PLD = zeros(T,H,Ω)

PLD[1,:,:] = PLD_excel[1:168,1:Ω]
PLD[2,:,:] = PLD_excel[169:336,1:Ω]
PLD[3,:,:] = PLD_excel[337:504,1:Ω]
PLD[4,:,:] = PLD_excel[505:672,1:Ω]

Random.seed!(13124152)

# ==================== 1 ==================
X = PLD[1,:,:]

Random.seed!(13124152)

# cluster X into 2 clusters using K-means
R = kmeans(X, 2)

@assert nclusters(R) == 2 # verify the number of clusters

a1 = assignments(R) # get the assignments of points to clusters
c = counts(R) # get the cluster sizes
M = R.centers # get the cluster centers

counter1=1
counter2=1

X3 = zeros(4,168,c[1])
X2 = zeros(4,168,c[2])
v2 = zeros(74)
v3 = zeros(126)

for i in 1:Ω
    if a1[i] == 2
        for t in 1:T
            X2[t,:,counter1] = PLD[t,:,i]
        end
        v2[counter1] = i
        counter1+=1
    else
        for t in 1:T
            X3[t,:,counter2] = PLD[t,:,i]
        end
        v3[counter2] = i
        counter2+=1
    end
end

# ==================== 2 ==================
R2 = kmeans(X2[2,:,:], 2)

@assert nclusters(R2) == 2 # verify the number of clusters

a2 = assignments(R2) # get the assignments of points to clusters
c = counts(R2) # get the cluster sizes
M = R2.centers # get the cluster centers

counter1=1
counter2=1
X4 = zeros(4,168,9)
X5 = zeros(4,168,65)
v4 = zeros(9)
v5 = zeros(65)

for i in 1:74
    if a2[i] == 2
        for t in 1:T
            X4[t,:,counter1] = X2[t,:,i]
        end
        v4[counter1] = v2[i]
        counter1+=1
    else
        for t in 1:T
            X5[t,:,counter2] = X2[t,:,i]
        end
        v5[counter2] = v2[i]
        counter2+=1
    end
end

# ==================== 3 ==================
R3 = kmeans(X3[2,:,:], 2)

@assert nclusters(R3) == 2 # verify the number of clusters

a3 = assignments(R3) # get the assignments of points to clusters
c = counts(R3) # get the cluster sizes
M = R3.centers # get the cluster centers

counter1=1
counter2=1
X6 = zeros(4,168,55)
X7 = zeros(4,168,71)
v6 = zeros(55)
v7 = zeros(71)

for i in 1:126
    if a3[i] == 2
        for t in 1:T
            X6[t,:,counter1] = X3[t,:,i]
        end
        v6[counter1] = v3[i]
        counter1+=1
    else
        for t in 1:T
            X7[t,:,counter2] = X3[t,:,i]
        end
        v7[counter2] = v3[i]
        counter2+=1
    end
end

# ==================== 4 ==================
R4 = kmeans(X4[3,:,:], 2)

@assert nclusters(R4) == 2 # verify the number of clusters

a4 = assignments(R4) # get the assignments of points to clusters
c = counts(R4) # get the cluster sizes
M = R4.centers # get the cluster centers

counter1=1
counter2=1
X8 = zeros(4,168,7)
X9 = zeros(4,168,2)
v8 = zeros(7)
v9 = zeros(2)

for i in 1:9
    if a4[i] == 2
        for t in 1:T
            X8[t,:,counter1] = X4[t,:,i]
        end
        v8[counter1] = v4[i]
        counter1+=1
    else
        for t in 1:T
            X9[t,:,counter2] = X4[t,:,i]
        end
        v9[counter2] = v4[i]
        counter2+=1
    end
end

# ==================== 5 ==================
R5 = kmeans(X5[3,:,:], 2)

@assert nclusters(R5) == 2 # verify the number of clusters

a5 = assignments(R5) # get the assignments of points to clusters
c = counts(R5) # get the cluster sizes
M = R5.centers # get the cluster centers

counter1=1
counter2=1
X10 = zeros(4,168,18)
X11 = zeros(4,168,47)
v10 = zeros(18)
v11 = zeros(47)

for i in 1:65
    if a5[i] == 1
        for t in 1:T
            X10[t,:,counter1] = X5[t,:,i]
        end
        v10[counter1] = v5[i]
        counter1+=1
    else
        for t in 1:T
            X11[t,:,counter2] = X5[t,:,i]
        end
        v11[counter2] = v5[i]
        counter2+=1
    end
end

# ==================== 6 ==================
R6 = kmeans(X6[3,:,:], 2)

@assert nclusters(R6) == 2 # verify the number of clusters

a6 = assignments(R6) # get the assignments of points to clusters
c = counts(R6) # get the cluster sizes
M = R6.centers # get the cluster centers

counter1=1
counter2=1
X12 = zeros(4,168,33)
X13 = zeros(4,168,22)
v12 = zeros(33)
v13 = zeros(22)

for i in 1:55
    if a6[i] == 1
        for t in 1:T
            X12[t,:,counter1] = X6[t,:,i]
        end
        v12[counter1] = v6[i]
        counter1+=1
    else
        for t in 1:T
            X13[t,:,counter2] = X6[t,:,i]
        end
        v13[counter2] = v6[i]
        counter2+=1
    end
end

# ==================== 7 ==================
R7 = kmeans(X7[3,:,:], 2)

@assert nclusters(R7) == 2 # verify the number of clusters

a7 = assignments(R7) # get the assignments of points to clusters
c = counts(R7) # get the cluster sizes
M = R7.centers # get the cluster centers

counter1=1
counter2=1
X14 = zeros(4,168,19)
X15 = zeros(4,168,52)
v14 = zeros(19)
v15 = zeros(52)

for i in 1:71
    if a7[i] == 2
        for t in 1:T
            X14[t,:,counter1] = X7[t,:,i]
        end
        v14[counter1] = v7[i]
        counter1+=1
    else
        for t in 1:T
            X15[t,:,counter2] = X7[t,:,i]
        end
        v15[counter2] = v7[i]
        counter2+=1
    end
end