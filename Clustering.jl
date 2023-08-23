using Pkg
using Random
using Clustering

PLD_excel = Matrix(DataFrame(XLSX.readtable("Spot_price_scenarios.xlsx", "S_insample (2)")))
PLD = zeros(T,H,Ω)

PLD[1,:,:] = PLD_excel[1:168,1:Ω]
PLD[2,:,:] = PLD_excel[169:336,1:Ω]
PLD[3,:,:] = PLD_excel[337:504,1:Ω]
PLD[4,:,:] = PLD_excel[505:672,1:Ω]

Y = mean(PLD[1,:,:],dims=1)

Random.seed!(13124152)

# cluster X into 2 clusters using K-means
R = kmeans(Y, 3)

@assert nclusters(R) == 3 # verify the number of clusters

a = assignments(R) # get the assignments of points to clusters
c = counts(R) # get the cluster sizes
M = R.centers # get the cluster centers


# make a random dataset with 1000 random 5-dimensional points
X = PLD[1,:,:]

Random.seed!(13124152)

# cluster X into 2 clusters using K-means
R = kmeans(X, 3)

@assert nclusters(R) == 3 # verify the number of clusters

a = assignments(R) # get the assignments of points to clusters
c = counts(R) # get the cluster sizes
M = R.centers # get the cluster centers

counter1=1
counter2=1

for i in 1:Ω
    if a[i] == 2
        println(i)
    end
end


X3 = zeros(4,168,c[1])
X2 = zeros(4,168,c[2])

for i in 1:Ω
    if a[i] == 2
        for t in 1:T
            X2[t,:,counter1] = PLD[t,:,i]
        end
        counter1+=1
    else
        for t in 1:T
            X3[t,:,counter2] = PLD[t,:,i]
        end
        counter2+=1
    end
end

R2 = kmeans(X2[2,:,:], 2)

@assert nclusters(R2) == 2 # verify the number of clusters

a = assignments(R2) # get the assignments of points to clusters
c = counts(R2) # get the cluster sizes
M = R2.centers # get the cluster centers

counter1=1
counter2=1
X4 = zeros(4,168,c[1])
X5 = zeros(4,168,c[2])

for i in 1:Ω
    if a[i] == 1
        for t in 1:T
            X4[t,:,counter1] = PLD[t,:,i]
        end
        counter1+=1
    else
        for t in 1:T
            X5[t,:,counter2] = PLD[t,:,i]
        end
        counter2+=1
    end
end
