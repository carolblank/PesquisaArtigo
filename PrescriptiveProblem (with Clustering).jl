using Gurobi
using JuMP
using DataFrames
using XLSX
using Statistics

function n(t,w)
    if t == 1
        n = 1
    elseif t == 2
        if w in v2
            n = 2
        else w in v3
            n = 3
        end
    elseif t == 3
        if w in v4
            n = 4
        elseif w in v5
            n = 5
        elseif w in v6
            n = 6
        elseif w in v7
            n = 7
        end
    else
        if w in v8
            n = 8
        elseif w in v9
            n = 9
        elseif w in v10
            n = 10
        elseif w in v11
            n = 11
        elseif w in v12
            n = 12
        elseif w in v13
            n = 13
        elseif w in v14
            n = 14
        elseif w in v15
            n = 15
        end
    end
return n
end


# Parametros do caso
Ω = 200
T = 4
H = 168
GF = 100

geracao = Matrix(DataFrame(XLSX.readtable("Wind_plant_scenarios_and_parms.xlsx", "G_insample (2)")))
g = zeros(T,H,Ω)

g[1,:,:] = geracao[1:168,1:Ω]./76.4
g[2,:,:] = geracao[169:336,1:Ω]./76.4
g[3,:,:] = geracao[337:504,1:Ω]./76.4
g[4,:,:] = geracao[505:672,1:Ω]./76.4

PLD_excel = Matrix(DataFrame(XLSX.readtable("Spot_price_scenarios.xlsx", "S_insample (2)")))
PLD = zeros(T,H,Ω)

PLD[1,:,:] = PLD_excel[1:168,1:Ω]
PLD[2,:,:] = PLD_excel[169:336,1:Ω]
PLD[3,:,:] = PLD_excel[337:504,1:Ω]
PLD[4,:,:] = PLD_excel[505:672,1:Ω]


# Parametros de contrato/quantidade
xmin = ones(15)*(-1)
xmax = ones(15)
ymin = 0
ymax = 1
q = 100         # Quantidade do contrato semenal
media_PLD = [
    mean(PLD[1,:,:]),

    mean(X2[2,:,:]), mean(X3[2,:,:]),

    mean(X4[3,:,:]), mean(X5[3,:,:]), 
    mean(X6[3,:,:]), mean(X7[3,:,:]),

    mean(X8[4,:,:]), mean(X9[4,:,:]), 
    mean(X10[4,:,:]), mean(X11[4,:,:]), 
    mean(X12[4,:,:]), mean(X13[4,:,:]),
    mean(X14[4,:,:]), mean(X15[4,:,:])
]
λ = 0.93.*media_PLD .+15.31     # Preço contrato semanal
Q = 100                         # Quantidade do contrato mensal
P = 0.93*mean(PLD) +20          # Preço contrato mensal

# Parametros de risco
p = 1 / Ω .* ones(Ω);
α = 0.95
ξ = 0.5

# Otimização

model = Model(Gurobi.Optimizer)

@variable(model, R[1:T,1:Ω])
@variable(model, z)
@variable(model, δ[1:Ω] >= 0)
@variable(model, x[1:15])
@variable(model, y)

@objective(
    model,
    Max,
    ξ * (z - sum(δ[w] * p[w] / (1 - α) for w in 1:Ω)) +
    (1 - ξ) * sum(p[w] *sum(R[t,w] for t in 1:T) for w in 1:Ω)
);

@constraint(
    model,
    [t ∈ 1:T, w ∈ 1:Ω],
    R[t,w] ==
    sum((g[t,h,w]*GF*PLD[t,h,w]) + (P - PLD[t,h,w])*Q*y + (λ[n(t,w)] - PLD[t,h,w])*q*x[n(t,w)] for h in 1:H)
);

@constraint(model, [t ∈ 1:T, w ∈ 1:Ω], Q*y + q*x[n(t,w)] <= GF)

@constraint(model, [n ∈ 1:15], xmin[n] <= x[n])

@constraint(model, [n ∈ 1:15], x[n] <= xmax[n])

@constraint(model, ymin <= y)

@constraint(model, y <= ymax)

@constraint(model, z - 1/(1-α)*sum(p[w]*δ[w] for w in 1:Ω) >= 0)

@constraint(model, [w ∈ 1:Ω], δ[w] >= z - sum(R[t,w] for t in 1:T))

status = optimize!(model)
x_otimo = JuMP.value.(x)
y_otimo = JuMP.value(y)
R_otimo = JuMP.value.(R)
resp = JuMP.objective_value(model)