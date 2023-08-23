using Gurobi
using JuMP
using DataFrames
using XLSX
using Statistics




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
xmin = ones(T)*(-1)
xmax = ones(T)
ymin = 0
ymax = 1
q = 100         # Quantidade do contrato semenal
λ = [
    [mean(PLD[1,:,:])],

    [mean(X2[2,:,:]), mean(X3[2,:,:])],

    [mean(PLD[3,:,1:23]), mean(PLD[3,:,24:45]), mean(PLD[3,:,46:67]),
     mean(PLD[3,:,68:90]), mean(PLD[3,:,91:112]), mean(PLD[3,:,113:134]),
     mean(PLD[3,:,135:156]), mean(PLD[3,:,157:178]), mean(PLD[3,:,179:200])],

    [mean(PLD[4,:,1:7]), mean(PLD[4,:,8:15]), mean(PLD[4,:,16:23]),
     mean(PLD[4,:,24:30]), mean(PLD[4,:,31:37]), mean(PLD[4,:,38:45]),
     mean(PLD[4,:,46:52]), mean(PLD[4,:,53:59]), mean(PLD[4,:,60:67]),
     mean(PLD[4,:,68:74]), mean(PLD[4,:,75:82]), mean(PLD[4,:,83:90]),
     mean(PLD[4,:,91:97]), mean(PLD[4,:,98:104]), mean(PLD[4,:,105:112]),
     mean(PLD[4,:,113:119]), mean(PLD[4,:,120:126]), mean(PLD[4,:,127:134]),
     mean(PLD[4,:,135:141]), mean(PLD[4,:,142:148]), mean(PLD[4,:,149:156]),
     mean(PLD[4,:,157:163]), mean(PLD[4,:,164:170]), mean(PLD[4,:,171:178]),
     mean(PLD[4,:,179:185]), mean(PLD[4,:,186:192]), mean(PLD[4,:,193:200])]
]               # Preço contrato semanal
Q = 100         # Quantidade do contral mensal
P = 100         # Preço contrato mensal

# Parametros de risco
p = 1 / Ω .* ones(Ω);
α = 0.95

# Otimização

model = Model(Gurobi.Optimizer)

@variable(model, R[1:T,1:Ω])
@variable(model, z)
@variable(model, δ[1:Ω] >= 0)
@variable(model, x[1:T])
@variable(model, y)

@objective(
    model,
    Max,
    sum(p[w] *sum(R[t,w] for t in 1:T) for w in 1:Ω)
);

@constraint(
    model,
    [t ∈ 1:T, w ∈ 1:Ω],
    R[t,w] ==
    sum((g[t,h,w]*GF*PLD[t,h,w]) + (P - PLD[t,h,w])*Q*y + (λ[t][n(t,w)] - PLD[t,h,w])*q*x[t] for h in 1:H)
);

@constraint(model, [t ∈ 1:T], xmin[t] <= x[t])

@constraint(model, [t ∈ 1:T], x[t] <= xmax[t])

@constraint(model, ymin <= y)

@constraint(model, y <= ymax)

@constraint(model, z - 1/(1-α)*sum(p[w]*δ[w] for w in 1:Ω) >= 0)

@constraint(model, [w ∈ 1:Ω], δ[w] >= z - sum(R[t,w] for t in 1:T))

status = optimize!(model)
x_otimo = JuMP.value.(x)
y_otimo = JuMP.value(y)
R_otimo = JuMP.value.(R)
resp = JuMP.objective_value(model)