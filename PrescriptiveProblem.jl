using Gurobi
using JuMP
using DataFrames
using XLSX
using Statistics

function n(t,w)
    if t == 1
        n = 1
    elseif t == 2
        if w in 1:67
            n = 1
        elseif w in 68:134
            n = 2
        else
            n = 3
        end
    elseif t == 3
        if w in 1:23
            n = 1
        elseif w in 24:45
            n = 2
        elseif w in 46:67
            n = 3
        elseif w in 68:90
            n = 4
        elseif w in 91:112
            n = 5
        elseif w in 113:134
            n = 6
        elseif w in 135:156
            n = 7
        elseif w in 157:178
            n = 8
        else
            n = 9
        end
    else
        if w in 1:7
            n = 1
        elseif w in 8:15
            n = 2
        elseif w in 16:23
            n = 3
        elseif w in 24:30
            n = 4
        elseif w in 31:37
            n = 5
        elseif w in 38:45
            n = 6
        elseif w in 46:52
            n = 7
        elseif w in 53:59
            n = 8
        elseif w in 60:67
            n = 9
        elseif w in 68:74
            n = 10
        elseif w in 75:82
            n = 11
        elseif w in 83:90
            n = 12
        elseif w in 91:97
            n = 13
        elseif w in 98:104
            n = 14
        elseif w in 105:112
            n = 15
        elseif w in 113:119
            n = 16
        elseif w in 120:126
            n = 17
        elseif w in 127:134
            n = 18
        elseif w in 135:141
            n = 19
        elseif w in 142:148
            n = 20
        elseif w in 149:156
            n = 21
        elseif w in 157:163
            n = 22
        elseif w in 164:170
            n = 23
        elseif w in 171:178
            n = 24
        elseif w in 179:185
            n = 25
        elseif w in 186:192
            n = 26
        else
            n = 27
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
xmin = ones(T)*(-1)
xmax = ones(T)
ymin = 0
ymax = 1
q = 100         # Quantidade do contrato semenal
λ = [
    [mean(PLD[1,:,:])],

    [mean(PLD[2,:,1:67]), mean(PLD[2,:,68:134]), mean(PLD[2,:,135:200])],

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