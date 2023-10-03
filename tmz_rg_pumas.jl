using ModelingToolkit
using OrdinaryDiffEq, Plots, DifferentialEquations
using PumasQSP

@parameters gamma_1 psi C0 D0 r K BW IC50_1 Imax_1 IC50_2 gamma_2 Imax_2 xi VD1 Cl1 k23 ka1 k32 Cl2 ka2 Vpla Q Vtis
@variables t C(t) D(t) AbsTMZ(t) PlaTMZ(t) CSFTMZ(t) AbsRG(t) PlaRG(t) TisRG(t)# cAUC(t) pRGauc(t) #cCSFTMZ(t) cPlaRG(t) pi1(t) exp1(t) exp2(t) exp12(t) Imax_sum(t) E_num(t) E_den(t) E(t) t1(t) t2(t) log_C0K(t) fun(t) delta(t)

Diff = Differential(t)

function erelu(x)
    return ifelse(x >= 0 ? x : 0)
end

eqns = [
    Diff(C) ~ C * r * log(K / C) - (((Imax_1 * (((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + Imax_2 * (((xi * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + (Imax_1 + Imax_2) * ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) - Imax_1 * Imax_2 * ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2))) / ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + (((xi * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) + 1)) * (K * exp((log(C0/K)) * exp(-r * ((-log(log(C/K) / log(C0/K)) / r) + 72)))) / (72 * C)) * C,
    Diff(D) ~ (((Imax_1 * (((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + Imax_2 * (((xi * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + (Imax_1 + Imax_2) * ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) - Imax_1 * Imax_2 * ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2))) / ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + (((xi * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) + 1)) * (K * exp((log(C0/K)) * exp(-r * ((-log(log(C/K) / log(C0/K)) / r) + 72)))) / (72 * C)) * C,
    # Diff(C) ~ C * r * log(K / C) - ((Imax_1*((CSFTMZ / 140)/pi1)^gamma_1 + Imax_2*((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2 + (Imax_1 + Imax_2) * (((CSFTMZ / 140)/pi1)^gamma_1 * ((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2) - Imax_1 * Imax_2 * (((CSFTMZ / 140)/pi1)^gamma_1 * ((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2)) / (((CSFTMZ / 140)/pi1)^gamma_1 + ((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2 + (((CSFTMZ / 140)/pi1)^gamma_1 * ((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2) + 1) * K * exp(log(C0/K) * exp(-r * (-log(log(C/K) / log(C0/K)) / r + t1 + 72))) / (72 * C) * C),
    # Diff(D) ~ ((Imax_1*((CSFTMZ / 140)/pi1)^gamma_1 + Imax_2*((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2 + (Imax_1 + Imax_2) * (((CSFTMZ / 140)/pi1)^gamma_1 * ((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2) - Imax_1 * Imax_2 * (((CSFTMZ / 140)/pi1)^gamma_1 * ((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2)) / (((CSFTMZ / 140)/pi1)^gamma_1 + ((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2 + (((CSFTMZ / 140)/pi1)^gamma_1 * ((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2) + 1) * K * exp(log(C0/K) * exp(-r * (-log(log(C/K) / log(C0/K)) / r + t1 + 72))) / (72 * C) * C),
    # cCSFTMZ ~ CSFTMZ / 140,
    # cPlaRG ~ PlaRG / (1000 * Vpla),

    # # #domain error
    # # cCSFTMZ = erelu(cCSFTMZ)
    # # cPlaRG = erelu(cPlaRG)
    # # C = erelu(C)

    # pi1 ~ psi * IC50_1,
    # exp1 ~ (cCSFTMZ/pi1)^gamma_1,
    # exp2 ~ ((xi * cPlaRG)/pi1)^gamma_2,
    # exp12 ~ exp1 * exp2,

    # Imax_sum ~ Imax_1 + Imax_2,
    # E_num ~ Imax_1 * exp1 + Imax_2 * exp2 + Imax_sum * exp12 - Imax_1 * Imax_2 * exp12,
    # E_den ~ exp1 + exp2 + exp12 + 1,
    # E ~ E_num / E_den,

    # t1 ~ -log(log(C/K) / log(C0/K)) / r,
    # t2 ~ t1 + 72,
    # log_C0K ~ log(C0/K),
    # fun ~ K * exp(log_C0K * exp(-r * t2)),
    # delta ~ E * fun / (72 * C),

    # D(C) ~ C * r * log(K / C) - delta * C,
    # D(D) ~ delta * C,
    
    Diff(AbsTMZ) ~ -ka1 * AbsTMZ,
    Diff(PlaTMZ) ~ ka1 * AbsTMZ - (Cl1 / VD1) * PlaTMZ - k23 * PlaTMZ + k32 * CSFTMZ,
    Diff(CSFTMZ) ~ k23 * PlaTMZ - k32 * CSFTMZ,
    Diff(AbsRG) ~ -ka2 * AbsRG,
    Diff(PlaRG) ~ ka2 * AbsRG - (Cl2 / Vpla) * PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * PlaRG,
    Diff(TisRG) ~ -(Q / Vtis) * TisRG + (Q / Vpla) * PlaRG,
    # Diff(cAUC) ~ C,
    # Diff(pRGauc) ~ PlaRG
]

p = [
    gamma_1 => 0.42676, #0.5689
    psi => 1,          #0.8
    C0 => 17.7,
    D0 => 0.0,
    r => 0.007545/24,
    K => 158.04,
    BW => 70.0,
    IC50_1 => 15.5936*0.000001*194.151, #5.807*0.000001*194.151
    Imax_1 => 1.1026,   #0.905
    IC50_2 => 0.00924,
    gamma_2 => 2.712,
    Imax_2 => 1.0,
    xi => IC50_1/IC50_2,
    VD1 => 30.3,
    Cl1 => 10.5705229706946, #V
    k23 => 0.000823753578728557, #V
    ka1 => 9.75543575220511,     #V
    k32 => 0.76,
    Cl2 => 32.6682,  #V
    ka2 => 0.0385233, #V
    Vpla => 0.934662,
    Q => 0.0302696,
    Vtis => 0.00299745
]
# C, D, AbsTMZ, PlaTMZ, CSFTMZ, AbsRG, PlaRG, TisRG = 17.7,0.0,0.0,0.0,0.0,0.0,0.0,0.0

@named de = ODESystem(eqns)
#, t, [C, D, AbsTMZ, PlaTMZ, CSFTMZ, AbsRG, PlaRG, TisRG], ode_params, tspan=(0.0,1000.0))

sys = structural_simplify(de)
# sys=de
u0 = [
    C => 17.7,
    D => 0.0,
    AbsTMZ => 0.0,
    PlaTMZ => 0.0,
    CSFTMZ => 0.0,
    AbsRG => 0.0,
    PlaRG => 0.0,
    TisRG => 0.0
]
tspan = (0.0, 50.0) .* 24.0
prob_jac = ODEProblem(sys, u0, tspan, p, jac=true)
sol = solve(prob_jac)

plot(sol)

#add dosing

@unpack AbsTMZ = sys
tmz_dose = Bolus(AbsTMZ, )

trial = Trial(sys; tspan=tspan)
prob = InverseProblem(trial, de, [])
sol = simulate(trial, prob)
plot(sol)
