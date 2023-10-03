include("utils.jl");
using ModelingToolkit
using OrdinaryDiffEq, Plots, DifferentialEquations
using BenchmarkTools
using ForwardDiff

@parameters gamma_1 psi C0 D0 r K BW IC50_1 Imax_1 IC50_2 gamma_2 Imax_2 xi VD1 Cl1 k23 ka1 k32 Cl2 ka2 Vpla Q Vtis dose_vector
@variables t _CSFTMZ(t) _PlaRG(t) _C(t) C(t) D(t) AbsTMZ(t) PlaTMZ(t) CSFTMZ(t) AbsRG(t) PlaRG(t) TisRG(t) cAUC(t) pRGauc(t) 

Diff = Differential(t)

function erelu(x)
    return ifelse(x >= 0, x, eps())
end
@register erelu(x)


eqns = [
    _CSFTMZ ~ erelu(CSFTMZ),
    _PlaRG ~ erelu(PlaRG),
    _C ~ erelu(C),

    Diff(C) ~ _C * r * log(K / _C) - (((Imax_1 * (((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + Imax_2 * (((xi * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + (Imax_1 + Imax_2) * ((((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * (((xi * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) - Imax_1 * Imax_2 * ((((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * (((xi * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2))) / ((((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + (((xi * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + ((((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * (((xi * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) + 1)) * (K * exp((log(C0/K)) * exp(-r * ((-log(log(C/K) / log(C0/K)) / r) + 72)))) / (72 * _C)) * _C,
    Diff(D) ~ (((Imax_1 * (((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + Imax_2 * (((xi * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + (Imax_1 + Imax_2) * ((((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * (((xi * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) - Imax_1 * Imax_2 * ((((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * (((xi * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2))) / ((((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + (((xi * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + ((((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * (((xi * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) + 1)) * (K * exp((log(C0/K)) * exp(-r * ((-log(log(C/K) / log(C0/K)) / r) + 72)))) / (72 * _C)) * _C,

    Diff(AbsTMZ) ~ -ka1 * AbsTMZ,
    Diff(PlaTMZ) ~ ka1 * AbsTMZ - (Cl1 / VD1) * PlaTMZ - k23 * PlaTMZ + k32 * _CSFTMZ,
    Diff(CSFTMZ) ~ k23 * PlaTMZ - k32 * _CSFTMZ,
    Diff(AbsRG) ~ -ka2 * AbsRG,
    Diff(PlaRG) ~ ka2 * AbsRG - (Cl2 / Vpla) * _PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * _PlaRG,
    Diff(TisRG) ~ -(Q / Vtis) * TisRG + (Q / Vpla) * _PlaRG,
    
    Diff(cAUC) ~ _C,
    Diff(pRGauc) ~ _PlaRG
]

hours = 24.0
end_time = 28.0*5.0
end_treat = 42.0
tspan = (0.0, end_time+7.0) .* hours

#dose amounts
avg_human_surface_area = 1.7 #m^2
tmz_treat_dose = 75.0*avg_human_surface_area
tmz_adjuv_dose = 150.0*avg_human_surface_area
dose_amount = 1800.0*avg_human_surface_area

tmz_treat_dosetimes = spaced_list(end_treat,1.0,0.0,0.0).*hours
tmz_adjuv_dosetimes = spaced_list(end_time,5.0,23.0,end_treat+28.0).*hours
rg_dosetimes = spaced_list(end_time-1.0,18.0,10.0,0.0).*hours

function tmz_dose!(integ, u, p, ctx)
    SciMLBase.set_proposed_dt!(integ, 0.001)
    integ.u[u.AbsTMZ] += ctx
end

tmz_treat = tmz_treat_dosetimes => (tmz_dose!, [AbsTMZ], [], tmz_treat_dose)
tmz_adjuv = tmz_adjuv_dosetimes => (tmz_dose!, [AbsTMZ], [], tmz_adjuv_dose)

function rg_dose!(integ, u, p, ctx)
    SciMLBase.set_proposed_dt!(integ, 0.01)
    dose_idx = findfirst(==(integ.t), rg_dosetimes)
    dose_to_add = integ.p[p[end][2]][dose_idx]
    integ.u[u.AbsRG] += dose_to_add
end

doses = ones(length(rg_dosetimes)).*dose_amount
rg_dose = rg_dosetimes => (rg_dose!, [AbsRG], [Q], nothing)

@named osys = ODESystem(eqns;
                        discrete_events = [rg_dose,tmz_treat, tmz_adjuv])

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
    # dose_vector => doses
]

u0 = [
    C => 17.7,
    D => 0.0,
    AbsTMZ => 0.0,
    PlaTMZ => 0.0,
    CSFTMZ => 0.0,
    AbsRG => 0.0,
    PlaRG => 0.0,
    TisRG => 0.0,
    cAUC => 0.0,
    pRGauc => 0.0
]

oprob = ODEProblem(structural_simplify(osys), u0, tspan, p, jac=true)

sol = solve(oprob)
plotter(sol)

@btime sol = solve(oprob)

#get gradients


function loss(θ, ode_params)
    prob_jac = remake(oprob, p=[ode_params;θ])
    int_sol = solve(prob_jac)

    # cell = (int_sol[9,end]-ode_params[end-1][1])/(ode_params[end-1][2]-ode_params[end-1][1])
    # drug = (int_sol[10,end]-ode_params[end-1][3])/(ode_params[end-1][4]-ode_params[end-1][3])
    return  int_sol[9,end] + int_sol[10,end]#cell + drug
end

ForwardDiff.gradient(θ -> loss(θ, p), doses)