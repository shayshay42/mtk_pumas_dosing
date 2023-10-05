using ModelingToolkit
using OrdinaryDiffEq, Plots, DifferentialEquations
using BenchmarkTools
using ForwardDiff
using Serialization

# utils includes a function for creating dosing schedules (spaced_list) and a function for plotting the simulation (plotter)
include("assets/utils.jl");
# rg_vp.jls contains the population parameters for the RG model
population = deserialize("assets/rg_vp.jls")
# rg_min_tumor_auc_max_dose.jls contains the minimum tumor auc when each patient is treated with maximum dose
tumor_min_scaling = deserialize("assets/rg_min_tumor_auc_max_dose.jls")
# rg_max_tumor_auc_min_dose.jls contains the maximum tumor auc when each patient is treated with minimum dose
tumor_max_scaling = deserialize("assets/rg_max_tumor_auc_min_dose.jls")
# rg_min_dose_auc.jls contains the minimum drug auc when each patient is treated with minimum dose
drug_min_scaling = deserialize("assets/rg_min_dose_auc.jls")
# rg_max_dose_auc.jls contains the maximum drug auc when each patient is treated with maximum dose
drug_max_scaling = deserialize("assets/rg_max_dose_auc.jls")
#make a voctor of tuples of the form (x,y)
scaler = [collect(t) for t in zip(tumor_min_scaling, tumor_max_scaling, drug_min_scaling, drug_max_scaling)]

par = @parameters gamma_1 psi C0 D0 r K BW IC50_1 Imax_1 IC50_2 gamma_2 Imax_2 xi VD1 Cl1 k23 ka1 k32 Cl2 ka2 Vpla Q Vtis scaling dose_vector
sta = @variables t _CSFTMZ(t) _PlaRG(t) _C(t) C(t) D(t) AbsTMZ(t) PlaTMZ(t) CSFTMZ(t) AbsRG(t) PlaRG(t) TisRG(t) cAUC(t) pRGauc(t) 

Diff = Differential(t)

# this function allows us to circumvent Domain errors that occur when dosing events are encountered
function erelu(x)
    return ifelse(x >= 0, x, eps())
end
@register erelu(x)

#stiff awareness stuff and complexity
eqns = [
    _CSFTMZ ~ erelu(CSFTMZ),
    _PlaRG ~ erelu(PlaRG),
    _C ~ erelu(C),  #check is because fractional power in denominator so see complex part

    Diff(C) ~ _C * r * log(K / _C) - (((Imax_1 * (((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + Imax_2 * ((((IC50_1/IC50_2) * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + (Imax_1 + Imax_2) * ((((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) - Imax_1 * Imax_2 * ((((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2))) / ((((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + ((((IC50_1/IC50_2) * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + ((((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) + 1)) * (K * exp((log(C0/K)) * exp(-r * ((-log(log(C/K) / log(C0/K)) / r) + 72)))) / (72 * _C)) * _C,
    Diff(D) ~ (((Imax_1 * (((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + Imax_2 * ((((IC50_1/IC50_2) * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + (Imax_1 + Imax_2) * ((((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) - Imax_1 * Imax_2 * ((((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2))) / ((((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + ((((IC50_1/IC50_2) * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + ((((_CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (_PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) + 1)) * (K * exp((log(C0/K)) * exp(-r * ((-log(log(C/K) / log(C0/K)) / r) + 72)))) / (72 * _C)) * _C,

    Diff(AbsTMZ) ~ -ka1 * AbsTMZ,
    Diff(PlaTMZ) ~ ka1 * AbsTMZ - (Cl1 / VD1) * PlaTMZ - k23 * PlaTMZ + k32 * _CSFTMZ,
    Diff(CSFTMZ) ~ k23 * PlaTMZ - k32 * _CSFTMZ,
    Diff(AbsRG) ~ -ka2 * AbsRG,
    Diff(PlaRG) ~ ka2 * AbsRG - (Cl2 / Vpla) * _PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * _PlaRG,
    Diff(TisRG) ~ -(Q / Vtis) * TisRG + (Q / Vpla) * _PlaRG,
    
    Diff(cAUC) ~ _C, #AUC loss with fewer bits of precision
    Diff(pRGauc) ~ _PlaRG
]

# time frames of the trial
hours = 24.0
end_time = 28.0*5.0
end_treat = 42.0
tspan = (0.0, end_time+7.0) .* hours

tmz_treat_dosetimes = spaced_list(end_treat,1.0,0.0,0.0).*hours
tmz_adjuv_dosetimes = spaced_list(end_time,5.0,23.0,end_treat+28.0).*hours
rg_dosetimes = spaced_list(end_time-1.0,18.0,10.0,0.0).*hours

#dose amounts
avg_human_surface_area = 1.7 #m^2
tmz_treat_dose = 75.0*avg_human_surface_area
tmz_adjuv_dose = 150.0*avg_human_surface_area
dose_amount = 1800.0*avg_human_surface_area

#affect function for the discontinuous events
function tmz_dose!(integ, u, p, ctx)
    SciMLBase.set_proposed_dt!(integ, 0.01)
    integ.u[u.AbsTMZ] += ctx
end

tmz_treat = tmz_treat_dosetimes => (tmz_dose!, [AbsTMZ], [], tmz_treat_dose)
tmz_adjuv = tmz_adjuv_dosetimes => (tmz_dose!, [AbsTMZ], [], tmz_adjuv_dose)

function rg_dose!(integ, u, p, ctx)
    SciMLBase.set_proposed_dt!(integ, 0.01)
    integ.u[u.AbsRG] += integ.p[end][end][ctx[integ.t]]
end

doses = ones(length(rg_dosetimes)).*dose_amount
time_dose_map = Dict(zip(rg_dosetimes, Int64.(1:length(rg_dosetimes))))
rg_treat = rg_dosetimes => (rg_dose!, [AbsRG], [dose_vector], time_dose_map)

@named osys = ODESystem(eqns, t, sta[2:end], par, tspan=tspan, #excludes the independent variable t from states
                        discrete_events = [tmz_treat, tmz_adjuv, rg_treat])
# parameters(osys)
#trial set of parameters
p = [
    gamma_1 => 0.42676, #0.5689
    psi => 1.0,          #0.8
    C0 => 17.7,
    D0 => 0.0,   #UNUSED
    r => 0.007545/24,
    K => 158.04,
    BW => 70.0,   #UNUSED
    IC50_1 => 15.5936*0.000001*194.151, #5.807*0.000001*194.151
    Imax_1 => 1.1026,   #0.905
    IC50_2 => 0.00924,
    gamma_2 => 2.712,
    Imax_2 => 1.0,
    xi => IC50_1/IC50_2,  #UNUSED
    VD1 => 30.3,
    Cl1 => 10.5705229706946, #V
    k23 => 0.000823753578728557, #V
    ka1 => 9.75543575220511,     #V
    k32 => 0.76,
    Cl2 => 32.6682,  #V
    ka2 => 0.0385233, #V
    Vpla => 0.934662,
    Q => 0.0302696,
    Vtis => 0.00299745,
    scaling => scaler[1],
    dose_vector => doses
]
# p = (ode_params, dose_vector=?doses)
#initial conditions 
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

#simplification removes intermediate/auxilary equations from the system
simp=structural_simplify(osys)

ode_params = population[:,1]
scale = scaler[1]
doseit = ones(90).*1500.0
part = [ode_params;[scale];[doseit]]
prob=ODEProblem(simp, u0, tspan, part) #jac=true errors
# sol = solve(prob)
# plotter(sol)

# loss function balances the tumor auc and drug auc
function loss(θ,params)
    int_prob = remake(prob, p=(params[1:end-1],[params[end],θ]))
    int_sol = solve(int_prob)
    cell = (int_sol[9,end]- params[end][1])/( params[end][2]- params[end][1])
    drug = (int_sol[10,end]- params[end][3])/( params[end][4]- params[end][3])
    return cell + drug
end

ode_params = [population[:,2];[scaler[2]]]
# fill(1495.0,90)
loss(doses, ode_params)
@btime loss(doses, ode_params)

ForwardDiff.gradient(θ -> loss(θ, ode_params), doses)
@btime ForwardDiff.gradient(θ -> loss(θ, ode_params), doses)


# indexof(sym, syms) = findfirst(isequal(sym), syms)
# indexof(dose_vector, parameters(simp))
# indexof(scaling, parameters(simp))
# pnew = ModelingToolkit.varmap_to_vars(p, parameters(simp))

# oprob = ODEProblem(simp, u0, tspan, p) #errors with jac=true

# sol = solve(oprob)
# plotter(sol)

# @btime sol = solve(oprob)

# init = ModelingToolkit.varmap_to_vars(u0, sta[5:end])
# opar = ModelingToolkit.varmap_to_vars(p, par)
# prob = ODEProblem(simp, init, tspan, opar)
# sol = solve(prob)

# prob3 = remake(prob2, p=([population[:,1];[scaler[1]];[doseit].*0.999]))
# sol3 = solve(prob3)

# # Extracting the values and arrays from prob3.p
# values_from_prob3 = prob3.p[1:23]
# array1_from_prob3 = prob3.p[24]
# array2_from_prob3 = prob3.p[25]

# # Updating the nested array values
# # Assuming you want to use values of 1500.0 like in prob2.p
# updated_array2_from_prob3 = fill(1500.0, length(array2_from_prob3))

# # Constructing the new prob3.p with tuple structure
# prob3_par = (values_from_prob3, [array1_from_prob3, updated_array2_from_prob3])

# prob4 = remake(prob3, p=prob3_par)
# sol4 = solve(prob4)
# plotter(sol4)

# #get gradients
# indexof(sym, syms) = findfirst(isequal(sym), syms)
# indexof(dose_vector, parameters(simp))
# pnew = ModelingToolkit.varmap_to_vars(p, par)

# function tt()
#     prob = remake(oprob, p=p)
#     int_sol = solve(prob)
#     return int_sol[9,end] + int_sol[10,end]
# end
# tt()
# prob1 = ODEProblem(simp, u0, tspan, [population[1];[dose_vector=>ones(90).*5000.0]])
# prob = remake(prob1, p=[p[1:end-1];[dose_vector=>ones(90).*5000.0]])
# prob.p
# int_sol = solve(prob)
# int_sol[9,end] + int_sol[10,end]

# function loss(θ,params)
#     # prob = remake(oprob, p=[params;θ])
#     prob = ODEProblem(simp, u0, tspan, [p[1:end-1];[dose_vector=>ones(90).*2000.0]])
#     println(prob.p)
#     int_sol = solve(prob)

#     scale = params[end][end]
#     # println(scale)
#     cell = (int_sol[9,end]-scale[1])/(scale[2]-scale[1])
#     drug = (int_sol[10,end]-scale[3])/(scale[4]-scale[3])
#     return cell + drug
# end

# loss([dose_vector=>[1:90...].*5000.0], p[1:end-1])

# ForwardDiff.gradient(θ -> loss(θ, p), doses)


