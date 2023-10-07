using ModelingToolkit
using OrdinaryDiffEq, Plots, DifferentialEquations
using BenchmarkTools
using ForwardDiff
using Serialization
using LinearAlgebra

# utils includes a function for creating dosing schedules (spaced_list) and a function for plotting the simulation (plotter)
include("utils.jl");

par = [
    @parameters gamma_1 = 0.42676 [tunable = true, bounds = (0, Inf)]
    @parameters psi = 1.0 [tunable = true, bounds = (0, Inf)]
    @parameters C0 = 17.7 [tunable = true, bounds = (0, Inf)]
    @parameters D0 = 0.0 [tunable = true, bounds = (0, Inf)]
    @parameters r = 0.007545/24 [tunable = true, bounds = (0, Inf)]
    @parameters K = 158.04 [tunable = true, bounds = (0, Inf)]
    @parameters BW = 70.0 [tunable = true, bounds = (0, Inf)]
    @parameters IC50_1 = 15.5936*0.000001*194.151 [tunable = true, bounds = (0, Inf)]
    @parameters Imax_1 = 1.1026 [tunable = true, bounds = (0, Inf)]
    @parameters IC50_2 = 0.00924 [tunable = true, bounds = (0, Inf)]
    @parameters gamma_2 = 2.712 [tunable = true, bounds = (0, Inf)]
    @parameters Imax_2 = 1.0 [tunable = true, bounds = (0, Inf)]
    @parameters xi = IC50_1/IC50_2 [tunable = true, bounds = (0, Inf)]
    @parameters VD1 = 30.3 [tunable = true, bounds = (0, Inf)]
    @parameters Cl1 = 10.5705229706946 [tunable = true, bounds = (0, Inf)]
    @parameters k23 = 0.000823753578728557 [tunable = true, bounds = (0, Inf)]
    @parameters ka1 = 9.75543575220511 [tunable = true, bounds = (0, Inf)]
    @parameters k32 = 0.76 [tunable = true, bounds = (0, Inf)]
    @parameters Cl2 = 32.6682 [tunable = true, bounds = (0, Inf)]
    @parameters ka2 = 0.0385233 [tunable = true, bounds = (0, Inf)]
    @parameters Vpla = 0.934662 [tunable = true, bounds = (0, Inf)]
    @parameters Q = 0.0302696 [tunable = true, bounds = (0, Inf)]
    @parameters Vtis = 0.00299745 [tunable = true, bounds = (0, Inf)]
    @parameters scaling = [0.1,1.0,0.1,1.0] [description = "min and max scaling values for the loss function terms (arbitrary here)"]
    @parameters dose_vector = fill(1000.0,90)
]

sta = [
    @variables t = 0.0 [bounds=(0,Inf)]
    @variables C(t) = 17.7 [bounds=(0,Inf)]
    @variables D(t) = 0.0 [bounds=(0,Inf)]
    @variables AbsTMZ(t) = 0.0 [bounds=(0,Inf)]
    @variables PlaTMZ(t) = 0.0 [bounds=(0,Inf)]
    @variables CSFTMZ(t) = 0.0 [bounds=(0,Inf)]
    @variables AbsRG(t) = 0.0 [bounds=(0,Inf)]
    @variables PlaRG(t) = 0.0 [bounds=(0,Inf)]
    @variables TisRG(t) = 0.0 [bounds=(0,Inf)]
]

Diff = Differential(t)

# this function allows us to circumvent Domain errors that occur when dosing events are encountered
# function erelu(x)
#     return ifelse(x >= 0, x, eps())
# end
# @register erelu(x)

eqns = [
    # CSFTMZ ~ erelu(CSFTMZ),
    # _PlaRG ~ erelu(PlaRG),
    # _C ~ erelu(C),  #fractional power in denominator of Hill(see complex part)

    Diff(C) ~ C * r * log(K / C) - (((Imax_1 * (((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + Imax_2 * ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + (Imax_1 + Imax_2) * ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) - Imax_1 * Imax_2 * ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2))) / ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) + 1)) * (K * exp((log(C0/K)) * exp(-r * ((-log(log(C/K) / log(C0/K)) / r) + 72)))) / (72 * C)) * C,
    Diff(D) ~ (((Imax_1 * (((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + Imax_2 * ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + (Imax_1 + Imax_2) * ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) - Imax_1 * Imax_2 * ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2))) / ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) + 1)) * (K * exp((log(C0/K)) * exp(-r * ((-log(log(C/K) / log(C0/K)) / r) + 72)))) / (72 * C)) * C,

    Diff(AbsTMZ) ~ -ka1 * AbsTMZ,
    Diff(PlaTMZ) ~ ka1 * AbsTMZ - (Cl1 / VD1) * PlaTMZ - k23 * PlaTMZ + k32 * CSFTMZ,
    Diff(CSFTMZ) ~ k23 * PlaTMZ - k32 * CSFTMZ,
    Diff(AbsRG) ~ -ka2 * AbsRG,
    Diff(PlaRG) ~ ka2 * AbsRG - (Cl2 / Vpla) * PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * PlaRG,
    Diff(TisRG) ~ -(Q / Vtis) * TisRG + (Q / Vpla) * PlaRG,
    
    #Alternative way to compute the AUC
    # Diff(cAUC) ~ C, 
    # Diff(pRGauc) ~ PlaRG
]

# time frames of the trial
hours = 24.0
end_time = 28.0*5.0
end_treat = 42.0
tspan = (0.0, end_time+7.0) .* hours

#uses function from utils.jl to create staggered/intermittent periodic schedule
tmz_treat_dosetimes = spaced_list(end_treat,1.0,0.0,0.0).*hours
tmz_adjuv_dosetimes = spaced_list(end_time,5.0,23.0,end_treat+28.0).*hours
rg_dosetimes = spaced_list(end_time-1.0,18.0,10.0,0.0).*hours

#dose amounts
avg_human_surface_area = 1.7 #m^2

#Standard of care and trialed dose_amounts
tmz_treat_dose = 75.0*avg_human_surface_area
tmz_adjuv_dose = 150.0*avg_human_surface_area
dose_amount = 1800.0*avg_human_surface_area

#affect function for the discontinuous events
#DRUG 1
function tmz_dose!(integ, u, p, ctx)
    SciMLBase.set_proposed_dt!(integ, 0.001)
    integ.u[u.AbsTMZ] += ctx
end

tmz_treat = tmz_treat_dosetimes => (tmz_dose!, [AbsTMZ], [], tmz_treat_dose)
tmz_adjuv = tmz_adjuv_dosetimes => (tmz_dose!, [AbsTMZ], [], tmz_adjuv_dose)

#DRUG 2
function rg_dose!(integ, u, p, ctx)
    SciMLBase.set_proposed_dt!(integ, 0.001)
    integ.u[u.AbsRG] += integ.p[end][end][ctx[integ.t]]
end

#initial dose vector set to maximum dose for all time points
doses = ones(length(rg_dosetimes)).*dose_amount
#mapping time point to index of the dose in the dose_vector
time_dose_map = Dict(zip(rg_dosetimes, Int64.(1:length(rg_dosetimes))))
rg_treat = rg_dosetimes => (rg_dose!, [AbsRG], [dose_vector], time_dose_map)

@named osys = ODESystem(eqns, t, sta[2:end], par; tspan=tspan, #excludes the independent variable t from states
                        discrete_events = [tmz_treat, tmz_adjuv, rg_treat])
#simplification removes intermediate/auxilary equations from the system
simp=structural_simplify(osys)

#trial set of parameters
ode_params = [
    0.42676, 
    1.0,          
    17.7,
    0.0,   
    0.007545/24,
    158.04,
    70.0,   
    15.5936*0.000001*194.151,
    1.1026,
    0.00924,
    2.712,
    1.0,
    15.5936*0.000001*194.151 / 0.00924,
    30.3,
    10.5705229706946,
    0.000823753578728557,
    9.75543575220511,     
    0.76,
    32.6682,
    0.0385233,
    0.934662,
    0.0302696,
    0.00299745
]
# #initial conditions 
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

# specifying a parameter structure to consistently use in the loss function
scale = [0.1,1.0,0.1,1.0]
doseit = fill(1500.0,90)
part = [ode_params;[scale];[doseit]]
prob=ODEProblem(simp, u0, tspan, part)#, jac=true) errors with jacobian

#inject times holds all the points where discontinuities exist
inject_times = sort(unique([rg_dosetimes;tmz_treat_dosetimes;tmz_adjuv_dosetimes]))

# loss function balances the tumor auc and drug auc using trapez from utils.jl
function loss(θ,params)
    int_prob = remake(prob, p=(params[1:end-1],[params[end],θ]))
    int_sol = solve(int_prob, dense=false, d_discontinuities=inject_times, dtmax=2)
    cell = (trapez(int_sol.t,int_sol[1,:])- params[end][1])/( params[end][2]- params[end][1])
    drug = (trapez(int_sol.t,int_sol[7,:])- params[end][3])/( params[end][4]- params[end][3])
    return cell + drug
end

pop_p = [ode_params;[scale]]
loss(doses, ode_params)
# @btime loss(doses, ode_params)

#getting gradient for the dose vector initializing at the maximum dose
ForwardDiff.gradient(θ -> loss(θ, pop_p), doses)
# @btime ForwardDiff.gradient(θ -> loss(θ, ode_params), doses)

#@info "Starting optimization"
# adtype = Optimization.AutoForwardDiff()
# optf = Optimization.OptimizationFunction((x, p) -> loss(x, p), adtype)
# optprob = Optimization.OptimizationProblem(optf, doses, ode_params, lb=zeros(length(doses)), ub=doses)
# t_iter=100
# callback, losses, iter, anim = create_callback(t_iter,verbose=false, animate=true, progress_bar=true, saving=false, early_stopping=true)
# opt = Optim.LBFGS()
# res = Optimization.solve(optprob, opt, maxiters=t_iter, callback=callback)