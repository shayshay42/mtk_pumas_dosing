using ModelingToolkit
using OrdinaryDiffEq, Plots, DifferentialEquations
using BenchmarkTools
using DiffEqCallbacks

# utils includes a function for creating dosing schedules (spaced_list) and a function for plotting the simulation (plotter)
include("assets/utils.jl");

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

@variables t = 0.0 [bounds=(0,Inf)]
@variables C(t) = 17.7 [bounds=(0,Inf)]
@variables D(t) = 0.0 [bounds=(0,Inf)]
@variables AbsTMZ(t) = 0.0 [bounds=(0,Inf)]
@variables PlaTMZ(t) = 0.0 [bounds=(0,Inf)]
@variables CSFTMZ(t) = 0.0 [bounds=(0,Inf)]
@variables AbsRG(t) = 0.0 [bounds=(0,Inf)]
@variables PlaRG(t) = 0.0 [bounds=(0,Inf)]
@variables TisRG(t) = 0.0 [bounds=(0,Inf)]

Diff = Differential(t)

# this function allows us to circumvent Domain errors that occur when dosing events are encountered
# function erelu(x)
#     return ifelse(x >= 0, x, eps())
# end
# @register erelu(x)

#stiff awareness stuff and complexity
eqns = [
    # CSFTMZ ~ erelu(CSFTMZ),
    # _PlaRG ~ erelu(PlaRG),
    # _C ~ erelu(C),  #check is because fractional power in denominator so see complex part

    Diff(C) ~ C * r * log(K / C) - (((Imax_1 * (((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + Imax_2 * ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + (Imax_1 + Imax_2) * ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) - Imax_1 * Imax_2 * ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2))) / ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) + 1)) * (K * exp((log(C0/K)) * exp(-r * ((-log(log(C/K) / log(C0/K)) / r) + 72)))) / (72 * C)) * C,
    Diff(D) ~ (((Imax_1 * (((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + Imax_2 * ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + (Imax_1 + Imax_2) * ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) - Imax_1 * Imax_2 * ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2))) / ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) + ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2) + ((((CSFTMZ / 140)/(psi * IC50_1))^gamma_1) * ((((IC50_1/IC50_2) * (PlaRG / (1000 * Vpla)))/(psi * IC50_1))^gamma_2)) + 1)) * (K * exp((log(C0/K)) * exp(-r * ((-log(log(C/K) / log(C0/K)) / r) + 72)))) / (72 * C)) * C,

    Diff(AbsTMZ) ~ -ka1 * AbsTMZ,
    Diff(PlaTMZ) ~ ka1 * AbsTMZ - (Cl1 / VD1) * PlaTMZ - k23 * PlaTMZ + k32 * CSFTMZ,
    Diff(CSFTMZ) ~ k23 * PlaTMZ - k32 * CSFTMZ,
    Diff(AbsRG) ~ -ka2 * AbsRG,
    Diff(PlaRG) ~ ka2 * AbsRG - (Cl2 / Vpla) * PlaRG + (Q / Vtis) * TisRG - (Q / Vpla) * PlaRG,
    Diff(TisRG) ~ -(Q / Vtis) * TisRG + (Q / Vpla) * PlaRG,
    
    # Diff(cAUC) ~ C, #AUC loss with fewer bits of precision
    # Diff(pRGauc) ~ PlaRG
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
    SciMLBase.set_proposed_dt!(integ, 0.1)
    # integ.opts.abstol = 1e-9
    integ.u[u.AbsTMZ] += ctx
end

tmz_treat = tmz_treat_dosetimes => (tmz_dose!, [AbsTMZ], [], tmz_treat_dose)
tmz_adjuv = tmz_adjuv_dosetimes => (tmz_dose!, [AbsTMZ], [], tmz_adjuv_dose)

function rg_dose!(integ, u, p, ctx)
    # SciMLBase.set_proposed_dt!(integ, 0.01)
    integ.u[u.AbsRG] += integ.p[end][end][ctx[integ.t]]
end

doses = ones(length(rg_dosetimes)).*dose_amount
time_dose_map = Dict(zip(rg_dosetimes, Int64.(1:length(rg_dosetimes))))
rg_treat = rg_dosetimes => (rg_dose!, [AbsRG], [], time_dose_map)

@named osys = ODESystem(eqns; tspan=tspan, #excludes the independent variable t from states
                        discrete_events = [tmz_treat])#, tmz_adjuv, rg_treat])
# parameters(osys)
#simplification removes intermediate/auxilary equations from the system
simp=structural_simplify(osys)

prob=ODEProblem(simp) #jac=true errors

choice_function(integrator) = (integrator.t in tmz_treat_dosetimes) ? 1 : 2
alg_switch = CompositeAlgorithm((Tsit5(), Rosenbrock23()), choice_function)
sol = solve(prob, AutoTsit5(Rodas4P()), callback=PositiveDomain())
sol = solve(prob, AutoTsit5(KenCarp4()), d_discontinuities=tmz_treat_dosetimes, dtmax=1)
@btime sol = solve(prob, AutoTsit5(Rodas4P()), callback=PositiveDomain())
@btime sol = solve(prob, callback=PositiveDomain())
@btime sol = solve(prob, AutoTsit5(KenCarp4()), d_discontinuities=tmz_treat_dosetimes, dtmax=4, saveat=1)
#13.071 ms (3710 allocations: 784.08 KiB)

# sol = solve(prob) - ERROR: DomainError with -169977.60499323878:
# Exponentiation yielding a complex result requires a complex argument.
# Replace x^y with (x+0im)^y, Complex(x)^y, or similar.

sol = solve(prob, alg_hints=[:stiff]) #, saveat=0.1)
# sol = solve(prob, alg_hints=[:stiff])
#same as above

sol = solve(prob, isoutofdomain=(u,p,t)->any(x->x<0,u), abstol=1e-9)

   plotter(sol)
