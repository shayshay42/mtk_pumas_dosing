using ModelingToolkit, OrdinaryDiffEq, Test

@parameters t α β γ δ
@variables x(t) y(t)
D = Differential(t)

eqs = [D(x) ~ α*x - β*x*y,
       D(y) ~ -δ*y + γ*x*y]

@named sys = ODESystem(eqs)
@named sys2 = liouville_transform(sys)
@variables trJ

u0 = [x => 1.0,
      y => 1.0,
      trJ => 1.0]

prob = ODEProblem(sys,u0,tspan,p)
sol = solve(prob,Tsit5())