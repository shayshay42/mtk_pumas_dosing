using PumasQSP
using Catalyst
using Plots

rn = @reaction_network begin
    # Define model
    @species skin(t) = 0.0 body(t) = 0.0
    @parameters kabs = 2.0 kelim = 1.0
    kabs, skin --> body
    kelim, body --> ∅
end
model  = convert(ODESystem, rn)

  # Unpack parameters and define doses
@unpack skin = model
bolus_dose_1 = Bolus(skin, 1.0, 0.25)
bolus_dose_2 = Bolus(skin, 0.5, 0.75)
bolus_doses = [bolus_dose_1, bolus_dose_2]

trial = Trial(model; doses = bolus_doses, saveat = 0.01,tspan = (0.0, 1.5))
prob = InverseProblem(trial, model, [])
sol = simulate(trial, prob)
plot(sol)
