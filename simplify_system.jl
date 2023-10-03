    cCSFTMZ ~ CSFTMZ / 140,
    cPlaRG ~ PlaRG / (1000 * Vpla),

    # #domain error
    # cCSFTMZ = erelu(cCSFTMZ)
    # cPlaRG = erelu(cPlaRG)
    # C = erelu(C)

    pi1 ~ psi * IC50_1,
    exp1 ~ ((CSFTMZ / 140)/pi1)^gamma_1,
    exp2 ~ ((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2,
    exp12 ~ (((CSFTMZ / 140)/pi1)^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2),

    Imax_sum ~ Imax_1 + Imax_2,
    E_num ~ Imax_1 * (((CSFTMZ / 140)/pi1)^gamma_1) + Imax_2 * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2) + (Imax_1 + Imax_2) * ((((CSFTMZ / 140)/pi1)^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2)) - Imax_1 * Imax_2 * ((((CSFTMZ / 140)/pi1)^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2)),
    E_den ~ (((CSFTMZ / 140)/pi1)^gamma_1) + (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2) + ((((CSFTMZ / 140)/pi1)^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2)) + 1,
    E ~ (Imax_1 * (((CSFTMZ / 140)/pi1)^gamma_1) + Imax_2 * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2) + (Imax_1 + Imax_2) * ((((CSFTMZ / 140)/pi1)^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2)) - Imax_1 * Imax_2 * ((((CSFTMZ / 140)/pi1)^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2))) / ((((CSFTMZ / 140)/pi1)^gamma_1) + (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2) + ((((CSFTMZ / 140)/pi1)^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2)) + 1),

    t1 ~ -log(log(C/K) / log(C0/K)) / r,
    t2 ~ (-log(log(C/K) / log(C0/K)) / r) + 72,
    log_C0K ~ log(C0/K),
    fun ~ K * exp((log(C0/K)) * exp(-r * ((-log(log(C/K) / log(C0/K)) / r) + 72))),
    delta ~ ((Imax_1 * (((CSFTMZ / 140)/pi1)^gamma_1) + Imax_2 * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2) + (Imax_1 + Imax_2) * ((((CSFTMZ / 140)/pi1)^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2)) - Imax_1 * Imax_2 * ((((CSFTMZ / 140)/pi1)^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2))) / ((((CSFTMZ / 140)/pi1)^gamma_1) + (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2) + ((((CSFTMZ / 140)/pi1)^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2)) + 1)) * (K * exp((log(C0/K)) * exp(-r * ((-log(log(C/K) / log(C0/K)) / r) + 72)))) / (72 * C),

    D(C) ~ C * r * log(K / C) - (((Imax_1 * (((CSFTMZ / 140)/pi1)^gamma_1) + Imax_2 * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2) + (Imax_1 + Imax_2) * ((((CSFTMZ / 140)/pi1)^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2)) - Imax_1 * Imax_2 * ((((CSFTMZ / 140)/pi1)^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2))) / ((((CSFTMZ / 140)/pi1)^gamma_1) + (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2) + ((((CSFTMZ / 140)/pi1)^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2)) + 1)) * (K * exp((log(C0/K)) * exp(-r * ((-log(log(C/K) / log(C0/K)) / r) + 72)))) / (72 * C)) * C,
    D(D) ~ (((Imax_1 * (((CSFTMZ / 140)/pi1)^gamma_1) + Imax_2 * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2) + (Imax_1 + Imax_2) * ((((CSFTMZ / 140)/pi1)^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2)) - Imax_1 * Imax_2 * ((((CSFTMZ / 140)/pi1)^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2))) / ((((CSFTMZ / 140)/pi1)^gamma_1) + (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2) + ((((CSFTMZ / 140)/pi1)^gamma_1) * (((xi * (PlaRG / (1000 * Vpla)))/pi1)^gamma_2)) + 1)) * (K * exp((log(C0/K)) * exp(-r * ((-log(log(C/K) / log(C0/K)) / r) + 72)))) / (72 * C)) * C,