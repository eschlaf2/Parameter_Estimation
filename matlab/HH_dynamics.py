def HH_dynamics(state):

    import numpy as np

    V, n, h, B, gB, EB, VBth, SB, tauB, Icur = state

    F = np.zeros_like(state)

    # Set parameters

    mInf = 1 / (1 + np.exp((-V - 34.5) / 10))  # HH options
    nInf = 1 / (1 + np.exp((-V - 29.5) / 10))
    hInf = 1 / (1 + np.exp((V + 59.4) / 10.7))
    tauN = 0.25 + 4.35 * np.exp(-abs(V + 10) / 10)
    tauH = 0.15 + 1.15 / (1 + np.exp((V + 33.5) / 15))
    BInf = 1. / (1 + np.exp(-(V - VBth) / SB))
    tauB = tauB
    C = 0.9
    EK = -95
    ENa = 50
    EL = -70
    gNa = 100
    gK = 7
    gL = 0.25

    # Calculate changes
    Vdot = (Icur -  # drive current
            gK * n**4 * (V - EK) -  # Potassium
            gNa * mInf ** 3 * h * (V - ENa) -  # Sodium
            gB * B * (V - EB) -  # mystery
            gL * (V - EL)) / C  # Leak

    ndot = (nInf - n) / tauN  # F2(n)

    hdot = (hInf - h) / tauH  # F3(h)

    Bdot = (BInf - B) / tauB  # F4(B)

    F[1: 4, :] = [Vdot, ndot, hdot, Bdot]  # state change
