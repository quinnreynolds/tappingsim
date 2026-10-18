"""Models of the pressure drop through the packed burden at a tap-hole.

The furnace models are configured with one of these functions to describe the
resistance to flow as the molten phases pass through the burden and enter the
tap-hole. Each takes the dimensions of the tap-hole and of the burden, the
properties of the burden particles, and the viscosity and density of the fluid,
and returns the coefficients of the pressure drop dP = A*u^2 + B*u, where u is
the fluid velocity.

"""


def bedmodel_kozenycarman(tapholediameter, bedmindiameter, bedmaxdiameter,
                          bedparticlediameter, bedparticlesphericity,
                          bedporosity, viscosity, density):
    """Function to calculate pressure drop through packed burden as flow enters
    the tap-hole.

    Parameters
    ----------
    tapholediameter : float
        Diameter of tap-hole channel, m.
    bedmindiameter : float
        Diameter of cavity in burden in front of tap-hole, m (set to zero
        for no cavity).
    bedmaxdiameter : float
        Extent of burden from tap-hole, m (set to a high value for entire
        furnace).
    bedparticlediameter : float
        Diameter of the constituent particles of the burden layer, m.
    bedparticlesphericity : float
        Sphericity of the constituent particles of the burden layer.
    bedporosity : float
        Porosity of the burden layer.
    viscosity : float
        Viscosity of the fluid, Pa.s.
    density : float
        Density of the fluid, kg/m3.

    Returns
    -------
    tuple of float
        A and B constants in pressure drop expression.

    Note
    ----
    dP = A*u^2 + B*u, where u is the fluid velocity. The Kozeny-Carman
    correlation is valid for laminar flow conditions only. [Kozeny1927]_
    [Carman1937]_

    References
    ----------
    .. [Kozeny1927] J. Kozeny. Ueber kapillare Leitung des Wassers im Boden.
       Sitzungsber Akad. Wiss., Wien, 136(2a): 271-306, 1927.
    .. [Carman1937] P.C. Carman. Fluid flow through granular beds. Transactions,
       Institution of Chemical Engineers, London, 15: 150-166, 1937.

    """
    rt, rmin = 0.5 * tapholediameter, 0.5 * bedmindiameter
    eff_d = bedparticlediameter * bedparticlesphericity
    rrmax = rt / (0.5 * bedmaxdiameter)
    if rmin > rt:
        rrmin = rt / rmin
        bconst = (180 * viscosity * rt * (1 - bedporosity) * (1 - bedporosity)
                  / (2 * eff_d * eff_d * bedporosity * bedporosity * bedporosity)
                  * (rrmin - rrmax))
    else:
        rrmin = rmin / rt
        bconst = (180 * viscosity * rt * (1 - bedporosity) * (1 - bedporosity)
                  / (2 * eff_d * eff_d * bedporosity * bedporosity * bedporosity)
                  * (2 - rrmin - rrmax))
    return 0, bconst


def bedmodel_ergun(tapholediameter, bedmindiameter, bedmaxdiameter,
                   bedparticlediameter, bedparticlesphericity, bedporosity,
                   viscosity, density):
    """Function to calculate pressure drop through packed burden as flow enters
    the tap-hole.

    Parameters
    ----------
    tapholediameter : float
        Diameter of tap-hole channel, m.
    bedmindiameter : float
        Diameter of cavity in burden in front of tap-hole, m (set to zero
        for no cavity).
    bedmaxdiameter : float
        Extent of burden from tap-hole, m (set to a high value for entire
        furnace).
    bedparticlediameter : float
        Diameter of the constituent particles of the burden layer, m.
    bedparticlesphericity : float
        Sphericity of the constituent particles of the burden layer.
    bedporosity : float
        Porosity of the burden layer.
    viscosity : float
        Viscosity of the fluid, Pa.s.
    density : float
        Density of the fluid, kg/m3.

    Returns
    -------
    tuple of float
        A and B constants in pressure drop expression.

    Note
    ----
    dP = A*u^2 + B*u, where u is the fluid velocity. The Ergun correlation is
    valid for both laminar and turbulent flow regimes. [Ergun1952]_

    References
    ----------
    .. [Ergun1952] S. Ergun. Fluid flow through packed columns. Chem. Eng.
       Prog. 48: 89-94, 1952.
    """
    rt, rmin = 0.5 * tapholediameter, 0.5 * bedmindiameter
    eff_d = bedparticlediameter * bedparticlesphericity
    rrmax = rt / (0.5 * bedmaxdiameter)
    cmult = rt * (1 - bedporosity) / (eff_d * eff_d *
                                      bedporosity * bedporosity * bedporosity)
    amult = 0.145833333 * density * eff_d * cmult
    bmult = 75 * viscosity * (1 - bedporosity) * cmult
    if rmin > rt:
        rrmin = rt / rmin
        aconst = amult * (rrmin * rrmin * rrmin - rrmax * rrmax * rrmax)
        bconst = bmult * (rrmin - rrmax)
    else:
        rrmin = rmin / rt
        aconst = amult * (4 - 3 * rrmin * rrmin *
                          rrmin - rrmax * rrmax * rrmax)
        bconst = bmult * (2 - rrmin - rrmax)
    return aconst, bconst
