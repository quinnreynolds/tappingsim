"""Classes to create transfer launder objects.

Launders are typically open-channel gravity flow conduits which redirect the
flow of tapped material from the furnace into storage ladles or other
post-treatment units.

"""

import math

import numpy
from scipy.constants import g


def _sqrt(x):
    """Scalar square root, returning nan for negative input as numpy.sqrt does, but
    several times faster than numpy.sqrt for scalar arguments."""
    return math.sqrt(x) if x >= 0 else math.nan


class SimpleSisoLaunder:
    """Tapping launder class. This version is a simple single-in-single-out
    launder channel with no accumulation or spillage effects.

    Attributes
    ----------
    vdotmetal_out : float
        The output metal flowrate from the launder in the current state,
        m3/s.
    vdotslag_out : float
        The output slag flowrate from the launder in the current state,
        m3/s.

    """

    def __init__(self):
        self.vdotmetal_out = 0
        self.vdotslag_out = 0

    def calc_dt(self, dt, vdotmetal_in, vdotslag_in):
        """
        Integrate object state forward over a single time step, and update
        output flowrates.

        Parameters
        ----------
        dt : float
            Length of time step, in s.
        vdotmetal_in : float
            Flowrate of metal into the launder, m3/s.
        vdotslag_in : float
            Flowrate of slag into the launder, m3/s.

        Returns
        -------
        None.

        """
        self.vdotmetal_out = vdotmetal_in
        self.vdotslag_out = vdotslag_in


class SkimmerLaunder:
    """Tapping launder class. This version is a two-tank model of a launder feeding into a
    pre-weir slag trough, through an underflow weir, and then into a post-weir metal trough.
    Open-channel gravity overflow outlets exit from the top of each trough section.

    The trough sections are assumed to have a trapezoidal cross section, and are defined by
    length, height, and lower and upper trough widths.

    Parameters
    ----------
    slag_trough_length : float
        The internal length of the pre-weir (slag) trough, m.
    slag_trough_upper_width : float
        The internal width of the slag trough at the top of the trough, m.
    slag_trough_lower_width : float
        The internal width of the slag trough at the floor of the trough, m.
    slag_trough_height : float
        The internal height of the slag trough, m.
    slag_trough_outlet_height : float
        The height of the slag trough overflow outlet crest above the trough floor, m.
    slag_trough_outlet_width : float
        The width of the slag trough overflow outlet, m.
    weir_underflow_height : float
        The height of the underflow channel under the weir, measured from the trough floor, m.
    weir_underflow_width : float
        The width of the underflow channel under the weir, m.
    metal_trough_length : float
        The internal length of the post-weir (metal) trough, m.
    metal_trough_upper_width : float
        The internal width of the metal trough at the top of the trough, m.
    metal_trough_lower_width : float
        The internal width of the metal trough at the floor of the trough, m.
    metal_trough_height : float
        The internal height of the metal trough, m.
    metal_trough_outlet_height : float
        The height of the metal trough overflow outlet crest above the trough floor, m.
    metal_trough_outlet_width : float
        The width of the metal trough overflow outlet, m.
    slag_trough_outlet_slag_Cd : float
        Discharge coefficient for the slag phase leaving through the slag trough outlet,
        dimensionless.
    slag_trough_outlet_metal_Cd : float
        Discharge coefficient for the metal phase leaving through the slag trough outlet,
        dimensionless.
    metal_trough_outlet_slag_Cd : float
        Discharge coefficient for the slag phase leaving through the metal trough outlet,
        dimensionless.
    metal_trough_outlet_metal_Cd : float
        Discharge coefficient for the metal phase leaving through the metal trough outlet,
        dimensionless.
    density_metal : float
        The density of the metal phase, kg/m3.
    density_slag : float
        The density of the slag phase, kg/m3.
    max_substep : float, optional
        The longest internal time step used by calc_dt, in s. Longer time steps are divided
        into equal sub-steps no longer than this. The default is 1 s.

    Attributes
    ----------
    slag_trough_length : float
        The internal length of the pre-weir (slag) trough, m.
    slag_trough_upper_width : float
        The internal width of the slag trough at the top, m.
    slag_trough_lower_width : float
        The internal width of the slag trough at the floor, m.
    slag_trough_height : float
        The internal height of the slag trough, m.
    slag_trough_outlet_height : float
        The height of the slag trough outlet crest above the floor, m.
    slag_trough_outlet_width : float
        The width of the slag trough overflow outlet, m.
    weir_underflow_height : float
        The height of the underflow channel under the weir, m.
    weir_underflow_width : float
        The width of the underflow channel under the weir, m.
    metal_trough_length : float
        The internal length of the post-weir (metal) trough, m.
    metal_trough_upper_width : float
        The internal width of the metal trough at the top, m.
    metal_trough_lower_width : float
        The internal width of the metal trough at the floor, m.
    metal_trough_height : float
        The internal height of the metal trough, m.
    metal_trough_outlet_height : float
        The height of the metal trough outlet crest above the floor, m.
    metal_trough_outlet_width : float
        The width of the metal trough overflow outlet, m.
    slag_trough_outlet_slag_Cd : float
        Slag discharge coefficient of the slag trough outlet.
    slag_trough_outlet_metal_Cd : float
        Metal discharge coefficient of the slag trough outlet.
    metal_trough_outlet_slag_Cd : float
        Slag discharge coefficient of the metal trough outlet.
    metal_trough_outlet_metal_Cd : float
        Metal discharge coefficient of the metal trough outlet.
    density_metal : float
        The density of the metal phase, kg/m3.
    density_slag : float
        The density of the slag phase, kg/m3.
    max_substep : float
        The longest internal time step used by calc_dt, in s.
    slag_trough_hmetal : float
        The current level of metal in the slag trough, m.
    slag_trough_hslag : float
        The current level of the top of the slag layer in the slag trough, m. Levels are
        measured from the trough floor, so this is the total bath depth and not the slag
        layer thickness.
    metal_trough_hmetal : float
        The current level of metal in the metal trough, m.
    metal_trough_hslag : float
        The current level of the top of the slag layer in the metal trough, m. Measured from
        the trough floor, as for the slag trough.
    slag_trough_vdotmetal_out : float
        The current outlet volume flowrate of metal from the slag trough overflow, m3/s.
    slag_trough_vdotslag_out : float
        The current outlet volume flowrate of slag from the slag trough overflow, m3/s.
    metal_trough_vdotmetal_out : float
        The current outlet volume flowrate of metal from the metal trough overflow, m3/s.
    metal_trough_vdotslag_out : float
        The current outlet volume flowrate of slag from the metal trough overflow, m3/s.
    outlet_premult : float
        Constant pre-multiplier, 2/3*sqrt(g), in the overflow outlet flow expressions.

    Note
    ----
    Both troughs are initially empty. Material enters the slag trough, and flows to the metal
    trough through the underflow channel under the weir, driven by the difference in
    hydrostatic pressure between the two troughs. Each trough drains through its own
    open-channel overflow outlet, which is modelled as a broad-crested weir with flowrate
    proportional to the overflow head to the power 3/2. Metal entrainment in the overflow
    streams is not yet modelled, so the slag outflows are currently free of metal.

    The state is integrated in time with a linearly-implicit (Rosenbrock-Euler) scheme, which
    remains stable for time steps far beyond the limit of explicit schemes. That limit is
    proportional to the throughput and tends to zero as the flow does, because the flow
    through the weir underflow channel varies as the square root of the pressure difference.
    The step size is therefore limited by accuracy, not stability, and the scheme is first
    order. See calc_dt.

    """

    # Mass balance coefficients linking the flows to the rates of change of the phase
    # volumes. Rows: metal and slag volume in the slag trough, then in the metal trough.
    # Columns: weir metal flow, weir slag flow, then metal and slag outflow from the slag
    # trough, then metal and slag outflow from the metal trough.
    _FLOWS_TO_VOLUME_RATES = numpy.array(
        [
            [-1, 0, -1, 0, 0, 0],
            [0, -1, 0, -1, 0, 0],
            [1, 0, 0, 0, -1, 0],
            [0, 1, 0, 0, 0, -1],
        ],
        dtype=float,
    )

    # Volume perturbation used to estimate the sensitivity of the flows to the state, m3.
    _JACOBIAN_EPS = 1e-7

    def __init__(
        self,
        # pre-weir
        slag_trough_length,
        slag_trough_upper_width,
        slag_trough_lower_width,
        slag_trough_height,
        slag_trough_outlet_height,
        slag_trough_outlet_width,
        # weir
        weir_underflow_height,
        weir_underflow_width,
        # post-weir
        metal_trough_length,
        metal_trough_upper_width,
        metal_trough_lower_width,
        metal_trough_height,
        metal_trough_outlet_height,
        metal_trough_outlet_width,
        # overflows
        slag_trough_outlet_slag_Cd,
        slag_trough_outlet_metal_Cd,
        metal_trough_outlet_slag_Cd,
        metal_trough_outlet_metal_Cd,
        # properties
        density_metal,
        density_slag,
        # numerics
        max_substep=1.0,
    ):
        if max_substep <= 0:
            raise ValueError("max_substep must be positive")
        self.slag_trough_length = slag_trough_length
        self.slag_trough_upper_width = slag_trough_upper_width
        self.slag_trough_lower_width = slag_trough_lower_width
        self.slag_trough_height = slag_trough_height
        self.slag_trough_outlet_height = slag_trough_outlet_height
        self.slag_trough_outlet_width = slag_trough_outlet_width
        self.weir_underflow_height = weir_underflow_height
        self.weir_underflow_width = weir_underflow_width
        self.metal_trough_length = metal_trough_length
        self.metal_trough_upper_width = metal_trough_upper_width
        self.metal_trough_lower_width = metal_trough_lower_width
        self.metal_trough_height = metal_trough_height
        self.metal_trough_outlet_height = metal_trough_outlet_height
        self.metal_trough_outlet_width = metal_trough_outlet_width
        self.slag_trough_outlet_slag_Cd = slag_trough_outlet_slag_Cd
        self.slag_trough_outlet_metal_Cd = slag_trough_outlet_metal_Cd
        self.metal_trough_outlet_slag_Cd = metal_trough_outlet_slag_Cd
        self.metal_trough_outlet_metal_Cd = metal_trough_outlet_metal_Cd
        self.density_metal = density_metal
        self.density_slag = density_slag
        self.metal_trough_vdotmetal_out = 0
        self.metal_trough_vdotslag_out = 0
        self.slag_trough_vdotmetal_out = 0
        self.slag_trough_vdotslag_out = 0
        self.metal_trough_hmetal = 0
        self.metal_trough_hslag = 0
        self.slag_trough_hmetal = 0
        self.slag_trough_hslag = 0
        self.max_substep = max_substep
        self.outlet_premult = 2 / 3 * numpy.sqrt(g)

    def _levels_to_volumes(self, hm_ss, hs_ss, hm_ms, hs_ms):
        """Convert phase levels in the two troughs to phase volumes.

        The trough width varies linearly with height, so the cross-sectional area up to a
        level h is h*(W0 + 0.5*dW*h/H), where W0 is the lower width, dW is the difference
        between upper and lower widths, and H is the trough height.

        Parameters
        ----------
        hm_ss : float
            Level of metal in the slag trough, m.
        hs_ss : float
            Level of the top of the slag layer in the slag trough, m.
        hm_ms : float
            Level of metal in the metal trough, m.
        hs_ms : float
            Level of the top of the slag layer in the metal trough, m.

        Returns
        -------
        tuple of float
            Volume of metal in the slag trough, volume of slag in the slag trough, volume of
            metal in the metal trough, and volume of slag in the metal trough, m3.

        """
        # slag section
        H, L = self.slag_trough_height, self.slag_trough_length
        W0 = self.slag_trough_lower_width
        dW = self.slag_trough_upper_width - self.slag_trough_lower_width
        Vm_ss = L * hm_ss * (2 * H * W0 + hm_ss * dW) / (2 * H)
        Vs_ss = L * (hs_ss - hm_ss) * (2 * H * W0 + (hs_ss + hm_ss) * dW) / (2 * H)
        # metal section
        H, L = self.metal_trough_height, self.metal_trough_length
        W0 = self.metal_trough_lower_width
        dW = self.metal_trough_upper_width - self.metal_trough_lower_width
        Vm_ms = L * hm_ms * (2 * H * W0 + hm_ms * dW) / (2 * H)
        Vs_ms = L * (hs_ms - hm_ms) * (2 * H * W0 + (hs_ms + hm_ms) * dW) / (2 * H)
        return Vm_ss, Vs_ss, Vm_ms, Vs_ms

    def _volumes_to_levels(self, Vm_ss, Vs_ss, Vm_ms, Vs_ms):
        """Convert phase volumes in the two troughs to phase levels.

        This is the inverse of _levels_to_volumes, obtained by solving the quadratic
        relationship between trough volume and level. A straight-walled trough (upper and
        lower widths equal) is handled separately to avoid division by zero.

        Parameters
        ----------
        Vm_ss : float
            Volume of metal in the slag trough, m3.
        Vs_ss : float
            Volume of slag in the slag trough, m3.
        Vm_ms : float
            Volume of metal in the metal trough, m3.
        Vs_ms : float
            Volume of slag in the metal trough, m3.

        Returns
        -------
        tuple of float
            Level of metal in the slag trough, level of the top of the slag layer in the slag
            trough, level of metal in the metal trough, and level of the top of the slag
            layer in the metal trough, m.

        """
        # slag section
        H, L = self.slag_trough_height, self.slag_trough_length
        W0 = self.slag_trough_lower_width
        dW = self.slag_trough_upper_width - self.slag_trough_lower_width
        # (abs(dW) <= 1e-8 is numpy.isclose(dW, 0) without the large scalar overhead)
        if abs(dW) <= 1e-8:
            hm_ss = Vm_ss / (L * W0)
            hs_ss = hm_ss + Vs_ss / (L * W0)
        else:
            hm_ss = (_sqrt(H * (H * L * W0**2 + 2 * Vm_ss * dW) / L) - H * W0) / dW
            hs_ss = (
                _sqrt(H * (H * L * W0**2 + 2 * Vm_ss * dW + 2 * Vs_ss * dW) / L) - H * W0
            ) / dW
        # metal section
        H, L = self.metal_trough_height, self.metal_trough_length
        W0 = self.metal_trough_lower_width
        dW = self.metal_trough_upper_width - self.metal_trough_lower_width
        if abs(dW) <= 1e-8:
            hm_ms = Vm_ms / (L * W0)
            hs_ms = hm_ms + Vs_ms / (L * W0)
        else:
            hm_ms = (_sqrt(H * (H * L * W0**2 + 2 * Vm_ms * dW) / L) - H * W0) / dW
            hs_ms = (
                _sqrt(H * (H * L * W0**2 + 2 * Vm_ms * dW + 2 * Vs_ms * dW) / L) - H * W0
            ) / dW
        return hm_ss, hs_ss, hm_ms, hs_ms

    def _calc_vdot_out(self, dt):
        """Calculate volumetric outflow rates of metal and slag from the slag and metal
        trough outlets for the current state, and update the slag_trough_vdotmetal_out,
        slag_trough_vdotslag_out, metal_trough_vdotmetal_out and metal_trough_vdotslag_out
        attributes.

        Overflow from a trough occurs when the top of the slag layer exceeds the outlet
        crest height. Each phase flows over the outlet as over a broad-crested weir, with
        flowrate proportional to Cd*width*head**1.5. If the metal level is also above the
        crest, the slag layer head is the slag layer thickness, and the metal head is
        expressed as an equivalent metal head by scaling the slag layer thickness by the
        density ratio of slag to metal.
        TODO: Metal entrainment in the slag outflow is currently set to zero, pending
        implementation of overflow models as used in the ladle classes.

        Parameters
        ----------
        dt : float
            Length of time step, in s.

        Returns
        -------
        None.

        """
        # slag section
        if self.slag_trough_hslag < self.slag_trough_outlet_height:
            # slag and metal are both below outlet
            self.slag_trough_vdotmetal_out = 0
            self.slag_trough_vdotslag_out = 0
        else:
            # TODO: Implement overflow entrainment models similar to ladle models - the hardcoded
            #       vfrac = 0 here is a placehonder, a simple zero-entrainment model.
            #       Replace with e.g. vfrac = self.overflowmodel(h_overflow - h_metal)
            vfrac = 0
            if self.slag_trough_hmetal < self.slag_trough_outlet_height:
                # slag is above outlet, metal is below outlet
                # TODO need to correct slag phase density if entrained metal is present?
                h1 = self.slag_trough_hslag - self.slag_trough_outlet_height
                qslag = (
                    self.outlet_premult
                    * self.slag_trough_outlet_slag_Cd
                    * self.slag_trough_outlet_width
                    * h1**1.5
                )
                self.slag_trough_vdotmetal_out = vfrac * qslag
                self.slag_trough_vdotslag_out = (1 - vfrac) * qslag
            else:
                # slag and metal are both above outlet
                # TODO need to correct slag phase density if entrained metal is present?
                h1 = self.slag_trough_hslag - self.slag_trough_hmetal
                h2 = self.slag_trough_hmetal - self.slag_trough_outlet_height
                rho_ratio = self.density_slag / self.density_metal
                qslag = (
                    self.outlet_premult
                    * self.slag_trough_outlet_slag_Cd
                    * self.slag_trough_outlet_width
                    * h1**1.5
                )
                qmetal = (
                    self.outlet_premult
                    * self.slag_trough_outlet_metal_Cd
                    * self.slag_trough_outlet_width
                    * ((rho_ratio * h1 + h2) ** 1.5 - (rho_ratio * h1) ** 1.5)
                )
                self.slag_trough_vdotmetal_out = qmetal + vfrac * qslag
                self.slag_trough_vdotslag_out = (1 - vfrac) * qslag

        # metal section
        if self.metal_trough_hslag < self.metal_trough_outlet_height:
            # slag and metal are both below outlet
            self.metal_trough_vdotmetal_out = 0
            self.metal_trough_vdotslag_out = 0
        else:
            # TODO: Implement overflow entrainment models similar to ladle models - the hardcoded
            #       vfrac = 0 here is a placehonder, a simple zero-entrainment model.
            #       Replace with e.g. vfrac = self.overflowmodel(h_overflow - h_metal)
            vfrac = 0
            if self.metal_trough_hmetal < self.metal_trough_outlet_height:
                # slag is above outlet, metal is below outlet
                # TODO need to correct slag phase density if entrained metal is present?
                h1 = self.metal_trough_hslag - self.metal_trough_outlet_height
                qslag = (
                    self.outlet_premult
                    * self.metal_trough_outlet_slag_Cd
                    * self.metal_trough_outlet_width
                    * h1**1.5
                )
                self.metal_trough_vdotmetal_out = vfrac * qslag
                self.metal_trough_vdotslag_out = (1 - vfrac) * qslag
            else:
                # slag and metal are both above outlet
                # TODO need to correct slag phase density if entrained metal is present?
                h1 = self.metal_trough_hslag - self.metal_trough_hmetal
                h2 = self.metal_trough_hmetal - self.metal_trough_outlet_height
                rho_ratio = self.density_slag / self.density_metal
                qslag = (
                    self.outlet_premult
                    * self.metal_trough_outlet_slag_Cd
                    * self.metal_trough_outlet_width
                    * h1**1.5
                )
                qmetal = (
                    self.outlet_premult
                    * self.metal_trough_outlet_metal_Cd
                    * self.metal_trough_outlet_width
                    * ((rho_ratio * h1 + h2) ** 1.5 - (rho_ratio * h1) ** 1.5)
                )
                self.metal_trough_vdotmetal_out = qmetal + vfrac * qslag
                self.metal_trough_vdotslag_out = (1 - vfrac) * qslag

    def _calc_vdot_weir(self):
        """Calculate the volumetric flowrates of metal and slag through the weir underflow
        channel for the current state.

        The flowrates are calculated from the hydrostatic pressure difference between the
        troughs, using the Bernoulli equation for the velocity of each phase. Depending on
        the position of the slag-metal interface and the free surface relative to the top of
        the underflow channel, the channel is filled by metal, by metal and slag, or only
        partially filled.

        Returns
        -------
        tuple of float
            Metal and slag flowrates from the slag trough to the metal trough (negative if
            reversed), m3/s.

        """
        # Assumptions and conventions:
        #  - Positive flow direction is from slag section to metal section
        #  - No deformation of slag-metal interface at weir inflow
        A_w = self.weir_underflow_height * self.weir_underflow_width
        P_ss = (
            self.slag_trough_hmetal * self.density_metal * g
            + (self.slag_trough_hslag - self.slag_trough_hmetal) * self.density_slag * g
        )
        P_ms = (
            self.metal_trough_hmetal * self.density_metal * g
            + (self.metal_trough_hslag - self.metal_trough_hmetal) * self.density_slag * g
        )
        deltaP = P_ss - P_ms
        u_mw = numpy.sign(deltaP) * numpy.sqrt(2 * numpy.abs(deltaP) / self.density_metal)
        u_sw = numpy.sign(deltaP) * numpy.sqrt(2 * numpy.abs(deltaP) / self.density_slag)
        # slag-metal interface in slag trough is above weir underflow channel
        if self.slag_trough_hmetal > self.weir_underflow_height:
            vdot_metal_12 = A_w * u_mw
            vdot_slag_12 = 0
        # slag-metal interface in slag trough is below weir underflow channel
        elif self.slag_trough_hslag > self.weir_underflow_height:
            vdot_metal_12 = self.slag_trough_hmetal / self.weir_underflow_height * A_w * u_mw
            vdot_slag_12 = (1 - self.slag_trough_hmetal / self.weir_underflow_height) * A_w * u_sw
        # free surface in slag trough is below weir underflow channel
        else:
            vdot_metal_12 = self.slag_trough_hmetal / self.weir_underflow_height * A_w * u_mw
            vdot_slag_12 = (
                ((self.slag_trough_hslag - self.slag_trough_hmetal) / self.weir_underflow_height)
                * A_w
                * u_sw
            )
        return vdot_metal_12, vdot_slag_12

    def _calc_flows(self, V):
        """Calculate all internal and outlet flowrates for a given set of phase volumes.

        The level attributes and the outlet flowrate attributes are overwritten to match
        the given volumes.

        Parameters
        ----------
        V : array_like of float
            Volume of metal in the slag trough, slag in the slag trough, metal in the metal
            trough, and slag in the metal trough, m3.

        Returns
        -------
        numpy.ndarray
            Flowrates, m3/s: metal and slag through the weir, metal and slag out of the slag
            trough overflow, and metal and slag out of the metal trough overflow.

        """
        (
            self.slag_trough_hmetal,
            self.slag_trough_hslag,
            self.metal_trough_hmetal,
            self.metal_trough_hslag,
        ) = self._volumes_to_levels(*V)
        vdot_metal_12, vdot_slag_12 = self._calc_vdot_weir()
        self._calc_vdot_out(0.0)  # dt is not used by the outlet calculation
        return numpy.array(
            [
                vdot_metal_12,
                vdot_slag_12,
                self.slag_trough_vdotmetal_out,
                self.slag_trough_vdotslag_out,
                self.metal_trough_vdotmetal_out,
                self.metal_trough_vdotslag_out,
            ]
        )

    def _step(self, dt, vdotmetal_in, vdotslag_in):
        """Advance the state over a single linearly-implicit (Rosenbrock-Euler) step.

        With the phase volumes V, the flows q(V) and the mass balance dV/dt = F*q + inflow
        (where F is a constant matrix), the step is
        dV = dt*(I - dt*J)^-1 * dV/dt, with the Jacobian J = F*dq/dV estimated by forward
        differences. The flowrates over the step that exactly satisfy the volume balance
        are q + dq/dV*dV.

        Parameters
        ----------
        dt : float
            Length of time step, in s.
        vdotmetal_in : float
            Flowrate of metal into the slag trough, m3/s.
        vdotslag_in : float
            Flowrate of slag into the slag trough, m3/s.

        Returns
        -------
        numpy.ndarray
            Effective flowrates over the step, in the order given by _calc_flows, m3/s.

        """
        V = numpy.array(
            self._levels_to_volumes(
                self.slag_trough_hmetal,
                self.slag_trough_hslag,
                self.metal_trough_hmetal,
                self.metal_trough_hslag,
            )
        )
        flows = self._calc_flows(V)
        dflows_dV = numpy.empty((len(flows), len(V)))
        for j in range(len(V)):
            V_perturbed = V.copy()
            V_perturbed[j] += self._JACOBIAN_EPS
            dflows_dV[:, j] = (self._calc_flows(V_perturbed) - flows) / self._JACOBIAN_EPS
        dV_dt = self._FLOWS_TO_VOLUME_RATES @ flows + numpy.array(
            [vdotmetal_in, vdotslag_in, 0, 0]
        )
        jacobian = self._FLOWS_TO_VOLUME_RATES @ dflows_dV
        V_delta = dt * numpy.linalg.solve(numpy.eye(len(V)) - dt * jacobian, dV_dt)
        (
            self.slag_trough_hmetal,
            self.slag_trough_hslag,
            self.metal_trough_hmetal,
            self.metal_trough_hslag,
        ) = self._volumes_to_levels(*(V + V_delta))
        return flows + dflows_dV @ V_delta

    def calc_dt(self, dt, vdotmetal_in, vdotslag_in):
        """
        Integrate object state forward over a single time step, and update
        output flowrates.

        The time step is divided into equal sub-steps no longer than max_substep, each
        integrated with a linearly-implicit (Rosenbrock-Euler) scheme. The inflow rates are
        held constant over the time step. The outlet flowrates are set to the average of the
        flowrates over the sub-steps, so that the change in the volumes held in the troughs
        is exactly the volume added by the inflows less the volume leaving by the outlets.

        Parameters
        ----------
        dt : float
            Length of time step, in s.
        vdotmetal_in : float
            Flowrate of metal into the slag trough, m3/s.
        vdotslag_in : float
            Flowrate of slag into the slag trough, m3/s.

        Returns
        -------
        None.

        Note
        ----
        The scheme is first order in time. It has no stability limit for the weir flow
        coupling between the troughs, which limits explicit schemes to time steps that are
        proportional to the throughput (about 0.2 s at 0.01 m3/s for the troughs in the unit
        tests). The step size is limited by accuracy instead: for time steps up to about 2 s
        the error in the phase volumes is of the order of 1% of the trough volume or
        better, and it increases with the step size. Very large steps are unreliable, for
        example filling an empty launder in a single step of 10 s or more can produce
        negative volumes.

        """
        n_substeps = max(1, int(numpy.ceil(dt / self.max_substep)))
        flows = sum(
            self._step(dt / n_substeps, vdotmetal_in, vdotslag_in) for _ in range(n_substeps)
        ) / n_substeps
        (
            self.slag_trough_vdotmetal_out,
            self.slag_trough_vdotslag_out,
            self.metal_trough_vdotmetal_out,
            self.metal_trough_vdotslag_out,
        ) = flows[2:]
