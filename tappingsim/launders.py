"""Classes to create transfer launder objects.

Launders are typically open-channel gravity flow conduits which redirect the
flow of tapped material from the furnace into storage ladles or other
post-treatment units.

"""

import numpy
from scipy.constants import g


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

    Attributes
    ----------
    TODO

    """

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
    ):
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
        self.outlet_premult = 2 / 3 * numpy.sqrt(g)

    def _levels_to_volumes(self, hm_ss, hs_ss, hm_ms, hs_ms):
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
        # slag section
        H, L = self.slag_trough_height, self.slag_trough_length
        W0 = self.slag_trough_lower_width
        dW = self.slag_trough_upper_width - self.slag_trough_lower_width
        if numpy.isclose(dW, 0):
            hm_ss = Vm_ss / (L * W0)
            hs_ss = hm_ss + Vs_ss / (L * W0)
        else:
            hm_ss = (numpy.sqrt(H * (H * L * W0**2 + 2 * Vm_ss * dW) / L) - H * W0) / dW
            hs_ss = (
                numpy.sqrt(H * (H * L * W0**2 + 2 * Vm_ss * dW + 2 * Vs_ss * dW) / L) - H * W0
            ) / dW
        # metal section
        H, L = self.metal_trough_height, self.metal_trough_length
        W0 = self.metal_trough_lower_width
        dW = self.metal_trough_upper_width - self.metal_trough_lower_width
        if numpy.isclose(dW, 0):
            hm_ms = Vm_ms / (L * W0)
            hs_ms = hm_ms + Vs_ms / (L * W0)
        else:
            hm_ms = (numpy.sqrt(H * (H * L * W0**2 + 2 * Vm_ms * dW) / L) - H * W0) / dW
            hs_ms = (
                numpy.sqrt(H * (H * L * W0**2 + 2 * Vm_ms * dW + 2 * Vs_ms * dW) / L) - H * W0
            ) / dW
        return hm_ss, hs_ss, hm_ms, hs_ms

    def _calc_vdot_out(self, dt):
        """Calculate volumetric outflow rates of metal and slag from the slag and metal
        section outlets for the current state, and update the vdotmetal_out and
        vdotslag_out attributes.

        Overflow occurs when the total bath level exceeds the ladle depth.
        TODO: Metal entrainment in the slag outflow is calculated using the
        configured overflow model.

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
        # Assumptions and conventions:
        #  - Time step dt is "small", explicit forward integration is acceptable
        #  - Positive flow direction is from slag section to metal section
        #  - No deformation of slag-metal interface at weir inflow
        A_w = self.weir_underflow_height * self.weir_underflow_width
        P_ss = (
            self.slag_trough_hmetal * self.density_metal * g
            + (self.slag_trough_hslag - self.slag_trough_hmetal) * self.density_slag * g
        )
        P_ms = (
            self.metal_trough_hmetal * self.density_metal * g
            + (self.metal_trough_hmetal - self.metal_trough_hmetal) * self.density_slag * g
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
        # update trough outlet flowrates
        self._calc_vdot_out(dt)
        # calculate mass balance and new levels
        V_deltas = (
            numpy.array(
                [
                    vdotmetal_in - vdot_metal_12 - self.slag_trough_vdotmetal_out,
                    vdotslag_in - vdot_slag_12 - self.slag_trough_vdotslag_out,
                    vdot_metal_12 - self.metal_trough_vdotmetal_out,
                    vdot_slag_12 - self.metal_trough_vdotslag_out,
                ]
            )
            * dt
        )
        V_0 = numpy.array(
            self._levels_to_volumes(
                self.slag_trough_hmetal,
                self.slag_trough_hslag,
                self.metal_trough_hmetal,
                self.metal_trough_hslag,
            )
        )
        V_new = V_0 + V_deltas
        (
            self.slag_trough_hmetal,
            self.slag_trough_hslag,
            self.metal_trough_hmetal,
            self.metal_trough_hslag,
        ) = self._volumes_to_levels(*V_new)
