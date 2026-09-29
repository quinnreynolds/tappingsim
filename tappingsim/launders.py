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
    """Tapping launder class. This version is a two-tank model of a launder feeding an
    underflow weir to separate slag and metal.

    Attributes
    ----------
    vdotmetal_out : float
        The output metal flowrate from the launder in the current state,
        m3/s.
    vdotslag_out : float
        The output slag flowrate from the launder in the current state,
        m3/s.

    """

    def __init__(
        self,
        # post-weir
        metal_section_length,
        metal_section_width,
        metal_section_shape,
        metal_outlet_height,
        # pre-weir
        slag_section_length,
        slag_section_width,
        slag_section_shape,
        slag_outlet_height,
        # weir
        weir_underflow_height,
        weir_underflow_width,
        # properties
        density_metal,
        density_slag,
    ):
        self.metal_section_length = metal_section_length
        self.metal_section_width = metal_section_width
        self.metal_section_shape = metal_section_shape  # 'box' or 'v'
        self.metal_outlet_height = metal_outlet_height
        self.slag_section_length = slag_section_length
        self.slag_section_width = slag_section_width
        self.slag_section_shape = slag_section_shape  # 'box' or 'v'
        self.slag_outlet_height = slag_outlet_height
        self.weir_underflow_height = weir_underflow_height
        self.weir_underflow_width = weir_underflow_width
        self.density_metal = density_metal
        self.density_slag = density_slag
        self.metal_section_vdotmetal_out = 0
        self.metal_section_vdotslag_out = 0
        self.slag_section_vdotmetal_out = 0
        self.slag_section_vdotslag_out = 0
        self.metal_section_hmetal = 0
        self.metal_section_hslag = 0
        self.slag_section_hmetal = 0
        self.slag_section_hslag = 0

    def _levels_to_volumes(self, hm_ss, hs_ss, hm_ms, hs_ms):
        if self.slag_section_shape == "box":
            dVdh = self.slag_section_length * self.slag_section_width
            Vm_ss = dVdh * hm_ss
            Vs_ss = dVdh * (hs_ss - hm_ss)
        elif self.slag_section_shape == "v":
            raise NotImplementedError("V profile launders not implemented yet.")
        else:
            raise NotImplementedError(f"Unknown slag launder profile '{self.slag_section_shape}'")

        if self.metal_section_shape == "box":
            dVdh = self.metal_section_length * self.metal_section_width
            Vm_ms = dVdh * hm_ms
            Vs_ms = dVdh * (hs_ms - hm_ms)
        elif self.slag_section_shape == "v":
            raise NotImplementedError("V profile launders not implemented yet.")
        else:
            raise NotImplementedError(
                f"Unknown metal launder profile '{self.metal_section_shape}'"
            )

        return Vm_ss, Vs_ss, Vm_ms, Vs_ms

    def _volumes_to_levels(self, Vm_ss, Vs_ss, Vm_ms, Vs_ms):
        if self.slag_section_shape == "box":
            dhdV = 1 / (self.slag_section_length * self.slag_section_width)
            hm_ss = dhdV * Vm_ss
            hs_ss = hm_ss + dhdV * Vs_ss
        elif self.slag_section_shape == "v":
            raise NotImplementedError("V profile launders not implemented yet.")
        else:
            raise NotImplementedError(f"Unknown slag launder profile '{self.slag_section_shape}'")

        if self.metal_section_shape == "box":
            dhdV = 1 / (self.metal_section_length * self.metal_section_width)
            hm_ms = dhdV * Vm_ms
            hs_ms = hm_ms + dhdV * Vs_ms
        elif self.slag_section_shape == "v":
            raise NotImplementedError("V profile launders not implemented yet.")
        else:
            raise NotImplementedError(
                f"Unknown metal launder profile '{self.metal_section_shape}'"
            )

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
        xarea_ss = self.slag_section_length * self.slag_section_width
        xarea_ms = self.metal_section_length * self.metal_section_width

        # slag section
        if self.slag_section_hslag < self.slag_outlet_height:
            self.slag_section_vdotmetal_out = 0
            self.slag_section_vdotslag_out = 0
        else:
            # TODO: implement overflow entrainment models similar to ladle models -
            #       vfrac = 0 is a simple no-entrainment model
            # vfrac = self.overflowmodel(self.depth - self.hmetal)
            vfrac = 0
            if self.slag_section_hmetal < self.slag_outlet_height:
                vslagout = (self.slag_section_hslag - self.slag_outlet_height) * xarea_ss
                self.slag_section_vdotmetal_out = vslagout * vfrac / dt
                self.slag_section_vdotslag_out = vslagout * (1 - vfrac) / dt
            else:
                vslagout = (self.slag_section_hslag - self.slag_section_hmetal) * xarea_ss
                self.slag_section_vdotmetal_out = (
                    (self.slag_section_hmetal - self.slag_outlet_height) * xarea_ss
                    + vslagout * vfrac
                ) / dt
                self.slag_section_vdotslag_out = vslagout * (1 - vfrac) / dt

        # metal section
        if self.metal_section_hslag < self.metal_outlet_height:
            self.metal_section_vdotmetal_out = 0
            self.metal_section_vdotslag_out = 0
        else:
            # TODO: implement overflow entrainment models similar to ladle models -
            #       vfrac = 0 is a simple no-entrainment model
            # vfrac = self.overflowmodel(self.depth - self.hmetal)
            vfrac = 0
            if self.metal_section_hmetal < self.metal_outlet_height:
                vslagout = (self.metal_section_hslag - self.metal_outlet_height) * xarea_ms
                self.metal_section_vdotmetal_out = vslagout * vfrac / dt
                self.metal_section_vdotslag_out = vslagout * (1 - vfrac) / dt
            else:
                vslagout = (self.metal_section_hslag - self.metal_section_hmetal) * xarea_ms
                self.metal_section_vdotmetal_out = (
                    (self.metal_section_hmetal - self.metal_outlet_height) * xarea_ms
                    + vslagout * vfrac
                ) / dt
                self.metal_section_vdotslag_out = vslagout * (1 - vfrac) / dt

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
        # FRACTIONAL TIME STEP - FILLING

        # Assumptions and conventions:
        #  - Positive flow direction is from slag section to metal section
        #  - No deformation of slag-metal interface at weir inflow

        A_w = self.weir_underflow_height * self.weir_underflow_width
        P_ss = (
            self.slag_section_hmetal * self.density_metal * g
            + (self.slag_section_hslag - self.slag_section_hmetal) * self.density_slag * g
        )
        P_ms = (
            self.metal_section_hmetal * self.density_metal * g
            + (self.metal_section_hmetal - self.metal_section_hmetal) * self.density_slag * g
        )
        deltaP = P_ss - P_ms
        u_mw = numpy.sign(deltaP) * numpy.sqrt(2 * numpy.abs(deltaP) / self.density_metal)
        u_sw = numpy.sign(deltaP) * numpy.sqrt(2 * numpy.abs(deltaP) / self.density_slag)

        # slag-metal interface in slag section is above weir underflow channel
        if self.slag_section_hmetal > self.weir_underflow_height:
            vdot_metal_12 = A_w * u_mw
            vdot_slag_12 = 0
        # slag-metal interface in slag section is below weir underflow channel
        elif self.slag_section_hslag > self.weir_underflow_height:
            vdot_metal_12 = self.slag_section_hmetal / self.weir_underflow_height * A_w * u_mw
            vdot_slag_12 = (1 - self.slag_section_hmetal / self.weir_underflow_height) * A_w * u_sw
        # free surface in slag section is below weir underflow channel
        else:
            vdot_metal_12 = self.slag_section_hmetal / self.weir_underflow_height * A_w * u_mw
            vdot_slag_12 = (
                ((self.slag_section_hslag - self.slag_section_hmetal) / self.weir_underflow_height)
                * A_w
                * u_sw
            )

        V_deltas = (
            numpy.array(
                [
                    vdotmetal_in - vdot_metal_12,
                    vdotslag_in - vdot_slag_12,
                    vdot_metal_12,
                    vdot_slag_12,
                ]
            )
            * dt
        )
        V_0 = numpy.array(
            self._levels_to_volumes(
                self.slag_section_hmetal,
                self.slag_section_hslag,
                self.metal_section_hmetal,
                self.metal_section_hslag,
            )
        )
        V_new = V_0 + V_deltas
        (
            self.slag_section_hmetal,
            self.slag_section_hslag,
            self.metal_section_hmetal,
            self.metal_section_hslag,
        ) = self._volumes_to_levels(*V_new)

        # FRACTIONAL TIME STEP - EMPTYING
        self._calc_vdot_out(dt)
        V_deltas = (
            numpy.array(
                [
                    self.slag_section_vdotmetal_out,
                    self.slag_section_vdotslag_out,
                    self.metal_section_vdotmetal_out,
                    self.metal_section_vdotslag_out,
                ]
            )
            * dt
        )
        V_0 = numpy.array(
            self._levels_to_volumes(
                self.slag_section_hmetal,
                self.slag_section_hslag,
                self.metal_section_hmetal,
                self.metal_section_hslag,
            )
        )
        V_new = V_0 + V_deltas
        (
            self.slag_section_hmetal,
            self.slag_section_hslag,
            self.metal_section_hmetal,
            self.metal_section_hslag,
        ) = self._volumes_to_levels(*V_new)
        
        
