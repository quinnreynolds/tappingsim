"""Classes to create tapping system simulator objects.

It is often of interest to consider some or all components of a furnace tapping
system as a single integrated simulation, in order to investigate interactions
between different components and the effect of over-arching control and
operation philosophies.

"""

import math
from numbers import Integral


class SAF():
    """Wrapper class for creating a submerged-arc furnace tapping system
    simulator object. This system contains only the furnace.

    Parameters
    ----------
    furnace : object
        An instance of the SubmergedArcFurnace class.

    Attributes
    ----------
    furnace : object
        A reference to the provided furnace object.
    timetotaliser : float
        A running counter for time, in s.
    powertotaliserkWh : float
        A running counter for energy input to the furnace, in kWh.
    metalmasstotaliser : float
        A running counter of metal tapped from the furnace, in kg.
    slagmasstotaliser : float
        A running counter of slag tapped from the furnace, in kg.

    """

    def __init__(self, furnace):
        self.furnace = furnace
        self.timetotaliser = 0
        self.powertotaliserkWh = 0
        self.metalmasstotaliser = 0
        self.slagmasstotaliser = 0

    def open_taphole(self):
        """
        Set furnace tap-hole to open.

        """
        self.furnace.tapholeopen_yn = True

    def close_taphole(self):
        """
        Set furnace tap-hole to closed.

        """
        self.furnace.tapholeopen_yn = False

    def reset_time_totaliser(self):
        """
        Reset time totaliser counter to zero.

        """
        self.timetotaliser = 0

    def reset_power_totaliser(self):
        """
        Reset power totaliser counter to zero.

        """
        self.powertotaliserkWh = 0

    def reset_mass_totaliser(self):
        """
        Reset tap mass totaliser counters to zero.

        """
        self.metalmasstotaliser = 0
        self.slagmasstotaliser = 0

    def reset_all_totalisers(self):
        """
        Reset all totaliser counters to zero.

        """
        self.reset_time_totaliser()
        self.reset_power_totaliser()
        self.reset_mass_totaliser()

    def calc_dt(self, dt):
        """
        Integrate simulator state forward over a single time step.

        Parameters
        ----------
        dt : float
            Length of time step, in s.

        Returns
        -------
        None.

        """
        self.furnace.calc_dt(dt)
        self.timetotaliser += dt
        self.powertotaliserkWh += dt * (self.furnace.powerMVA
                                        * self.furnace.powerfactor / 3.6)
        self.metalmasstotaliser += dt * (self.furnace.vdotmetal_out
                                         * self.furnace.densitymetal)
        self.slagmasstotaliser += dt * (self.furnace.vdotslag_out
                                        * self.furnace.densityslag)


class SAFWithLadles():
    """Wrapper class for creating a submerged-arc furnace tapping system
    simulator object. This system contains the furnace, a transfer launder,
    and an arbitrary number of ladles in series.

    Parameters
    ----------
    furnace : object
        An instance of the SubmergedArcFurnace class or equivalent.
    launder : object
        An instance of the SimpleSisoLaunder class or equivalent. The launder
        object must implement the vdotmetal_out and vdotslag_out attributes,
        the metal and slag flowrates out of the launder into the first ladle,
        as these are read by the calc_dt method.
    ladles: list of object
        A list of instances of the CylindricalLadle class or equivalent. The
        ladle objects must have a diameter attribute, the internal diameter of
        the ladle, as this is read by the ladle_masses method, so ladles of
        other shapes, such as the PanLadle class, cannot be used here.

    Attributes
    ----------
    furnace : object
        A reference to the provided furnace object.
    launder : object
        A reference to the provided launder object.
    ladles: list of object
        A reference to the provided list of ladle objects.
    timetotaliser : float
        A running counter for time, in s.
    powertotaliserkWh : float
        A running counter for energy input to the furnace, in kWh.
    metalmasstotaliser : float
        A running counter of metal tapped from the furnace, in kg.
    slagmasstotaliser : float
        A running counter of slag tapped from the furnace, in kg.

    """

    def __init__(self, furnace, launder, ladles):
        self.furnace = furnace
        self.launder = launder
        self.ladles = ladles
        self.timetotaliser = 0
        self.powertotaliserkWh = 0
        self.metalmasstotaliser = 0
        self.slagmasstotaliser = 0

    def open_taphole(self):
        """
        Set furnace tap-hole to open.

        """
        self.furnace.tapholeopen_yn = True

    def close_taphole(self):
        """
        Set furnace tap-hole to closed.

        """
        self.furnace.tapholeopen_yn = False

    def empty_ladles(self):
        """
        Set slag and metal levels in all ladles to zero.

        """
        for ldl in self.ladles:
            ldl.hmetal, ldl.hslag = 0, 0

    def reset_time_totaliser(self):
        """
        Reset time totaliser counter to zero.

        """
        self.timetotaliser = 0

    def reset_power_totaliser(self):
        """
        Reset power totaliser counter to zero.

        """
        self.powertotaliserkWh = 0

    def reset_mass_totaliser(self):
        """
        Reset tap mass totaliser counters to zero.

        """
        self.metalmasstotaliser = 0
        self.slagmasstotaliser = 0

    def reset_all_totalisers(self):
        """
        Reset all totaliser counters to zero.

        """
        self.reset_time_totaliser()
        self.reset_power_totaliser()
        self.reset_mass_totaliser()

    def ladle_masses(self):
        """
        Calculate masses of metal and slag in each ladle.

        Parameters
        ----------
        None.

        Returns
        -------
        metalmasses : list of float
            Mass of metal in each ladle, in kg.
        slagmasses : list of float
            Mass of slag in each ladle, in kg.

        """
        metalmasses, slagmasses = [], []
        for ldl in self.ladles:
            area = math.pi * 0.25 * ldl.diameter**2
            metalmasses.append(area * ldl.hmetal * self.furnace.densitymetal)
            slagmasses.append(area * (ldl.hslag - ldl.hmetal)
                              * self.furnace.densityslag)
        return metalmasses, slagmasses

    def calc_dt(self, dt):
        """
        Integrate simulator state forward over a single time step.

        Parameters
        ----------
        dt : float
            Length of time step, in s.

        Returns
        -------
        None.

        """
        self.furnace.calc_dt(dt)
        if self.furnace.tapholeopen_yn:
            self.launder.calc_dt(dt, self.furnace.vdotmetal_out,
                                 self.furnace.vdotslag_out)
            self.ladles[0].calc_dt(dt, self.launder.vdotmetal_out,
                                   self.launder.vdotslag_out)
            for lprev, lnext in zip(self.ladles[:-1], self.ladles[1:]):
                lnext.calc_dt(dt, lprev.vdotmetal_out, lprev.vdotslag_out)

        self.timetotaliser += dt
        self.powertotaliserkWh += dt * (self.furnace.powerMVA
                                        * self.furnace.powerfactor) / 3.6
        self.metalmasstotaliser += dt * (self.furnace.vdotmetal_out
                                         * self.furnace.densitymetal)
        self.slagmasstotaliser += dt * (self.furnace.vdotslag_out
                                        * self.furnace.densityslag)


def _check_hole_index(furnace, index):
    """Check that an index refers to one of the tap-holes of a furnace.

    Raises
    ------
    IndexError
        If the index is not an integer from 0 to the number of tap-holes less
        one.

    """
    if (isinstance(index, bool) or not isinstance(index, Integral)
            or not 0 <= index < furnace.ntapholes):
        raise IndexError(
            f"tap-hole index must be an integer from 0 to "
            f"{furnace.ntapholes - 1}, not {index!r}")


class DCAF():
    """Wrapper class for creating a direct-current arc furnace tapping system
    simulator object. This system contains only the furnace.

    The tap-holes are opened and closed independently, by index. Each tap-hole
    can pass both metal and slag, depending on the levels in the furnace, so
    the masses tapped through each tap-hole are totalised separately for both
    phases.

    Parameters
    ----------
    furnace : object
        An instance of the DCArcFurnace class.

    Attributes
    ----------
    furnace : object
        A reference to the provided furnace object.
    timetotaliser : float
        A running counter for time, in s.
    powertotaliserkWh : float
        A running counter for energy input to the furnace, in kWh.
    metalmasstotalisers : list of float
        A running counter of metal tapped through each tap-hole, in kg.
    slagmasstotalisers : list of float
        A running counter of slag tapped through each tap-hole, in kg.

    Note
    ----
    The mass totalisers use the flowrates averaged over each time step that
    the furnace reports, so they are consistent with the change in the furnace
    inventory.

    """

    def __init__(self, furnace):
        self.furnace = furnace
        self.timetotaliser = 0
        self.powertotaliserkWh = 0
        self._zero_mass_totalisers()

    def _zero_mass_totalisers(self):
        n = self.furnace.ntapholes
        self.metalmasstotalisers = [0 for _ in range(n)]
        self.slagmasstotalisers = [0 for _ in range(n)]

    @property
    def metalmasstotaliser(self):
        """A running counter of metal tapped through all of the tap-holes, in
        kg."""
        return sum(self.metalmasstotalisers)

    @property
    def slagmasstotaliser(self):
        """A running counter of slag tapped through all of the tap-holes, in
        kg."""
        return sum(self.slagmasstotalisers)

    def open_taphole(self, index):
        """
        Set a furnace tap-hole to open.

        Parameters
        ----------
        index : int
            The position of the tap-hole in the furnace tap-hole lists.

        """
        _check_hole_index(self.furnace, index)
        self.furnace.tapholesopen_yn[index] = True

    def close_taphole(self, index):
        """
        Set a furnace tap-hole to closed.

        Parameters
        ----------
        index : int
            The position of the tap-hole in the furnace tap-hole lists.

        """
        _check_hole_index(self.furnace, index)
        self.furnace.tapholesopen_yn[index] = False

    def reset_time_totaliser(self):
        """
        Reset time totaliser counter to zero.

        """
        self.timetotaliser = 0

    def reset_power_totaliser(self):
        """
        Reset power totaliser counter to zero.

        """
        self.powertotaliserkWh = 0

    def reset_mass_totaliser(self):
        """
        Reset tap mass totaliser counters to zero.

        """
        self._zero_mass_totalisers()

    def reset_all_totalisers(self):
        """
        Reset all totaliser counters to zero.

        """
        self.reset_time_totaliser()
        self.reset_power_totaliser()
        self.reset_mass_totaliser()

    def _update_totalisers(self, dt):
        furnace = self.furnace
        self.timetotaliser += dt
        self.powertotaliserkWh += dt * furnace.powerMW / 3.6
        for i in range(furnace.ntapholes):
            self.metalmasstotalisers[i] += dt * (furnace.vdotmetal_outs[i]
                                                 * furnace.densitymetal)
            self.slagmasstotalisers[i] += dt * (furnace.vdotslag_outs[i]
                                                * furnace.densityslag)

    def calc_dt(self, dt):
        """
        Integrate simulator state forward over a single time step.

        Parameters
        ----------
        dt : float
            Length of time step, in s.

        Returns
        -------
        None.

        """
        self.furnace.calc_dt(dt)
        self._update_totalisers(dt)


class DCAFWithLadles(DCAF):
    """Wrapper class for creating a direct-current arc furnace tapping system
    simulator object. This system contains the furnace and, for each
    tap-hole, an optional transfer launder and an arbitrary number of ladles
    in series.

    Parameters
    ----------
    furnace : object
        An instance of the DCArcFurnace class.
    launders : list of object
        The launder receiving the flow from each tap-hole, in the order of the
        furnace tap-holes. Each is an instance of the SimpleSisoLaunder class
        or equivalent, which must implement the vdotmetal_out and
        vdotslag_out attributes, the metal and slag flowrates out of the
        launder into the first ladle, as these are read by the calc_dt
        method. An entry of None means that the flow from the tap-hole is
        only totalised, and the tap-hole must then have no ladles.
    ladles : list of list of object
        The ladles receiving the flow from each tap-hole, in the order of the
        furnace tap-holes. Each entry is a list of instances of the
        CylindricalLadle class or equivalent, filled in series from the
        launder, or an empty list for a tap-hole with no ladles. A launder or
        a ladle can belong to only one tap-hole. The ladle objects must have a
        diameter attribute, the internal diameter of the ladle, and an outlet
        attribute with the phase densities, as these are read by the
        ladle_masses method, so ladles of other shapes, such as the PanLadle
        class, cannot be used here.

    Attributes
    ----------
    furnace : object
        A reference to the provided furnace object.
    launders : list of object
        A reference to the provided list of launders.
    ladles : list of list of object
        A reference to the provided list of lists of ladles.
    timetotaliser : float
        A running counter for time, in s.
    powertotaliserkWh : float
        A running counter for energy input to the furnace, in kWh.
    metalmasstotalisers : list of float
        A running counter of metal tapped through each tap-hole, in kg.
    slagmasstotalisers : list of float
        A running counter of slag tapped through each tap-hole, in kg.

    Note
    ----
    As for the SAFWithLadles system, the launder and ladles for a tap-hole are
    advanced only during time steps in which that tap-hole is open.

    """

    def __init__(self, furnace, launders, ladles):
        n = furnace.ntapholes
        launders, ladles = list(launders), list(ladles)
        if len(launders) != n or len(ladles) != n:
            raise ValueError(
                f"launders and ladles must each have an entry for each of the "
                f"{n} furnace tap-holes, but have {len(launders)} and "
                f"{len(ladles)}")
        ladles = [[] if chain is None else list(chain) for chain in ladles]
        for i, (launder, chain) in enumerate(zip(launders, ladles)):
            if launder is None and chain:
                raise ValueError(
                    f"tap-hole {i} has ladles but no launder to feed them")
        units = [u for u in launders if u is not None]
        units += [ldl for chain in ladles for ldl in chain]
        if len({id(u) for u in units}) != len(units):
            raise ValueError(
                "each launder and ladle can belong to only one tap-hole")
        super().__init__(furnace)
        self.launders = launders
        self.ladles = ladles

    def empty_ladles(self):
        """
        Set slag and metal levels in all ladles to zero.

        """
        for chain in self.ladles:
            for ldl in chain:
                ldl.hmetal, ldl.hslag = 0, 0

    def ladle_masses(self):
        """
        Calculate masses of metal and slag in each ladle, using the densities
        of the ladle.

        Parameters
        ----------
        None.

        Returns
        -------
        metalmasses : list of list of float
            Mass of metal in each ladle for each tap-hole, in kg.
        slagmasses : list of list of float
            Mass of slag in each ladle for each tap-hole, in kg.

        """
        metalmasses, slagmasses = [], []
        for chain in self.ladles:
            metalmasses.append([])
            slagmasses.append([])
            for ldl in chain:
                area = math.pi * 0.25 * ldl.diameter**2
                metalmasses[-1].append(area * ldl.hmetal
                                       * ldl.outlet.density_metal)
                slagmasses[-1].append(area * (ldl.hslag - ldl.hmetal)
                                      * ldl.outlet.density_slag)
        return metalmasses, slagmasses

    def calc_dt(self, dt):
        """
        Integrate simulator state forward over a single time step.

        Parameters
        ----------
        dt : float
            Length of time step, in s.

        Returns
        -------
        None.

        """
        furnace = self.furnace
        furnace.calc_dt(dt)
        for i, (launder, chain) in enumerate(zip(self.launders, self.ladles)):
            if launder is None or not furnace.tapholesopen_yn[i]:
                continue
            launder.calc_dt(dt, furnace.vdotmetal_outs[i],
                            furnace.vdotslag_outs[i])
            vdotmetal, vdotslag = launder.vdotmetal_out, launder.vdotslag_out
            for ldl in chain:
                ldl.calc_dt(dt, vdotmetal, vdotslag)
                vdotmetal, vdotslag = ldl.vdotmetal_out, ldl.vdotslag_out
        self._update_totalisers(dt)
