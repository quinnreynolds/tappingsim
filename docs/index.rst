tappingsim
==========

A Python package for reduced-order modelling of metallurgical furnace tapping systems.

.. image:: https://zenodo.org/badge/256487822.svg
   :target: https://zenodo.org/doi/10.5281/zenodo.10129773

Overview
--------

**tappingsim** provides tools for simulating the tapping process on submerged-arc furnaces (SAFs)
and similar metallurgical reactors. The package uses reduced-order models to track the
time-evolution of metal and slag levels within the furnace hearth, taphole flow rates, and
downstream ladle filling.

The core models are:

- **Furnace hearth**: mass balances for metal and slag phases, driven by smelting power and
  tapped out through a configurable taphole model.
- **Taphole flow**: packed-bed pressure drop (Kozeny-Carman or Ergun) and pipe friction factor
  (Bellos, Cheng, or Serghides) correlations to compute volumetric outflow.
- **Launders**: pass-through and two-tank skimmer models routing flow from the furnace taphole
  to ladles.
- **Ladles**: cylindrical vessels that accumulate metal and slag, and overflow through an
  open-channel outlet.
- **Components** (``tappingsim.components``): the shared parts from which the models are built
  and configured: the burden pressure-drop and friction-factor models, the open-channel flow over
  overflow outlets with models for metal entrained in the slag overflow, and the time integration
  schemes: a linearly-implicit one for the stiff overflow flows of the ladle and launder models,
  and an explicit fourth-order Runge-Kutta one for the furnace.

Getting Started
---------------

Install from the repository root::

    pip install -e .

For development (includes pytest)::

    pip install -e ".[dev]"

See the :doc:`examples/index` section for worked examples using Jupyter notebooks.

.. toctree::
   :maxdepth: 2
   :caption: Contents

   api/index
   examples/index

Indices and Tables
------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
