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
- **Launders**: simple pass-through models routing flow from the furnace taphole to ladles.
- **Ladles**: cylindrical vessels that accumulate metal and slag, with overflow models.

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
