.. pyger documentation master file, created by
   sphinx-quickstart on Tue May  3 13:57:42 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Pyger Documentation
===================

Overview
********

Pyger is a tool for calculating allowed hot dwell times and associated cooldown times given spacecraft thermal constraints. More specifically Pyger:

* Uses a Monte-Carlo approach to sample a realistic ensemble of simulation starting temperatures.

* Generates both typically achievable (50%) and best-case (90%) hot dwell times.

* Generates typically achievable cooldown times for individual MSIDs.

* Includes all xija-based thermal models and the pline mode for available hot time (not cooldown time).

* Hot constraints implemented as Python classes derived from a base ConstraintModel class.

* Cooldown constraints implemented as Python Numpy Arrays

* Outputs data to pickle files for later use.

Due to the large number of constraints and the potential complexity of their interaction, it is suggested that dwell time characteristics be calculated on an individual model by model basis and superimposed while post processing the data rather than including multiple models in the initial dwell time simulations. Although the latter is possible in the current Pyger framework, it limits one's ability to draw insights from the data.

Please see the Jupyter notebook for examples on how to use Pyger `in the repository hosted on Github here`_.

.. _in the repository hosted on Github here: https://github.com/sot/pyger/blob/1.0/doc/Pyger_Demo.ipynb

Contents:

.. toctree::
   :maxdepth: 4

   pyger


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

