temperature_model
=================

Data-driven model that maps pyrolysis temperature (and optionally feedstock)
to biochar composition targets (H/C ratio, O/C ratio, aromaticity) and a set
of reference properties.  Fit from the **UC Davis Biochar Database**
(https://biochar.ucdavis.edu/) with aromaticity calibrated from NMR data
collected by Wood, Mašek & Erastova (2024).

Module-level helpers
---------------------

.. autofunction:: biochar.temperature_model.properties

.. autodata:: biochar.temperature_model.VALID_FEEDSTOCKS

``TemperatureModel``
---------------------

.. autoclass:: biochar.temperature_model.TemperatureModel
   :members: composition, predict, predict_all, valid_feedstocks, provenance
   :show-inheritance:

Build helper
------------

.. autofunction:: biochar.temperature_model.build_model

.. autofunction:: biochar.temperature_model.compare_models
