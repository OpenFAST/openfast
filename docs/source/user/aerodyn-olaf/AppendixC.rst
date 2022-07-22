.. _OLAF-List-of-Output-Channels:

Appendix C: OLAF List of Output Channels
========================================

This is a list of all possible output parameters from the OLAF module.
The names are grouped by meaning, but can be ordered in the OUTPUTS
section of the *AeroDyn15* primary input file, as the user sees fit.
:math:`N\beta` refers to output node, :math:`\beta`, where :math:`\beta`
is a number in the range [1,9], corresponding to entry, :math:`\beta`,
in the **OutNd** list. :math:`B\alpha` is prefixed to each output name,
where :math:`\alpha` is a number in the range [1,3], corresponding to
the blade number.


.. list-table:: Available OLAF Output Channels
   :widths: 25 15 50
   :header-rows: 1
   :align: center
   :name: Tab:OLAFoutputs

   *  - Channel Name(s)
      - Units
      - Description
   *  - :math:`B \alpha N \beta Gam`
      - :math:`m^2/s`
      - Circulation along the blade


..
      ============================ ============= ===========================
      Channel Name(s)              Units         Description
      ============================ ============= ===========================
      :math:`B \alpha N \beta Gam` :math:`m^2/s` Circulation along the blade
      ============================ ============= ===========================
