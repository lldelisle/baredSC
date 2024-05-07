Usage
=====

.. contents:: 
    :local:

General usage
-------------

Here is a description of all possible parameters.
The tool take around 1 minute for 1d, 2000 cells, default parameters, independently of the number of gaussian and around 30 seconds for 300 cells.
For the 2d 300 cells is around 3 minutes, 2000 cells around 15 minutes.
Increasing the number of samples in the MCMC or in the burning phase will increase the time.

baredSC_1d
-----------

.. argparse::
   :ref: baredSC.common.parse_arguments('baredSC_1d')
   :prog: baredSC_1d

combineMultipleModels_1d
------------------------

.. argparse::
   :ref: baredSC.common.parse_arguments('combineMultipleModels_1d')
   :prog: combineMultipleModels_1d

baredSC_2d
-----------

.. argparse::
   :ref: baredSC.common.parse_arguments('baredSC_2d')
   :prog: baredSC_2d

combineMultipleModels_2d
------------------------

.. argparse::
   :ref: baredSC.common.parse_arguments('combineMultipleModels_2d')
   :prog: combineMultipleModels_2d
