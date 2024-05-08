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
   :module: baredSC.baredSC_1d
   :func: parse_arguments
   :prog: baredSC_1d

combineMultipleModels_1d
------------------------

.. argparse::
   :module: baredSC.combineMultipleModels_1d
   :func: parse_arguments
   :prog: combineMultipleModels_1d

baredSC_2d
-----------

.. argparse::
   :module: baredSC.baredSC_2d
   :func: parse_arguments
   :prog: baredSC_2d

combineMultipleModels_2d
------------------------

.. argparse::
   :module: baredSC.combineMultipleModels_2d
   :func: parse_arguments
   :prog: combineMultipleModels_2d
