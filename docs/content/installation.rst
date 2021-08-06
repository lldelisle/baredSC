Installation
============

.. contents:: 
    :local:
    
Requirements
------------

It as only been tested on linux but should work on MacOS.

It requires python >= 3.7 (tested 3.7.3 and 3.9.1)

Dependencies of classical python packages:

* numpy (tested 1.16.4 and 1.19.5)
* matplotlib (tested 3.1.1 and 3.3.4)
* pandas (tested 0.25.0 and 1.2.1)
* scipy (tested 1.3.0 and 1.6.0)

Dependencies of a python package from Jean-Baptiste Delisle dedicated to mcmc:

* `samsam <https://obswww.unige.ch/~delisle/samsam/doc/>`_ (above 0.1.2)

Installation
------------

You can install it with pip:

.. code:: bash

    $ pip install baredSC


Or with conda:

.. code:: bash

    $ conda create -n baredSC -c bioconda -c conda-forge baredsc

