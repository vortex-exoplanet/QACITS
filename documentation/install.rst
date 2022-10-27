.. _install:

Installing and Configuring
####################################

Obtaining QACITS
==========================================
``QACITS`` is hosted on github at https://github.com/vortex-exoplanet/QACITS.  
The master branch is a development branch and may be frequently changed. For 
publishable/reproducible results, we strongly recommend using the latest tagged 
release from:

https://github.com/vortex-exoplanet/QACITS/releases

Installing QACITS
==========================================

You can download ``QACITS`` from its GitHub repository as a zip file. A ``setup.py`` file (setuptools) is included in the root folder of ``QACITS``.
Enter the package's root folder and run:

.. code-block:: bash

  $ python setup.py install

Python Package Dependencies
==========================================

``QACITS`` requires Python 3.6+ and the following packages:

* photutils
* astropy
* numpy
* scipy
* matplotlib (for visualization of results only)

Run-time Environment and Deployment 
==========================================

While in closed-loop observation mode, QACITS should run on a real-time clock workstation, featuring high core CPU count and sufficient main memory. 