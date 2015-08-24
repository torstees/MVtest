.. MV-Test documentation master file, created by
   sphinx-quickstart on Fri Aug 21 14:11:20 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

MVTest -- GWAS Analysis
===================================

Contents:

.. toctree::
   :maxdepth: 2

   pygwas/modules


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

What is MVTest?
===============

*TODO: Write some background information about the application and it's
scientific basis.*


Installation
++++++++++++
MVTest requires python 2.7 or later (but 3.0 and later) as well as the
following packages:

* NumPy (version 1.7.2 or later)   www.numpy.org
* SciPY (version 0.13.2 or later)  www.scipy.org

If these aren't already installed, and you don't have root access to the
machine, please see the section, Miniconda for easy instructions on
installing a compact version of python suitable for running MVTest. MVTest's
installer can install these for you, however, it assumes that you have
write access to your python library, which will not be the case by default
on shared systems.

Download the package at: TODO: URL

To install the software, run the setup script as shown below:
> `python setup.py install`

If no errors are reported, it should be installed and ready to use.

System Requirements
+++++++++++++++++++
Aside from the library dependencies, MVTest's requirements depend largely on
the number of SNPs and individuals being analyzed as well as the data format
being used. In general, GWAS sized datasets will require several gigabytes of
memory when using the traditional pedigree format, however, even 10s of
thousands of subjects can be analyzed with less than 1 gigabyte of RAM when
the data is formatted as transposed pedigree or PLINK's default bed format.

Otherwise, it is recommended that the system be run on a unix-like system
such as Linux or OS X, but it should work under windows as well (this is
not a fully supported platform).

## Running Unit Tests
MVTest comes with a unit test suite which can be run prior to installation.
To run the tests, simply run the following command from within the root
directory of the extracted archive's contents:
python setup.py test

If no errors are reported, then mvtest should run correctly on your system.

Miniconda
+++++++++
Miniconda is a minimal version of the package manager used by the Anaconda
python distribution. It makes it easy to create local installations of python
with the latest versions of the common scientific libraries for users who don't
have root access to their target machines.

Firstly, download and install the appropriate version of miniconda at the
project website. Please be sure to choose the Python 2 version:
http://conda.pydata.org/miniconda.html

While it is doing the installation, please allow it to update your PATH
information. Also, be sure to follow directions such as starting a new
shell to allow those changes to take effect.

Once those changes have taken effect, install setuptools and scipy:
conda install setuptools scipy

Installing SciPy will also force the installation of NumPy, whish is
also required for running mvtest. (setuptools includes easy_install).

Once that has been completed successfully, you should be ready to follow
the standard instructions for installing mvtest.

Documentation
+++++++++++++
Documentation for mvtest is still under construction. However, the application
provides reasonable inline help using standard unix help arguments:

> `mvtest.py -h`

or

> `mvtest.py --help`

In general, overlapping functionality should mimic that of PLINK.

Command-Line Arguments
++++++++++++++++++++++
Command line arguments used by mvtest often mimick those used by PLINK, except
where there is no matching functionality (or the functionality differs
significantly.)

For the parameters listed below, when a parameter requires a value, the value
must follow the argument with a single space separating the two (no '=' signs.)
For flags with no specified value, passing the flag indicates that condition
is to be "activated".