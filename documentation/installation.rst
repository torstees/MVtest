Installation
============
MVTest requires python 2.7 or later (but 3.0 and later) as well as the
following packages:

* NumPy (version 1.7.2 or later)   www.numpy.org
* SciPY (version 0.13.2 or later)  www.scipy.org

If these aren't already installed, and you don't have root access to the
machine, please see the section, :ref:`Miniconda` or :ref:`virtual-env`
for easy instructions on different ways of installing tools as a
restricted user. MVTest's installer can install these for you, however,
it assumes that you have write access to your python library, which
will not be the case by default on shared systems.

Download the package at: TODO: URL

To install the software, run the setup script as shown below:
$ `python setup.py install`

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

Running Unit Tests
++++++++++++++++++
MVTest comes with a unit test suite which can be run prior to installation.
To run the tests, simply run the following command from within the root
directory of the extracted archive's contents:

$ `python setup.py test`

If no errors are reported, then mvtest should run correctly on your system.

.. _virtual-env:

Virtual Env
+++++++++++
Virtual ENV is a powerful too for python programers and users alike, as it
allows for users to deploy different versions of python applications without
the need for root access to the machine.

Because MVTest requires version 2.7, you'll need to ensure that your machine's
python version is in compliance. Virtual Env basically uses the the system
version of python, but creates a user owned environment wrapper allowing
users to install libraries easily without administrative rights to the
machine.

For a helpful introduction to VirtualEnv, please have a look at the
tutorial: http://www.simononsoftware.com/virtualenv-tutorial/

.. _miniconda:

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
$ `conda install pip scipy`

Installing SciPy will also force the installation of NumPy, which is
also required for running mvtest. (setuptools includes easy_install).

Once that has been completed successfully, you should be ready to follow
the standard instructions for installing mvtest.