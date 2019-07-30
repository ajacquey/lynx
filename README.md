<h1 align="center">
  <br>
  <a href="https://gitext.gfz-potsdam.de/ajacquey/lynx">LYNX</a>
  <br>
  Lithosphere dYnamic Numerical toolboX
  <br>
  A MOOSE-based application
  <br>
</h1>

<h4 align="center">A numerical simulator for modelling deformation of the lithosphere, based on <a href="http://mooseframework.org/" target="blank">MOOSE</a>.</h4>

<p align="center">
  <a href="LICENSE">
    <img src="https://img.shields.io/badge/license-GPLv3-blue.svg"
         alt="GPL License">
  </a>
</p>

## About
LYNX (Lithosphere dYnamic Numerical toolboX) is a numerical simulator for modelling coupled Thermo-Hydro-Mechanical processes of porous rocks.
The simulator is developed by [Antoine Jacquey](http://www.gfz-potsdam.de/en/staff/antoine-jacquey/) <a href="https://orcid.org/0000-0002-6259-4305" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a> and [Mauro Cacace](http://www.gfz-potsdam.de/en/section/basin-modeling/staff/profil/mauro-cacace/) <a href="https://orcid.org/0000-0001-6101-9918" target="orcid.widget" rel="noopener noreferrer" style="vertical-align:top;"><img src="https://orcid.org/sites/default/files/images/orcid_16x16.png" style="width:1em;margin-right:.5em;" alt="ORCID iD icon"></a> at the [GFZ Potsdam, German Research Centre for Geosciences](http://www.gfz-potsdam.de/en/home/) from the section [Basin Modelling](http://www.gfz-potsdam.de/en/section/basin-modeling/).


LYNX is a MOOSE-based application. Visit the [MOOSE framework](http://mooseframework.org) page for more information.

## Licence
LYNX is distributed under the [GNU GENERAL PUBLIC LICENSE v3](https://gitext.gfz-potsdam.de/ajacquey/lynx/blob/master/LICENSE).


## Getting Started

#### Minimum System Requirements
The following system requirements are from the MOOSE framework (see [Getting Started](http://mooseframework.org/getting-started/) for more information):
* Compiler: C++11 Compliant GCC 4.8.4, Clang 3.4.0, Intel20130607
* Python 2.7+
* Memory: 16 GBs (debug builds)
* Processor: 64-bit x86
* Disk: 30 GBs
* OS: UNIX compatible (OS X, most flavors of Linux)

#### 1. Setting Up a MOOSE Installation
To install LYNX, you need first to have a working and up-to-date installation of the MOOSE framework.  
To do so, please visit the [Getting Started](http://mooseframework.org/getting-started/) page of the MOOSE framework and follow the instructions. If you encounter difficulties at this step, you can ask for help on the [MOOSE-users Google group](https://groups.google.com/forum/#!forum/moose-users).

#### 2. Clone LYNX
LYNX can be cloned directly from [GitLab](https://gitext.gfz-potsdam.de/ajacquey/lynx) using [Git](https://git-scm.com/). In the following, we refer to the directory `projects` which you created during the MOOSE installation (by default `~/projects`):  

    cd ~/projects
    git clone https://gitext.gfz-potsdam.de/ajacquey/lynx.git
    cd ~/projects/lynx
    git checkout master

*Note: the "master" branch of LYNX is the "stable" branch which is updated only if all tests are passing.*

#### 3. Compile LYNX
You can compile LYNX by following these instructions:

    cd ~/projects/lynx
    make -j4

#### 4. Test LYNX
To make sure that everything was installed properly, you can run the tests suite of LYNX:

    cd ~/projects/lynx
    ./run_tests -j2

If all the tests passed, then your installation is working properly. You can now use the LYNX simulator!

## Usage
To run LYNX from the command line with multiple processors, use the following command:

    mpiexec -n <nprocs> ~/projects/lynx/lynx-opt -i <input-file>

Where `<nprocs>` is the number of processors you want to use and `<input-file>` is the path to your input file (extension `.i`).  

Information about the structure of the LYNX input files can be found in the documentation (link to follow).

## Cite

If you use LYNX for your work please cite:
* This repository:  
Antoine B. Jacquey, & Mauro Cacace. (2019, July 30). LYNX: Lithosphere dYnamic Numerical toolboX, a MOOSE-based application v1.0.


Please read the [CITATION](https://gitext.gfz-potsdam.de/ajacquey/lynx//blob/master/CITATION) file for more information.

## Publications using LYNX

More to come...
