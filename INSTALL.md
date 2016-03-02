# Installation Notes

The Monitor uses the LSST DM Stack's Butler, and also the SNCosmo package (among other things). Here's how to get everything working.

## 1. Get Anaconda Python

According to the [LSST community](https://community.lsst.org/t/up-and-running-with-sims-maf-contrib/383), and
various peopls around the DESC, it seems best to install the DM stack using `conda` - or in fact, `miniconda`.
Click [here](http://conda.pydata.org/miniconda.html) and download python 2.7, if you don't have Anaconda python already.
```
bash ~/Downloads/Miniconda-latest-MacOSX-x86_64.sh
```
This installs miniconda into `${HOME}/miniconda2`. You should now pre-pend your `PATH` with `${HOME}/minconda2/bin` so that
this is becomes your default version of python.

## 2. Install the LSST DM Stack

First you need to add the LSST "channel:"
```
conda config --add channels http://research.majuric.org/conda/stable
```
Now do:
```
conda install lsst-distrib
```
and find something else to do for half an hour. This will install the LSST packages in `${HOME}/miniconda2/pkgs/`.

## 3. Get set up to use the LSST DM Stack

See the project's notes [here](https://confluence.lsstcorp.org/display/LSWUG/Using+the+LSST+Stack).
LSST packages have to be added to your PATH, etc. before you can use them. This is handled by `eups`.
For analyzing simulated images, we need to do:
```
    source eups-setups.sh
    setup -T v11_0 obs_lsstSim --keep
```
This needs to be done every time you start a new shell - so these lines could be worth adding to your `.bashrc` file or equivalent.

## 4. Install python packages

We use several useful python packages that don't come with the DM stack:
```
pip install sncosmo
pip install pandas
pip install nose
```
