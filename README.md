# The DESC Monitor
<!-- ![Travis CI build status](https://travis-ci.org/DarkEnergyScienceCollaboration/Monitor.svg) Travis will struggle with DM stack dependencies... -->

Extracts the light curves of ("monitors") all time-variable objects of cosmological interest. The basic concept is a general purpose light curve extraction tool, that can be sub-classed to specialize in supernova monitoring or lensed quasar monitoring. In an ideal world, this code would merely take care of the book-keeping involved in assembling light curves in the appropriate format for modeling codes like `SNCosmo` and `SLTimer` by simply querying the various DM Level 2 databases. However, we are aware of the possibility that we may need to implement multi-object extended `MultiFit` (a.k.a. "SuperFit") routines into the `Monitor`, if the Level 2 light curves do not provide sufficient accuracy. (This is a particular worry for the `SLMonitor`, which will need to deal with highly blended objects). The `Monitor` is therefore being developed as "Level 3" code, following DM standards and module structure.

## Setup and testing:
From bash:
```
$ source setup/setup.sh
$ nosetests
```
or use the c-shell alternative. The `Monitor` uses some DM stack code, notably the `Butler`; see the [installation notes](https://github.com/DarkEnergyScienceCollaboration/Monitor/blob/master/INSTALL.md) for help getting set up.

## Demo

See our demo ipython notebook [here](https://github.com/DarkEnergyScienceCollaboration/Monitor/blob/master/examples/lightcurve_example.ipynb).


## People

* [Bryce Kalmbach](https://github.com/DarkEnergyScienceCollaboration/Monitor/issues/new?body=@jbkalmbach) (UW)
* [Phil Marshall](https://github.com/DarkEnergyScienceCollaboration/Monitor/issues/new?body=@drphilmarshall) (SLAC)

## License etc

This is open source software, available under the BSD license. If you are interested in this project, please do drop us a line via the hyperlinked contact names above, or by [writing us an issue](https://github.com/DarkEnergyScienceCollaboration/Monitor/issues/new).
