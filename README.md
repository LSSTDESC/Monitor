# The DESC Monitor

Extracts the light curves of ("monitors") all time-variable objects of cosmological interest. The basic concept is a general purpose light curve extraction tool, that can be sub-classed to specialize in supernova monitoring or lensed quasar monitoring. In an ideal world, this code would merely take care of the book-keeping involved in assembling light curves in the appropriate format for modeling codes like `SNCosmo` and `SLTimer` by simply querying the various DM Level 2 databases. However, we are aware of the possibility that we may need to implement multi-object "SuperFit" routines into the Monitor, if the L2 light curves do not provide sufficient accuracy. (This is a particular worry for the `SLMonitor`, which will need to deal with highly blended objects).

## Setup and testing:
From bash
```
$ source <Monitor install directory>/setup/setup.sh
$ python <Monitor install directory>/tests/test_Monitor.py
```

## Demo

Placeholder:
```
    python Monitor/monitor.py
```

## People

* [Bryce Kalmbach](https://github.com/DarkEnergyScienceCollaboration/Monitor/issues/new?body=@jbkalmbach) (UW)
* [Phil Marshall](https://github.com/DarkEnergyScienceCollaboration/Monitor/issues/new?body=@drphilmarshall) (SLAC)

## License etc

This is open source software, available under the BSD license. If you are interested in this project, please do drop us a line via the hyperlinked contact names above, or by [writing us an issue](https://github.com/DarkEnergyScienceCollaboration/Monitor/issues/new).
