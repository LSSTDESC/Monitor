# The DESC Monitor
![Travis CI build status](https://travis-ci.org/LSSTDESC/Monitor.svg)

Extracts the light curves of ("monitors") all time-variable objects of cosmological interest. The basic concept is a general purpose light curve extraction tool, that can be sub-classed to specialize in supernova monitoring or lensed quasar monitoring. In an ideal world, this code would merely take care of the book-keeping involved in assembling light curves in the appropriate format for modeling codes like `SNCosmo` and `SLTimer` by simply querying the various DM Level 2 databases. However, we are aware of the possibility that we may need to implement multi-object extended `MultiFit` (a.k.a. "SuperFit") routines into the `Monitor`, if the Level 2 light curves do not provide sufficient accuracy. (This is a particular worry for the `SLMonitor`, which will need to deal with highly blended objects). The `Monitor` is therefore being developed as "Level 3" code, following DM standards and module structure.

## Setup and testing:
From bash:
```
$ source setup/setup.sh
$ nosetests
```
or use the c-shell alternative. The `Monitor` uses some DM stack code, notably the `Butler`; see the [installation notes](https://github.com/DarkEnergyScienceCollaboration/Monitor/blob/master/INSTALL.md) for help getting set up.

##### Setting up ssh tunnel for database access:

In order to produce light curves with output from the [Twinkles](https://github.com/DarkEnergyScienceCollaboration/Twinkles/tree/master/python/desc/twinkles) project one needs to be able to access the SQL database where DM processed output is stored. To do that one needs to setup an ssh tunnel for access. We use the same tools as the connection to the UW LSST CATSIM Database with instructions [here](https://confluence.lsstcorp.org/display/SIM/Accessing+the+UW+CATSIM+Database). Follow the step at the beginning to install the necessary tools and then replace step 1 code with the following command line entry with your NERSC username in the proper place:

```
ssh -L 3307:scidb1.nersc.gov:3306 your_username@cori.nersc.gov
```

Do not worry about step 2 and continue to step 3 where you are instructed to create a db-auth.paf file, but replace the code on the website with the following parameters:
```
database: {
  authInfo: {
    host: '127.0.0.1'
    port: 3307
    user: $db_username
    password: $db_password
  }
}
```
If you do not have the $db_username or $db_password and are interested in access please contact a member of the Monitor team for more information.

## Demo

We have a couple example notebooks available for users to get started using the `Monitor`:

* [depth_curve_example.ipynb](examples/depth_curve_example.ipynb)
  * Provides an introduction on the `Monitor`'s light curve construction tools and shows how to open a database connection to get data from NERSC. This is the best starting point if you want to try using the `Monitor` for the first time.

* [simple_error_model.ipynb](examples/simple_error_model.ipynb)
  * An introduction to the analysis possible with the `Monitor`. This is where we are showing results for our latest analysis on the [Twinkles](https://github.com/DarkEnergyScienceCollaboration/Twinkles/tree/master/python/desc/twinkles) project. It is available for other users to use as a jumping off point for their own analysis.


## People

* [Bryce Kalmbach](https://github.com/DarkEnergyScienceCollaboration/Monitor/issues/new?body=@jbkalmbach) (UW)
* [Phil Marshall](https://github.com/DarkEnergyScienceCollaboration/Monitor/issues/new?body=@drphilmarshall) (SLAC)
* [Jim Chiang](https://github.com/DarkEnergyScienceCollaboration/Monitor/issues/new?body=@jchiang87) (SLAC)
* [Rahul Biswas](https://github.com/DarkEnergyScienceCollaboration/Monitor/issues/new?body=@rbiswas4) (UW)

## License etc

This is open source software, available under the BSD license. If you are interested in this project, please do drop us a line via the hyperlinked contact names above, or by [writing us an issue](https://github.com/DarkEnergyScienceCollaboration/Monitor/issues/new).
