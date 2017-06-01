# Example Jupyter Notebooks

## Intro notebooks

We have a couple example notebooks available for users to get started using the `Monitor`:

* [depth_curve_example.ipynb](depth_curve_example.ipynb)
  * Provides an introduction on the `Monitor`'s light curve construction tools and shows how to open a database connection to get data from NERSC. This is the best starting point if you want to try using the `Monitor` for the first time.

* [simple_error_model.ipynb](simple_error_model.ipynb)
  * An introduction to the analysis possible with the `Monitor`. This is where we are showing results for our latest analysis on the [Twinkles](https://github.com/DarkEnergyScienceCollaboration/Twinkles/tree/master/python/desc/twinkles) project. It is available for other users to use as a jumping off point for their own analysis.

## Additional examples

* [create_truth_db.ipynb](create_truth_db.ipynb)
  * Shows how we create a "truth" database with the simulation inputs for a given set of visits using an LSST Opsim database. This "truth" database can then be used to compare to the results from processing simulated images with the LSST DM stack.