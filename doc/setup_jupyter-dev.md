## Steps for setting up jupyter-dev for Monitor

1. ####Setup Kernel Spec

  Setup a `kernelspec` (instructions come from [this issue](https://github.com/jupyterhub/jupyterhub/issues/847#issuecomment-267044166)) by running:
  
  `ipython kernel install --user --name lsst`
  
  and you'll get a directory under:
  
  `~/.local/share/jupyter/kernels/lsst`
  
  with a file named `kernel.json`
  
2. ####Copy the script [lsst-kernel.sh](lsst-kernel.sh).

 Replace $HOME_DIR with where you have pserv and Monitor repos cloned.

3. ####Modify `~/.local/share/jupyter/kernels/lsst/kernel.json`.

  Change:
  
  ```
 "argv": [
  "/usr/bin/python",
  "-m",
  "ipykernel",
  "-f",
  "{connection_file}"
 ],
 ```
 
 to:
 
 ```
 "argv": [   
   "/path/to/lsst-kernel.sh",   
   "-f",   
   "{connection_file}" 
 ],
 ```
 
4. ####Test it out!

 You should be good to go! Try running the first four cells of the [example notebook](../examples/lightcurve_example.ipynb) **making sure the your notebook is running in the `lsst` kernel**. To change kernels use the `Change Kernel` option under `Kernel` in the jupyter notebook menu bar.
 
 The rest will not work yet because we are still working on establishing database connections from jupyter-dev, but once that is done we will amend this.
