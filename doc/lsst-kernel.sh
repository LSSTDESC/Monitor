#!/usr/bin/env bash

source /global/common/cori/contrib/lsst/lsstDM/setupStack-12_1-Run3.1-a.sh
setup lsst_sims

module load python/2.7-anaconda

#Setup Monitor packages here
source $HOME_DIR/pserv/setup/setup.sh
source $HOME_DIR/Monitor/setup/setup.sh

exec python -m ipykernel $@
