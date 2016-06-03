from __future__ import absolute_import, division, print_function

import os

import pymssql
from lsst.utils import getPackageDir
import lsst.sims.catUtils.baseCatalogModels as bcm
from lsst.daf.persistence import DbAuth

config = bcm.BaseCatalogConfig()
config.load(os.path.join(getPackageDir("sims_catUtils"), "config", "db.py"))
DBConnection = pymssql.connect(user=DbAuth.username(config.host, config.port),
                               password=DbAuth.password(config.host, config.port),
                               database=config.database, port=config.port)
db = DBConnection.cursor()


def _build_columns(columnVars):
    return tuple(xx.strip() for xx in columnVars.split(','))

def get_snParams(db, snidSeq):
    """
    """
    columnVars='snid, redshift, snra, sndec, x0, x1, c'
    config = bcm.BaseCatalogConfig()
    config.load(os.path.join(getPackageDir("sims_catUtils"), "config", "db.py"))
    query = """SELECT """ 
    query += columnVars
    query += """ FROM twinkSN WHERE snid in {}""".format(snidSeq)
    print(query)
    db.execute(query)
    return db.fetchall()


if __name__ == '__main__':
    
    x = get_snParams(db, ((6001163623700, 6000324908000)))
    print(x)
