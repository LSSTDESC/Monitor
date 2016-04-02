import pandas as pd
import lsst.daf.persistence as dp
from sqlalchemy import create_engine


def forcedSourcesFromRepo(tract, visit, filterName, butler, dataFrame=None):
    dataID = dict(tract=0, visit=visit, filter=filterName)
    t = butler.get('forced_src', dataId=dataId, immediate=True)
    df  = pd.dataFrame()
    for col in t.schema.getNames():
        df[col] = t[t.schema[col].asKey()]
    for key in dataID:
        df[key] = dataID[key]
    if not dataFrame is None:
	return pd.concat(df, dataFrame)
    else:
	return df

def ingestDataFrame(dataBase, exists='append'):


    if exists == 'append':
        df.to_sql(dataBase, con=engine, if_exists='append')
    else:
        df.to_sql(dataBase, con=engine, if_exists='replace')

class ForcedSources(object):

    def __init__(self, repo, database, outputdir):

        self.repo = repo
        self.database = database
        self.butler = dp.Butler(repo, outputRoot=output) 
        engine = create_engine('sqlite:///' + database)


