from __future__ import absolute_import
import numpy as np
from lsst.sims.catUtils.exampleCatalogDefinitions import PhoSimCatalogPoint

__all__= ["TruthCatalogPoint"]

class TruthCatalogPoint(PhoSimCatalogPoint):

    catalog_type = 'truth_catalog_POINT'

    column_outputs = ['prefix', 'uniqueId', 'raJ2000', 'decJ2000', 'phoSimMagNorm', 'sedFilepath',
                      'redshift', 'shear1', 'shear2', 'kappa', 'raOffset', 'decOffset',
                      'spatialmodel', 'galacticExtinctionModel', 'galacticAv', 'galacticRv',
                      'internalExtinctionModel']

    default_columns = [('redshift', 0., float), ('shear1', 0., float), ('shear2', 0., float),
                       ('kappa', 0., float), ('raOffset', 0., float), ('decOffset', 0., float),
                       ('galacticExtinctionModel', 'CCM', (str, 3)), ('galacticRv', 3.1, float),
                       ('internalExtinctionModel', 'none', (str, 4))]

    default_formats = {'S': '%s', 'f': '%.9g', 'i': '%i'}

    spatialModel = "point"

    transformations = {'raJ2000': np.degrees, 'decJ2000': np.degrees}
