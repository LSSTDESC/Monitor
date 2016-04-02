"""
Utilities that we will often be using
"""
from future.builtins import dict
__all__ = ['aliasDictionary', 'mapSeq2Standard']

def aliasDictionary(sequence, aliases):
    """
    Given a dictionary of
    aliases of the form {standard_string_names: list of possible aliases}
    this finds occurances of aliases (ignoring capitalization) in the sequence
    of strings called sequence and returns a dictionary
    Dict(alias: standard_string_name)

    Parameters
    ----------
    sequence : a sequence of strings
    aliases : Dictionary with list-valued values
        dictionary with keys = standard desired names, value = possible aliases
        It is assumed that a single alias only works for a single key

    Returns
    -------
    Dictionary with keys in the sequence of strings that have been identified
    as a possible alias of standard_string_names with
    value =standard_string_names

    Examples
    --------
    >>> aliases = dict(time=['mjd','expmjd'], flux=['counts'], fluxerr=['flux_err', 'fluxerror'], zpsys=['magsys'])
    >>> testSeq = ['mJd', 'band', 'zp', 'Flux', 'fluxError', 'zpsys']
    >>> aliasDictionary(testSeq, aliases) == {'Flux': 'flux', 'fluxError': 'fluxerr', 'mJd': 'time'}
    True

    Notes:
    -----
    1. It is obviously imperative that the list of possible aliases combined
    for all keys should not have duplicates
    2. The dictionary can be used to change non-standard column names in a
    `pd.DataFrame` through the use of `dataFrame.rename(columns=aliasDictionary,
    inplace=True)`
    """

    # Check that there is no string which might be the alias of more than one
    # variable
    _valueList = list(x for l in aliases.values() for x in l)
    assert len(set(_valueList)) == len(_valueList)

    # insert the standard name in the list of values
    # (necesary for the case where variable names are std names with caps)
    _ = [aliases.update({key: [key] + aliases[key]}) for key in aliases.keys()]

    # invert aliases dictionary to get
    inverse_dict = dict()
    _ = [inverse_dict.update({alias:aliastuple[0]})
         for aliastuple in aliases.iteritems() for alias in aliastuple[1]]

    # Account for capitalization issues
    testDict = dict((s.lower(), s) for s in sequence)
    aliasedSet = set(inverse_dict.keys()) & set(testDict.keys())


    # return the dictionary
    aliasDict = dict((testDict[key], inverse_dict[key]) for key in aliasedSet
                     if testDict[key] != inverse_dict[key])

    return aliasDict

def mapSeq2Standard(sequence, aliasDict):
    """
    given a sequence of strings return a list of the strings, replacing those
    strings in the input sequence which are keys in the dictionary by their
    dictionary values

    Parameters
    ---------
    sequence : sequence of strings

    aliasDict : Dictionary
        Dictionary of aliases where the keys are strings in the sequence and
        the values are the standard names associated with those strings

    Returns
    -------
    list of strings

    Examples
    --------
    >>> aliases = dict(time=['mjd','expmjd'], flux=['counts'], fluxerr=['flux_err', 'fluxerror']) 
    >>> testSeq = ['mJd', 'band', 'zp', 'Flux', 'fluxError']
    >>> aliasDict = aliasDictionary(testSeq, aliases)
    >>> mapSeq2Standard(testSeq, aliasDict) == ['time', 'band', 'zp', 'flux', 'fluxerr']
    True

    Notes
    -----

    """
    std = [aliasDict[x] if x in aliasDict else x for x in sequence]
    return std
