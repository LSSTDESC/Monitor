import desc.monitor as monitor
def test_aliasDictionary():

    testSeq = ['mJd', 'band', 'zp']
    aliases = dict(time=['date', 'mjd', 'expMJD'], filter=['band', 'throughput'])
    aliasDict = monitor.aliasDictionary(testSeq, aliases)
    assert aliasDict == {'mJd': 'time', 'band': 'filter'}

def test_map2StandardSeq():

    
    testSeq = ['mJd', 'band', 'zp']
    aliases = dict(time=['date', 'mjd', 'expMJD'], filter=['band', 'throughput'])
    aliasDict = monitor.aliasDictionary(testSeq, aliases)
    assert map

def test_aliasDictionary_withcaps():

    testSeq = ['mJd', 'band', 'zp', 'Flux']
    aliases = dict(time=['date', 'mjd', 'expMJD'],
                   filter=['band', 'throughput'], flux=['counts'])
    aliasDict = monitor.aliasDictionary(testSeq, aliases)
    assert aliasDict == {'mJd': 'time', 'band': 'filter', 'Flux':'flux'}
