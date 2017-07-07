import os
script_path = os.path.dirname(os.path.realpath(__file__))


import rpy2.robjects as ro
ro.r('''
	source('%s/callDE.R')
''' % script_path)

_runLimma = ro.globalenv['runLimma']
_runVoomLimma = ro.globalenv['runVoomLimma']

def runLimma(data, rids, sample_class):
    res = _runLimma(data, rids, sample_class)
    res = dict(zip(res.names, map(list,list(res))))
    return res

def runVoomLimma(data, rids, sample_class):
    res = _runVoomLimma(data, rids, sample_class)
    res = dict(zip(res.names, map(list,list(res))))
    return res

import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate() # to convert numpy matrix to R matrix

def list2factor(l):
	## convert a python list of int 
	## to an R factor 
	fac = ro.IntVector(l)
	fac = ro.FactorVector(fac)
	return fac

def numpy_to_pandas_df(mat, rids, cids, rid_name='genes'):
	## convert mat, rids, cids to a pandas df 
	d = dict(zip(cids, map(lambda x : x.tolist(), mat.T)))
	d[rid_name] = rids
	df = pd.DataFrame(d)
	ordered_cols = [rid_name] + cids
	df = df[ordered_cols]
	return df



