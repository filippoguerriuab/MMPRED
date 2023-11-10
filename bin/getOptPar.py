from pathlib import Path


def getOptPar(iargs, parName, parDefValue, parType=None, optArgOptArgs=None):

	if parName in iargs:
		par_value_idx = iargs.index(parName)+1
		parValue = iargs[par_value_idx]
		if parType == 'float':
			parValue = float(parValue)	

		elif parType == 'int':
			parValue = int(parValue)

		elif parType == 'list':
			sep = getOptPar([], 'sep', ',')	
			parValue = sep.split(parValue)

		elif parType == 'path':
			parValue = Path(parValue)

		elif parType == 'str':
			parValue = str(parValue)
	

	else:
		parValue = parDefValue
		

	return parValue
    

