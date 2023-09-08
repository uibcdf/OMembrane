from pathlib import PurePath

parent = PurePath(__file__).parent

demo = {}

# POPC membrane

demo['POPC'] = {}
demo['POPC']['POPC.msmpk'] = parent.joinpath('POPC/POPC.msmpk').__str__()
demo['POPC']['POPC.psf'] = parent.joinpath('POPC/POPC.psf').__str__()
demo['POPC']['POPC.dcd'] = parent.joinpath('POPC/POPC.dcd').__str__()



del (PurePath, parent)
