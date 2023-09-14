"""
OMembrane
This must be a short description of the project
"""

from ._version import __version__


__documentation_web__ = 'https://www.uibcdf.org/OMembrane'
__github_web__ = 'https://github.com/uibcdf/OMembrane'
__github_issues_web__ = __github_web__ + '/issues'


from . import lipid
from . import build
from . import analysis
#from . import permeability
from .demo import demo

