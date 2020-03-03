"""
minimal OMUSE example
"""

from omuse.units import units
from omuse.community.qgmodel.interface import QGmodel

code=QGmodel(redirection="none")

code.parameters.dt=0.5 | units.hour

code.evolve_model(1.| units.day)

print(code.grid.psi.max().in_(units.Sv/units.km))
