"""
script to calculate growth of wave as function of fetch using SWAN
"""

import numpy
import pickle

from omuse.units import units

from omuse.community.swan.interface import Swan

from amuse.io import write_set_to_file

def rectangle(**kwargs):
    print("calculating...")
    
    L=kwargs.get("L", 100 | units.km)
    N=kwargs.get("N",100)
    flow=kwargs.get("flow", 0.0521 | units.rev/units.s)
    fhigh=kwargs.get("fhigh",1. | units.rev/units.s)
    u10=kwargs.get("u10", 10. | units.m/units.s)
    aspect=kwargs.get("aspect",0.2)
    depth=kwargs.get("depth",100 | units.m)
    msc=kwargs.get("msc",32)
    mdc=kwargs.get("mdc",36)
    
    s=Swan()

    s.parameters.grid_origin_x=0. | units.m
    s.parameters.grid_origin_y=0. | units.m
    s.parameters.grid_orientation=0. | units.deg
    s.parameters.grid_length_x=aspect*L
    s.parameters.grid_length_y=L
    s.parameters.grid_nmesh_x=int(aspect*N)
    s.parameters.grid_nmesh_y=N
    s.parameters.number_of_directions=mdc
    s.parameters.number_of_frequencies=msc
    s.parameters.lowest_frequency=flow
    s.parameters.highest_frequency=fhigh

    s.parameters.input_grid_origin_x=0. | units.m
    s.parameters.input_grid_origin_y=0. | units.m
    s.parameters.input_grid_orientation=0. | units.deg
    s.parameters.input_grid_dx=L/N
    s.parameters.input_grid_dy=L/N
    s.parameters.input_grid_nmesh_x=int(aspect*N)
    s.parameters.input_grid_nmesh_y=N
  
    s.parameters.uniform_wind_velocity=u10
    s.parameters.uniform_wind_direction=90. | units.deg
  
    s.parameters.wrap_x_coordinate=True
        
    s.parameters.use_gen3_parameters=True
    s.parameters.use_breaking_parameters=True
    s.parameters.use_triads_parameters=True
    s.parameters.use_friction_parameters=True

    s.forcings.depth=depth
  
    s.evolve_model(0. | units.s)
  
    label=str(abs(hash(str(kwargs))))
    with open(label+".args","w") as f:
        pickle.dump(kwargs,f)
    write_set_to_file(s.grid,label+".amuse","amuse",append_to_file=False)
    print("done")

def new_option_parser(usage=None):
    from amuse.units.optparse import OptionParser
    result = OptionParser(usage=usage)
    result.add_option("--msc", 
        dest="msc", default = 32 , type=int, 
        help="number of freq. bins [%default]")
    result.add_option("--mdc", 
        dest="mdc", default = 36 , type=int, 
        help="number of directional bins [%default]")
    result.add_option("--wind", 
        dest="u10", default = 10. | units.m/units.s, unit=units.m/units.s, type=float, 
        help="wind velocity at 10m in %unit [%default]")
    result.add_option("--depth", 
        dest="depth", default = 100. | units.m, unit=units.m, type=float, 
        help="depth in %unit [%default]")
    result.add_option("--box", 
        dest="L", default = 100. | units.km, unit=units.km, type=float, 
        help="boxsize in %unit [%default]")
    result.add_option("--flow", 
        dest="flow", default = 0.0521 | units.rev/units.s, unit=units.rev/units.s, type=float, 
        help="lowest freq in %unit [%default]")
    result.add_option("--fhigh", 
        dest="fhigh", default = 1. | units.rev/units.s, unit=units.rev/units.s, type=float, 
        help="high freq in %unit [%default]")
    result.add_option("--N", 
        dest="N", default = 100, type=int, 
        help="resolution [%default]")        
    result.add_option("--aspect", 
        dest="aspect", default = 0.2, type=float, 
        help="aspect ratio of domain [%default]")        
    return result
        
if __name__=="__main__":
    p=new_option_parser(__doc__)
    o, arguments  = p.parse_args()
    test=rectangle(**vars(o))
