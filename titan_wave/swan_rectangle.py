import numpy
import pickle

from omuse.community.swan.interface import SwanInterface,Swan
from omuse.units import units

from matplotlib import pyplot

from amuse.test.amusetest import TestWithMPI

from amuse.io import write_set_to_file

#~ import logging
#~ logging.basicConfig(level=logging.DEBUG)
#~ logging.getLogger("code").setLevel(logging.DEBUG)

def runlabel( **kwargs):
    keys=sorted(kwargs)
    args=[(x,kwargs[x]) for x in keys]
    return str(abs(hash(str(args))))

def rectangle(**kwargs):
        print("kwargs:", kwargs)
        label=runlabel(**kwargs)
        print("label:",label)

        grav=kwargs["grav"]
        rho_air=kwargs["rho_air"]
        rho_water=kwargs["rho_water"]
        L=kwargs["L"]
        N=kwargs["N"]
        flow=kwargs["flow"]
        fhigh=kwargs["fhigh"]
        u10=kwargs["u10"]
        aspect=kwargs["aspect"]
        depth=kwargs["depth"]
        under_relaxation_factor=kwargs["under_relaxation_factor"]
        msc=kwargs["msc"]
        mdc=kwargs["mdc"]
        visc=kwargs["viscosity"]
        
        s=Swan(redirection="none")

        s.parameters.gravitational_acceleration=grav
        s.parameters.air_density=rho_air
        s.parameters.water_density=rho_water
        #~ s.parameters.maximum_error_level=2
        s.parameters.under_relaxation_factor=under_relaxation_factor

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
        s.parameters.minimum_wind_speed=0.1 | units.m/units.s


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

        if visc.number>0: 
            s.parameters.use_input_turbulent_visc=True
            s.parameters.turbulent_viscosity_factor=1.
            
        s.parameters.use_gen3_parameters=True
        s.parameters.use_breaking_parameters=True
        s.parameters.use_triads_parameters=True
        s.parameters.use_friction_parameters=True

        s.forcings.depth=depth
        if visc.number>0: 
            s.forcings.visc=visc
            
        s.evolve_model(0. | units.s)
                
        with open(label+".args","w") as f:
            pickle.dump(kwargs,f)
        write_set_to_file(s.grid.copy(),label+".amuse","amuse",append_to_file=False)

def new_option_parser(usage=None):
    from amuse.units.optparse import OptionParser
    result = OptionParser(usage=usage)
    result.add_option("--wind", 
        dest="u10", default = 10. | units.m/units.s, unit=units.m/units.s, type=float, 
        help="wind velocity at 10m [%default]")
    result.add_option("--depth", 
        dest="depth", default = 100. | units.m, unit=units.m, type=float, 
        help="depth [%default]")
    result.add_option("--box", 
        dest="L", default = 100. | units.km, unit=units.km, type=float, 
        help="boxsize [%default]")
    result.add_option("--flow", 
        dest="flow", default = 0.0521 | units.rev/units.s, unit=units.rev/units.s, type=float, 
        help="lowest freq [%default]")
    result.add_option("--fhigh", 
        dest="fhigh", default = 1. | units.rev/units.s, unit=units.rev/units.s, type=float, 
        help="high freq [%default]")
    result.add_option("--grav", 
        dest="grav", default = 9.81 | units.m/units.s**2, unit=units.m/units.s**2, type=float, 
        help="gravitational acceleration [%default]")
    result.add_option("--air", 
        dest="rho_air", default = 1.28 | units.kg/units.m**3, unit=units.kg/units.m**3, type=float, 
        help="density of air [%default]")
    result.add_option("--water", 
        dest="rho_water", default = 1025. | units.kg/units.m**3, unit=units.kg/units.m**3, type=float, 
        help="density of water [%default]")
    result.add_option("--N", 
        dest="N", default = 100, type=int, 
        help="resolution [%default]")        
    result.add_option("--aspect", 
        dest="aspect", default = 0.2, type=float, 
        help="aspect ratio of domain [%default]")
    result.add_option("--frelax", 
        dest="under_relaxation_factor", default = 0.001, type=float, 
        help="solver under relaxation factor [%default]")
    result.add_option("--nfreq", 
        dest="msc", default = 32, type=int, 
        help="number of frequency bins [%default]")
    result.add_option("--ndir", 
        dest="mdc", default = 36, type=int, 
        help="number of directional bins [%default]")
    result.add_option("--visc", 
        dest="viscosity", default = 0. | units.cm**2/units.s, unit=units.cm**2/units.s, type=float, 
        help="kinematic/turbulent viscosity in %unit [%default]")
    return result
        
if __name__=="__main__":
    p=new_option_parser()
    o, arguments  = p.parse_args()
    test=rectangle(**vars(o))
