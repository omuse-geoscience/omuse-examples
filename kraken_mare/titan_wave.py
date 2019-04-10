import numpy
import cPickle

from omuse.community.swan.interface import SwanInterface,Swan
from omuse.units import units
from amuse.units.trigo import cos,sin

from matplotlib import pyplot

from amuse.test.amusetest import TestWithMPI

from amuse.io import write_set_to_file,read_set_from_file

#~ import logging
#~ logging.basicConfig(level=logging.DEBUG)
#~ logging.getLogger("code").setLevel(logging.DEBUG)
  
def runlabel( **kwargs):
    keys=sorted(kwargs)
    args=[(x,kwargs[x]) for x in keys]
    return str(abs(hash(str(args))))

def swan_eq(**kwargs):

        grav=kwargs["grav"]
        rho_air=kwargs["rho_air"]
        rho_water=kwargs["rho_water"]
        flow=kwargs["flow"]
        fhigh=kwargs["fhigh"]
        u10=kwargs["u10"]
        udir=kwargs["udir"]
        gridname=kwargs["gridname"]
        coord=kwargs["coord"]
        depth=kwargs["depth"]
        visc=kwargs["viscosity"]
        planetary_radius=kwargs["planetary_radius"]
        pixel_scale=kwargs["pixel_scale"]
        under_relaxation_factor=kwargs["under_relaxation_factor"]
        msc=kwargs["msc"]
        mdc=kwargs["mdc"]

        print kwargs

        nodes=read_set_from_file("kraken_nodes","amuse",close_file=True)
        elements=read_set_from_file("kraken_elements","amuse",close_file=True)

        s=Swan(grid_type="unstructured", input_grid_type="unstructured", 
            redirection="none",
            coordinates=coord)

        s.parameters.gravitational_acceleration=grav
        s.parameters.air_density=rho_air
        s.parameters.water_density=rho_water
        s.parameters.planetary_radius=planetary_radius
        s.parameters.maximum_error_level=2
        s.parameters.under_relaxation_factor=under_relaxation_factor


        ncells=len(elements)
        nverts=len(nodes)

        s.parameters.number_of_cells=ncells
        s.parameters.number_of_vertices=nverts
        s.parameters.number_of_directions=mdc
        s.parameters.number_of_frequencies=msc
        s.parameters.lowest_frequency=flow
        s.parameters.highest_frequency=fhigh
        s.parameters.minimum_wind_speed=0.1 | units.m/units.s
        s.parameters.max_iterations_stationary=50
        s.parameters.use_gen3_parameters=True
        s.parameters.use_breaking_parameters=True
        s.parameters.use_triads_parameters=True
        s.parameters.use_friction_parameters=True
        s.parameters.use_input_wind=True

        if visc.number>0: 
            s.parameters.use_input_turbulent_visc=True
            s.parameters.turbulent_viscosity_factor=1.

        #~ s.parameters.uniform_wind_velocity=u10
        #~ s.parameters.uniform_wind_direction=udir

        print s.parameters
        #~ raise

        channel=nodes.new_channel_to(s.nodes)
        if coord=="spherical":
          channel.copy_attributes(["lon","lat","vmark"])
        else:
          channel.transform(["x","y","vmark"],lambda x,y,z: (x*pixel_scale,y*pixel_scale,z),["x","y","vmark"])
          #~ channel.copy_attributes(["x","y","vmark"])          
        channel=elements.new_channel_to(s.elements)
        channel.copy_attributes(["n1","n2","n3"])

        #~ print abs(s.nodes.lon-nodes.lon).max()
        #~ print abs(s.nodes.lat-nodes.lat).max()
    
        forcings=s.forcings.empty_copy()

        if depth=="default":
          forcings.depth=nodes.depth
        elif depth=="x2":
          forcings.depth=nodes.depth*2
        elif depth=="x10":
          forcings.depth=nodes.depth*10
        elif depth=="sqrt":
          maxdepth=nodes.depth.max()
          forcings.depth=maxdepth*(nodes.depth/maxdepth)**0.5
        else:
          raise Exception("unknown depth option")

        wind_vx=u10*cos(udir)
        wind_vy=u10*sin(udir)

        if coord=="spherical":
          forcings.wind_vx=wind_vx
          forcings.wind_vy=wind_vy
        else:          
          forcings.wind_vx=wind_vx*sin(nodes.lon)+wind_vy*cos(nodes.lon)
          forcings.wind_vy=-wind_vx*cos(nodes.lon)+wind_vy*sin(nodes.lon)

          #~ pyplot.quiver(s.forcings.x.number,s.forcings.y.number,
             #~ forcings.wind_vx.number,forcings.wind_vy.number, scale=50)
          #~ pyplot.show()
          #~ raise

        if visc.number>0: 
            forcings.visc=visc
            forcings.new_channel_to(s.forcings).copy_attributes(["depth","wind_vx","wind_vy","visc"])
        else:
            forcings.new_channel_to(s.forcings).copy_attributes(["depth","wind_vx","wind_vy"])
    
        s.evolve_model(0. | units.s)
                        
        label=runlabel(**kwargs)
        print "label:",label
        with open("args_"+label+".pkl","w") as f:
            cPickle.dump(kwargs,f)
        write_set_to_file(s.nodes,gridname+"_nodes_eq_"+label+".amuse","amuse",append_to_file=False)

def new_option_parser(usage=__doc__):
    from amuse.units.optparse import OptionParser
    result = OptionParser(usage=usage)
    result.add_option("--wind", 
        dest="u10", default = 2. | units.m/units.s, unit=units.m/units.s, type=float, 
        help="wind velocity at 10m in %unit [%default]")
    result.add_option("--winddir", 
        dest="udir", default = 0. | units.deg, unit=units.deg, type=float, 
        help="wind direction at 10m in %unit [%default]")
    result.add_option("--flow", 
        dest="flow", default = 0.02 | units.rev/units.s, unit=units.rev/units.s, type=float, 
        help="lowest freq [%default]")
    result.add_option("--fhigh", 
        dest="fhigh", default = 0.5 | units.rev/units.s, unit=units.rev/units.s, type=float, 
        help="high freq [%default]")
    result.add_option("--grav", 
        dest="grav", default = 1.352 | units.m/units.s**2, unit=units.m/units.s**2, type=float, 
        help="gravitational acceleration [%default]")
    result.add_option("--air", 
        dest="rho_air", default = 5.4 | units.kg/units.m**3, unit=units.kg/units.m**3, type=float, 
        help="density of air [%default]")
    result.add_option("--water", 
        dest="rho_water", default = 600. | units.kg/units.m**3, unit=units.kg/units.m**3, type=float, 
        help="density of water [%default]")
    result.add_option("--gridname", 
        dest="gridname", default = "kraken", type=str, 
        help="name of input grid [%default]")
    result.add_option("--coord", 
        dest="coord", default = "spherical", type=str, 
        help="coordinate system (cartesian or spherical) [%default]")
    result.add_option("--depth", 
        dest="depth", default = "default", type=str, 
        help="depth option (default, x10, ..) [%default]")
    result.add_option("--visc", 
        dest="viscosity", default = 0. | units.cm**2/units.s, unit=units.cm**2/units.s, type=float, 
        help="kinematic/turbulent viscosity in %unit [%default]")        
    result.add_option("--Rplanet", 
        dest="planetary_radius", default = 2575.5 | units.km, unit=units.km, type=float, 
        help="Planetary radius (for spherical only) in %unit [%default]")
    result.add_option("--pixel", 
        dest="pixel_scale", default = 5.4 | units.km, unit=units.km, type=float, 
        help="conversion from projection map coordinates (cartesian only) %unit [%default]")
    result.add_option("--frelax", 
        dest="under_relaxation_factor", default = 0.001, type=float, 
        help="solver under relaxation factor [%default]")
    result.add_option("--nfreq", 
        dest="msc", default = 32, type=int, 
        help="number of frequency bins [%default]")
    result.add_option("--ndir", 
        dest="mdc", default = 36, type=int, 
        help="number of directional bins [%default]")
    return result
        
if __name__=="__main__":
    p=new_option_parser()
    o, arguments  = p.parse_args()
    test=swan_eq(**vars(o))
