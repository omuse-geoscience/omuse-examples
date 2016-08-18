"""
script running the Haringvliet test case for SWAN

see Zijlema 2010 (Coastal Engineering 57, 267-277)
"""
import numpy

from omuse.community.swan.interface import Swan
from omuse.units import units

from amuse.io import write_set_to_file

from read_triangle_mesh import read_triangle_mesh

from plot_moments import plot_moments

def read_bot_data(filename="f31hari.bot",n=88,m=117):
    f=open(filename,"r")
    lines=f.readlines()
    f.close()
    
    dat=[]
    for line in lines:
      for s in line.split():
        dat.append(float(s))
    
    return numpy.transpose(numpy.array(dat).reshape((m,n)))

def TestHaringvliet_Zijlema2010(flow=0.0521 | units.rev/units.s,
                                fhigh=1. | units.rev/units.s,
                                msc=32,
                                mdc=36 
                                ):        
    
    bathymetry=read_bot_data()
    bathymetry=numpy.array(bathymetry,dtype="float32")

    rt=read_triangle_mesh("f32hari") 
    rt.read_grid()
    nodes,elements=rt.get_sets() 

    ncells=len(elements)
    nverts=len(nodes)
      
    s=Swan(grid_type="unstructured",redirection="none",channel_type="sockets")
    
    s.parameters.number_of_cells=ncells
    s.parameters.number_of_vertices=nverts
    s.parameters.number_of_directions=mdc
    s.parameters.number_of_frequencies=msc
    s.parameters.lowest_frequency=flow
    s.parameters.highest_frequency=fhigh

    s.parameters.input_grid_origin_x=0. | units.m
    s.parameters.input_grid_origin_y=0. | units.m
    s.parameters.input_grid_orientation=0. | units.deg
    s.parameters.input_grid_dx=250. | units.m
    s.parameters.input_grid_dy=250. | units.m
    s.parameters.input_grid_nmesh_x=87
    s.parameters.input_grid_nmesh_y=116
  
    s.parameters.constant_water_level=1.7 | units.m
  
    s.parameters.uniform_wind_velocity=14. | units.m/units.s
    s.parameters.uniform_wind_direction=8.8 | units.deg
  
    s.parameters.unstructured_boundary_spec_file="f31har03.bnd"
    s.parameters.boundary_marker=2
  
    s.parameters.use_gen3_parameters=True
    s.parameters.use_breaking_parameters=True
    s.parameters.use_triads_parameters=True
    s.parameters.use_friction_parameters=True
    s.parameters.use_input_wind=False
                
    exc=s.get_exc_value(1)        
    bathymetry[bathymetry==-99.]=exc
    
    channel=nodes.new_channel_to(s.nodes)
    channel.copy_attributes(["x","y","vmark"])
    channel=elements.new_channel_to(s.elements)
    channel.copy_attributes(["n1","n2","n3"])

    s.forcings.depth=bathymetry | units.m
            
    s.evolve_model(0. | units.s)

    write_set_to_file(s.nodes,"nodes.amuse","amuse",append_to_file=False)
    write_set_to_file(s.elements,"elements.amuse","amuse",append_to_file=False)

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser(usage=__doc__)
    return result
        
if __name__=="__main__":
    o, arguments  = new_option_parser().parse_args()
  
    TestHaringvliet_Zijlema2010()
    plot_moments().plot()
