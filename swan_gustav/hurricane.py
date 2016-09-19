"""
script to to calculate wave field during hurricane evolution.

The computational domain is the north western atlantic (from the example 
grid at http://www.caseydietrich.com). SWAN is used to calculate the 
evolution of the action density describing the wave field spectrum under 
the forcing of a hurricane model (Holland model) with a storm track  in 
atcf (http://www.nrlmry.navy.mil/atcf_web/docs/database/new/abrdeck.html) 
format as input.
"""

import numpy

from omuse.community.swan.interface import Swan

from omuse.community.adcirc.read_grid import adcirc_grid_reader, adcirc_parameter_reader
from omuse.community.adcirc.write_grid import adcirc_grid_writer, adcirc_parameter_writer

from omuse.ext.hurricane_models import HollandHurricane

from omuse.units import units

from matplotlib import pyplot,tri

from amuse.io import write_set_to_file

from amuse.units.quantities import VectorQuantity
    
class SwanHurricane(object):

    def __init__(self, grid_file="grid.input", met_file="fort.22", dt=600 | units.s,
                  flow=2*numpy.pi*0.0521 | units.rad/units.s,
                  fhigh=2*numpy.pi | units.rad/units.s,msc=32,mdc=36):
        self._dt=dt
        self._msc=msc
        self._mdc=mdc
        self._flow=flow
        self._fhigh=fhigh
        self.read_grid(grid_file)
        self.initialize_met(met_file)
        self.initialize_swan()
        self.refresh_nodes()

    @property
    def parameters(self):
        return self.swan.parameters
    @property
    def model_time(self):
        return self.swan.model_time
    @property
    def swan_nodes(self):
        return self.swan.nodes
    @property
    def swan_forcings(self):
        return self.swan.forcings
    @property
    def met_nodes(self):
        return self.met.nodes
    @property
    def swan_elements(self):
        return self.swan.elements
    
    def read_grid(self, grid_file):
        gr=adcirc_grid_reader(grid_file,coordinates="spherical")
        gr.read_grid()
        nodes,elements,elev_boundaries,flow_boundaries=gr.get_sets()
        
        elements.n1=elements.nodes[:,0]
        elements.n2=elements.nodes[:,1]
        elements.n3=elements.nodes[:,2]
        
        i=0
        for e in elev_boundaries+flow_boundaries:
          i+=1
          nodes[e.nodes].vmark=i 

        self.nodes=nodes
        self.elements=elements
        self._elev_boundaries=elev_boundaries
        self._flow_boundaries=flow_boundaries
        self._neta=gr.parameters["NETA"]

    def initialize_met(self, met_file):
        self.met=HollandHurricane(self.nodes, file_or_track=met_file)

    def initialize_swan(self):
        swan=Swan(grid_type="unstructured", input_grid_type="unstructured",
                  mode="dynamic", coordinates="spherical")
        
        swan.parameters.number_of_cells=len(self.elements)
        swan.parameters.number_of_vertices=len(self.nodes)
        swan.parameters.number_of_directions=self._mdc
        swan.parameters.number_of_frequencies=self._msc
        swan.parameters.lowest_frequency=self._flow
        swan.parameters.highest_frequency=self._fhigh
      
        swan.parameters.use_gen3_parameters=True
        swan.parameters.use_breaking_parameters=True
        swan.parameters.use_triads_parameters=True
        swan.parameters.use_friction_parameters=True
        
        swan.parameters.use_input_wind=True
        swan.parameters.use_input_current=False
        swan.parameters.use_input_water_level=False

        swan.parameters.use_csigma_cfl_limiter=True             
        swan.parameters.use_ctheta_cfl_limiter=True             
                     
        swan.parameters.timestep=self._dt
            
        channel=self.nodes.new_channel_to(swan.nodes)
        channel.copy_attributes(["lon","lat","vmark"])
        channel=self.elements.new_channel_to(swan.elements)
        channel.copy_attributes(["n1","n2","n3"])
          
        forcings=swan.forcings.empty_copy()
      
        forcings.depth=self.nodes.depth
        forcings.wind_vx=0. | units.m/units.s
        forcings.wind_vy=0. | units.m/units.s
        
        channel=forcings.new_channel_to(swan.forcings)
        channel.copy_attributes(["depth","wind_vx","wind_vy"])
        self.swan=swan
      
    def evolve_model(self,tend):
        tnow=self.model_time
        timestep=self._dt
        met_to_swan=self.met.nodes.new_channel_to(self.swan.forcings)

        while tnow<tend-timestep/2: 
            self.met.evolve_model(tnow+timestep/2)
            met_to_swan.transform( ["wind_vx","wind_vy"], None, ["vx","vy"])
            self.swan.evolve_model(tnow+timestep)
            self.met.evolve_model(tnow+timestep)            
            tnow=self.model_time

        self.refresh_nodes()  

    def refresh_nodes(self):
        self.swan.nodes.new_channel_to(self.nodes).copy_attributes(["ac2", "wave_tau_x","wave_tau_y"])
        self.met.nodes.new_channel_to(self.nodes).transform( ["wind_vx","wind_vy","pressure"], None, ["vx","vy","pressure"])
        
class plot_swh(object):

    def moment(self,f,df,dtheta,ac2,n=0):
        return df*dtheta*(f[:,:]**(n+2)*ac2).sum(axis=-1).sum(axis=-1)/(2*numpy.pi)**(n)

    def __init__(self,nodes, elements, flow=2*numpy.pi*0.0521 | units.rad/units.s,
                  fhigh=2*numpy.pi | units.rad/units.s,msc=32,mdc=36):
        pyplot.ion()
        fig=pyplot.figure(figsize=(16,6))
        self.fig=fig
        pyplot.show()
          
        x=nodes.lon.number
        y=nodes.lat.number
        n1=elements.n1
        n2=elements.n2
        n3=elements.n3
  
        self.elements=elements
        self.nodes=nodes
  
        elements=numpy.column_stack((n1,n2,n3))
        self.triangulation=tri.Triangulation(x,y,elements)
          
        thetas=numpy.arange(mdc)*2*numpy.pi/mdc
        fac=(fhigh/flow)**(1./(msc-1))
        fs=flow*fac**numpy.arange(msc)
        fs=VectorQuantity.new_from_array(fs)
        
        f,theta=numpy.meshgrid(fs ,thetas)

        self.f=VectorQuantity.new_from_array(f)
  
        self.dlf=numpy.log(fac)
        self.dtheta=2*numpy.pi/mdc

    def redraw(self):
        self.fig.clf()
        self.f1=pyplot.subplot(121)
        self.f2=pyplot.subplot(122)
        self.f1.set_aspect('equal')
        self.f2.set_aspect('equal')

        val=self.nodes.pressure.value_in(units.mbar)
        im=self.f1.tripcolor(self.triangulation,val,shading='gouraud')
        self.f1.triplot(self.triangulation,lw=0.25)
        cb=pyplot.colorbar(im,ax=self.f1)
        cb.set_label("atmospheric pressure (mbar)")

        m0=self.moment(self.f,self.dlf,self.dtheta,self.nodes.ac2,0)
        swh=4*(m0)**0.5
                      
        self.f2.cla()      
        val=swh.value_in(units.m)
        im=self.f2.tripcolor(self.triangulation,val,shading='gouraud',cmap="bwr")
        self.f2.triplot(self.triangulation,lw=0.25)
        cb=pyplot.colorbar(im,ax=self.f2)
        cb.set_label("significant wave height (m)")
        
        pyplot.draw()
        
def evolve(code,tend=11. | units.day,timestep=0.5 | units.hour,p=None):
    i=0
    tnow=code.model_time

    write_set_to_file(code.nodes, "nodes-%6.6i"%i,"amuse",append_to_file=False)
    write_set_to_file(code.elements, "elements-%6.6i"%i,"amuse",append_to_file=False)

    print "starting main loop.."
    print "(this may take a while)"

    while tnow<tend: 
        i+=1
        code.evolve_model(tnow+timestep)
        tnow=code.model_time
        print i,tnow,tnow.in_(units.day)
        if p:
            p.redraw()

        write_set_to_file(code.nodes, "nodes-%6.6i"%i,"amuse",append_to_file=False)

def run(tend=11. | units.day, timestep=0.5 | units.hour, show_plot=True, storm_track_file="gustav.atcf"):
    code=SwanHurricane(met_file=storm_track_file)

    plot=None
    if show_plot:
        plot=plot_swh(code.nodes, code.elements)
        plot.redraw()

    evolve(code, tend, timestep, plot)

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser(usage=__doc__)
    result.add_option("--track", 
        dest="storm_track_file", default = "gustav.atcf", type=str,
        help="file with storm track [%default]")
    result.add_option("--tend", 
        dest="tend", default = 11. | units.day, unit=units.day, type=float, 
        help="end time of simulation in %unit [%default]")
    result.add_option("--timestep", 
        dest="timestep", default = 1. | units.hour, unit=units.hour, type=float, 
        help="time step between data dumps/ plots in %unit [%default]")
    result.add_option("-p", "--plot",
                  action="store_true", dest="show_plot", default=False,
                  help="show plot")
    return result
        
if __name__=="__main__":
    o, arguments  = new_option_parser().parse_args()
    run(**vars(o))
