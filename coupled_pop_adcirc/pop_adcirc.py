#!/usr/bin/env python

import os
import functools

import time
import threading

import numpy

from omuse.units import units

from matplotlib import pyplot,tri

from omuse.community.pop.interface import POP
from omuse.community.adcirc.interface import Adcirc

from omuse.community.adcirc.read_grid import adcirc_grid_reader, adcirc_parameter_reader
from omuse.community.adcirc.write_grid import get_default_parameter_set

from amuse.datamodel.staggeredgrid import StaggeredGrid
from omuse.ext.grid_remappers import interpolating_2D_remapper,conservative_spherical_remapper

from distributed_amuse import init_local_only

#~ import logging
#~ logging.basicConfig(level=logging.DEBUG)
#~ logging.getLogger("code").setLevel(logging.DEBUG)


#
# tbd:
# - ramping/ relaxation of adcirc 
# - selecting a good subregion for pop interpolation (maybe build this in remapper?)
#
# issues/ questions:
# - how to handle mismatch between grids
# 
#


class POP_Adcirc(object):
    def __init__(self, timestep=None, initialize_adcirc_state=False, adcirc_ramp=2. | units.day,
                 elev_boundary_sponge=False,layer_depth=500. | units.m, distributed=False, boundary_forcing="elevation"):
        self.timestep=timestep
        self.adcirc_ramp=adcirc_ramp
        self.elev_boundary_sponge=elev_boundary_sponge
        self.layer_depth=layer_depth
        self.boundary_forcing=boundary_forcing
        if distributed:
          self.initialize_distributed_amuse()
        self.initialize_pop()
        self.initialize_adcirc()
        if initialize_adcirc_state: 
            self.initialize_adcirc_state()
        self.initialize_channels()
        self.model_time=0. | units.s
        self.pop_time_offset=self.pop.model_time
        self.adcirc_time_offset=self.adcirc.model_time


    def initialize_distributed_amuse(self):
        da = init_local_only()
        da.use_for_all_workers()
        self.distributed_amuse=da

    def initialize_pop(self):
        p=POP(redirection="file", number_of_workers=8,redirect_stdout_file="pop.out")

        cwd=os.getcwd()
        p.change_directory(cwd)
        
        popdatadir="/home/inti/code/amuse/trunk/sandbox/pelupes/pop/"

        #set the grid we want to use
        p.set_horiz_grid_file(popdatadir+'data/input/grid/horiz_grid_20010402.ieeer8')
        p.set_vert_grid_file(popdatadir+'data/input/grid/in_depths.dat')
        p.set_topography_file(popdatadir+'data/input/grid/topography_20010702.ieeei4')
        
        #set the restart file
        p.set_ts_file(popdatadir+'data/input/restart/r.x1_SAMOC_control.00750101')
        
        #setup the forcing
        p.set_shf_monthly_file(popdatadir+'data/input/shf_monthly/shf.normal_year+flux.mon')
        p.set_sfwf_monthly_file(popdatadir+'data/input/sfwf/sfwf_phc0-50_ncarp_r46_flux.mon')
        p.set_ws_monthly_file(popdatadir+'data/input/ws_monthly/ws.1958-2000.mon')
        
        self.pop_grid=p.get_grid()        
        self.pop_forcings_grid=StaggeredGrid(p.elements, p.forcings, p._compute_cell_corners)
        self.pop=p
        
        self.timestep=self.timestep or self.pop.timestep/2
        
    def read_adcirc_grid(self, grid_file):
        gr=adcirc_grid_reader(grid_file,coordinates="spherical")
        gr.read_grid()
        nodes,elements,elev_boundaries,flow_boundaries=gr.get_sets()
        
        if self.boundary_forcing=="flow":
          # change elevation to flow boundary
          for b in elev_boundaries:
            b.type=22
          flow_boundaries=elev_boundaries+flow_boundaries
          elev_boundaries=[]
          
        
        elements.n1=elements.nodes[:,0]
        elements.n2=elements.nodes[:,1]
        elements.n3=elements.nodes[:,2]
        
        i=0
        for e in elev_boundaries+flow_boundaries:
          i+=1
          nodes[e.nodes].vmark=i 
          
        # put maximum on depth to simulate top layer  
        nodes.depth=numpy.minimum(nodes.depth ,self.layer_depth  )

        # fix south east corner of bathymetry to prevent step right on edge
        # this is actually counterproductive... 
        #~ lon=nodes.lon
        #~ lat=nodes.lat
        #~ a=numpy.where( (lat<11.5 |units.deg)*(lon>-61. | units.deg))[0]
        #~ nodes[a].depth=numpy.minimum(nodes[a].depth, 100. | units.m) 
        # maybe increasing depth??

        self.nodes=nodes
        self.elements=elements
        self._elev_boundaries=elev_boundaries
        self._flow_boundaries=flow_boundaries
        self._neta=gr.parameters["NETA"]
        self._nflux=gr.parameters["NFLUX"]


    def initialize_adcirc(self):

        self.read_adcirc_grid("grid.input")

        param=adcirc_parameter_reader("param.input")
        param.read_parameters(NETA=self._neta, NFLUX=self._nflux)
      
        #~ param.parameters['NBFR']=-1
        param.parameters['NWS']=0

        adcirc=Adcirc(coordinates="spherical", redirection="file",redirect_stdout_file="adcirc.out")
        adcirc.set_rootdir(os.getcwd())

        #~ adcirc._parameters=param.parameters
        
        adcirc._parameters=get_default_parameter_set()
        
        adcirc._parameters["NCOR"]=1 # set coriolis f internally  
        
        adcirc.assign_grid_and_boundary(self.nodes, self.elements, self._elev_boundaries, self._flow_boundaries)

        adcirc.parameters.use_interface_parameters=True
        adcirc.parameters.use_interface_grid=True
        if self.boundary_forcing=="elevation":
          adcirc.parameters.use_interface_elevation_boundary=True
        elif self.boundary_forcing=="flow":
          adcirc.parameters.use_interface_flow_boundary=True
        adcirc.parameters.use_interface_met_forcing=True
        
        adcirc.parameters.A_H=1000. | units.m**2/units.s
        #~ adcirc.parameters.A_H=param.parameters["ESLM"] | units.m**2/units.s
        timestep=self.timestep
        n=1
        #~ while timestep>abs(param.parameters["DTDP"]) | units.s: 
        while timestep>300. | units.s: 
          n+=1
          timestep=timestep/n
        adcirc.parameters.timestep=timestep
        #~ adcirc.parameters.bottom_friction_law=["linear","quadratic","hybrid"][param.parameters["NOLIBF"]]
        #~ adcirc.parameters.linear_bottom_friction_coeff=param.parameters["TAU"]| units.s**-1
        adcirc.parameters.bottom_friction_law="hybrid"
        adcirc.parameters.hybrid_bottom_friction_hbreak=50. | units.m
        #~ adcirc.parameters.bottom_friction_law="linear"
        #~ adcirc.parameters.linear_bottom_friction_coeff=0.00001| units.s**-1
        adcirc.parameters.quadratic_bottom_friction_coeff=0.003#param.parameters["CF"]
        adcirc.parameters.use_predictor_corrector=True#param.parameters["DTDP"]<0
        adcirc.parameters.central_longitude=param.parameters["SLAM0"] | units.deg
        adcirc.parameters.central_latitude=param.parameters["SFEA0"] | units.deg
        adcirc.parameters.GWCE_weighting_factor=-1.#param.parameters["TAU0"]
        
        #~ adcirc.parameters.GWCE_weighting_factor=0.005
        
        
        #~ adcirc.parameters.spatial_derivative_advective_term_parameter=0
        #~ adcirc.parameters.time_derivative_advective_term_parameter=0
        #~ adcirc.parameters.finite_amplitude_term_parameter=0
        
        adcirc.parameters.minimum_depth=5.| units.m

        if self.adcirc_ramp:
          adcirc.parameters.use_ramping=True
          adcirc.parameters.ramping_time=self.adcirc_ramp
        print adcirc.parameters

        self.adcirc=adcirc
        
        self.adcirc_grid=adcirc.get_grid()
        self.adcirc_forcings_grid=StaggeredGrid(self.adcirc_grid.elements, nodes=adcirc.forcings)

        elements=adcirc.elements.copy(filter_attributes=lambda x,y: y in ["lat","lon","n1","n2","n3"] )
        nodes=adcirc.nodes.copy(filter_attributes=lambda x,y: y in ["lat","lon"] )

        self.adcirc_memory_grid = StaggeredGrid( elements, nodes=nodes )
        
    def initialize_adcirc_boundary_channel(self):
        # two things remain to be fixed: 
        # - subgrid (as generated by pop nodes slicing) is not a structuredgrid yet
        # so remapping doesn't work directly (proper fix is to fix subgridding)
        # - the definitions of lat and lon are inconsistent between codes, leading
        # to problems with remapping. Various possible fixes: make def consistent,  or
        # define something in the remapper (coordinate trasnform keyword, or is the channel
        # transform (not implemented yet for remapping channel) something tobe used here?? 
        elev_boundaries=list(self.adcirc.elevation_boundaries())
        flow_boundaries=list(self.adcirc.flow_boundaries())
        
        remapper=functools.partial(interpolating_2D_remapper, axes_names=["lon","lat"])
        
        if self.boundary_forcing=="elevation":
          boundary_=elev_boundaries[0]
        elif self.boundary_forcing=="flow":
          boundary_=flow_boundaries[0]
        boundary=boundary_.empty_copy()
        boundary.node=boundary_.node 
        boundary.lat=self.adcirc.nodes[boundary.node].lat
        boundary.lon=self.adcirc.nodes[boundary.node].lon + (2*numpy.pi| units.rad)

        #~ from amuse.io import write_set_to_file
        #~ write_set_to_file(elev_boundaries[0],"elev_boundary","amuse")
        #~ raise

        # correction for grid mismatch (??) :
        boundary.lat+= (0.25 | units.deg)
        boundary.lon+= (0.25 | units.deg)

        # hardcoded subregion
        pop_region=self.pop.nodes[280:310,218:317]
        pop_region_copy=pop_region.copy()
        
        self._adcirc_boundary=boundary
        self.boundary_channel=pop_region_copy.new_remapping_channel_to(boundary, remapper)
        self._pop_boundary_region=pop_region
        self._pop_boundary_region_copy=pop_region_copy
        self._aux_boundary_channel1=pop_region.new_channel_to(self._pop_boundary_region_copy)
        self._aux_boundary_channel2=self._adcirc_boundary.new_channel_to(boundary_)
        
    def update_adcirc_boundary(self):
        if self.boundary_forcing=="elevation":
          self.update_elevation_boundary()
        elif self.boundary_forcing=="flow":
          self.update_flow_boundary()

    def update_elevation_boundary(self):
        self._aux_boundary_channel1.copy_attributes(["ssh"])        
        self.boundary_channel.copy_attributes(["ssh"])

        if self.elev_boundary_sponge:
            n=self.elev_nsponge
            sponge=numpy.arange(n)/(1.*n)
            self._adcirc_boundary[:n].ssh*=sponge
            self._adcirc_boundary[-n:].ssh*=sponge[::-1]

        self._aux_boundary_channel2.transform(["eta"],None,["ssh"])

    def update_flow_boundary(self):
        self._aux_boundary_channel1.copy_attributes(["vx_barotropic","vy_barotropic"])        
        self.boundary_channel.copy_attributes(["depth","vx_barotropic","vy_barotropic"])
        self._aux_boundary_channel2.transform(["flux_x","flux_y"],lambda x,y,z : (x*y, x*z), ["depth","vx_barotropic","vy_barotropic"])

    def initialize_channels(self):
        self.initialize_adcirc_boundary_channel()
        self.forcings_channel = self.pop_forcings_grid.new_remapping_channel_to(self.adcirc_forcings_grid, conservative_spherical_remapper)
        #~ self.memory_channel = self.pop_grid.new_remapping_channel_to(self.adcirc_memory_grid, conservative_spherical_remapper)

    def update_adcirc_forcings(self):
        self.forcings_channel.copy_attributes(["tau_x","tau_y"])
        #~ tau_x=self.adcirc_forcings_grid.nodes.tau_x
        #~ tau_y=self.adcirc_forcings_grid.nodes.tau_y
        #~ coriolis_f=self.adcirc_forcings_grid.nodes.coriolis_f
        #~ print "tau_x (min,max,mean):", tau_x.min(),tau_x.max(),tau_x.mean()
        #~ print "tau_y (min,max,mean):", tau_y.min(),tau_y.max(),tau_y.mean()
        #~ print "coriolis_f (min,max,mean):", coriolis_f.min(),coriolis_f.max(),coriolis_f.mean()


    def initialize_adcirc_state(self):
        channel1 = self.pop_grid.new_remapping_channel_to(self.adcirc_memory_grid, conservative_spherical_remapper)
        channel2 = self.adcirc_memory_grid.nodes.new_channel_to(self.adcirc_grid.nodes)
        channel1.copy_attributes(["vx_barotropic", "vy_barotropic","ssh"])
        
        def f(vx,vy,eta):
            return vx,vy,eta, 0. | units.m/units.s
            
        channel2.transform(["vx","vy","eta","deta_dt"], f, ["vx_barotropic", "vy_barotropic","ssh"])
        

    def evolve_model(self, tend, timestep=None):
        timestep=timestep or self.timestep or tend-self.model_time
        while self.model_time<tend-timestep/2:
            next_timestep=self.pop.timestep_next
            print "update boundary.."
            self.update_adcirc_boundary()
            print "update forcings.."
            self.update_adcirc_forcings()            
            print "evolve pop..."
            self.pop.evolve_model(self.model_time+next_timestep + self.pop_time_offset)
            print "evolve adcirc.."
            self.adcirc.evolve_model(self.model_time+next_timestep + self.adcirc_time_offset)
            print "done"
            self.model_time+=next_timestep
        

class plot_POP_Adcirc(object):
    def __init__(self,pop_grid,adcirc_grid):
        lat = adcirc_grid.nodes.lat.value_in(units.deg)
        lon = adcirc_grid.nodes.lon.value_in(units.deg)
        triangles = numpy.column_stack( ( adcirc_grid.elements.n1, adcirc_grid.elements.n2, 
                        adcirc_grid.elements.n3) )
        triangulation=tri.Triangulation(lon,lat,triangles)
        
        self._adcirc_lat=lat
        self._adcirc_lon=lon
        self._adcirc_triangulation=triangulation

        self.vv1x=numpy.zeros_like(lat)
        self.vv1y=numpy.zeros_like(lat)
        self.f1v=numpy.zeros_like(lat)

        lat = pop_grid.nodes.lat.value_in(units.deg)
        lon = pop_grid.nodes.lon.value_in(units.deg)
        #~ lat = numpy.swapaxes(lat,0,1)
        #~ lon = numpy.swapaxes(lon,0,1)
        #~ self._pop_subregion=slice(218,317),slice(266,320)
        self._pop_subregion=slice(266,320),slice(218,317)
        self._pop_lat=lat[self._pop_subregion]
        self._pop_lon=lon[self._pop_subregion]-360

        self.vv2x=numpy.zeros_like(self._pop_lat)
        self.vv2y=numpy.zeros_like(self._pop_lat)
        self.f2v=numpy.zeros_like(self._pop_lat)

        pyplot.ion()
        #~ pyplot.figure()

        f, (ax1, ax2) = pyplot.subplots(nrows=1, ncols=2, sharex=True, sharey=True,figsize=(16,6))
        ax1.set_adjustable('box-forced')
        ax2.set_adjustable('box-forced')
        self.fig=f
        self.ax1=ax1
        self.ax2=ax2
        self.firstcall=True


        self.update(pop_grid,adcirc_grid)
        pyplot.show()        

    def update(self,pop_grid, adcirc_grid):
        self.vv1x[...]=adcirc_grid.nodes.vx.value_in(units.m / units.s)
        self.vv1y[...]=adcirc_grid.nodes.vy.value_in(units.m / units.s)
        self.f1v[...]=adcirc_grid.nodes.eta.value_in(units.m)
        self.vv2x[...]=pop_grid.nodes.vx_barotropic.value_in(units.m / units.s)[self._pop_subregion]
        self.vv2y[...]=pop_grid.nodes.vy_barotropic.value_in(units.m / units.s)[self._pop_subregion]
        self.f2v[...]=pop_grid.nodes.ssh.value_in(units.m)[self._pop_subregion]

    def redraw(self):
      
        vmin=-1
        vmax=1
              
        if not self.firstcall:
          ax1=self.ax1.axis()
          ax2=self.ax2.axis()
        else:
          ax1=None
          ax2=None

        self.ax1.cla()
        self.ax2.cla()

        self.ax1.tripcolor(self._adcirc_triangulation, self.f1v,shading='flat',vmin=vmin,vmax=vmax)
        #ax1.triplot(self.triangulation, 'k-')
        #~ self.ax1.quiver(self._adcirc_lon, self._adcirc_lat, self.vv1x, self.vv1y, units='xy', width=0.01)

        self.ax2.pcolormesh(self._pop_lon, self._pop_lat, self.f2v, cmap=pyplot.cm.jet,vmin=vmin,vmax=vmax)
        #~ self.ax2.quiver(self._pop_lon, self._pop_lat, self.vv2x, self.vv2y, units='xy')

        #~ self.fig.tight_layout()
        #~ f.savefig("test_remapping.png", dpi=300)
        if ax1:
          self.ax1.axis(ax1)
          self.ax2.axis(ax2)
        self.firstcall=False

        pyplot.draw()

def evolve(code,p=None, timestep=None, outdir="output",tend=30. | units.day,output_step=10, write_3d=False):
    from amuse.io import write_set_to_file

    i=0
    tnow=code.model_time
    timestep=timestep or code.timestep 
    

    print "timestep=",timestep.in_(units.hour)
    t1=time.time()
    write_set_to_file(code.pop_grid.nodes, os.path.join(outdir,"pop_nodes-%6.6i"%i),"amuse",append_to_file=False)
    write_set_to_file(code.pop_forcings_grid.nodes, os.path.join(outdir,"pop_forcing-%6.6i"%i),"amuse",append_to_file=False)
    write_set_to_file(code.pop_grid.elements, os.path.join(outdir,"pop_elements-%6.6i"%i),"amuse",append_to_file=False)
    write_set_to_file(code.adcirc_grid.nodes, os.path.join(outdir,"adcirc_nodes-%6.6i"%i),"amuse",append_to_file=False)
    write_set_to_file(code.adcirc_forcings_grid.nodes, os.path.join(outdir,"adcirc_forcing-%6.6i"%i),"amuse",append_to_file=False)
    write_set_to_file(code.adcirc_grid.elements, os.path.join(outdir,"adcirc_elements-%6.6i"%i),"amuse",append_to_file=False)

    if write_3d:
      write_set_to_file(code.pop.nodes3d[266:320,218:317,:].copy(), os.path.join(outdir,"pop_nodes3d-%6.6i"%i),"amuse",append_to_file=False)
    t2=time.time()
    print "writing took:", t2-t1
    print "starting main loop"

    while tnow<tend-timestep/2: 
        i+=1
        code.evolve_model(tnow+timestep)
        tnow=code.model_time
        print i,tnow.in_(units.day),(code.pop.model_time-code.pop_time_offset).in_(units.day),
        print (code.adcirc.model_time-code.adcirc_time_offset).in_(units.day)
        print i,tnow/timestep,(code.pop.model_time-code.pop_time_offset)/timestep,
        print (code.adcirc.model_time-code.adcirc_time_offset)/timestep
        eta=code.adcirc.nodes.eta
        print "eta (min,max,mean):", eta.min(), eta.max(), eta.mean()
        if p:
            p.update(code.pop_grid,code.adcirc_grid)

        if i%output_step==0:
            write_set_to_file(code.pop_grid.nodes, os.path.join(outdir,"pop_nodes-%6.6i"%i),"amuse",append_to_file=False)
            write_set_to_file(code.pop_forcings_grid.nodes, os.path.join(outdir,"pop_forcing-%6.6i"%i),"amuse",append_to_file=False)
            write_set_to_file(code.pop_grid.elements, os.path.join(outdir,"pop_elements-%6.6i"%i),"amuse",append_to_file=False)
            write_set_to_file(code.adcirc_grid.nodes, os.path.join(outdir,"adcirc_nodes-%6.6i"%i),"amuse",append_to_file=False)
            write_set_to_file(code.adcirc_forcings_grid.nodes, os.path.join(outdir,"adcirc_forcing-%6.6i"%i),"amuse",append_to_file=False)

            if write_3d:
              write_set_to_file(code.pop.nodes3d[266:320,218:317,:].copy(), os.path.join(outdir,"pop_nodes3d-%6.6i"%i),"amuse",append_to_file=False)
            #~ write_set_to_file(code.adcirc_grid.elements, "adcirc_elements-%6.6i"%i,"amuse",append_to_file=False)

def evolve_w_plot(code,plot, **kwargs):    
    plot.redraw()
        
    timer = plot.fig.canvas.new_timer(interval = 200)
    timer.add_callback(plot.redraw)
    timer.start()

    t = threading.Thread(target=evolve, args=(code,plot), kwargs=kwargs)
    t.daemon = True
    t.start()

    raw_input()

def arguments():
    from amuse.units.optparse import OptionParser
    from shutil import copy2
    import cPickle
    
    result = OptionParser()
    result.add_option("--dir", 
        dest="directory", default = "./output", help="Input/output file directory [%default]")
    #~ result.add_option("-i", 
        #~ dest="file_index", default = 0, type="int", help="Input file index [%default]")
    result.add_option("--plot",
        dest="do_plot", default = False, help="show plot [%default]")
    result.add_option("--end", 
        dest="end_time", default = 10 | units.day, unit=units.day, type=float, 
        help="end time of simulation in %unit[%default]")
    result.add_option("--step", 
        dest="output_step", default = 10, type="int", help="steps between output [%default]")
    result.add_option("--force", 
        dest="forcing", default = "elevation", help="forcing used on boundary (elevation or flow) [%default]")

    o,a=result.parse_args()

    #~ label=str(abs(hash(str(vars(o)))))
    label="run"
    with open(os.path.join(o.directory, label+".args"),"w") as f:
      f.write(str(vars(o)))
      #~ cPickle.dump(str(vars(o)),f)

    copy2(__file__, o.directory)

    return o,a

if __name__=="__main__":
      o, arguments  = arguments()
      
      code=POP_Adcirc(layer_depth=500. | units.m, boundary_forcing=o.forcing)
      if o.do_plot:
        plot=plot_POP_Adcirc(code.pop_grid,code.adcirc_grid)
        evolve_w_plot(code, plot, tend=o.end_time,output_step=o.output_step, outdir=o.directory)
      else:
        evolve(code, tend=o.end_time,output_step=o.output_step, outdir=o.directory)
