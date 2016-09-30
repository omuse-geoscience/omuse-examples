"""

example script to calculate the dynamics of a mid-lattitude square ocean 
basin forced with a constant, lattitude dependend wind, replicating the
(single-gyre) setup of Viebahn 2014 (International Journal of Bifurcation 
and Chaos, Vol. 24, No. 2).

Plots the stream function upon completion.
"""

import numpy
 
from omuse.units import units

from amuse.io import write_set_to_file

from omuse.community.adcirc.interface import Adcirc
from omuse.community.qgmodel.interface import single_gyre_wind_model

from omuse.ext.simple_triangulations import unstructured_square_domain_sets

from matplotlib import pyplot,tri

def viebahn2014_qgmodel(N=50,  
                        H=4000. | units.m, 
                        reynolds_number=1.,
                        dm=0.04,
                        beta0=1.8616e-11 |(units.m * units.s)**-1,
                        L=1.e6 | units.m
                        ):

    rho=1000. | units.kg/units.m**3
    dx=L/N
    A_H=beta0*(dm*L)**3
    tau=reynolds_number*A_H*rho*beta0*H
    U=tau/beta0/rho/H/L
    delta_m=(A_H/beta0)**(1./3)/L
    delta_i=(U/beta0)**0.5/L
    timescale=1./(beta0*L)
  
    print "Viebahn 2014 setup"
    print "N=%i, Reynolds_number=%f"%(N,reynolds_number)
    print "dm (derived):", (A_H/beta0)**(1./3)/L
    print "tau:", tau.value_in(units.Pa)
    print "A:", A_H
    print "timescale:", timescale.in_(units.s)
    print "delta_m:", delta_m
    print "delta_i:", delta_i
  
    nodes,elements,elev_boundary,flow_boundary=unstructured_square_domain_sets(L=L,N=N)
    nodes.depth=H  
    
    code=Adcirc(redirection="none")
  
    code.assign_grid_and_boundary(nodes,elements,elev_boundary, flow_boundary)

    code.parameters.Lx=L
    code.parameters.use_interface_parameters=True
    code.parameters.use_interface_grid=True
    code.parameters.A_H=A_H
    code.parameters.timestep=100. | units.s
    code.parameters.use_predictor_corrector=True
    code.parameters.use_interface_met_forcing=True
    code.parameters.use_ramping=True
  
    print "parameters:"
    print code.parameters
    print 
  
    x=code.forcings.x
    y=code.forcings.y
  
    tau_x,tau_y=single_gyre_wind_model(x,y,L,tau)
    forcings=code.forcings.empty_copy()
    channel=forcings.new_channel_to(code.forcings)
    forcings.coriolis_f=(1.e-4| units.s**-1)+(beta0*(y-L/2))
    forcings.tau_x=tau_x
    forcings.tau_y=tau_y
    channel.copy_attributes(["tau_x","tau_y","coriolis_f"])
  
    return code

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser(usage=__doc__)
    result.add_option("--resolution", 
        dest="N", default = 25, type=int,
        help="grid resolution [%default]")
    result.add_option("--tend", 
        dest="tend", default = 10. | units.day, unit=units.day, type=float, 
        help="end time of simulation in %unit [%default]")
    result.add_option("--depth", 
        dest="H", default = 4000. | units.m, unit=units.m, type=float,
        help="ocean depth in %unit [%default]")
    result.add_option("--boxsize", 
        dest="L", default = 1000. | units.km, unit=units.km, type=float,
        help="box size in %unit [%default]")
    result.add_option("--reynolds_number", 
        dest="reynolds_number", default = 1., type=float,
        help="Reynolds number [%default]")
    return result

if __name__=="__main__":
    o, arguments  = new_option_parser().parse_args()
    
    tend=o.tend
    N=o.N
    H=o.H
    L=o.L
    reynolds_number=o.reynolds_number
    
    sys=viebahn2014_qgmodel(N=N,H=H, L=L, reynolds_number=reynolds_number)  
  
    print "Setup complete, starting evolve."
  
    sys.evolve_model(o.tend)
  
    x=sys.nodes.x
    y=sys.nodes.y
    eta=sys.nodes.eta
    n1=sys.elements.n1
    n2=sys.elements.n2
    n3=sys.elements.n3
  
    elements=numpy.column_stack((n1,n2,n3))
    triangulation=tri.Triangulation(x.value_in(units.km),y.value_in(units.km),elements)
  
    pyplot.clf()
    pyplot.gca().set_aspect('equal')
    pyplot.tripcolor(triangulation,eta.number,shading='gouraud')
    pyplot.colorbar(label=r"$\eta$  [ m ]"  )
    pyplot.xlabel("x (km)")
    pyplot.show()
  
