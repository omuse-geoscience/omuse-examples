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

from omuse.community.qgmodel.interface import QGmodel,single_gyre_wind_model

from matplotlib import pyplot

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
  
    print("Viebahn 2014 setup")
    print("N=%i, Reynolds_number=%f"%(N,reynolds_number))
    print("dm (derived):", (A_H/beta0)**(1./3)/L)
    print("tau:", tau.value_in(units.Pa))
    print("A:", A_H)
    print("timescale:", timescale.in_(units.s))
    print("delta_m:", delta_m)
    print("delta_i:", delta_i)
  
    qg=QGmodel(redirection="none")
  
    qg.parameters.Lx=L
    qg.parameters.Ly=L
    qg.parameters.dx=dx
    qg.parameters.dy=dx
    qg.parameters.dt=900. | units.s
    qg.parameters.A_H=A_H
    qg.parameters.interface_wind=True
    qg.parameters.rho=rho
    qg.parameters.beta0=beta0
    qg.parameters.ocean_depth=H  
    qg.parameters.tau=tau
  
    def wind_function(x,y):
      return single_gyre_wind_model(x,y,L,tau)
    
    qg.set_wind(wind_function)
  
    return qg

def new_option_parser():
    from amuse.units.optparse import OptionParser
    result = OptionParser(usage=__doc__)
    result.add_option("--resolution", 
        dest="N", default = 50, type=int,
        help="grid resolution [%default]")
    result.add_option("--tend", 
        dest="tend", default = 100 | units.day, unit=units.day, type=float, 
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
  
    print("Setup complete, starting evolve.")
  
    sys.evolve_model(o.tend)
  
    psi=(sys.parameters.ocean_depth*sys.grid[:,:,0].psi).value_in(units.Sv)
  
    pyplot.imshow(psi.transpose(), origin="lower", 
                    extent=[0,L.value_in(units.km),0,L.value_in(units.km)])
    pyplot.colorbar(label=r"$\psi \times H$  [ Sv ]"  )
    pyplot.xlabel("x (km)")
    pyplot.show()
  
