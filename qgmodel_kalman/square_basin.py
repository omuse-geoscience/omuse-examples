"""

example script to calculate the dynamics of a mid-lattitude square ocean 
basin forced with a constant, lattitude dependend wind, replicating the
(single-gyre) setup of Viebahn 2014 (International Journal of Bifurcation 
and Chaos, Vol. 24, No. 2).

Plots the stream function upon completion.
"""

import sys

import numpy
 
from omuse.units import units

from amuse.io import write_set_to_file

from omuse.community.qgmodel.interface import QGmodel,single_gyre_wind_model, dijkstra_wind_model, jans_wind_model

from matplotlib import pyplot

w_dict=dict(single=single_gyre_wind_model, dijkstra=dijkstra_wind_model,jan=jans_wind_model)

def viebahn2014_qgmodel(N=50,  
                        H=4000. | units.m, 
                        A_H=100 | units.m**2/units.s,
                        beta0=1.8616e-11 |(units.m * units.s)**-1,
                        L=1.e6 | units.m,
                        wind_model=single_gyre_wind_model,
                        tau=0.04 | units.Pa,
                        code=QGmodel,
                        timestep=900 | units.s,
                        timestep_method="leapfrog",
                        verbose=False
                        ):

    rho=1000. | units.kg/units.m**3
    dx=L/N
    dm=(A_H/beta0)**(1./3)/L
    reynolds_number=tau/(A_H*rho*beta0*H)

    U=tau/beta0/rho/H/L
    delta_m=(A_H/beta0)**(1./3)/L
    delta_i=(U/beta0)**0.5/L
    timescale=1./(beta0*L)
  
  
    if verbose:
      print("Viebahn 2014 setup")
      print("N=%i, Reynolds_number=%f"%(N,reynolds_number))
      print("dm (derived):", (A_H/beta0)**(1./3)/L)
      print("tau:", tau.value_in(units.Pa))
      print("A:", A_H)
      print("timescale:", timescale.in_(units.s))
      print("delta_m:", delta_m)
      print("delta_i:", delta_i)
  
    qg=code(redirection="none")
  
    qg.parameters.Lx=L
    qg.parameters.Ly=L
    qg.parameters.dx=dx
    qg.parameters.dy=dx
    qg.parameters.dt=timestep
    qg.parameters.A_H=A_H
    qg.parameters.interface_wind=True
    qg.parameters.rho=rho
    qg.parameters.beta0=beta0
    qg.parameters.ocean_depth=H  
    qg.parameters.tau=tau
    qg.parameters.timestep_method=timestep_method
    
    def wind_function(x,y):
      return wind_model(x,y,L,tau)
    
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
    result.add_option("--dt", 
        dest="dtplot", default = 10 | units.day, unit=units.day, type=float, 
        help="plot time step in %unit [%default]")
    result.add_option("--depth", 
        dest="H", default = 4000. | units.m, unit=units.m, type=float,
        help="ocean depth in %unit [%default]")
    result.add_option("--boxsize", 
        dest="L", default = 1000. | units.km, unit=units.km, type=float,
        help="box size in %unit [%default]")
    #~ result.add_option("--reynolds_number", 
        #~ dest="reynolds_number", default = 1., type=float,
        #~ help="Reynolds number [%default]")
    result.add_option("--wind_model", 
        dest="wind_model", default = "single", type=str,
        help="wind model: single, dijkstra or jan [%default]")
    result.add_option("--viscosity", 
        dest="A", default = 1000. | units.m**2/units.s, unit=units.m**2/units.s, type=float,
        help="lateral viscosity %unit [%default]")
    result.add_option("--wind_stress", 
        dest="tau", default = 0.04 | units.Pa, unit=units.Pa, type=float,
        help="wind_stress %unit [%default]")
    result.add_option("--dry_run", 
        dest="dry_run", default = False, action="store_true",
        help="dry run [%default]")


    return result

if __name__=="__main__":
    o, arguments  = new_option_parser().parse_args()
    
    tend=o.tend
    N=o.N
    H=o.H
    L=o.L
    tau=o.tau
    wind_model=w_dict[o.wind_model]
    A_H=o.A
    
    model=viebahn2014_qgmodel(N=N,H=H, L=L,  
        wind_model=wind_model, tau=tau, A_H=A_H)  
  
    if o.dry_run:
        print("stopping because of --dry_run")
        sys.exit(0)
  
    print("Setup complete, starting evolve.")

    pyplot.ion()
    
    tnow=0.*o.tend
    dt=o.dtplot
    
    f=pyplot.figure(figsize=(8,6))
  
    pyplot.show()
    while tnow<tend-dt/2:
  
      model.evolve_model(tnow+dt)
      tnow=model.model_time
      print("reached: ", tnow.in_(units.day), tnow/tend)
  
      psi=(model.parameters.ocean_depth*model.grid[:,:,0].psi).value_in(units.Sv)
    
      f.clf()
      pyplot.imshow(psi.transpose(), origin="lower", 
                      extent=[0,L.value_in(units.km),0,L.value_in(units.km)],vmin=-50,vmax=50.,
                      cmap="bwr")
      pyplot.colorbar(label=r"$\psi \times H$  [ Sv ]"  )
      pyplot.xlabel("x (km)")
      pyplot.draw()
      pyplot.pause(0.01)

    print("done")
    
    print(psi.shape)
    print(psi[0,:10])
    print(psi[-1,:10])
    print(psi[:10,0])
    print(psi[:10,-1])
    
    input()
