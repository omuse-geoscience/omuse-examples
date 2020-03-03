"""
plot swh as function of fetch (and a few theoretical curves)
(calculates data if not present on disk)
"""

import os

import numpy

from matplotlib import pyplot

from amuse.io import read_set_from_file

def moment(f,df,dtheta,ac2,n=0):    
    return df*dtheta*(f**(n+2)*ac2).sum(axis=-1).sum(axis=-1)/(2*numpy.pi)**(n)

from amuse.units import units

def significant_wave_height(f,ac2):
    dlf=numpy.log(f[0,1]/f[0,0])
    dtheta=2*units.pi/f.shape[0]
    return 4*moment(f,dlf,dtheta,ac2,0)**0.5

def jonswap_swh(f):
    return 0.0016*f**0.5

def wilson_swh(f):
    return 0.3*(1-(1+0.004*f**0.5)**-2)

def holthuijssen_2006_swh(f,hinf=0.24, k1=4.14e-4,m1=0.79,p=0.572):
    return hinf*numpy.tanh(k1*f**m1)**p

def tanh_swh(f,hinf=0.24, k1=2.88e-3/0.24,m1=0.45):
    return hinf*numpy.tanh(k1*f**m1)


def plot_swh_ac(**kwargs):

    L=kwargs["L"]
    N=kwargs["N"]
    flow=kwargs["flow"]
    fhigh=kwargs["fhigh"]
    u10=kwargs["u10"]
    aspect=kwargs["aspect"]
    depth=kwargs["depth"]
    grav=9.81| units.m/units.s**2


    label=str(abs(hash(str(kwargs))))
 
    directory="./"
    filename=label+".amuse"

    grid=read_set_from_file(os.path.join(directory,filename),"amuse")
        
    ac2=grid.ac2
    
    mxc,myc,mdc,msc=ac2.shape
    
    thetas=numpy.arange(mdc)*2*numpy.pi/mdc
    fac=(fhigh/flow)**(1./(msc-1))
    fs=flow.value_in(units.rev/units.s)*fac**numpy.arange(msc)
    
    f,theta=numpy.meshgrid(fs,thetas)

    fx=f*numpy.cos(theta)
    fy=f*numpy.sin(theta)
    
    f=f | units.rev/units.s

    swh=significant_wave_height(f,ac2)

    pyplot.ion()
    fig=pyplot.figure(figsize=(8,6))
    ax1=fig.add_axes([0.1,0.1,0.8,0.8])
    pyplot.show()
    
    print("swh:",swh.min().in_(units.m),swh.max().in_(units.m))
    
    x=grid.x[:,0].value_in(units.km)
    y=grid.y[0,:].value_in(units.km)
    f=grid.y[0,:]
    
    swh_=swh.mean(axis=0)
    ax1.plot(y,swh_.value_in(units.m),'k',lw=2.)
    gu2=grav/u10**2
    f=gu2*grid.y[0,:]
    sswh=swh_*gu2

    flim=numpy.array([0,10000.]) 
    whlim=numpy.array([0,0.24])
    xlim=flim/gu2
    ylim=whlim/gu2

    print("scaled fetch:", f.min(),f.max())
    print("scaled waveheight:", sswh.min(),sswh.max())
    ax1.plot(y,(jonswap_swh(f)/gu2).value_in(units.m),'r:',lw=1.5)
    ax1.plot(y,(wilson_swh(f)/gu2).value_in(units.m),'g:',lw=1.5)
    ax1.plot(y,(holthuijssen_2006_swh(f)/gu2).value_in(units.m),'b:',lw=1.5)
    ax1.plot(y,(tanh_swh(f)/gu2).value_in(units.m),'c:',lw=1.5)
    
    ax1.set_xlabel('fetch (km)')
    ax1.set_ylabel('Significant wave height (m)')
    ax1.set_xlim(*xlim.value_in(units.km))
    ax1.set_ylim(*ylim.value_in(units.m))

    ax2a=ax1.twiny()
    ax2=ax2a.twinx()

    ax2a.set_xlabel('(dimensionless)')
    ax2.set_ylabel('(dimensionless)')
    ax2.set_xlim(*flim)
    ax2.set_ylim(*whlim)

    pyplot.savefig("swh.png")
    
    input()

if __name__=="__main__":
    from wavegrowth import new_option_parser, rectangle
    p=new_option_parser(__doc__)
    o, arguments  = p.parse_args()
    
    kwargs=vars(o)
    label=str(abs(hash(str(kwargs))))
    if not os.path.isfile(label+'.args'):
        rectangle(**kwargs)

    plot_swh_ac(**kwargs)
  

