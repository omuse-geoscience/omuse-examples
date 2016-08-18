"""
plot swh and action density spectrum as function of fetch
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

def plot_swh_ac(**kwargs):

    L=kwargs["L"]
    N=kwargs["N"]
    flow=kwargs["flow"]
    fhigh=kwargs["fhigh"]
    u10=kwargs["u10"]
    aspect=kwargs["aspect"]
    depth=kwargs["depth"]
    grav=9.81 | units.m/units.s**2


    label=str(abs(hash(str(kwargs))))
 
    directory="./"
    filename=label+".amuse"

    grid=read_set_from_file(os.path.join(directory,filename),"amuse")
        
    ac2=grid.ac2
    
    mxc,myc,mdc,msc=ac2.shape
            
    thetas=numpy.arange(mdc+1)*2*numpy.pi/mdc
    fac=(fhigh/flow)**(1./(msc-1))
    fs=flow.value_in(units.rev/units.s)*fac**numpy.arange(msc)
    
    f,theta=numpy.meshgrid(fs,thetas)

    fx=f*numpy.cos(theta)
    fy=f*numpy.sin(theta)
    
    f=f | units.rev/units.s

    swh=significant_wave_height(f[:-1,:],ac2).value_in(units.m)

    pyplot.ion()
    fig=pyplot.figure(figsize=(16,6))
    pyplot.show()
    
    f1=pyplot.subplot(131)
    f2=pyplot.subplot(132)
    f3=pyplot.subplot(133)
    
    f2.set_aspect('equal')
    f3.set_aspect('equal')

    x=grid.x[:,0].value_in(units.km)
    y=grid.y[0,:].value_in(units.km)
    f=grid.y[0,:]
    
    f1.plot(y,swh.mean(axis=0))
    gu2=grav/u10**2
    f=gu2*grid.y[0,:]

    f1.plot(y,(jonswap_swh(f)/gu2).value_in(units.m),'r')
    f1.plot(y,(wilson_swh(f)/gu2).value_in(units.m),'g')
    
    f1.set_xlabel('fetch (km)')
    f1.set_ylabel('Significant wave height (m)')

        
    ax=f2.imshow(numpy.transpose(swh),origin="lower")
    f2.set_xlabel('pixel')
    f2.set_ylabel('pixel')
    pyplot.draw()


    def plot_ac2(ix,iy):
          f3.cla()
          vmin=0
          vmax=None
          f3.pcolormesh(fx,fy,ac2[ix,iy,:,:].number,vmin=vmin,vmax=vmax)
          f3.set_xlabel('f_x (rev/s)')
          f3.set_ylabel('f_y (rev/s)')
          pyplot.draw()
    
    def onclick(event):
        if event.xdata:
          plot_ac2(int(event.xdata),int(event.ydata))
    
    cid = fig.canvas.mpl_connect('motion_notify_event', onclick)
    
    raw_input()

if __name__=="__main__":
    from wavegrowth import new_option_parser, rectangle
    p=new_option_parser(__doc__)
    o, arguments  = p.parse_args()
    
    kwargs=vars(o)
    label=str(abs(hash(str(kwargs))))
    if not os.path.isfile(label+'.args'):
        rectangle(**kwargs)

    plot_swh_ac(**kwargs)
  

