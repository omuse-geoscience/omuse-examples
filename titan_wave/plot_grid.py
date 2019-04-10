"""
script to plot significant wave height, mean period as a function of fetch, 
as well as spectrum (hover with the mouse over the map of the domain).

"""

import os
import cPickle

import numpy

from matplotlib import pyplot

from amuse.io import read_set_from_file

def moment(f,df,dtheta,ac2,n=0):    
    return df*dtheta*(f[:-1,:]**(n+2)*ac2).sum(axis=-1).sum(axis=-1)/(2*numpy.pi)**(n)

from amuse.units import units

from swan_rectangle import new_option_parser, rectangle, runlabel


def significant_wave_height(f,ac2):
    dlf=numpy.log(f[0,1]/f[0,0])
    dtheta=2*units.pi/(f.shape[0]-1)
    return 4*moment(f,dlf,dtheta,ac2,0)**0.5

def mean_period(f,ac2):
    dlf=numpy.log(f[0,1]/f[0,0])
    dtheta=2*units.pi/(f.shape[0]-1)
    return moment(f,dlf,dtheta,ac2,0)/moment(f,dlf,dtheta,ac2,1)

def jonswap_swh(f):
    return 0.0016*f**0.5

def wilson_swh(f):
    return 0.3*(1-(1+0.004*f**0.5)**-2)

def plot_swh_ac(label,**kwargs):

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

    label=runlabel(**kwargs)
 
    directory="./"
    filename=label+".amuse"

    grid=read_set_from_file(os.path.join(directory,filename),"amuse")
        
    ac2=grid.ac2
    
    mxc,myc,mdc,msc=ac2.shape
    
    #~ print mxc,myc,mdc,msc
        
    thetas=numpy.arange(mdc+1)*2*numpy.pi/mdc
    fac=(fhigh/flow)**(1./(msc-1))
    fs=flow.value_in(units.rev/units.s)*fac**numpy.arange(msc)
    
    f,theta=numpy.meshgrid(fs,thetas)
        
    #~ print theta.shape,f.shape,ac2[0,0,:,:].shape
    
    fx=f*numpy.cos(theta)
    fy=f*numpy.sin(theta)
    
    f=f | units.rev/units.s

    swh=significant_wave_height(f,ac2).value_in(units.m)
    tmean=mean_period(f,ac2).value_in(units.s)

    pyplot.ion()
    fig=pyplot.figure(figsize=(16,6))
    pyplot.show()
    
    f1=pyplot.subplot(131)
    f2=pyplot.subplot(132)
    f3=pyplot.subplot(133)
    
    f2.set_aspect('equal')
    f3.set_aspect('equal')
    
    print "swh (min,max):",swh.min(),swh.max()
    print "tmean (min,max):",tmean.min(),tmean.max()
    
    x=grid.x[:,0].value_in(units.km)
    y=grid.y[0,:].value_in(units.km)
    f=grid.y[0,:]
    
    f1.plot(y,swh.mean(axis=0))
    gu2=grav/u10**2
    f=gu2*grid.y[0,:]

    print "scaled fetch:", f.min(),f.max()
    f1.plot(y,(jonswap_swh(f)/gu2).value_in(units.m),'r')
    f1.plot(y,(wilson_swh(f)/gu2).value_in(units.m),'g')
    
    f1.set_xlabel('fetch (km)')
    f1.set_ylabel('Significant wave height (m)')
    
    ax2=f1.twinx()
    ax2.plot(y,tmean.mean(axis=0), 'b:')
    ax2.set_ylabel('mean wave period (s)')
    
    ax=f2.imshow(numpy.transpose(swh),origin="lower",vmin=0,vmax=1.)
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
        #~ print 'button=%d, x=%d, y=%d, xdata=%f, ydata=%f'%(
            #~ event.button, event.x, event.y, event.xdata, event.ydata)
        if event.xdata:
          plot_ac2(int(event.xdata),int(event.ydata))
    
    cid = fig.canvas.mpl_connect('motion_notify_event', onclick)#button_press_event
    
    raw_input()

if __name__=="__main__":
    p=new_option_parser("usage: %prog <argfile> [options]")
    o, arguments  = p.parse_args()
    
    if arguments:
      argfile=arguments[0]
      if not os.path.isfile(argfile):
        raise Exception( argfile+" does not exists")
      f=open(argfile,"r")
      kwargs=cPickle.load(f)
      f.close()
    else:
      kwargs=vars(o)

    label=runlabel(**kwargs)

    if not os.path.isfile(label+'.args'):
        rectangle(**kwargs)

    plot_swh_ac(label,**kwargs)
  

