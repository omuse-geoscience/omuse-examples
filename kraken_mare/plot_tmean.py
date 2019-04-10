import os
import numpy
import cPickle

from mpl_toolkits.basemap import Basemap
from matplotlib import pyplot,tri

from amuse.io import read_set_from_file

from amuse.units import units

from titan_wave import new_option_parser,swan_eq,runlabel

from omuse.ext.dispersion_relations import DeepWaterWaves

from omuse.ext.spherical_geometry import distance, triangle_area

from amuse.ext.radial_profile import radial_density

def moment(f,df,dtheta,ac2,n=0):    
    return df*dtheta*(f[:,:]**(n+2)*ac2).sum(axis=-1).sum(axis=-1)/(2*numpy.pi)**(n)

def significant_wave_height(f,ac2):
    dlf=numpy.log(f[0,1]/f[0,0])
    dtheta=2*units.pi/(f.shape[0]-1)
    return 4*moment(f,dlf,dtheta,ac2,0)**0.5

def mean_period(f,ac2):
    dlf=numpy.log(f[0,1]/f[0,0])
    dtheta=2*units.pi/(f.shape[0]-1)
    return moment(f,dlf,dtheta,ac2,0)/moment(f,dlf,dtheta,ac2,1)

def plot_tmean(plabel,**kwargs):
    label=runlabel(**kwargs)
    if not os.path.isfile("args_"+label+".pkl"):
        swan_eq(**kwargs)

    nodes0=read_set_from_file(kwargs['gridname']+"_nodes","amuse")
    nodes=read_set_from_file(kwargs['gridname']+"_nodes_eq_"+label+".amuse","amuse")
    elements=read_set_from_file(kwargs['gridname']+"_elements","amuse")
    
    flow=kwargs["flow"]
    fhigh=kwargs["fhigh"]    
    mdc,msc=nodes[0].ac2.shape

    thetas=numpy.arange(mdc)*2*numpy.pi/mdc
    fac=(fhigh/flow)**(1./(msc-1))
    fs=flow.value_in(units.rev/units.s)*fac**numpy.arange(msc)
    
    f,theta=numpy.meshgrid(fs,thetas)
            
    fx=f*numpy.cos(theta)
    fy=f*numpy.sin(theta)
    
    f=f | units.rev/units.s

    tmean=mean_period(f,nodes.ac2)

    lon=nodes0.lon.value_in(units.deg)
    lat=nodes0.lat.value_in(units.deg)
          
    n1=elements.n1[:]
    n2=elements.n2[:]
    n3=elements.n3[:]

    elements=numpy.column_stack((n1,n2,n3))
            
    val=tmean.value_in(units.s)        
            
    f=pyplot.figure(figsize=(12,8))

    m = Basemap(projection='stere',lon_0=90.,lat_0=90.,rsphere=495.,
        llcrnrlat=52.,urcrnrlat=79.5,llcrnrlon=52,urcrnrlon=158)
    
    xm,ym=m(lon,lat)    
    triangulation=tri.Triangulation(xm,ym,elements)

    dm=pyplot.tripcolor(triangulation,val,shading='none',vmin=0,vmax=val.max()*1.1,cmap="jet") #Set1
    cbar = m.colorbar(dm,location='right',pad="10%")
    cbar.set_label('mean period (s)')
    
    #~ pyplot.triplot(triangulation,lw=0.8)
    #~ pyplot.xlim(0,500)
    
    m.drawmeridians(numpy.arange(0,360,15),labels=[1,0,0,1],latmax=88.,linewidth=1.5)
    m.drawparallels([55,60,65,70,75,80],labels=[0,1,1,0],linewidth=1.5)

    pyplot.savefig(plabel+"_tmean.png")
    #~ pyplot.show()
    pyplot.clf()
    
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
  
    plot_tmean("test",**kwargs)
