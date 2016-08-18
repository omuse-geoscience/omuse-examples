import numpy

from omuse.units import units,constants

from matplotlib import pyplot,tri

from amuse.io import write_set_to_file
from amuse.io import read_set_from_file

def moment(f,df,dtheta,ac2,n=0):
    return df*dtheta*(f[:,:]**(n+2)*ac2).sum(axis=-1).sum(axis=-1)/(2*numpy.pi)**(n)

class plot_moments(object):

    def __init__(self,flow=2*numpy.pi*0.0521 | units.rad/units.s,
                  fhigh=2*numpy.pi | units.rad/units.s,msc=32,mdc=36, nodes=None,elements=None):

        if not nodes or not elements:
            nodes=read_set_from_file("nodes.amuse","amuse")
            elements=read_set_from_file("elements.amuse","amuse")
  
        x=nodes.x
        y=nodes.y
        ac2=nodes.ac2
        depth=nodes.depth
        n1=elements.n1
        n2=elements.n2
        n3=elements.n3
  
        elements=numpy.column_stack((n1,n2,n3))
        self.triangulation=tri.Triangulation(x.value_in(units.km),y.value_in(units.km),elements)
          
        flow=flow.value_in(units.rev/units.s)  
        fhigh=fhigh.value_in(units.rev/units.s)  

        thetas=numpy.arange(mdc)*2*numpy.pi/mdc
        fac=(fhigh/flow)**(1./(msc-1))
        self.dlf=numpy.log(fac)
        self.dtheta=2*numpy.pi/mdc
        fs=flow*fac**numpy.arange(msc)
      
        f,theta=numpy.meshgrid(fs,thetas)
                    
        self.f=f | units.rev/units.s
        self.theta=theta | units.rad
        self.ac2=ac2
        self.depth=depth

    def plot(self):
        fig=pyplot.figure(figsize=(12,6))

        
        f1=pyplot.subplot(121)
        f2=pyplot.subplot(122)
        f1.set_aspect('equal')
        f2.set_aspect('equal')
      
        m0=moment(self.f,self.dlf,self.dtheta,self.ac2,0)
        m1=moment(self.f,self.dlf,self.dtheta,self.ac2,1)

        swh=4*(m0)**0.5
        tm01=m0/m1

        val=swh.value_in(units.m)
        vmin=0
        vmax=4
        im1=f1.tripcolor(self.triangulation,val,shading='gouraud',vmin=vmin,vmax=vmax)
        f1.set_xlabel("x (km)")
        f1.set_ylabel("y (km)")
        fig.colorbar(im1,ax=f1,label="Significant wave height Hm0 (m)")
        
        val=tm01.value_in(units.s)
        vmin=1.
        vmax=7.    
        im2=f2.tripcolor(self.triangulation,val,shading='gouraud',vmin=vmin,vmax=vmax)
        f2.set_xlabel("x (km)")
        fig.colorbar(im2,ax=f2, label="mean wave period Tm01 (s)")

        pyplot.show()

if __name__=="__main__":
    plot_moments().plot()
