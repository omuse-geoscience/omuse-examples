import numpy
import functools

from omuse.units import units
from omuse.units.quantities import to_quantity

from matplotlib import pyplot

def evensen_spectrum(k, theta, c=1| units.m**2, r=600. | units.km):
    sigma=2*numpy.pi/r
    return c * numpy.exp(-k**2/sigma**2)


class pseudoRandomField(object):
    def __init__(self, N, L, spectrum=None, amplitudes="exp", ascale=1| units.m**2, rscale=600. | units.km):
        if spectrum is None:
            spectrum=functools.partial(evensen_spectrum, c=ascale, r=rscale)
        self.spectrum=spectrum
        self.L=L
        self.N=N
        self.amplitude_distribution=amplitudes
        self.init_model()

    def kgrid(self,L=50. | units.m, N=100):
        k=numpy.fft.fftfreq(N, L/N)
        dk=k[1]
        kx,ky=numpy.meshgrid(k.number,k.number)
        return (2*numpy.pi*kx) | k.unit,(2*numpy.pi*ky) | k.unit, 2*numpy.pi*dk

    def init_model(self):
        kx,ky,dk=self.kgrid(self.L,self.N)
        self._dk=dk
        k=(kx**2+ky**2)**0.5

        kmax=k.max()
        mask=(k.number==0)
        k[mask]=kmax
        kmin=k.min()

        theta=numpy.arctan2(ky.number,kx.value_in(ky.unit))
        theta=theta*(theta>=0)+(2*numpy.pi+theta)*(theta<0)

        self._theta=theta
        self._mask=mask
        self._k=k
        self.calc_spectrum()

        self.amplitudes_phases()
    
    def amplitudes_phases(self):
        phase=numpy.random.uniform(0.,1.,self._k.shape)*2*numpy.pi
        if self.amplitude_distribution=="exp":
          amplitudes=-numpy.log(numpy.random.uniform(0.,1.,self._k.shape))
        elif self.amplitude_distribution=="constant":
          amplitudes=1
        else: raise Exception("amplitudes dist unknown")
        self._amplitudes=amplitudes
        self._phase=phase

    def calc_spectrum(self):
        kspectrum=self.spectrum(self._k,self._theta)
        kspectrum[self._mask]*=0
        self._kspectrum=kspectrum
            
    def field(self):
        amplitude=2**0.5*self._amplitudes*self._kspectrum*self._dk**2          
        phase=self._phase
        # cn=An/2i*exp(i*ph) formula, /2 taken care of taking real part
        # Re(ifft(f))= ifft(f/2+conj(f)/2)
        f=to_quantity(amplitude*numpy.exp(phase*1j)/1j) 
        return self.N**2*numpy.real( numpy.fft.ifft2(f.number)) | f.unit

class pseudoRandomField_freeslip(object):
    def __init__(self, N, L,spectrum=None, amplitudes="exp", ascale=1| units.m**2, rscale=600. | units.km):
        if spectrum is None:
            spectrum=functools.partial(evensen_spectrum, c=ascale, r=rscale)
        self.spectrum=spectrum
        self.L=2*L
        self.N=2*N
        self.amplitude_distribution=amplitudes
        self.init_model()

    def kgrid(self,L=50. | units.m, N=100):
        k=numpy.fft.fftfreq(N, L/N)
        dk=k[1]
        kx,ky=numpy.meshgrid(k.number,k.number)
        return (2*numpy.pi*kx) | k.unit,(2*numpy.pi*ky) | k.unit, 2*numpy.pi*dk

    def init_model(self):
        kx,ky,dk=self.kgrid(self.L,self.N)
        self._dk=dk
        k=(kx**2+ky**2)**0.5

        kmax=k.max()
        mask=(k.number==0)
        k[mask]=kmax
        kmin=k.min()

        theta=numpy.arctan2(ky.number,kx.value_in(ky.unit))
        theta=theta*(theta>=0)+(2*numpy.pi+theta)*(theta<0)

        self._theta=theta
        self._mask=mask
        self._k=k
        self.calc_spectrum()

        self.amplitudes_phases()
    
    def amplitudes_phases(self):
        phase=numpy.random.uniform(0.,1.,self._k.shape)*2*numpy.pi
        if self.amplitude_distribution=="exp":
          amplitudes=-numpy.log(numpy.random.uniform(0.,1.,self._k.shape))
        elif self.amplitude_distribution=="constant":
          amplitudes=1
        else: raise Exception("amplitudes dist unknown")
        self._amplitudes=amplitudes
        self._phase=phase

    def calc_spectrum(self):
        kspectrum=self.spectrum(self._k,self._theta)
        kspectrum[self._mask]*=0
        self._kspectrum=kspectrum
            
    def reinit(self):
        self.amplitudes_phases()

    def field(self):
        amplitude=2**0.5*self._amplitudes*self._kspectrum*self._dk**2          
        phase=self._phase
        # cn=An/2i*exp(i*ph) formula, /2 taken care of taking real part
        # Re(ifft(f))= ifft(f/2+conj(f)/2)
        f=to_quantity(amplitude*numpy.exp(phase*1j)/1j)
        
        N=self.N
        f[N/2+1:,:]=-f[N/2-1:0:-1,:]
        f[:,N/2+1:]=-f[:,N/2-1:0:-1]
        f[N/2+1:,N/2+1:]=f[N/2-1:0:-1,N/2-1:0:-1]
        
        f[N/2,:]=0*f[N/2,:]
        f[:,N/2]=0*f[:,N/2]
        f[0,:]=0*f[0,:]
        f[:,0]=0*f[:,0]
        
        result=2*self.N**2*numpy.real( numpy.fft.ifft2(f.number)) | f.unit
        return result[:N/2,:N/2]
    
    
  
if __name__=="__main__":
        
    # evensen spectrum:
    def spectrum(k, theta, c=1| units.m**2, r=600. | units.km):
        sigma=2*numpy.pi/r
        return c * numpy.exp(-k**2/sigma**2)
  
    L=3500 | units.km
    N=100
    dx=L/N
    L=N*dx
        
    prf=pseudoRandomField_freeslip(N,L,spectrum,amplitudes="constant")

    def expected_variance(k, dk, c=1| units.m**2, r=300. | units.km):
        sigma=2*numpy.pi/r
        return dk**4*c**2*(numpy.exp(-2*k**2/sigma**2)).sum()

    expected=expected_variance(prf._k,prf._dk)
    
    f2=[] | units.none**2
    for i in range(100):
      prf.amplitudes_phases()
      f=prf.field()
      f2.append(f.std()**2)
    print "var, expected:",f2.mean(),expected      

    
    field=prf.field()
    
    #~ print field[-1,:]

    f=pyplot.figure( figsize=(8,6))
    
    vmax=None
    
    l=L.value_in(units.km)
    pyplot.imshow(field.value_in(units.none),vmin=None,vmax=vmax,extent=[0,l,0,l])
    pyplot.axes().set_aspect("equal")
    ax=pyplot.colorbar()
    ax.set_label("field")
    pyplot.xlabel("x (km)")
    pyplot.ylabel("y (km)")
    pyplot.show()


    raise
    
    swh_expected=0.24*u10**2/g
    swh=4*sea.wavefield.std()
    print "wavefield 4*sigma:",swh.in_(units.m)
    
    # check ensemble
    #~ swh2=[] | units.m**2
    #~ for i in range(100):
      #~ sea.amplitudes_phases()
      #~ wavefield=sea.calc_seastate(0. | units.s)
      #~ swh2.append(wavefield.std()**2)
    #~ print 4*(swh2.mean())**0.5/swh_expected, swh2.std()      
    #~ raise

    #~ xx=numpy.arange(N)
    #~ yy=numpy.arange(N)
    #~ x,y=numpy.meshgrid(xx,yy)    
    #~ fig = pyplot.figure()
    #~ ax = fig.gca(projection='3d')
    #~ surface=sea.wavefield.value_in(units.m)
    #~ ax.plot_surface(x,y,surface,cmap="jet", linewidth=0)
    #~ ax.set_zlim(-50,50)
    #~ pyplot.show()
    #~ raise
    
    dt=0.05 | units.s
    
    print "dt,fps",dt,1/dt
    
