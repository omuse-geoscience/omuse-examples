import numpy

from matplotlib import pyplot

from omuse.units import units

from omuse.community.qgmodel.interface import QGmodel, jans_wind_model

from square_basin import viebahn2014_qgmodel
from pseudo_random import pseudoRandomField_freeslip

from kalman import EnsembleKalman_2

def test_noisemodel():
    def evensen_spectrum(k, theta, c=1| units.m**2, r=600. | units.km):
        sigma=2*numpy.pi/r
        return c * numpy.exp(-k**2/sigma**2)

    noisemodel=pseudoRandomField_freeslip(N=20,L=3500. | units.km, 
        spectrum=functools.partial(evensen_spectrum, c=1. | units.Sv*units.km**2, r=600. | units.km),
        amplitudes="constant")

    field=noisemodel.field()
    
    print field.std().in_(units.Sv)
    print ((2*2*noisemodel._dk**4*noisemodel._kspectrum**2).sum()**0.5).in_(units.Sv)

def test(N=10):
  
    class modelObservation(object):
        def __init__(self, model):
            self.model=model
            self.covariance_matrix=1.*numpy.eye(5)
            self.observation_points=[[500,500],[500,3000],[3000,500],[3000,3000],[1500,1500]] | units.km
        def model_observation(self, model):
            savepoints=model.grid.samplePoints(x=self.observation_points[:,0],
                                               y=self.observation_points[:,1])
            data=self.model.parameters.ocean_depth*savepoints.psi[:,0]
            return data.value_in(units.Sv)
        def observation(self, N):
            obs=self.model_observation(self.model)
            return numpy.random.multivariate_normal(obs, self.covariance_matrix, N).T
        def evolve_model(self,tend):
            self.model.evolve_model(tend+(20. | units.day))
        
    class member(QGmodel):
        def get_state(self):
            psi=self.parameters.ocean_depth*self.grid[:,:,0].psi
            return psi.value_in(units.Sv).flatten()
        def set_state(self, state):
            psi=state.reshape(self.grid.shape) | units.Sv
            self.grid.psi=psi/self.parameters.ocean_depth
      
    noisemodel=pseudoRandomField_freeslip(N=50,L=3500. | units.km, 
        amplitudes="constant", ascale=100000. | units.Sv*units.km**2, rscale=1200. | units.km )

    def factory(i):
        code= viebahn2014_qgmodel(N=50,  
                        H=4000. | units.m, 
                        A_H=10000 | units.m**2/units.s,
                        beta0=1.8616e-11 |(units.m * units.s)**-1,
                        L=3500. | units.km,
                        wind_model=jans_wind_model,
                        tau=0.04 | units.Pa,
                        code=member,
                        timestep_method="euler",
                        timestep=900. | units.s
                        )
        noisemodel.reinit()
        code.grid[:-1,:-1,0].psi=noisemodel.field()/code.parameters.ocean_depth
        return code
        
    kalman=EnsembleKalman_2(N, factory)

    reference=viebahn2014_qgmodel(N=50,  
                        H=4000. | units.m, 
                        A_H=10000 | units.m**2/units.s,
                        beta0=1.8616e-11 |(units.m * units.s)**-1,
                        L=3500. | units.km,
                        wind_model=jans_wind_model,
                        tau=0.04 | units.Pa,
                        code=QGmodel,
                        timestep_method="leapfrog",
                        timestep=900. | units.s
                        )

    observation=modelObservation(reference)

    print "evolve reference model"

    observation.evolve_model(0. | units.day)

    print "so far..."

    shape=[51,51]

    tnow=0. | units.day
    dt=1. | units.day
    tend=100 | units.day
    
    pyplot.ion()
    f, ax=pyplot.subplots(ncols=3, figsize=(12,4))
    pyplot.show()
    
    print kalman.state_variance().max()**0.5
    ax[0].cla()
    psi=reference.parameters.ocean_depth*reference.grid[:,:,0].psi
    psi=psi.value_in(units.Sv)
    vmin=psi.min()
    vmax=psi.max()
    #~ vmax=max(abs(vmin),abs(vmax))
    #~ vmin=-vmax
    
    ax[0].imshow(psi.T, origin="lower", vmin=vmin,vmax=vmax)
    
    ax[1].cla()
    state=kalman.state().reshape(shape)
    print "state:", state.min(), state.max()
    ax[1].imshow(state.T, origin="lower")#, vmin=vmin,vmax=vmax)


    ax[2].cla()
    state=kalman.state_variance().reshape(shape)**0.5
    ax[2].imshow(state.T, origin="lower")#, vmin=vmin/5,vmax=vmax/5)


    pyplot.draw()
    f.canvas.flush_events()         

    
    while tnow<tend-dt/2:
        observation.evolve_model(tnow)
        kalman.evolve_model(tnow)
        tnow+=dt
        
        kalman.assimilate(observation)
    
        print "max variance, diff:",kalman.state_variance().max()**0.5,
        ax[0].cla()
        psi=reference.parameters.ocean_depth*reference.grid[:,:,0].psi
        psi=psi.value_in(units.Sv)
        vmin=psi.min()
        vmax=psi.max()
        #~ vmax=max(abs(vmin),abs(vmax))
        #~ vmin=-vmax
        #~ print "ref:", vmin,vmax
        
        ax[0].imshow(psi.T, origin="lower", vmin=vmin,vmax=vmax)
        ax[0].set_title("reference")
        
        ax[1].cla()
        state=kalman.state().reshape(shape)
        print abs(psi-state).mean()
        ax[1].imshow(state.T, origin="lower")#, vmin=vmin,vmax=vmax)
        ax[1].set_title("ensemble mean")


        ax[2].cla()
        state=kalman.state_variance().reshape(shape)**0.5
        ax[2].imshow(state.T, origin="lower")#, vmin=vmin/5,vmax=vmax/5)
        ax[2].set_title("ensemble variance **0.5")


        pyplot.draw()
        f.canvas.flush_events()         

def test_response(N=10):
  
    class modelObservation(object):
        def __init__(self, model):
            self.model=model
            self.covariance_matrix=0.1*numpy.eye(5)
            self.observation_points=[[500,500],[500,3000],[3000,500],[3000,3000],[1500,1500]] | units.km
        def model_observation(self, model):
            savepoints=model.grid.samplePoints(x=self.observation_points[:,0],
                                               y=self.observation_points[:,1])
            data=self.model.parameters.ocean_depth*savepoints.psi[:,0]
            return data.value_in(units.Sv)
        def observation(self, N):
            obs=self.model_observation(self.model)
            return numpy.random.multivariate_normal(obs, self.covariance_matrix, N).T
        def evolve_model(self,tend):
            self.model.evolve_model(tend+(20. | units.day))
        
    class member(QGmodel):
        def get_state(self):
            psi=self.parameters.ocean_depth*self.grid[:,:,0].psi
            return psi.value_in(units.Sv).flatten()
        def set_state(self, state):
            psi=state.reshape(self.grid.shape) | units.Sv
            self.grid.psi=psi/self.parameters.ocean_depth
      
    noisemodel=pseudoRandomField_freeslip(N=50,L=3500. | units.km, 
        amplitudes="constant", ascale=100000. | units.Sv*units.km**2, rscale=1200. | units.km )

    def factory(i):
        code= viebahn2014_qgmodel(N=50,  
                        H=4000. | units.m, 
                        A_H=10000 | units.m**2/units.s,
                        beta0=1.8616e-11 |(units.m * units.s)**-1,
                        L=3500. | units.km,
                        wind_model=jans_wind_model,
                        tau=0.04 | units.Pa,
                        code=member,
                        timestep_method="euler",
                        timestep=900. | units.s
                        )
        noisemodel.reinit()
        code.grid[:-1,:-1,0].psi=noisemodel.field()/code.parameters.ocean_depth
        return code
        
    kalman=EnsembleKalman_2(N, factory)

    reference=viebahn2014_qgmodel(N=50,  
                        H=4000. | units.m, 
                        A_H=10000 | units.m**2/units.s,
                        beta0=1.8616e-11 |(units.m * units.s)**-1,
                        L=3500. | units.km,
                        wind_model=jans_wind_model,
                        tau=0.04 | units.Pa,
                        code=QGmodel,
                        timestep_method="leapfrog",
                        timestep=900. | units.s
                        )

    observation=modelObservation(reference)

    print "evolve reference model"

    observation.evolve_model(0. | units.day)

    print "so far..."

    shape=[51,51]

    tnow=0. | units.day
    dt=1. | units.day
    tend=100 | units.day
    
    pyplot.ion()
    f, ax=pyplot.subplots(ncols=3, figsize=(12,4))
    pyplot.show()
    
    print kalman.state_variance().max()**0.5
    ax[0].cla()
    psi=reference.parameters.ocean_depth*reference.grid[:,:,0].psi
    psi=psi.value_in(units.Sv)
    vmin=psi.min()
    vmax=psi.max()
    #~ vmax=max(abs(vmin),abs(vmax))
    #~ vmin=-vmax
    
    ax[0].imshow(psi.T, origin="lower", vmin=vmin,vmax=vmax)
    
    ax[1].cla()
    state=kalman.state().reshape(shape)
    print "state:", state.min(), state.max()
    ax[1].imshow(state.T, origin="lower")#, vmin=vmin,vmax=vmax)


    ax[2].cla()
    #~ state=kalman.state_variance().reshape(shape)**0.5
    #~ ax[2].imshow(state.T, origin="lower")#, vmin=vmin/5,vmax=vmax/5)


    pyplot.draw()
    f.canvas.flush_events()         

    
    while tnow<tend-dt/2:
        observation.evolve_model(tnow)
        kalman.evolve_model(tnow)
        tnow+=dt
        
        kalman.assimilate(observation)
    
        print "max variance, diff:",kalman.state_variance().max()**0.5,
        ax[0].cla()
        psi=reference.parameters.ocean_depth*reference.grid[:,:,0].psi
        psi=psi.value_in(units.Sv)
        vmin=psi.min()
        vmax=psi.max()
        #~ vmax=max(abs(vmin),abs(vmax))
        #~ vmin=-vmax
        #~ print "ref:", vmin,vmax
        
        ax[0].imshow(psi.T, origin="lower", vmin=vmin,vmax=vmax)
        ax[0].set_title("reference")
        
        ax[1].cla()
        state=kalman.state().reshape(shape)
        print abs(psi-state).mean()
        ax[1].imshow(state.T, origin="lower")#, vmin=vmin,vmax=vmax)
        ax[1].set_title("ensemble mean")


        ax[2].cla()
        state=kalman._K[:,0].reshape(shape)
        ax[2].imshow(state.T, origin="lower")#, vmin=vmin/5,vmax=vmax/5)
        ax[2].set_title("response 1")


        pyplot.draw()
        f.canvas.flush_events()         


if __name__=="__main__":
    #~ test_noisemodel()
    test_response(50)
