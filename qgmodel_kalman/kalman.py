import numpy

class modelObservation(object):
    def model_observation(self):
        raise Exception("not implemented")
    def observation(self, N):
        obs=self.model_observation(self.model)
        return numpy.random.multivariate_normal(obs, self.covariance_matrix, N).T
    def evolve_model(self,tend):
        self.model.evolve_model(tend)

    #~ def get_state(self, model):
        #~ # return numpy.array([model.x,model.v])
        #~ return numpy.array([model.x, model.v, model.omega])
      
    #~ def set_state(self, model, state):
        #~ model.x=state[0]
        #~ model.v=state[1]
        #~ model.omega=state[2]
        
class EnsembleKalman(object):
    def __init__(self, N, factory):
        self.N=N
        self.models=[factory(i) for i in range(N)]
        
    def evolve_model(self,tend):
        for m in self.models:
          m.evolve_model(tend)
          
    def state_matrix(self):
        return numpy.column_stack([x.get_state() for x in self.models])

    def gain(self,X, H, R):
        Q=numpy.cov(X)
        QHt=Q.dot(H.T)
        return QHt.dot(numpy.linalg.pinv(H.dot(QHt)+R))

    def assimilate(self, observation):
        X=self.state_matrix()
        D=observation.observation(len(self.models))
        H=observation.observation_matrix
        #~ R=numpy.cov(D) # either use observed covariance
        R=observation.covariance_matrix # or covariance matrix
        K=self.gain(X,H,R)
        Xa=X+K.dot(D-H.dot(X))
        for i,m in enumerate(self.models):
          m.set_state(Xa[:,i])
        
    def state(self):
        return numpy.mean(self.state_matrix(), axis=1)

    def state_variance(self):
        return numpy.var(self.state_matrix(), axis=1)

    def state_covariance(self):
        return numpy.cov(self.state_matrix())


# no observation matrix:
class EnsembleKalman_2(EnsembleKalman):
    def gain(self,X, HX, R):
        EX=numpy.mean(X,axis=1)
        A=(X.T-EX).T
        EHX=numpy.mean(HX,axis=1)
        HA=(HX.T-EHX).T
        P=1./(self.N-1)*HA.dot(HA.T)+R
        return 1./(self.N-1)*A.dot(HA.T).dot(numpy.linalg.pinv(P))

    def assimilate(self, observation):
        X=self.state_matrix()
        HX=numpy.column_stack([observation.model_observation(x) for x in self.models])
        D=observation.observation(len(self.models))
        #~ R=numpy.cov(D) # either use observed covariance
        R=observation.covariance_matrix # or covariance matrix
        K=self.gain(X,HX,R)        
        self._K=K
        Xa=X+K.dot(D-HX)
        for i,m in enumerate(self.models):
          m.set_state(Xa[:,i])


# easily invertable R, scales lin in #obs:
class EnsembleKalman_3(EnsembleKalman):
    def gain(self,X, HX, Rinv):
        EX=numpy.mean(X,axis=1)
        A=(X.T-EX).T
        EHX=numpy.mean(HX,axis=1)
        HA=(HX.T-EHX).T
        Q=numpy.identity(self.N)+1./(self.N-1)*(HA.T).dot(Rinv).dot(HA)
        Qinv=numpy.linalg.inv(Q)
        Pinv=Rinv.dot( numpy.identity(len(Rinv)) - 1./(self.N-1) * HA.dot(Qinv).dot(HA.T).dot(Rinv))
        return 1./(self.N-1)*A.dot(HA.T).dot(Pinv)

    def assimilate(self, observation):
        X=self.state_matrix()
        HX=numpy.column_stack([observation.model_observation(x) for x in self.models])
        D=observation.observation(len(self.models))
        #~ R=numpy.cov(D) # either use observed covariance
        R=observation.covariance_matrix # or covariance matrix
        Rinv=numpy.linalg.inv(R) # precompute if not so cheap
        K=self.gain(X,HX,Rinv)        
        Xa=X+K.dot(D-HX)
        for i,m in enumerate(self.models):
          m.set_state(Xa[:,i])

