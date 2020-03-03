#!/usr/bin/env python
import numpy
import os
from matplotlib import pyplot

from amuse.units import units
from omuse.community.pop.interface import POP
from omuse.ext.eddy_tracker.interface import EddyTracker

from pop_highres import init_pop_highres

class POP_Eddy_Tracker(object):

    def __init__(self):

        #application settings, days between tracking and lat/lon bounding box
        self.days_between = 1
        lonmin=0.
        lonmax=50.
        latmin=-45.
        latmax=-20.

        #setup POP
        self.pop = init_pop_highres()

        #start pop
        self.start_time = self.pop.get_model_time()

        #setup eddy tracker
        self.tracker = EddyTracker(grid=p.nodes, domain='Regional',
            lonmin=lonmin, lonmax=lonmax, latmin=latmin, latmax=latmax,
            days_between=self.days_between)

        #start tracker
        self.tracker.find_eddies(self.pop.nodes.ssh, rtime=self.start_time)

    def evolve(self, stop_time, plot=True):
        """evolve until stop_time, track eddies every days_between days, use plot=False to disable drawing plots"""
        tend = self.pop.get_model_time() + (self.days_between | units.day)

        while (tend < stop_time):
            self.pop.evolve_model(tend)

            self.tracker.find_eddies(ssh=self.pop.nodes.ssh, rtime=self.pop.get_model_time())
            if plot:
                self.tracker.plot_eddies(rtime=tend)
            tend = self.pop.get_model_time() + (self.days_between | units.day)

    def stop(self):
        """stop tracking, ensure output is written"""
        self.tracker.stop(tend)
        self.pop.stop()



if __name__=="__main__":
    pet = POP_Eddy_Tracker()

    stop_time = pet.start_time + (60.0 | units.day)

    pet.evolve(stop_time)

    pet.stop()

    print("online eddy tracking successfully completed!")
