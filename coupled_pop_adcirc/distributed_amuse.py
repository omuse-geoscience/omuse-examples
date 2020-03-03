#!/usr/bin/python
import sys
from amuse.units import units
from amuse.community.distributed.interface import DistributedAmuseInterface, DistributedAmuse
from amuse.community.distributed.interface import Resource, Resources, Pilot, Pilots
import atexit


def init_local_only():
  
    print("Setting up distributed code")
    instance = DistributedAmuse()
    instance.initialize_code()

    instance.parameters.webinterface_port = 4556

    print("url:", instance.get_webinterface_url())  
    #~ print "Resources:"
    #~ print instance.resources
      
    #Claim nodes on the resources. In this example simply the "local" machine
    pilot = Pilot()
    pilot.resource_name='local'
    pilot.node_count=1
    pilot.time= 99|units.hour
    pilot.slots_per_node=99
    pilot.node_label='local'
    instance.pilots.add_pilot(pilot)
      
    #~ print "Reservations:"
    #~ print instance.pilots
    
    print("Waiting for reservations")
    instance.wait_for_pilots()
#    atexit.register(instance.stop)
    return instance


class local_only(object):
  def __init__(self):
  
    print("Setting up distributed code")
    instance = DistributedAmuse()
    instance.initialize_code()

#    instance.parameters.debug = True
#    instance.parameters.webinterface_port = 4556
    instance.commit_parameters()
  
    print("url:", instance.get_webinterface_url())
  
    print("Resources:")
    print(instance.resources)
      
    #Claim nodes on the resources. In this example simply the "local" machine
    pilot = Pilot()
    pilot.resource_name='local' # label of the resource to be used
    pilot.node_count=1 # desired number of nodes
    pilot.time= 99|units.hour # wallclock that resource remains available (mainly for job queues)
    pilot.slots_per_node=99 # slots is accounting measure for a job 
    pilot.node_label='local' # label for subgroups of the resource 
    
    instance.pilots.add_pilot(pilot)

    print("Reservations:")
    print(instance.pilots)
    
    print("Waiting for reservations")
    instance.wait_for_pilots()
    self.instance=instance
  def __enter__(self):
    return self.instance
  def __exit__(self,exc_type, exc_val, exc_tb):
    self.instance.stop()
    return False  
      
if __name__=="__main__":
  script = sys.argv[1]
  with local_only():
    exec(compile(open(script, "rb").read(), script, 'exec'))
