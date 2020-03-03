import sys
from amuse.units import units
from amuse.community.distributed.interface import DistributedAmuseInterface, DistributedAmuse
from amuse.community.distributed.interface import Resource, Resources, Pilot, Pilots
import atexit


def init_das5_only(username, num_nodes, num_cores):

    print("Setting up distributed code")
    instance = DistributedAmuse()
    instance.parameters.debug = False
    instance.parameters.worker_queue_timeout=1 | units.hour

    instance.parameters.webinterface_port = 4556
    print("url:", instance.get_webinterface_url())
    instance.commit_parameters()

    #print "Resources:"
    resource = Resource()
    resource.name = "DAS-5"
    resource.location = username + "@fs0.das5.cs.vu.nl"
    resource.scheduler_type = "slurm"
    resource.amuse_dir = "/home/" + username + "/amuse/amuse"
    resource.tmp_dir = "/home/" + username + "/tmp"
    instance.resources.add_resource(resource)
    #print instance.resources

    pilot = Pilot()
    pilot.resource_name="DAS-5"
    pilot.queue_name="defq"
    pilot.node_count=num_nodes
    pilot.time= 24|units.hour
    pilot.slots_per_node=num_cores
    pilot.label="DAS-5-Pilot"

    instance.pilots.add_pilot(pilot)

    #~ print "Reservations:"
    #~ print instance.pilots
    print("Waiting for reservations")
    instance.wait_for_pilots()

    return instance

