#!/usr/bin/env python
import os

from amuse.units import units
from omuse.community.pop.interface import POP

from distributed_amuse import init_das5_only

def init_pop_highres():

    #machine settings
    popdatadir = '/var/scratch/bwn200/pop/input/'
    username = 'bwn200'

    # distribution files have been provided for the following settings:
    # 16 cores: [19, 23, 28, 37, 56] nodes
    # 24 cores: [15, 19, 25, 37] nodes
    # if you want a different number of cores and/or nodes generate
    # a distribution file yourself using: github.com/nlesc/esalsa-tools
    num_nodes = 37
    num_cores = 16
    num_workers = num_nodes * num_cores

    distributed_amuse = init_das5_only(username, num_nodes, num_cores)
    distributed_amuse.use_for_all_workers()

    p=POP(channel_type="distributed", redirection="none", mode='3600x2400x42',
                number_of_workers=num_workers, max_message_length=1000000)
    p.change_directory(os.getcwd())

    #setup distribution file
    p.parameters.distribution_file = 'distribution_files/dist-60x60-' + str(num_nodes) + '-'+ str(num_cores)

    #set grid info
    p.set_horiz_grid_file(popdatadir+'grid/grid.3600x2400.fob.da')
    p.set_vert_grid_file(popdatadir+'grid/in_depths.42.dat')
    p.set_topography_file(popdatadir+'grid/kmt_pbc.p1_tripole.s2.0-og.20060315.no_caspian_or_black')
    p.set_bottom_cell_file(popdatadir+'grid/dzbc_pbc.p1_tripole.s2.0-og.20060315.no_caspian_or_black')
    p.set_ts_file(popdatadir+'r.t0.1_42l_greenland.01150501')

    #setup forcing files
    p.set_shf_monthly_file(popdatadir+'forcing/shf.NY+H+f.mon')
    p.set_sfwf_monthly_file(popdatadir+'forcing/sfwf.C+r+g8+f.mon')
    p.set_ws_monthly_file(popdatadir+'forcing/ws.o_n_avg.mon')

    return p


if __name__=="__main__":
    #start pop
    p = init_pop_highres()
    p.get_model_time()
    print("POP started in High Resolution mode")
    print("You may interact with the code using the OMUSE interface, the POP object is called p")

    #go interactive
    import readline
    import rlcompleter
    readline.parse_and_bind("tab: complete")
    import code
    code.interact(local=dict(globals(), **locals()) )

