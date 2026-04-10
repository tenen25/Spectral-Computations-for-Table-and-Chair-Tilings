import math
import cmath
import numpy as np
from Sub_Rule import rule_defn_read,  iterate_sub
from Aux_bands import save_bands, read_bands, print_several_bands
from Aux_matrix import op_mat_no_phase, sample_band_edges
#MAIN FUNCTION

# --- Potential values ---
POT_SHIFT = 0   # shifts spectra of the operators
SCALE     = 1   # scaling factor for the potentials
val = SCALE * np.array([
     POT_SHIFT,                        # letter 'a'
     POT_SHIFT + 9,                    # letter 'b'
     POT_SHIFT + 18,                   # letter 'c'
     POT_SHIFT + 27,                   # letter 'd'
 ])

# --- Tiling configuration ---
SUB_DEG   = 2   #Stretch factor of the substitution
RULE_NAME = "Table_Rule" #"Chair_Rule" 
tile_rule = rule_defn_read(RULE_NAME)

start_tile_1    = [[0, 2], [2, 0]]  # 1st table-legal periodic config
start_tile_2    = [[1, 3], [3, 1]]  # 2nd table-legal periodic config
start_tile_3    = [[1,2,3], [0,1,2]] # chair-legal periodic config
start_tile_ileg = [[0, 0], [0, 0]]  # illegal config

start_tile      = start_tile_ileg    # select starting patch here

# --- Computation parameters ---
max_itera   = 4     #number of maximal substitution steps
res         = 9     #sampling resolution in the Floquet torus
init_itera  = 0     #which iterations to plot range start
fin_itera   = 4     #which iterations to plot range end

init_size = len(start_tile) * len(start_tile[0])
wid  = init_size * SUB_DEG  ** max_itera
size = wid ** 2


# --- Run ---
if __name__ == "__main__":
       
    # Workflow options — uncomment as needed:
    """ Next line generates and saves spectrum of itera-th iteration """
    for it in range(0, max_itera+1):
        if it in range(0,5):
            resol=res
        else:
            resol=2                
        patch = iterate_sub(it, SUB_DEG, start_tile, tile_rule)
        save_bands(it, resol, patch, val)
    """ Next line prints to a plot spectra of iterations """
    print_several_bands(
        init_itera,                 #Beginning of iterations
        fin_itera,                  #End of iterations
        val,                        #Potential values
        init_size,                  #The initial size of the starting patch
        vert_dist=0.02,             #Visual distance between iterations in plot
        plot_tiny=False,            #Boolean if small bands are plotted or not
        tiny_tol=1e-6,              #A band is treated as tiny if right - left < tiny_tol.
        tiny_height_factor=0.3,     #Vertical length for the height of the tiny bands
        color='black',              #Color spectral bands
        linewidth=4.0,              #Height of the lines for spectral bands
        fig_width=5.0,
        save_path="table_spectrum_ileg_Pot=0,9,18,27_res=9_tiny=False_HeightFacr=0.02_Linewidth=4_iter=6.pdf"
    )
    
    print("Finished computation.")

