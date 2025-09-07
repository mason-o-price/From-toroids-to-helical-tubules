# Import libraries
import numpy as np
import pandas as pd
import json
import random

def main(interactions_file: str, angles_file: str, output_folder: str, simulation_name: str, json_folder: str, eb_scan: list, kb_scan: list, mu: float, stretching_mod: float, l0: float, timesteps: float, num_tests: int, max_N: int, k_insertion: float, k_fusion2: float):
    # Read the interaction matrix of the toroid we want to generate
    int_matrix = pd.read_csv(interactions_file, header = None)
    interactions = int_matrix.to_numpy()

    # read the binding angles (in degrees)
    df_theta = pd.read_csv(angles_file, header = None)
    df_theta = df_theta*np.pi/180
    np_theta = df_theta.to_numpy()

    #generate a list of the subunits and their edge types
    subunits = []
    for i, e in enumerate(int_matrix.columns[::3]):
        subunits.append([i+1, e+1, e+2, e+3])
    subunits = np.array(subunits)

    # read a template .json input file to use most of the same parameters
    templatefile = 'template_input.json'
    with open(templatefile) as tf:
        template = json.load(tf)
    tf.close()

    #write the json files
    for e_b in eb_scan:
        for k_b in kb_scan:
            for n in range(num_tests):
                new_json = template #start with the template

                # add the type parameters between edges
                num_types = len(interactions)
                for i in range(num_types):
                    if 1 in list(interactions[:,i]): # check that there is an interaction here
                        binding_edge = list(interactions[:,i]).index(1)+1 # find the entry 1 in the interaction matrix
                        angle = np_theta[(binding_edge-1)//3, i//3] #find the binding angle of this interaction
                        new_json['parameters']['edges']['type'+str(i+1)] = {'bending_modulus': {'type'+str(binding_edge): k_b},
                                                                'e_b': {'type'+str(binding_edge): e_b},
                                                                'stretch_mod': stretching_mod,
                                                                'theta0': {'type'+str(binding_edge): angle},
                                                                'l0': l0}
                    else: # if there is no interaction here (i.e. this egde is passivated), 
                            # then make it self-complementary with a very strong repulsive force
                        binding_edge = i+1
                        angle = 0
                        new_json['parameters']['edges']['type'+str(i+1)] = {'bending_modulus': {'type'+str(binding_edge): k_b},
                                                                'e_b': {'type'+str(binding_edge): 10000},
                                                                'stretch_mod': stretching_mod,
                                                                'theta0': {'type'+str(binding_edge): angle},
                                                                'l0': l0}
                        
                # define the subunit parameters
                for i in range(len(subunits)):
                    new_json['parameters']['subunits']['sub'+str(i+1)] = {'edges': {'edge1': 'type'+str(subunits[i][1]),
                                                                                    'edge2': 'type'+str(subunits[i][2]), 
                                                                                    'edge3': 'type'+str(subunits[i][3])},
                                                                                    'k_insertion': k_insertion,
                                                                                    'mu': mu,
                                                                                    'number_of_rotational_configs': 3,
                                                                                    'R_exc': 0.15}
                # change some other parameters (R_add, d_max, l_fuse) depending on the values of the bending/stretching modulus.
                alp=2.0
                kT=1.0
                l_thermal_stretch = min( 2.0*(alp*kT/stretching_mod)**0.5, 0.5*l0 )
                l_thermal_bend = min( (3.0*l0**2*alp*kT/(2.0*k_b))**0.5, 0.5*l0 )
                #l_thermal = max(2.0*(alp*kT/eps_s)**0.5, (3.0*l0**2*alp*kT/(2.0*kappa_b))**0.5)
                l_thermal = max( l_thermal_stretch, l_thermal_bend  )
                new_json['parameters']['mesh_parameters']['R_add']=l_thermal #/3.5 /5.0
                new_json['parameters']['mesh_parameters']['d_max']=l_thermal / 3.5 /1.41
                new_json['parameters']['mesh_parameters']['l_fuse']=l_thermal
                new_json['parameters']['mesh_parameters']['k_fusion2']=k_fusion2 # If still doesn't work, try 1e-4. 
                ensemble = random.randint(100, 50000000) # Find a random value for the ensemble value
                species = str(random.randint(1, 12)) # Randomize the starting subunit
                new_json['parameters']['run_parameters']['init_file'] = 'subunits/subunit_' + species + '.om'
                new_json['parameters']['run_parameters']['ensemble'] = ensemble
                new_json['parameters']['run_parameters']['timesteps'] = timesteps
                new_json['parameters']['run_parameters']['convert_to_lammps'] = False # usually: False
                new_json['parameters']['run_parameters']['convert_to_lammps_trajectory'] = False # usually: False
                new_json['parameters']['run_parameters']['adaptive_rates'] = False # usually: False
                new_json['parameters']['run_parameters']['convert_to_vtk'] = False # usually: False
                new_json['parameters']['run_parameters']['convert_to_dat'] = False # usually: False
                new_json['parameters']['run_parameters']['Nmax'] = max_N
                new_json['parameters']['run_parameters']['dtsave'] = int(1e6)
                new_json['parameters']['run_parameters']['dtsave_snapshot'] = int(1e6)
                new_json['parameters']['run_parameters']['checkpoint_freq'] = int(1e6)
                new_json['parameters']['run_parameters']['data_folder'] = output_folder + "output_" + simulation_name + "__eb"+str(e_b).replace('.', '_')+"_kb"+str(k_b).replace('.', '_')+ "_ensemble"+str(ensemble)+"/"
                with open(json_folder + "input_" + simulation_name + "_eb"+str(e_b)+"_kb"+str(k_b)+ "_ensemble"+str(ensemble)+".json", "w") as outfile:
                    outfile.write(json.dumps(new_json, indent=4))

if __name__ == "__main__":
    # Define the file locations for the interaction matrix and the binding angle matrix
    interactions_file = "interactionMatrices/(4,2,2,1)_interactionMatrix.csv"
    angles_file = "bindingAngles/(4,2,2,1)_bindingAngles.csv"

    # Define the data folder where the simulations will output the results
    output_folder = "/home/masonprice/toroids_new_tests/data_toroid_mu4_5/"
    simulation_name = "4221"

    # Define the folder where the jsons will be dumped
    json_folder = "new_jsons/jsons_toroid_mu4_5/"

    # Define the array of binding energy and bending modulus values
    eb_scan = [-5, -6, -7, -8, -9, -10, -11]
    kb_scan = [1e-3,  1e-2, 1e-1, 1e0, 1e1, 1e2, 1e3]

    # Set the simulation conditions
    mu = -4.5 # Chemical potential (usually: -4.5)
    stretching_mod = 100 # stretching modulus (usually: 100)
    l0 = 1.0 # preferred edge length (usually: 1.0)
    timesteps = 1e9 # maximum number of timesteps (usually: 1e8)
    num_tests = 50 # number of tests per data point (usually: 50)
    max_N = 200 # maximum number of triangles (usually: 200)
    k_insertion = 0.0001 # rate of insertion (was previously: 0.005. new value: 0.0001) --> must preserve detailed balance (i.e., P < 1)
    k_fusion2 = 0.0005 # rate of edge fusion (was previously: 0.001. new value: 0.0005) --> must preserve detailed balance (i.e., P < 1)

    # Call the main function
    main(interactions_file, angles_file, output_folder, simulation_name, json_folder, eb_scan, kb_scan, mu, stretching_mod, l0, timesteps, num_tests, max_N, k_insertion, k_fusion2)