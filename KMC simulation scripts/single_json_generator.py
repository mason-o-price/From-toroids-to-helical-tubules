# import libraries
import numpy as np
import pylab as pl
import pandas as pd
import json
import bson
import copy
import random

def main(intMatrixFile: str, angleMatrixFile: str, outputDataFolder: str, jsonFileName: str, mu: float, e_b: float, k_b: float, timesteps=1e8):
    # Read the interaction matrix of the toroid we want to generate
    int_matrix = pd.read_csv(intMatrixFile, header = None)
    print(int_matrix)

    # Set the stretching modulus and default edge-length
    stretching_mod = 100
    l0 = 1.0
    df_eb=int_matrix.copy()
    df_eb[df_eb!=1]=np.nan
    df_eb[df_eb==1]=e_b

    # Create the matrix of binding angles
    df_kb = int_matrix.copy()
    df_kb[df_kb != 1] = np.nan
    df_kb[df_kb == 1] = k_b

    # read the binding angles (in degrees)
    df_theta = pd.read_csv(angleMatrixFile, header = None)
    print(df_theta)
    df_theta = df_theta*np.pi/180
    np_theta = df_theta.to_numpy()

    #generate a list of the subunits and their edge types
    subunits = []
    for i, e in enumerate(int_matrix.columns[::3]):
        subunits.append([i+1, e+1, e+2, e+3])
    subunits = np.array(subunits)

    templatefile = 'template_input.json'
    with open(templatefile) as tf:
        template = json.load(tf)
    tf.close()
    new_json = template

    interactions = int_matrix.to_numpy()
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
    
    for i in range(len(subunits)):
        new_json['parameters']['subunits']['sub'+str(i+1)] = {'edges': {'edge1': 'type'+str(subunits[i][1]),
                                                                        'edge2': 'type'+str(subunits[i][2]), 
                                                                        'edge3': 'type'+str(subunits[i][3])},
                                                                        'k_insertion': 0.00025,
                                                                        'mu': mu,
                                                                        'number_of_rotational_configs': 3,
                                                                        'R_exc': 0.15}
        
    # change some other parameters (R_add, d_max, l_fuse) depending on the values of the bending/stretching modulus.
    alp=2.0
    kT=1.0
    l_thermal_stretch = min( 2.0*(alp*kT/stretching_mod)**0.5, 0.5*l0 )
    l_thermal_bend = min( (3.0*l0**2*alp*kT/(2.0*k_b))**0.5, 0.5*l0 )
    #l_thermal = max(2.0*(alp*kT/eps_s)**0.5, (3.0*l0**2*alp*kT/(2.0*kappa_b))**0.5)
    l_thermal = max(l_thermal_stretch, l_thermal_bend)
    new_json['parameters']['mesh_parameters']['R_add']=l_thermal #/3.5 /5.0
    new_json['parameters']['mesh_parameters']['d_max']=l_thermal / 3.5 /1.41
    new_json['parameters']['mesh_parameters']['l_fuse']=l_thermal
    ensemble = random.randint(1000, 500000)
    new_json['parameters']['run_parameters']['ensemble'] = ensemble
    new_json['parameters']['run_parameters']['timesteps'] = timesteps # usually: 1e7
    new_json['parameters']['run_parameters']['convert_to_lammps'] = True
    new_json['parameters']['run_parameters']['convert_to_vtk'] = True
    new_json['parameters']['run_parameters']['convert_to_dat'] = True
    new_json['parameters']['run_parameters']['Nmax'] = 200 # usually: 200
    new_json['parameters']['run_parameters']['dtsave'] = 100000
    new_json['parameters']['run_parameters']['dtsave_snapshot'] = 100000
    new_json['parameters']['run_parameters']['checkpoint_freq'] = 100000
    new_json['parameters']['run_parameters']['data_folder'] = outputDataFolder
    # print(new_json)

    with open(jsonFileName, "w") as outfile:
        outfile.write(json.dumps(new_json, indent=4))
        print(f"Wrote file {jsonFileName}")

if __name__ == "__main__":
    # Set the file/folder names
    intMatrixFile = "interactionMatrices/trefoil_interactionMatrix.csv" # interaction matrix
    angleMatrixFile = "bindingAngles/trefoil_angleMatrix.csv" # binding angles matrix
    outputDataFolder = "/home/mason/openmesh_output/data_trefoil_kb100_eb-8/" # folder to dump data from the output of the simulation
    jsonFileName = "input_trefoil_kb100_eb-8.json" # name of the json file

    # set the chemical potential, binding energy and bending modulus
    mu = -4.5
    e_b = -8
    k_b = 100

    # Set the number of timesteps
    timesteps = 1e8

    main(intMatrixFile, angleMatrixFile, outputDataFolder, jsonFileName, mu, e_b, k_b, timesteps)