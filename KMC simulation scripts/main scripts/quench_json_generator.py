# Import libraries
import os
import json
import csv
import warnings

def main(json_directory: str, data_log_directory: str, output_directory: str):
    for file_name in os.listdir(data_log_directory):
        #read the data
        with open(os.path.join(data_log_directory, file_name)) as f:
            reader = csv.reader(f, delimiter="\t")
            data_log = [line for line in reader]

        # Find the final binary key in the simulation
        final_bin_key = int(data_log[-1][1])
        
        # Close the file
        f.close() 

        # Take the last part of the file name, e.g. output_4221__eb-6_kb1000_0_ensemble32077217.dat --> eb-6_kb1000_0_ensemble32077217.dat
        if "__" not in file_name:
            raise Warning(f"The data_log file name must have a double-underscore '__' to separate the overall test vs specific ensemble, \n e.g., output_4221__eb-6_kb1000_0_ensemble32077217.dat --> eb-6_kb1000_0_ensemble32077217.dat \n Your data_log file name is: {file_name}")
        sim_id = file_name.split("__", 1)[1] # Note: this expects that the data_folder names have a double-underscore separating the overall name from the specific conditions
        sim_id = sim_id.replace(".dat", "") # Get rid of the .dat extension of the data_log file.

        # Now try to find the corresponding .json file
        matched = False # Assume by default we have not yet matched the data_log file to the json file. 
        for json_name in os.listdir(json_directory):
            if sim_id in json_name.replace(".", "_"): # Check if the simulation id string is in the json file name.
                matched = True # Record the fact that we found the corresponding json file

                # Read the json file
                with open(os.path.join(json_directory, json_name)) as jf:
                    json_file = json.load(jf)
                jf.close()
                
                # Create a copy of the json file to change for the quenching procedure
                quench_json = json_file

                # Find the location where the data was output
                output_data_folder = json_file['parameters']['run_parameters']['data_folder']
                if not isinstance(output_data_folder, str):
                    raise TypeError("output_data_folder must be a string")
                    
                # Write new time parameters for the json file
                quench_json['parameters']['run_parameters']['timesteps'] = 10000
                quench_json['parameters']['run_parameters']['dtsave'] = 100
                quench_json['parameters']['run_parameters']['nsteps_umbrella'] = 10000
                quench_json['parameters']['run_parameters']['checkpoint_freq'] = 100
                quench_json['parameters']['run_parameters']['dt_burnin'] = 100
                quench_json['parameters']['run_parameters']['dtsave_snapshot'] = 100
                quench_json['parameters']['run_parameters']['convert_to_lammps'] = True
                quench_json['parameters']['run_parameters']['convert_to_vtk'] = True
                
                # Set the initial structure to be the output of the simulation
                quench_json['parameters']['run_parameters']['init_file'] = output_data_folder + 'snapshots.om'

                # Set the initial file position to be the final binary key from the original simulation
                quench_json['parameters']['run_parameters']['init_file_pos'] = final_bin_key

                # Change the location of where the quench data is output
                original_data_folder = quench_json['parameters']['run_parameters']['data_folder']
                if "data" not in original_data_folder:
                    raise Warning(f"The name of the data folder must include 'data' in the parent folder for this test, \n e.g., /home/masonprice/data_toroid_new_sim/output_4221__eb-6_kb0_1_ensemble5838983 \n Your data folder name is: {original_data_folder}") 
                quench_json['parameters']['run_parameters']['data_folder'] = original_data_folder.replace("data", "quench_data")
                
                # Now save the json file
                with open(output_directory + "/quench_" + json_name, "w") as outfile:
                    outfile.write(json.dumps(quench_json, indent=4))
                outfile.close()

        if matched == False:
            warnings.warn(f"Failed to match {file_name}. \n Simulation ID: {sim_id}")
    
    # Signal when the function is done.
    print("Complete.")

if __name__ == "__main__":
    # Specify where you saved all of the original json input files to run the simulations in the first place
    json_directory = r''

    # Specify the location of the data_log.dat files that result from the simulations
    data_log_directory = r''

    # Specify the name of the output folder
    output_directory = ""

    # Test whether the output folder exists
    if not os.path.isdir(output_directory):
        raise Warning(f"The output directory {output_directory} does not exist. \n The directory must exist before running.")

    # Call the main function
    main(json_directory, data_log_directory, output_directory)