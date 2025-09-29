To run kinetic Monte Carlo (KMC) simulations under many different conditions, we prepare simulations and analyze the results as follows. 
1. Run multiple_json_generator.py to generate all of the JSON files used as inputs to the simulation. 
2. Transfer the JSON files to the high performance computing cluster (HPCC), using WINSCP or another tool. 
3. Run submitArray.sh on the HPCC to run the array of simulations in parallel, being sure to specify the number of jobs. One can count the number of JSON files using, for example, countFiles.sh
4. When the simulations finish, run moveData.sh to grab the essential files from each simulation output. 
5. Transfer the summarized data to a local machine.
6. Run quench_json_generator.py to generate JSON files to be used as inputs for the quenching procedure. 
7. Transfer the quench JSON files to the HPCC. 
8. Run submitQuenchArray.sh to run the array of quench simulations. 
9. Run moveData.sh to summarize the quench output data.
10. Move the quench data to a local machine. 
10. Run parameter_scan_analysis.py to analyze the simulation results, using both the primary simulation data and quench data. 


To run a single simulation (for example to create snapshot visual data), we do the following. 
1. Run single_json_generator.py with the desired parameters. 
2. Either run the simulation on a local machine or run submit.sh on the HPCC.
3. Run quench_json_generator.py as before. 
4. Run the quench simulation either on a local machine or on the HPCC using submitQuench.sh. 