# Import libraries
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import os
import csv
from matplotlib.colors import LogNorm
import warnings

def main(data_location: str, quench_location: str, num_runs: int):
    # Define the array of binding energy and bending modulus values that we are scanning over
    eb_scan = [f'eb-{x}' for x in [5,6,7,8,9,10,11]] # Binding energy
    kb_scan = ['kb0_001', 'kb0_01', 'kb0_1', 'kb1_0', 'kb10_0', 'kb100_0', 'kb1000_0'] # Bending modulus

    num_rows = 7 # Number of rows in the heatmap
    num_cols = 7 # Number of columns in the heatmap

    # Create a dictionary to store the categories of assembly (key = location (i,j), value = list of assembly types)
    categories = {}
    for row in range(num_rows):
        for col in range(num_cols):
            categories.update({(row, col): []})

    # Define matrices to store the other types of data we want to store
    MQS_matrix = np.zeros([num_rows, num_cols]) # Matrix for the mean quadratic strain values (MQS)
    N_total_matrix = np.zeros([num_rows, num_cols]) # Matrix for the total number of triangles
    P_closed_matrix = np.zeros([num_rows, num_cols]) # Matrix for the probability of closure
    yield_matrix = np.zeros([num_rows, num_cols]) # Matrix for the yield
    P_5fold = np.zeros([num_rows, num_cols]) # Matrix for the frequency of 5-fold toroids (specifically for the (4,2,2,1) toroid)

    for eb in eb_scan: # For each binding energy ...
        for kb in kb_scan: # For each bending modulus ...
            dir_list = os.listdir(data_location) # Find the directory
            for file_name in dir_list: # For each file name in the directory ...
                if ( eb in file_name ) and ( kb in file_name ) and ( file_name[0] == 'o' ): # If the file name includes the corresponding binding energy and bending modulus ...
                    # Read the data from the original simulation
                    with open(os.path.join(data_location, file_name)) as f:
                        reader = csv.reader(f, delimiter="\t")
                        file = [line for line in reader]

                    # Read the data from the quenching simulation
                    with open(os.path.join(quench_location, file_name)) as qf: 
                        quench_reader = csv.reader(qf, delimiter="\t")
                        quench_file = [line for line in quench_reader]
                    
                    # Find the corresponding row and column in the matrix (for the current values of eb and kb)
                    row = eb_scan.index(eb)
                    col = kb_scan.index(kb)

                    # ---- Read the data ----
                    N_total = int(file[-1][2]) # Find the total number of triangles at the end of the simulation

                    # Test if the structure is closed
                    is_closed = False # Assume by default the structure is open
                    if ( int(file[-1][-1]) == 0 ):
                        is_closed = True # If there are no unbound edges -> closed

                    MQS = float(quench_file[-1][-1]) # Read the Mean Quadratic Strain from the result of the quenching (MQS)

                    # ===== Start storing the read values in the matrices =====
                    MQS_matrix[row, col] += MQS  # Add the current value of MQS to the corresponding entry so that we can take the average
                    N_total_matrix[row, col] += N_total # Add the current value of N_total to the corresponding entry so that we can take the average

                    # Find the probability of being closed
                    if is_closed:
                        P_closed_matrix[row, col] += 1 # Increment the count of closed assemblies so that we can normalize to find the probability

                    # ---- Classify the assembly ----
                    assembly_type = None # Set the default assembly type to be None

                    # If the assembly has very low strain and is not undernucleated, classify it as on-target (1). 
                    if MQS < 1e-4 and N_total >= 6:
                        assembly_type = 1
                        yield_matrix[row, col] += 1 # Update the yield matrix
                    # If the assembly has less than 6 triangles, classify it as undernucleated (2)
                    if N_total < 6:
                        assembly_type = 2
                    # If the assembly has any strain (and is not undernucleated), classify it as overgrown (3)
                    if MQS >= 1e-4 and N_total >= 6:
                        assembly_type = 3
                    # If the assembly is closed but with high strain, classify it as off-target closed (4)
                    if ( is_closed ) and ( MQS > 1e-4 ):
                        assembly_type = 4

                    # Check if the assembly is a 5-fold
                    if ( is_closed ) and ( N_total == 120):
                        P_5fold[row, col] += 1
                        print(f"Found 5-fold at {file_name}")

                    # Throw a warning if we fail to classify the assembly
                    if assembly_type == None:
                        warnings.warn(f"Failed to classify assembly at {eb}, {kb} \n N_total={N_total}, is_closed={is_closed}, strain={MQS}")
                    
                    # Append the assembly type to the corresponding data point
                    categories[(row, col)].append(assembly_type) # add the value from this file to the current number in the matrix, so that in the end we divide by the number of runs to get the average

    # Normalize the matrices
    yield_matrix /= num_runs # Yield matrix
    N_total_matrix /= num_runs # Total number of triangles
    final_MQS = MQS_matrix/num_runs # Mean quadratic strain
    P_closed_matrix /= num_runs # Probability of being closed
    P_5fold /= num_runs # Probability of a 5-fold

    # --------------------- Most probable assembly type ---------------------
    # Find the most probable assembly type for each data point:
    most_probable = np.zeros([num_rows,num_cols]) # Create a matrix of zeros to store the most probable assembly type. 

    for row in range(num_rows):
        for col in range(num_cols):
            assemblies = categories[(row,col)]
            most_probable[row, col] = max(assemblies, key=assemblies.count)
        
    type_data = pd.DataFrame({'$10^{%d}$' % int(x-3): most_probable[:,x] for x in range(7)})
    type_data.index = [-5,-6,-7,-8,-9,-10, -11]
    print(f"Categories: {type_data}")
    # Create a figure
    plt.figure()
    ax = plt.axes()
    plt.title('Most probable assemblies') # Choose the title --> assembly type

    # Define the colors for the different types of assembly
    assembly_colors = [[0.09019608, 0.71372549, 0.7372549], # 1 <-> on-target
                        [0.99215686, 0.96470588, 0.83529412], # 2 <->  undernucleated
                        [0.98039216, 0.74117647, 0.63137255], # 3 <-> overgrown
                        [0.94117647, 0.44313725, 0.40392157]] # 4 <-> off-target, closed

    # Plot the heat map
    sns.heatmap(type_data, cmap = assembly_colors, linewidths = 0.5, linecolor = 'gray')

    # Make the axis scales equal
    ax.set_aspect('equal', adjustable='box') 
    plt.xlabel('$\kappa_b$ $(k_B T)$') # Label the x-axis
    plt.ylabel('$\epsilon_b$ $(k_B T)$') # Label the y-axis

    # Show the figure
    plt.show()

    # ---------- Plot a circle-probability figure of the assembly types ---------------
    # Print the yield matrix first
    print(f"Yield: {yield_matrix}")

    # Create a figure
    plt.figure()
    ax = plt.axes()
    plt.title('Most probable assemblies') # Choose the title --> assembly type (circles)
    plt.ylim((-11.5, -4.5)) # Set the upper- and lower-bounds for x
    plt.xlim((-3.5, 3.5)) # Set the upper- and lower-bounds for y 

    # Set the labels of the ticks on the x-axis
    # x_labels = [f"10^{x}" for x in range(-3, 4)] # Define the labels
    # ax.set_xticklabels(x_labels) # Set the property of the axes to be the labels we want

    # Define the values of x and y (for the plot)
    x_vals = [-3, -2,- 1, 0, 1, 2, 3]
    y_vals = [-5, -6, -7, -8, -9, -10, -11]

    # Plot all of the circles
    for i in range(7):
        for j in range(7):
            # Find the center of the circle
            x0 = x_vals[j]
            y0 = y_vals[i]

            # Retrieve the corresponding list of assembly results
            assembly_array = categories[(i, j)]

            # Retreive the most probable assembly type
            assembly_mode = int(most_probable[i, j])
            mode_probability = assembly_array.count(assembly_mode)/len(assembly_array)

            # Choose the color of the circl e based on the most probable type of assembly
            circle_color = assembly_colors[assembly_mode-1]

            # Calculate the radius of the circle
            P_max = 1 # Define the maximum probability
            P_min = 0.5 # Define the minimum probability
            r_max = 0.45 # Define the maximum circle radius
            r_min = r_max/np.sqrt(P_max/P_min) # Calcualte the minimum circle radius so that the area scales with the probability
            interp_param = (mode_probability - P_min)/P_min # Define a parameter for a linear interpolation (i.e., t such that 0<t<1)
            circle_r = np.sqrt(r_max**2*interp_param + (1-interp_param)*r_min**2)

            # Plot the corresponding circle
            circle = plt.Circle((x0, y0), circle_r, color = circle_color)
            plt.gca().add_patch(circle)

    # Make the axis scales equal
    ax.set_aspect('equal', adjustable='box')
    plt.xlabel('$\kappa_b$ $(k_B T)$') # Label the x-axis
    plt.ylabel('$\epsilon_b$ $(k_B T)$') # Label the y-axis

    # Show the figure
    plt.show()

    # --------------------- Plot the total number of triangles ---------------------
    # Turn the N_total data into a dataframe so that we can plot it as a heatmap
    N_total_data = pd.DataFrame({'$10^{%d}$' % int(x-3): N_total_matrix[:,x] for x in range(7)})
    N_total_data.index = [-5,-6,-7,-8,-9,-10,-11]
    print(f"Total number of triangles: {N_total_data}")

    # Create a figure
    plt.figure()
    ax = plt.axes()
    plt.title('$\mathcal{N}_{total}$') # Choose the title --> total number of triangles

    # Plot the heat map
    sns.heatmap(N_total_data, cmap = 'viridis', linewidths = 0.5, linecolor = 'gray')

    # Make the axis scales equal
    ax.set_aspect('equal', adjustable='box')
    plt.xlabel('$\kappa_b$ $(k_B T)$')
    plt.ylabel('$\epsilon_b$ $(k_B T)$')

    # Show the figure
    plt.show()

    # --------------------- Plot the strain (MQS) ---------------------
    # Turn the MQS data into a dataframe
    MQS_data = pd.DataFrame({'$10^{%d}$' % int(x-3): final_MQS[:,x] for x in range(7)})
    # pd.DataFrame({'$10^{-2}$': final_numbers[:,0], '$10^{-1}$':  final_numbers[:,1], '$10^{0}$':  final_numbers[:,2], '$10^{1}$':  final_numbers[:,3], '$10^{2}$':  final_numbers[:,4], '$10^{3}$':  final_numbers[:,5]})
    MQS_data.index = [-5,-6,-7,-8,-9,-10,-11]
    print(f"Mean Quadratic Strain (MQS): {MQS_data}")

    # Create a figure
    plt.figure() 
    ax = plt.axes()
    plt.title('$\Xi$') # Choose the title --> strain

    # Plot the heat map
    sns.heatmap(MQS_data, cmap = 'viridis', linewidths = 0.5, linecolor = 'gray', norm = LogNorm())

    # Make the axis scales equal
    ax.set_aspect('equal', adjustable='box')
    plt.xlabel('$\kappa_b$ $(k_B T)$') # Label the x-axis
    plt.ylabel('$\epsilon_b$ $(k_B T)$') # Label the x-axis

    # Show the figure
    plt.show()

    # --------------------- Plot the probability of closure ---------------------
    closed_data = pd.DataFrame({'$10^{%d}$' % int(x-3): P_closed_matrix[:,x] for x in range(7)})
    closed_data.index = [-5,-6,-7,-8,-9,-10,-11]
    print(f"Probability of closure: {closed_data}")

    # Create a figure
    plt.figure()
    ax = plt.axes()
    plt.title('$\mathbf{P}(closed)$') # Choose the title --> probability of closure

    # Plot the heatmap
    sns.heatmap(closed_data, cmap = 'mako', linewidths = 0.5, linecolor = 'gray', vmax = 1)

    # Make the axis scales equal
    ax.set_aspect('equal', adjustable='box')
    plt.xlabel('$\kappa_b$ $(k_B T)$') # Label the x-axis
    plt.ylabel('$\epsilon_b$ $(k_B T)$') # Label the y-axis

    # Show the figure
    plt.show()

    # --------------------- Plot the probability of forming a 5-fold ---------------------
    data_5fold = pd.DataFrame({'$10^{%d}$' % int(x-3): P_5fold[:,x] for x in range(7)})
    data_5fold.index = [-5,-6,-7,-8,-9,-10,-11]
    print(f"Probability of 5-folds: {data_5fold}")

    # Create a figure
    plt.figure() #dpi=200
    ax = plt.axes()
    plt.title('P(5-fold)')

    # Plot the heatmap
    sns.heatmap(data_5fold, cmap = 'rocket', linewidths = 0.5, linecolor = 'gray', vmin = 0)

    # Make the axis scales equal
    ax.set_aspect('equal', adjustable='box')
    plt.xlabel('$\kappa_b$ $(k_B T)$') # Label the x-axis
    plt.ylabel('$\epsilon_b$ $(k_B T)$') # Label the y-axis

    # Show the figure
    plt.show()

if __name__ == "__main__":
    # Define the locations of the data that we need 
    data_location = r'' # Location of the original simulation data files
    quench_location = r'' # Location of the quenching data files
    num_runs = 50 # Specify the total number of simulations ran per data point

    # Types of assembly:
    # 1: on-target
    # 2: undernucleated 
    # 3: over-grown
    # 4: off-target closed

    # Call the main function
    main(data_location, quench_location, num_runs)