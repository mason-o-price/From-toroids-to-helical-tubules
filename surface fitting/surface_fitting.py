import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import pandas as pd

## Define functions

# Define rotation matrices
def Rx(theta: float, X, Y, Z) -> list:
    # Note: X remains unchanged
    # Rotate Y & Z
    Ynew = np.cos(theta)*Y - np.sin(theta)*Z
    Znew = np.sin(theta)*Y + np.cos(theta)*Z
    return [X, Ynew, Znew]
  
def Ry(theta: float, X, Y, Z) -> list:
    # Note: Y remains unchanged
    # Rotate X & Z
    Xnew = np.cos(theta)*X + np.sin(theta)*Z
    Znew = -np.sin(theta)*X + np.cos(theta)*Z
    return [Xnew, Y, Znew]
  
def Rz(theta: float, X, Y, Z) -> list:
    # Note: Z remains unchanged
    # Rotate X & Y
    Xnew = np.cos(theta)*X - np.sin(theta)*Y
    Ynew = np.sin(theta)*X + np.cos(theta)*Y
    return [Xnew, Ynew, Z]

# Define a parametrized circle embedded in R3
def circle(R: float) -> list:
    resolution = 50

    # Define an array of the independent parameter
    t = np.linspace(0, 2*np.pi, resolution)

    # Define the parametrization
    x = R*np.cos(t)
    y = R*np.sin(t)
    z = t*0

    # Return the coordinates
    return [x,y,z]

# Define a parametrized line (along the z-axis)
def line(tmin: float, tmax: float) -> list:
    # Set the number of points to create
    resolution = 100
    
    # Create a linear space for the parameter
    t = np.linspace(tmin, tmax, resolution)

    # Define the parametrization
    x = 0
    y = 0
    z = t

    # Return the coordinates
    return [x,y,z]

# Define a parametrized helix in R3
def helix(R: float, p: float, tmin: float, tmax: float) -> list:
    # Set the number of points to create
    resolution = 100
    
    # Create a linear space for the parameter
    t = np.linspace(tmin, tmax, resolution)

    # Define the parametrization
    x = R*np.cos(t)
    y = R*np.sin(t)
    z = p*t/(2*np.pi)

    # Return the coordinates
    return [x,y,z]

# Define a parametrized wave curve
def wave_curve(R: float, p: float, tmin: float, tmax: float) -> list:
    # Set the number of points to create
    resolution = 100
    
    # Create a linear space for the parameter
    t = np.linspace(tmin, tmax, resolution)

    # Define the parametrization
    x = R*np.cos(t)
    y = 0*t
    z = p*t/(2*np.pi)

    # Return the coordinates
    return [x,y,z]

# Define a parametrized torus in terms of the minor- and major-radii
def torus(r: float, R: float) -> list:
    resolution = 100 # number of sample points for the parametrized surface in terms of the variables u,v
    
    # create a mesh of u,v
    u = np.linspace(0, 2*np.pi, resolution)
    v = np.linspace(0, 2*np.pi, resolution)
    u,v = np.meshgrid(u, v) 

    # Define the parametrization
    x = (R + r*np.cos(u))*np.cos(v)
    y = (R + r*np.cos(u))*np.sin(v)
    z = r*np.sin(u)

    # Return the coordinates
    return [x,y,z]

# Define a parametrized helical tubule in terms of the minor- and major-radii and pitch
def helical_tubule(r: float, R: float, p: float, umin: float, umax: float) -> list:
    resolution = 200 # number of sample points for the parametrized surface in terms of the variables u,v
    
    # create a mesh of u,v
    u = np.linspace(umin, umax, resolution)
    v = np.linspace(0, 2*np.pi, resolution)
    u,v = np.meshgrid(u, v)

    # Define the parametrization
    x = R*np.cos(u) - r*np.cos(u)*np.cos(v) + (p / (2*np.pi))*r*np.sin(u)*np.sin(v) / np.sqrt(R**2 + (p / (2*np.pi))**2)
    y = R*np.sin(u) - r*np.sin(u)*np.cos(v) - (p / (2*np.pi))*r*np.cos(u)*np.sin(v) / np.sqrt(R**2 + (p / (2*np.pi))**2)
    z = (p / (2*np.pi))*u + R*r*np.sin(v)/np.sqrt(R**2 + (p / (2*np.pi))**2)

    # Return the coordinates
    return [x,y,z]

# Define a function to transform a set of coordinates (rotate (theta0, alpha0, gamma0) -> translate + (x0, y0, z0))
def transform(params: list, x, y, z) -> list:
    # Unpack the parameters
    [x0, y0, z0, theta0, alpha0, gamma0] = params

    # Rotate the coordinates by angles theta0, alpha0, gamma0 about the x-, y- and z-axes respectively. 
    [x,y,z] = Rx(theta0, x, y, z)
    [x,y,z] = Ry(alpha0, x, y, z)
    [x,y,z] = Rz(gamma0, x, y, z)
    
    # Offset the coordinates by a translation vector (x0, y0, z0)
    [x,y,z] = [x+x0, y+y0, z+z0]

    return [x,y,z]

# Define a reverse function to transform a set of coordinates (translate -(x0, y0, z0) -> rotate -(theta0, alpha0, gamma0))
def transform_reverse(params: list, x, y, z) -> list:
    # Unpack the parameters
    [x0, y0, z0, theta0, alpha0, gamma0] = params

    # Offset the coordinates by a translation vector (x0, y0, z0)
    [x,y,z] = [x-x0, y-y0, z-z0]

    # Rotate the coordinates by angles theta0, alpha0, gamma0 about the x-, y- and z-axes respectively. 
    [x,y,z] = Rz(-gamma0, x, y, z)
    [x,y,z] = Ry(-alpha0, x, y, z)
    [x,y,z] = Rx(-theta0, x, y, z)
    
    return [x,y,z]

# Define a function to calculate the distance from a surface to the set of vertices
def distance(vertex: list, surface: list) -> float:
    # Unpack the coordinates of the surface
    [x,y,z] = surface

    # Return the 
    return np.min((vertex[0]-x)**2 + (vertex[1]-y)**2 + (vertex[2]-z)**2)

# Define a loss function for a 1D circle
def loss_circle(params: list, vertices: list) -> float:    
    # Unpack the parameters
    [R, x0, y0, z0, theta0, alpha0, gamma0] = params
    
    # Create the coordinates for the parametrized torus
    [x,y,z] = circle(R)
    
    # Offset the curve by a translation and rotation
    [x,y,z] = transform([x0, y0, z0, theta0, alpha0, gamma0], x, y, z)
    
    # Define the loss function L = sum(min(distance from vertices to surface)^2)
    return np.sum([distance(vertex, [x,y,z]) for vertex in vertices])

# Define a loss function for a line
def loss_line(params: list, vertices: list) -> float:
    # Unpack the parameters
    [tmin, tmax, x0, y0, z0, theta0, alpha0, gamma0] = params
    
    # Create the coordinates for a parametrized helical tubule
    [x,y,z] = line(tmin, tmax)

    # Offset the curve by a translation and rotation
    [x,y,z] = transform([x0, y0, z0, theta0, alpha0, gamma0], x, y, z)
    
    # Define the loss function L = sum(min(distance from vertices to surface)^2)
    return np.sum([distance(vertex, [x,y,z]) for vertex in vertices])

# Define a loss function for a 1D helix
def loss_helix(params: list, vertices: list) -> float:
    # Unpack the parameters
    [R, p, tmin, tmax, x0, y0, z0, theta0, alpha0, gamma0] = params
    
    # Create the coordinates for a parametrized helical tubule
    [x,y,z] = helix(R, p, tmin, tmax)

    # Offset the curve by a translation and rotation
    [x,y,z] = transform([x0, y0, z0, theta0, alpha0, gamma0], x, y, z)
    
    # Define the loss function L = sum(min(distance from vertices to surface)^2)
    return np.sum([distance(vertex, [x,y,z]) for vertex in vertices])

# Define a loss function for a wave curve
def loss_wave1D(params: list, vertices: list) -> float:
    # Unpack the parameters
    [R, p, tmin, tmax, x0, y0, z0, theta0, alpha0, gamma0] = params
    
    # Create the coordinates for a parametrized helical tubule
    [x,y,z] = wave_curve(R, p, tmin, tmax)

    # Offset the curve by a translation and rotation
    [x,y,z] = transform([x0, y0, z0, theta0, alpha0, gamma0], x, y, z)
    
    # Define the loss function L = sum(min(distance from vertices to surface)^2)
    return np.sum([distance(vertex, [x,y,z]) for vertex in vertices])

# Define a loss function for a 2D torus surface
def loss_torus(params: list, vertices: list) -> float:
    # Unpack the parameters
    [r, R, x0, y0, z0, theta0, alpha0, gamma0] = params
    
    # Create the coordinates for the parametrized torus
    [x,y,z] = torus(r, R)
        
    # Offset the surface by a translation and rotation
    [x,y,z] = transform([x0, y0, z0, theta0, alpha0, gamma0], x, y, z)
    
    # Define the loss function L = sum(min(distance from vertices to surface)^2)
    return np.sum([distance(vertex, [x,y,z]) for vertex in vertices])

# Define a loss function for a 2D helical tubule surface
def loss_helical(params: list, vertices: list) -> float:
    # Unpack the parameters
    [r, R, p, tmin, tmax, x0, y0, z0, theta0, alpha0, gamma0] = params
    
    # Create the coordinates for a parametrized helical tubule
    [x,y,z] = helical_tubule(r, R, p, tmin, tmax)

    # Offset the surface by a translation and rotation
    [x,y,z] = transform([x0, y0, z0, theta0, alpha0, gamma0], x, y, z)
    
    # Define the loss function L = sum(min(distance from vertices to surface)^2)
    return np.sum([distance(vertex, [x,y,z]) for vertex in vertices])

def meanXYZ(coordinates: list) -> float:
    XYZavg = np.array([0., 0., 0.])
    num_vertices = len(coordinates)
    for i in range(num_vertices):
        XYZavg += coordinates[i,:]
    XYZavg /= num_vertices
    return XYZavg

def main(file_name: str, surface_type: str, c: float, L: float, N: float, p: float):
    ## For a helical tubule
    # Import the coordinates of the vertices
    vertices = pd.read_csv(file_name, header = None)
    vertices = vertices.to_numpy()

    # Find the average coordinates
    [x0, y0, z0] = meanXYZ(vertices)

    # Perform optimization
    if surface_type == 'T':
        R0 = L/2/np.tan(2*np.pi/N)
        p0 = [R0, x0, y0, z0, 0, 0, 0]

        # Perform optimization
        initial_fit = minimize(loss_circle, p0, args=(vertices,))
        initial_params = initial_fit.x

        r0 = c/(2*np.pi)
        p1 = np.insert(initial_params, 0, r0)

        # Perform optimization
        result = minimize(loss_torus, p1, args=(vertices,))

        # Optimized parameters
        optimized_params = result.x

        print("Optimized params:")
        print("[r, R, x0, y0, z0, theta0, alpha0, gamma0]")
        print(optimized_params)

        [r, R, x0, y0, z0, theta0, alpha0, gamma0] = optimized_params
        [x,y,z] = torus(r, R)
        [x,y,z] = transform([x0, y0, z0, theta0, alpha0, gamma0], x,y,z) 

        # Plot the torus and the vertices
        fig = plt.figure() # dpi=200
        ax = fig.add_subplot(111, projection='3d')

        ax.plot_surface(x, y, z, color='b', alpha = 0.2, label='Surface fit')

        # ax.scatter3D(xyz_0[0], xyz_0[1], xyz_0[2], 'bo')

        ax.scatter3D(vertices[:,0], vertices[:,1], vertices[:,2], color = 'k', label = 'Vertices');
        plt.axis('equal')

        ax.set_axis_off()

        # Plot the circle inside
        [R, x0, y0, z0, theta0, alpha0, gamma0] = initial_params
        [x,y,z] = circle(R)
        [x,y,z] = transform([x0, y0, z0, theta0, alpha0, gamma0], x,y,z) 
        ax.view_init(elev=30, azim=120)
        ax.plot3D(x, y, z, 'g', label = 'Co-planar circular axis')

        plt.legend(frameon=False)
        plt.show()
        
    elif surface_type == 'H':
        # Choose an initial parameter tuple
        p0 = [-1, 1, x0, y0, z0, 0, 0, 0]

        # Fit a line to find the axis of the helix
        axis_fit = minimize(loss_line, p0, args=(vertices,))
        axis_params = axis_fit.x
        [tmin, tmax, x0, y0, z0, theta0, alpha0, gamma0] = axis_params
        R0 = L/2/np.tan(2*np.pi/N)

    # Choose a new initial parameter tuple
        p1 = [R0, p, tmin, tmax, x0, y0, z0, theta0, alpha0, gamma0]

        # Perform optimization
        helix_fit = minimize(loss_helix, p1, args=(vertices,))
        helix_params = helix_fit.x
        [R, p, tmin, tmax, x0, y0, z0, theta0, alpha0, gamma0] = helix_params

        r0 = c/(2*np.pi)
        p2 = [r0, R, p, tmin, tmax, x0, y0, z0, theta0, alpha0, gamma0]

        # Do the final optimization
        result = minimize(loss_helical, p2, args=(vertices,))
        # Optimized parameters
        optimized_params = result.x

        print("Optimized params:")
        print("[r, R, p, umin, umax, x0, y0, z0, theta0, alpha0, gamma0]")
        print(optimized_params)

        [r, R, p, umin, umax, x0, y0, z0, theta0, alpha0, gamma0] = optimized_params
        [x,y,z] = helical_tubule(r, R, p, umin, umax)
        [x,y,z] = transform([x0, y0, z0, theta0, alpha0, gamma0], x, y, z) 

        # Prepare a figure
        fig = plt.figure() # dpi=200
        ax = fig.add_subplot(111, projection='3d')

        # Plot the surface
        ax.plot_surface(x, y, z, color='b', alpha = 0.2, label='Surface fit')

        # Plot the vertices
        ax.scatter3D(vertices[:,0], vertices[:,1], vertices[:,2], color = 'k', label = 'Vertices')
        plt.axis('equal') # Make the axes scales the same

        ax.set_axis_off() # Turn off the axes

        # Plot the axis of the helix
        [tmin, tmax, x0, y0, z0, theta0, alpha0, gamma0] = axis_params
        [x,y,z] = line(tmin, tmax)
        [x,y,z] = transform([x0, y0, z0, theta0, alpha0, gamma0], x, y, z) 
        ax.plot3D(x, y, z, 'r', label = 'Axis')

        ## If you want to re-center the axis:
        # # Offset the axis to be centered at the origin
        # [x_new, y_new, z_new] = np.array(transform_reverse([x0, y0, z0, theta0, alpha0, gamma0], x, y, z))
        # ax.plot3D(x_new, y_new, z_new, 'r--')

        # Plot the helix inside
        [R, p, tmin, tmax, x0, y0, z0, theta0, alpha0, gamma0] = helix_params
        [x,y,z] = helix(R, p, tmin, tmax)
        [x,y,z] = transform([x0, y0, z0, theta0, alpha0, gamma0], x, y, z) 
        ax.plot3D(x, y, z, 'g', label = 'Co-axial helix')

        # If you want to re-center the helix and surface fit:
        # # Offset the helix to be centered at the origin
        # [x_new, y_new, z_new] = np.array(transform_reverse([x0, y0, z0, theta0, alpha0, gamma0], x, y, z))
        # ax.plot3D(x_new, y_new, z_new, 'g--')

        # # Offset the vertices to be centered at the origin
        # [x_new, y_new, z_new] = np.array(transform_reverse([x0, y0, z0, theta0, alpha0, gamma0], vertices[:,0], vertices[:,1], vertices[:,2]))
        # ax.scatter3D(x_new, y_new, z_new, color = 'k', facecolors='none', label = 'Offset vertices') # Plot the offset vertices

        # # Save the offset vertex coordinates
        # np.savetxt('offset_coordinates.csv', (x_new, y_new, z_new), delimiter=',')

        # # Plot the axes centered at the origin  
        # ax.scatter3D(0, 0, 0, 'b')
        # ax.plot3D([0, 1], [0, 0], [0, 0], 'b-') # x-axis
        # ax.plot3D([0, 0], [0, 1], [0, 0], 'g-') # y-axis
        # ax.plot3D([0, 0], [0, 0], [0, 1], 'r-') # z-axis

        ax.view_init(elev=10, azim=90) # Set the viewing angle
        plt.legend(frameon=False) # Make a legend
        plt.show() # Show the plot
    elif surface_type == 'S':
        # Choose an initial parameter tuple
        p0 = [-1, 1, x0, y0, z0, 0, 0, 0]

        # Fit a line to find the axis of the helix
        axis_fit = minimize(loss_line, p0, args=(vertices,))
        axis_params = axis_fit.x
        [tmin, tmax, x0, y0, z0, theta0, alpha0, gamma0] = axis_params

        R0 = 1
        p = 1
        # R0 = L/2/np.tan(2*np.pi/N) + r0

        # Choose a new initial parameter tuple
        p1 = [R0, p, tmin, tmax, x0, y0, z0, theta0, alpha0, gamma0]

        # Perform optimization for the co-axial helix
        wave_fit1D = minimize(loss_helix, p1, args=(vertices,))
        wave_params1D = wave_fit1D.x

        print("Optimized params:")
        print("[R, p, umin, umax, x0, y0, z0, theta0, alpha0, gamma0]")
        print(wave_params1D)

        # Prepare a figure
        fig = plt.figure() # dpi=200
        ax = fig.add_subplot(111, projection='3d')

        # Plot the vertices
        ax.scatter3D(vertices[:,0], vertices[:,1], vertices[:,2], color = 'k', label = 'Vertices')
        plt.axis('equal') # Make the axes scales the same

        ax.set_axis_off() # Turn off the axes

        # Plot the axis of the helix
        [tmin, tmax, x0, y0, z0, theta0, alpha0, gamma0] = axis_params
        [x,y,z] = line(tmin, tmax)
        [x,y,z] = transform([x0, y0, z0, theta0, alpha0, gamma0], x, y, z) 
        ax.plot3D(x, y, z, 'r', label = 'Axis')

        # Plot the helix inside
        [R, p, tmin, tmax, x0, y0, z0, theta0, alpha0, gamma0] = wave_params1D
        [x,y,z] = wave_curve(R, p, tmin, tmax)
        [x,y,z] = transform([x0, y0, z0, theta0, alpha0, gamma0], x, y, z) 
        ax.plot3D(x, y, z, 'g', label = 'Co-axial wave')

        ax.view_init(elev = 0, azim = 100) # Set the viewing angle
        plt.legend(frameon=False) # Make a legend
        plt.show() # Show the plot
    else:
        Warning(f"Unrecognized surface type. \n -> Choose 'T' for toroid, 'H' for helical tubule, or 'S' for snake. <-")

if __name__ == "__main__":
    # Specify the name of the file storing the coordinates of the vertices
    file_name = "(3,1,1,1)_coordinates.csv"

    surface_type = 'H' # 'T' = toroid, 'H' = helical tubule, 'S' = snake

    # Input parameters (T,L,D,R) with the corresponding axis symmetry N and best guess for pitch p.
    T = 3
    L = 1
    N = 4
    p = -1 # Positive for right-handed, negative for left-handed

    # Call the main function
    main(file_name, surface_type, T, L, N, p)