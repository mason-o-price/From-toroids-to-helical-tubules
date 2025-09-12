# Additional Information for the article "From toroids to helical tubules: Kirigami-inspired programmable assembly of two-periodic curved crystals from DNA origami"

Authors: Mason Price, Daichi Hayakawa, Thomas E. Videbæk, Rupam Saha, Botond Tyukodi, Michael F. Hagan, Seth Fraden, Gregory M. Grason, W. Benjamin Rogers

This respository includes:
- [`generating geometries`](https://github.com/mason-o-price/From-toroids-to-helical-tubules/tree/main/generating%20geometries): MATLAB code used to generate triangulated geometries for curved tubules. The primary script is [`toroidSimulation`](https://github.com/mason-o-price/From-toroids-to-helical-tubules/blob/main/generating%20geometries/toroidSimulation.m)
- [`KMC simulation scripts`](https://github.com/mason-o-price/From-toroids-to-helical-tubules/tree/main/KMC%20simulation%20scripts): Python scripts to generate input json files for KMC simulations, and Python/Bash scripts to analyze the simulation results.
- [`surface fitting`](https://github.com/mason-o-price/From-toroids-to-helical-tubules/tree/main/surface%20fitting): Code and data used to fit a parametrized surface (of a torus or helical tube) to a set of vertices. The python script is [`surface_fitting.py`](https://github.com/mason-o-price/From-toroids-to-helical-tubules/blob/main/surface%20fitting/surface_fitting.py) and the data for each structure is stored as a CSV with the notation (T,L,D,R).
