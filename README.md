# Bridge Testing Simulation
The simulation and optimization of a crane boom truss in Matlab was the final project for MTE 219. For our design, we chose to optimize a warren bridge.

# Description
The bridge was to be made from a 4" by 12" basswood sheet and span over a 30 cm gap. One side of the bridge was fastened down while the other was under the load. The bridge and its members were modelled. The method used is called the direct method and was used to create a stiffness matrix, which includes the stress of each member. By Looking at these failure modes and strengthening members under a critical load we were able to increase the bridge's strength to weight ratio (PV value). All values that were used in the stress calculations were obtained from a lab analyzing the material properties of basswood.

The report (.pdf) shows 3 optimized designed for the bridges. The initial design had a PV value at around 100 without optimization meaning that the bridge's performance was effectively doubled after optimization. Drawing of each of the bridges and their respective cutout patterns can be found in the report as well as the .SDLPRT files. The FEM_SampleTrussFiles_# show how the maximum mass for each iteration of the bridge was calculated. The MTE219_Prj_Grp_43.m is the final and best PV valued bridge and contains the most refined code and comments.

# Authors and acknowledgment
Authors: Luca Cristiano, David Feldt, Joshua Kurien
