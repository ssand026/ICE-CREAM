ICE-CREAM
(Integrated Codes to Efficiently Control Real-time Excitations in Architectured Materials)
---------

ICE-CREAM is an open-source software package designed for quantum control in 
real-space systems, utilizing a finite-element basis for flexible handling of complex 
nanostructures and other effective-mass systems.

The software package relies on two external packages:
1. MATLAB's PDE toolbox
	- required to generate meshes from external files or builtin geometry desciptions
	in the FE_mesh.load function.
2. GIBBON (https://www.gibboncode.org)
	- a constrained 3D Delaunay triangulation is required to generate meshes with 
	face constraints in the enclosing_sphere function, required for setting proper 
	boundary conditions for the Poisson equation.

All other functions require only the base MATLAB functionality. Written in MATLAB 2025b, 
backwards compatibility may be limited.


The general workflow is as follows:
1. Generate a "FE_mesh" object from a STL/STEP file, a builtin geometry description, or
	a set of node coordinates and their connectivities. Optionally, one can improve the
	numeric/mesh conditioning using "FE_mesh.reorder" and/or "FE_mesh.smooth" respectively.
2. Define the PDE coefficients, such as effective mass, confining potentials, etc. using 
	the generated mesh and the "FE_coeff" class.
3. Use FEmat and/or electric_dipole/magnetic_dipole to generate the matrix representation
	of the weak-form Hamiltonian. Custom dipole moments can be generated using "FEmat" and
	"FE_coeff" as well.
4. Set the boundary-conditions (defaults to Dirichlet), and choose whether to perform 
	calculations in the sparse/orthonormal basis; apply the basis/boundary-conditions
	to the Hamiltonian components using the "to_basis" function.
5. Compute eigenstates of the Hamiltonian, using the "eigenstates" function. The resulting
	eigenstates can be visualized using "meshplot2D" or "meshplot3D". Make sure to convert
	back to the "full" basis using "to_basis" before plotting.
6. Choose the initial/final states for the control problem. Selection rules for the system
	can be inferred by inputting the eigenstates and dipole-moments into "get_TDM".
7. Choose the control-duration and timestep sizes, the "get_QSL" function can be used to
	determine a lower-bound for the required control duration.
8. Set an intial guess for the control-field. A zero-valued intitial guess is often ideal.
9. (Optional) set a pulse-envelope for the control field via a custom array or "pulse_envelope".
10. Run the "optimize_field" function to generate an optimized control field.

Examples of this workflow can be found in the "examples" folder.