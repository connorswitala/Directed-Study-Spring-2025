


Different grids and boundary conditions are detailed below (I just swap between these two). 
	
	 When declaring BoundaryConditions BCs, the order is (Left Boundary, Right Boundary, Bottom Boundary, Top Boundary)
	
	 Boundary Conditions are: Inlet, Outlet, Symmetry, Adiabatic Wall (No heat transfer), Isothermal Wall (constant wall temperature)
	
	 There are currently four possible grids. The cylinder grid is a little nonsensical at the moment. 



	 (1) Ramp Grid - First input: Nx, Second input: Ny, Third input: Length of entrance section, Fourth input: Length of ramp section, Fifth input: Length of exit section
	 Sixth input: Height of domain, Seventh input: Angle of ramp in degrees.
	 
	 
	BoundaryConditions BCs(BoundaryCondition::Inlet, BoundaryCondition::Outlet, BoundaryCondition::Symmetry, BoundaryCondition::Symmetry);   
	RampGrid grid(Nx, Ny, 10, 10, 10, 6, 15);
	

	 (2) Cylinder Grid - First input: Nx, Second input: Ny, Third input: Cylinder radius, Fourth input: Radius of fluid domain at theta = pi, 
	 Fifth input: Thickness of first layer of cells around cylinder.  Sixth Input: Use pi / 2, Seventh input: Use 3 * pi / 2 
	
	 
	BoundaryConditions BCs(BoundaryCondition::Outlet, BoundaryCondition::Outlet, BoundaryCondition::IsothermalWall, BoundaryCondition::Inlet);       
	CylinderGrid grid(Nx, Ny, 0.1, 0.3, 0.45, 0.00001, pi / 2, 3 * pi / 2); 


	 (3) Flate Plate - First input: Nx, Second input: Ny, Third input: Length in x-direction, Fourth input: Length in y-diretion, Fifth input: First layer cell thickness
	
	
	BoundaryConditions BCs(BoundaryCondition::Inlet, BoundaryCondition::Outlet, BoundaryCondition::IsothermalWall, BoundaryCondition::Symmetry);   
	FlatPlateGrid grid(Nx, Ny, 1e-3, 1e-3, 5e-6);  

	
	 (4) Double Cone (Essentially just double ramp) - First input: Nx, Second input: Ny, Third Input, Length of entrance section, Fourth input: Length of first angled ramp section,
	 Fifth input: Length of second angled ramp section, Sixth input: Length of exit section, Seventh input: Angle in degrees of first ramp, Eigth input: Angle in degrees of second ramp,
	 Ninth input: domain height. 
	
	
	BoundaryConditions BCs(BoundaryCondition::Inlet, BoundaryCondition::Outlet, BoundaryCondition::Symmetry, BoundaryCondition::Symmetry);  
	DoubleConeGrid grid(Nx, Ny, 1, 1, 1, 1, 25, 50, 2.5);  