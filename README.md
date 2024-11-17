This code is an exercise for an aerodynamics course, where it uses the method of collocation to solve the Prandtl-Lanchester lifting-line equation. The goal is to compute the lift coefficient (Cl), drag coefficient (Cd), and efficiency (Cl/Cd) of a wing based on various design parameters, including taper ratio (TR), aspect ratio (AR), and twist. The code systematically varies these parameters, solves the governing equations using a linear solver, and outputs the results to a text file.

### Key components:
1. **Build Functions**: These generate necessary vectors (e.g., theta, ClAlpha) and matrices (e.g., A, C vector) based on the wing's parameters.
2. **Linear Solver**: The system of equations is solved using LU decomposition.
3. **Results Calculation**: For each combination of parameters (TR, AR, twist), the lift and drag coefficients are computed, and the efficiency is derived.
4. **Output**: The results are saved in a well-structured table in an output file.

This code serves as a tool to explore the aerodynamic efficiency of different wing configurations, providing insights into the relationship between geometric parameters and aerodynamic performance.