Code for reproducing the results in [1]

File description:
- Trimming_1D_bases						Plots the 1D bases.
- Trimming_1D_Laplace_cond_B_splines_Jacobi			Conditioning of the Jacobi preconditioned matrices for the B-spline basis.
- Trimming_1D_Laplace_cond_B_splines_no_prec			Conditioning of the non-preconditioned matrices for the B-spline basis.
- Trimming_1D_Laplace_cond_B_splines				Conditioning of the system matrices for the B-spline basis for various preconditioning strategies.
- Trimming_1D_Laplace_cond_Lagrange_deflation			Conditioning of the deflation preconditioned matrices for the Lagrange basis and various values of gamma.
- Trimming_1D_Laplace_cond_Lagrange_Jacobi			Conditioning of the Jacobi preconditioned matrices for the Lagrange basis.
- Trimming_1D_Laplace_cond_Lagrange				Conditioning of the system matrices for the Lagrange basis and various preconditioning strategies.	
- Trimming_1D_Laplace_conv_B_splines				Conditioning of the system matrices and convergence of iterative solvers for the B-spline basis for various preconditioning strategies.
- Trimming_1D_Laplace_conv_Lagrange				Conditioning of the system matrices and convergence of iterative solvers for the Lagrange basis for various preconditioning strategies.
- Trimming_2D_Laplace_cond_B_splines_counter_ex_Schz1_Schz2	Counter-examples for the additive Schwarz preconditioner with smooth spline bases.
- Trimming_2D_Laplace_cond_B_splines_counter_ex_SIPIC		Counter-examples for the SIPIC preconditioner with spline bases.
- Trimming_2D_Laplace_cond_B_splines_deflation			Comparison of the deflation-based preconditioner with and without rank-reduction.
- Trimming_2D_Laplace_cond_B_splines_Jacobi			Verification of the scaling relations for the Jacobi preconditioner for spline bases and various trimming configurations.
- Trimming_2D_Laplace_cond_B_splines				Conditioning of the system matrices for the B-spline basis for various trimmed geometries and preconditioning strategies.
- Trimming_2D_Laplace_cond_Lagrange_counter_ex_SIPIC1		First counter-example for SIPIC with the Lagrange basis.
- Trimming_2D_Laplace_cond_Lagrange_counter_ex_SIPIC2		Second counter-example for SIPIC with the Lagrange basis.
- Trimming_2D_Laplace_cond_Lagrange_Jacobi			Verification of the scaling relations for the Jacobi preconditioner for Lagrange bases and various trimming configurations.
- Trimming_2D_Laplace_cond_Lagrange				Conditioning of the system matrices for the Lagrange basis for various trimmed geometries and preconditioning strategies.
- Trimming_2D_Laplace_conv_B_splines				Conditioning of the system matrices and convergence of iterative solvers for the B-spline basis for various preconditioning strategies and various trimmed geometries.
- Trimming_2D_Laplace_Example1_B_splines			L^2 projection problem on a rotated lattice structure. Discretization with smooth B-spline bases of degree 2 or 3. Comparison of the deflation-based preconditioner with and without rank-reduction.
- Trimming_2D_Laplace_Example2_Lagrange				Poisson problem on a extruded domain discretized with quadratic Lagrange bases. 
- Trimming_2D_Laplace_Example3_Bsplines				Wave equation on a waveguide-inspired spiky structure. Simulation of the first time step of an implicit solver. Performance comparison of various preconditioning strategies.
- Trimming_2D_Laplace_Example3_full_sim				Wave equation on a waveguide-inspired spiky structure. Complete simulation.
- Trimming_2D_ring						Comparison between a boundary-fitted and an unfitted discretization.
- Trimming_2D_rotated_square_plot_partition			Rotated square example: visualization of the active mesh and partitioning into cut and uncut regions.


Reference:
Y. Voet, M. Möller, P. Antolin and C. Vuik. Deflation-based preconditioning for immersed finite element methods and immersogeometric analysis.
