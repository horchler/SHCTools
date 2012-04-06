%SHCTOOLS
%   Version 1.0, 4-6-12
%
%   SHC network creation, modification, and visualization, and conversion.
%       buildrho             - Create RHO matrix from SHC network structure.
%       loadnet              - Load SHC network structure from XML or MAT file.
%       plotactivity         - Plot simulated SHC network activity.
%       plotnet              - Visualization of SHC network topology.
%       savenet              - Save SHC network structure as XML or MAT file.
%       shc_addnetwork       - Insert a sub-network into SHC network structure.
%       shc_create           - RHO matrix from parameters or transition matrix.
%       shc_deletenetwork    - Remove a sub-network and its children.
%       shc_initialize    	 - Create required fields for SHC network structure.
%       shc_reset            - Reset SHC network structure to default values.
%       shc_validatenetwork	 - Validate SHC network structure.
%
%   Lotka-Volterra SHC network tools.
%       shc_lv_eigs          - Eigenvalues of Lotka-Volterra system.
%       shc_lv_integrate     - Solve stochastic Lotka-Volterra ODEs.
%       shc_lv_jacobian      - Jacobian of Lotka-Volterra system.
%       shc_lv_ode         	 - ODEs for N-dimensional Lotka-Volterra system.
%       shc_lv_symequilibria - Symbolic equlibrium points and eigenvalues.
%
%   Stone-Holmes distribution.
%       stoneholmescdf    	 - Cummulative distribution function.
%       stoneholmesfit       - Parameter estimates from data.
%       stoneholmesinv       - Inverse cummulative distribution.
%       stoneholmeslike   	 - Negative log-likelihood.
%       stoneholmesmedian 	 - Median of distribution.
%       stoneholmesmode      - Mode (maximum value) of distribution.
%       stoneholmespdf       - Probability density function.
%       stoneholmesrnd   	 - Random number generator.
%
%   XML validation files.
%       shc.dtd              - Document Type Definition.
%       shc.xsd              - XML Schema.

%   Tested with Matlab 7.13.0.564 (R2011b)
%   Mac OS X 10.6.8 (Build: 10K549), Java 1.6.0_31-b04-413-10M3623

%   Andrew D. Horchler, adh9@case.edu, Created 1-4-12
%   Revision: 1.0, 4-6-12