%SHCTools
%  Version 1.3, 4-5-14
%
%  Lotka-Volterra SHC network tools.
%      shc_lv_createcycle     - Create Lotka-Volterra SHC cycle from parameters.
%      shc_lv_eigs            - Eigenvalues/vectors of Lotka-Volterra system.
%      shc_lv_epsilonfit      - Nonlinear-compensated (fitted) noise magnitudes.
%      shc_lv_equilibria      - Solve for all equlibrium points symbolically.
%      shc_lv_ic              - Find initial conditions close to SHC manifold.
%      shc_lv_integrate       - Solve Lotka-Volterra SDEs numerically.
%      shc_lv_iscycle         - Check if Lotka-Volterra system is an SHC cycle.
%      shc_lv_isstable        - Check stability of Lotka-Volterra system nodes.
%      shc_lv_isuniform       - Check if Lotka-Volterra system is uniform SHC.
%      shc_lv_jacobian        - Jacobian of N-dimensional Lotka-Volterra system.
%      shc_lv_lambda_us       - Dominant unstable and stable eigenvalues.
%      shc_lv_meanperiod      - Simulate system to find average mean period.
%      shc_lv_mintransitiontime - Estimate minimum transition times of network.
%      shc_lv_neighborhood    - Scaled neighborhood-size for Lotka-Volterra SHC.
%      shc_lv_ode         	  - ODEs for N-dimensional Lotka-Volterra system.
%      shc_lv_params          - Connection matrix parameters from specification.
%      shc_lv_passagetime     - Mean first passage times of SHC network.
%      shc_lv_stability       - Saddle values of Lotka-Volterra system nodes.
%      shc_lv_symparams       - Create symbolic SHC parameters with assumptions.
%      shc_lv_taufit          - Nonlinear-compensated (fitted) mean sub-periods.
%
%  Stone-Holmes distribution.
%      stoneholmescdf    	  - Cummulative distribution function.
%      stoneholmeschi2gof     - Chi-squared goodness-of-fit test.
%      stoneholmescv          - Coefficient of variation of samples.
%      stoneholmesdemo        - Demonstrate Stone-Holmes distribution functions.
%      stoneholmesfit         - Parameter estimates from data.
%      stoneholmesinv         - Inverse cummulative distribution.
%      stoneholmesinvpassagetime - Noise as a function of mean passage time.
%      stoneholmesksdensity   - Kernel density estimation for distribution.
%      stoneholmeskstest      - Kolmogorov-Smirnov goodness-of-fit test.
%      stoneholmeslike   	  - Negative log-likelihood.
%      stoneholmesmedian 	  - Median of distribution.
%      stoneholmesmode        - Mode (maximum value) of distribution.
%      stoneholmespassagetime - Mean passage time.
%      stoneholmespdf         - Probability density function.
%      stoneholmesrnd   	  - Random number generator.
%      stoneholmesstat        - Mean and variance of distribution.
%      stoneholmesvar         - Variance of Stone-Holmes distribution samples.

%  This version tested with Matlab 8.5.0.173394 (R2015a)
%  Mac OS X 10.10.2 (Build: 14C1514), Java 1.7.0_60-b19
%  Compatibility maintained back through Matlab 7.4 (R2007a)

%  Andrew D. Horchler, adh9 @ case . edu
%  Created: 1-4-12, Version: 1.3, 4-5-14


%  Acknowledgment of support: This material is based upon work supported by the
%  National Science Foundation under Grant No. 1065489. Disclaimer: Any
%  opinions, findings, and conclusions or recommendations expressed in this
%  material are those of the author(s) and do not necessarily reflect the views
%  of the National Science Foundation.


%  Copyright © 2012-2015, Andrew D. Horchler
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are met:
%   * Redistributions of source code must retain the above copyright notice,
%     this list of conditions and the following disclaimer.
%   * Redistributions in binary form must reproduce the above copyright notice,
%     this list of conditions and the following disclaimer in the documentation
%     and/or other materials provided with the distribution.
%   * Neither the name of Case Western Reserve University nor the names of its
%     contributors may be used to endorse or promote products derived from this
%     software without specific prior written permission.
%
%  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
%  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
%  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
%  ARE DISCLAIMED. IN NO EVENT SHALL ANDREW D. HORCHLER BE LIABLE FOR ANY
%  DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
%  (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
%  LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
%  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
%  (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
%  SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.