%SHCTools
%  Version 1.1, 8-30-13
%
%  SHC network creation, modification, and visualization, and conversion.
%      buildrho               - Create RHO matrix from SHC network structure.
%      loadnet                - Load SHC network structure from XML or MAT file.
%      plotactivity           - Plot simulated SHC network activity.
%      plotnet                - Visualization of SHC network topology.
%      savenet                - Save SHC network structure as XML or MAT file.
%      shc_addnetwork         - Insert a sub-network into SHC network structure.
%      shc_create             - Create SHC network structure from parameters.
%      shc_deletenetwork      - Remove a sub-network and its children.
%      shc_initialize    	  - Create required SHC network structure fields.
%      shc_reset              - Reset SHC network structure to default values.
%      shc_validatenetwork	  - Validate SHC network structure.
%
%  Lotka-Volterra SHC network tools.
%      shc_lv_eigs            - Eigenvalues/vectors of Lotka-Volterra system.
%      shc_lv_epsilonfit      - Nonlinear-compensated (fitted) noise magnitude.
%      shc_lv_globalpassagetime - Find passage times from SHC network structure.
%      shc_lv_ic              - Find initial conditions close to SHC manifold.
%      shc_lv_integrate       - Solve Lotka-Volterra SDEs numerically.
%      shc_lv_iscycle         - Check if Lotka-Volterra system is an SHC cycle.
%      shc_lv_isstable        - Check stability of Lotka-Volterra system nodes.
%      shc_lv_jacobian        - Jacobian of N-dimensional Lotka-Volterra system.
%      shc_lv_lambda_us       - Dominant unstable and stable eigenvalues.
%      shc_lv_meanperiod      - Simulate system to find average mean period.
%      shc_lv_mintransitiontime - Estimate minimum transition times of network.
%      shc_lv_ode         	  - ODEs for N-dimensional Lotka-Volterra system.
%      shc_lv_params          - Find RHO matrix parameters from description.
%      shc_lv_passagetime     - Mean first passage times of SHC network.
%      shc_lv_stability       - Saddle values of Lotka-Volterra system nodes.
%      shc_lv_symequilibria   - Solve for all equlibrium points symbolically.
%      shc_lv_symparams       - Create symbolic SHC parameters with assumptions.
%      shc_lv_taufit          - Nonlinear-compensated (fitted) mean period.
%      shc_lv_transitiontime  - Inter-passage transition times of SHC network.
%
%  Stone-Holmes distribution.
%      stoneholmescdf    	  - Cummulative distribution function.
%      stoneholmeschi2gof     - Chi-squared goodness-of-fit test.
%      stoneholmescv          - Coefficient of variation of samples.
%      stoneholmesdemo        - Demonstrate Stone-Holmes distribution functions.
%      stoneholmesfit         - Parameter estimates from data.
%      stoneholmesinv         - Inverse cummulative distribution.
%      stoneholmesinvpassagetime - Noise as a function of mean passage time.
%      stoneholmeskstest      - Kolmogorov-Smirnov goodness-of-fit test.
%      stoneholmeslike   	  - Negative log-likelihood.
%      stoneholmesmedian 	  - Median of distribution.
%      stoneholmesmode        - Mode (maximum value) of distribution.
%      stoneholmespassagetime - Mean passage time.
%      stoneholmespdf         - Probability density function.
%      stoneholmesrnd   	  - Random number generator.
%      stoneholmesstat        - Mean and variance of distribution.
%      stoneholmesvar         - Variance of Stone-Holmes distribution samples.
%
%  XML validation files.
%      shc.dtd                - Document Type Definition.
%      shc.xsd                - XML Schema.

%  This version tested with Matlab 8.1.0.604 (R2013a)
%  Mac OS X 10.8.4 (Build: 12E55), Java 1.6.0_51-b11-457-11M4509
%  Compatibility maintained back through Matlab 7.4 (R2007a)

%  Andrew D. Horchler, adh9 @ case . edu
%  Created: 1-4-12, Revision: 1.1, 8-30-13


%  Acknowledgment of support: This material is based upon work supported by the
%  National Science Foundation under Grant No. 1065489. Disclaimer: Any
%  opinions, findings, and conclusions or recommendations expressed in this
%  material are those of the author(s) and do not necessarily reflect the views
%  of the National Science Foundation.


%  Copyright � 2012-2013, Andrew D. Horchler
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