%SHCTools
%  Version 1.0, 5-25-12
%
%  SHC network creation, modification, and visualization, and conversion.
%      buildrho               - Create RHO matrix from SHC network structure.
%      loadnet                - Load SHC network structure from XML or MAT file.
%      plotactivity           - Plot simulated SHC network activity.
%      plotnet                - Visualization of SHC network topology.
%      savenet                - Save SHC network structure as XML or MAT file.
%      shc_addnetwork         - Insert a sub-network into SHC network structure.
%      shc_create             - RHO matrix from parameters or transition matrix.
%      shc_deletenetwork      - Remove a sub-network and its children.
%      shc_initialize    	  - Create required SHC network structure fields.
%      shc_reset              - Reset SHC network structure to default values.
%      shc_validatenetwork	  - Validate SHC network structure.
%
%  Lotka-Volterra SHC network tools.
%      shc_lv_eigs            - Eigenvalues/vectors of Lotka-Volterra system.
%      shc_lv_ic              - Generate initial conditions on the SHC manifold.
%      shc_lv_integrate       - Solve Lotka-Volterra SDEs numerically.
%      shc_lv_jacobian        - Jacobian of N-dimensional Lotka-Volterra system.
%      shc_lv_ode         	  - ODEs for N-dimensional Lotka-Volterra system.
%      shc_lv_params          - Find RHO matrix parameters from description.
%      shc_lv_symequilibria   - Solve for all equlibrium points symbolically.
%
%  Stone-Holmes distribution.
%      stoneholmescdf    	  - Cummulative distribution function.
%      stoneholmesfit         - Parameter estimates from data.
%      stoneholmesinv         - Inverse cummulative distribution.
%      stoneholmeslike   	  - Negative log-likelihood.
%      stoneholmesmedian 	  - Median of distribution.
%      stoneholmesmode        - Mode (maximum value) of distribution.
%      stoneholmespassagetime - Mean passage time.
%      stoneholmespdf         - Probability density function.
%      stoneholmesrnd   	  - Random number generator.
%
%  XML validation files.
%      shc.dtd                - Document Type Definition.
%      shc.xsd                - XML Schema.

%  Tested with Matlab 7.14.0.739 (R2012a)
%  Mac OS X 10.6.8 (Build: 10K549), Java 1.6.0_31-b04-415-10M3635
%  Compatibility maintained back through Matlab 7.4 (R2007a)

%  Andrew D. Horchler, adh9 @ case . edu
%  Created: 1-4-12, Revision: 1.0, 5-25-12


%  Acknowledgment of support: This material is based upon work supported by the
%  National Science Foundation under Grant No. 1065489. Disclaimer: Any
%  opinions, findings, and conclusions or recommendations expressed in this
%  material are those of the author(s) and do not necessarily reflect the views
%  of the National Science Foundation.


%  Copyright � 2012, Andrew D. Horchler
%  All rights reserved.
%
%  Redistribution and use in source and binary forms, with or without
%  modification, are permitted provided that the following conditions are met:
%   * Redistributions of source code must retain the above copyright
%     notice, this list of conditions and the following disclaimer.
%   * Redistributions in binary form must reproduce the above copyright
%     notice, this list of conditions and the following disclaimer in the
%     documentation and/or other materials provided with the distribution.
%   * Neither the name of Case Western Reserve University nor the
%     names of its contributors may be used to endorse or promote products
%     derived from this software without specific prior written permission.
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