SHCTools
========
######Version 1.0, 6-6-12######
#####A Matlab Toolbox for the Simulation, Analysis, and Design of Continuous Dynamical Systems  
Based on Networks of Stable Heteroclinic Channels (SHCs).#####


SHC network creation, modification, and visualization, and conversion.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```buildrho```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Create RHO matrix from SHC network structure.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```loadnet```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Load SHC network structure from XML or MAT file.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```plotactivity```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Plot simulated SHC network activity.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```plotnet```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Visualization of SHC network topology.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```savenet```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Save SHC network structure as XML or MAT file.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```shc_addnetwork```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Insert a sub-network into SHC network structure.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```shc_create```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- RHO matrix from parameters or transition matrix.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```shc_deletenetwork```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Remove a sub-network and its children.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```shc_initialize```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Create required SHC network structure fields.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```shc_reset```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Reset SHC network structure to default values.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```shc_validatenetwork```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Validate SHC network structure.

Lotka-Volterra SHC network tools.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```shc_lv_eigs```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Eigenvalues/vectors of Lotka-Volterra system.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```shc_lv_ic```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Generate initial conditions on the SHC manifold.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```shc_lv_integrate```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Solve Lotka-Volterra SDEs numerically.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```shc_lv_jacobian```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Jacobian of N-dimensional Lotka-Volterra system.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```shc_lv_ode```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- ODEs for N-dimensional Lotka-Volterra system.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```shc_lv_params```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Find RHO matrix parameters from description.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```shc_lv_symequilibria```&nbsp;&nbsp;&nbsp;&nbsp;- Solve for all equlibrium points symbolically.

Stone-Holmes distribution.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```stoneholmescdf```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Cummulative distribution function.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```stoneholmesfit```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Parameter estimates from data.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```stoneholmesinv```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Inverse cummulative distribution.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```stoneholmeslike```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Negative log-likelihood.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```stoneholmesmedian```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Median of distribution.   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```stoneholmesmode```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Mode (maximum value) of distribution.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```stoneholmespassagetime```&nbsp;- Mean passage time.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```stoneholmespdf```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Probability density function.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```stoneholmesrnd```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Random number generator.

XML validation files.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```shc.dtd```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Document Type Definition.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;```shc.xsd```&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- XML Schema.

--------

Andrew D. Horchler, adh9 @ case . edu  
Created: 1-4-12, Revision: 1.0, 6-6-12  

This version tested with Matlab 7.14.0.739 (R2012a)  
Mac OS X 10.6.8 (Build: 10K549), Java 1.6.0_31-b04-415-10M3635  
Compatibility maintained back through Matlab 7.4 (R2007a)  

Acknowledgment of support: This material is based upon work supported by the National Science Foundation under Grant No. 1065489. Disclaimer: Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.

Copyright &copy; 2012, Andrew D. Horchler  
All rights reserved.  

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * Neither the name of Case Western Reserve University nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ANDREW D. HORCHLER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.