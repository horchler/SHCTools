SHCTools
========
#####A Matlab Toolbox for the Simulation, Analysis, and Design of Continuous Dynamical Systems Based&nbsp;on&nbsp;Networks&nbsp;of&nbsp;Stable&nbsp;Heteroclinic&nbsp;Channels&nbsp;(SHCs).#####
######Version 1.1, 6-8-13######

SHC network creation, modification, and visualization, and conversion.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```buildrho```](https://github.com/horchler/SHCTools/blob/master/SHCTools/buildrho.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Create RHO matrix from SHC network structure.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```loadnet```](https://github.com/horchler/SHCTools/blob/master/SHCTools/loadnet.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Load SHC network structure from XML or MAT file.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```plotactivity```](https://github.com/horchler/SHCTools/blob/master/SHCTools/plotactivity.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Plot simulated SHC network activity.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```plotnet```](https://github.com/horchler/SHCTools/blob/master/SHCTools/plotnet.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Visualization of SHC network topology.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```savenet```](https://github.com/horchler/SHCTools/blob/master/SHCTools/savenet.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Save SHC network structure as XML or MAT file.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_addnetwork```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_addnetwork.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Insert a sub-network into SHC network structure.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_create```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_create.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Create SHC network structure from parameters.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_deletenetwork```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_deletenetwork.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Remove a sub-network and its children.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_initialize```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_initialize.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Create required SHC network structure fields.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_reset```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_reset.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Reset SHC network structure to default values.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_validatenetwork```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_validatenetwork.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Validate SHC network structure.

Lotka-Volterra SHC network tools.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_eigs```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_eigs.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Eigenvalues/vectors of Lotka-Volterra system.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_epsilonfit```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_epsilonfit.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Nonlinear-compensated (fitted) noise magnitude.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_globalpassagetime```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_globalpassagetime.m)&nbsp;&nbsp;- Find passage times from SHC network structure.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_ic```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_ic.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Generate initial conditions on the SHC manifold.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_integrate```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_integrate.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Solve Lotka-Volterra SDEs numerically.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_iscycle```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_iscycle.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Check if Lotka-Volterra system is an SHC cycle.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_isstable```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_isstable.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Check stability of Lotka-Volterra system nodes.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_jacobian```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_jacobian.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Jacobian of N-dimensional Lotka-Volterra system.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_lambda_us```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_lambda_us.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Dominant unstable and stable eigenvalues.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_meanperiod```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_meanperiod.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Simulate system to find average mean period.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_mintransitiontime```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_mintransitiontime.m)&nbsp;&nbsp;- Estimate minimum transition times of network.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_ode```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_ode.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- ODEs for N-dimensional Lotka-Volterra system.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_params```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_params.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Find RHO matrix parameters from description.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_passagetime```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_passagetime.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Mean first passage times of SHC network.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_stability```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_stability.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Saddle values of Lotka-Volterra system nodes.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_symequilibria```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_symequilibria.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Solve for all equlibrium points symbolically.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_symparams```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_symparams.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Create symbolic SHC parameters with assumptions.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_taufit```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_taufit.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Nonlinear-compensated (fitted) mean period.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_transitiontime```](https://github.com/horchler/SHCTools/blob/master/SHCTools/shc_lv_transitiontime.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Inter-passage transition times of SHC network.

Stone-Holmes distribution.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```stoneholmescdf```](https://github.com/horchler/SHCTools/blob/master/SHCTools/stoneholmes/stoneholmescdf.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Cummulative distribution function.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```stoneholmeschi2gof```](https://github.com/horchler/SHCTools/blob/master/SHCTools/stoneholmes/stoneholmeschi2gof.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Chi-squared goodness-of-fit test.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```stoneholmesdemo```](https://github.com/horchler/SHCTools/blob/master/SHCTools/stoneholmes/stoneholmesdemo.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Demonstrate Stone-Holmes distribution functions.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```stoneholmesfit```](https://github.com/horchler/SHCTools/blob/master/SHCTools/stoneholmes/stoneholmesfit.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Parameter estimates from data.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```stoneholmesinv```](https://github.com/horchler/SHCTools/blob/master/SHCTools/stoneholmes/stoneholmesinv.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Inverse cummulative distribution.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```stoneholmesinvpassagetime```](https://github.com/horchler/SHCTools/blob/master/SHCTools/stoneholmes/stoneholmesinvpassagetime.m)&nbsp;- Noise as a function of mean passage time.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```stoneholmeskstest```](https://github.com/horchler/SHCTools/blob/master/SHCTools/stoneholmes/stoneholmeskstest.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Kolmogorov-Smirnov goodness-of-fit test.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```stoneholmeslike```](https://github.com/horchler/SHCTools/blob/master/SHCTools/stoneholmes/stoneholmeslike.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Negative log-likelihood.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```stoneholmesmedian```](https://github.com/horchler/SHCTools/blob/master/SHCTools/stoneholmes/stoneholmesmedian.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Median of distribution.   
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```stoneholmesmode```](https://github.com/horchler/SHCTools/blob/master/SHCTools/stoneholmes/stoneholmesmode.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Mode (maximum value) of distribution.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```stoneholmespassagetime```](https://github.com/horchler/SHCTools/blob/master/SHCTools/stoneholmes/stoneholmespassagetime.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Mean passage time.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```stoneholmespdf```](https://github.com/horchler/SHCTools/blob/master/SHCTools/stoneholmes/stoneholmespdf.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Probability density function.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```stoneholmesrnd```](https://github.com/horchler/SHCTools/blob/master/SHCTools/stoneholmes/stoneholmesrnd.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Random number generator.

XML validation files.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc.dtd```](https://github.com/horchler/SHCTools/blob/master/SHCTools/xml/shc.dtd)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Document Type Definition.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc.xsd```](https://github.com/horchler/SHCTools/blob/master/SHCTools/xml/shc.xsd)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- XML Schema.  
&nbsp;  

--------

References  
&nbsp;  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;V.S. Afraimovich, M.I. Rabinovich, and P. Varona, &#8220;Heteroclinic contours in neural ensembles and the winnerless competition,&#8221;  
&nbsp;&nbsp;&nbsp;&nbsp; *Int. J. Bifurcation Chaos*, Vol. 14, 2004, pp. 1195&ndash;1208. [[PDF Preprint](http://arxiv.org/abs/nlin/0304016)] [[http://dx.doi.org/10.1142/S0218127404009806](http://dx.doi.org/10.1142/S0218127404009806)]  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;M.I. Rabinovich, R. Huerta, P. Varona, and V.S. Afraimovich, &#8220;Transient Cognitive Dynamics, Metastability, and Decision  
&nbsp;&nbsp;&nbsp;&nbsp;Making,&#8221; *PLoS Comp. Biol.*, Vol. 4, No. 5, 2008, e1000072. [[http://dx.doi.org/10.1371/journal.pcbi.1000072](http://dx.doi.org/10.1371/journal.pcbi.1000072)]  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;J.W. Reyn, &#8220;A Stability Criterion for Separatrix Polygons in the Phase Plane,&#8221; *Nieuw Archief Voor Wiskunde (3)*,  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Vol. 27, 1979, pp. 238&ndash;254. [[PDF](https://github.com/horchler/SHCTools/blob/master/references/reyn_criterion_1979.pdf?raw=true)]  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;K.M. Shaw, Y.-M. Park, H.J. Chiel, and P.J. Thomas, &#8220;Phase Resetting in an Asymptotically Phaseless System:  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;On the Phase Response of Limit Cycles Verging on a Heteroclinic Orbit,&#8221; *SIAM J. Appl. Dyn. Syst.*,  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Vol. 11, No. 1, 2012, pp. 350&ndash;391. [[http://dx.doi.org/10.1137/110828976](http://dx.doi.org/10.1137/110828976)]  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;E. Stone and P. Holmes, &#8220;Random Perturbations of Heteroclinic Attractors,&#8221; *SIAM J. Appl. Math.*,  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;Vol. 50, No. 3, 1990, pp. 726&ndash;743. [[http://jstor.org/stable/2101884](http://jstor.org/stable/2101884)]  
&nbsp;  

--------

Andrew D. Horchler, *adh9 @ case . edu*, [biorobots.case.edu](http://biorobots.case.edu/)  
Created: 1-4-12, Revision: 1.1, 6-8-13  

This version tested with Matlab 8.0.0.783 (R2012b)  
Mac OS X 10.8.4 (Build: 12E55), Java 1.6.0_45-b06-451-11M4406  
Compatibility maintained back through Matlab 7.4 (R2007a)  
&nbsp;  

--------

Acknowledgment of support: This material is based upon work supported by the [National Science Foundation](http://www.nsf.gov/) under  
[Grant No.&nbsp;1065489](http://www.nsf.gov/awardsearch/showAward.do?AwardNumber=1065489). Disclaimer: Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.  
&nbsp;  

Copyright &copy; 2012&ndash;2013, Andrew D. Horchler  
All rights reserved.  

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * Neither the name of Case Western Reserve University nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ANDREW D. HORCHLER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.