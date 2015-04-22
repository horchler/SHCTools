SHCTools &ndash; Horchler, et al. 2015
========
#####Matlab Tools for the Simulation, Analysis, and Design of Stable Heteroclinic Channel (SHC) cycles.#####
######Version 1.0, 4-5-15######

This branch of the SHCTools repository maintains the code for the paper *Designing Responsive Pattern Generators: Stable Heteroclinic Channel Cycles for Modeling and Control* (Horchler, et al. 2015). The code in this branch is meant for readability and is specific to the subject the paper: mean first passage time and Lotka-Volterra SHC cycles. Many optimizations and features found in the [master SHCTools branch](https://github.com/horchler/SHCTools/tree/master) are not present. The functions here will remain static and will only be updated to fix bugs.

&nbsp;  
&nbsp;  
How to install:  
 1. Download and expand the *[SHCTools-SHC-Design-Paper.zip](https://github.com/horchler/SHCTools/archive/SHC-Design-Paper.zip)* ZIP archive of the repository.  
 2. Move the resultant *SHCTools-SHC-Design-Paper* folder to the desired permanent location.  
 3. In Matlab, navigate to *SHCTools-SHC-Design-Paper/SHCTools-Horchler_etal_2015/*. Add this folder to your Matlab search ```path``` if desired.  
 4. Code to produce the basic form of the figures in the paper can be found in *SHCTools-SHC-Design-Paper/Figures-Horchler_etal_2015/*. These M-files use the functions in the *SHCTools-Horchler_etal_2015/* folder. Either move these M-files or set the search path as described in step&nbsp;3.  
 5. Minor edits and bug reports and fixes can be submitted by [filing an issue](https://github.com/horchler/SHCTools/issues) or via email (please be sure to to refer to the *SHC-Design-Paper* branch). To add new functionality or make propose major changes, please [fork the repository](https://help.github.com/articles/fork-a-repo) and work from the main branch. Any new features should be accompanied by some means of testing. Email or file an issue if you have any questions.  
&nbsp;  

--------

Lotka-Volterra SHC cycle tools.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_createcycle```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_lv_createcycle.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Create Lotka-Volterra SHC cycle connection matrix.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_demo```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_lv_demo.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Eigenvalues/vectors of Lotka-Volterra system.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_eigs```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_lv_eigs.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Basic demonstration of SHC cycle design and simulation.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_epsilonfit```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_lv_epsilonfit.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Nonlinear-compensated (fitted) noise magnitudes.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_ic```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_lv_ic.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Find initial conditions on the SHC manifold.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_integrate```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_lv_integrate.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Solve Lotka-Volterra SDEs numerically.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_jacobian```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_lv_jacobian.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Jacobian of N-dimensional Lotka-Volterra system.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_lambda_us```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_lv_lambda_us.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Dominant unstable and stable eigenvalues.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_meanperiod```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_lv_meanperiod.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Simulate system to measure mean sub-periods.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_mintransitiontime```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_lv_mintransitiontime.m)&nbsp;&nbsp;&nbsp;&nbsp;- Minimum transition times between nodes.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_neighborhood```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_lv_neighborhood.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Get/set scaled neighborhood size for SHC cycle.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_ode```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_lv_ode.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- ODEs for N-dimensional Lotka-Volterra system.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_params```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_lv_params.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Connection matrix parameters from specification.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_passagetime```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_lv_passagetime.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Mean first passage times of SHC cycle.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_lv_taufit```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_lv_taufit.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Nonlinear-compensated (fitted) mean sub-periods.  

Passage time.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_epsilonfit```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_epsilonfit.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Nonlinear-compensated (fitted) noise magnitude.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_passagetime```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_passagetime.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Mean first passage time of a linearized system.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```shc_taufit```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/shc_taufit.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Nonlinear-compensated (fitted) mean period.  

Utilities.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```coefvar```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/coefvar.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Coefficient of variation for data samples.  
&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;[```wrightOmegaq```](https://github.com/horchler/SHCTools/blob/SHC-Design-Paper/SHCTools-Horchler_etal_2014/wrightOmegaq.m)&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;- Wright *&omega;* function, a solution to *&omega;*+ln(*&omega;*) = *Z*.  
&nbsp;  

--------

Key references.  
 1. V.S. Afraimovich, M.I. Rabinovich, and P. Varona, &#8220;Heteroclinic contours in neural ensembles and the winnerless competition principle,&#8221; *Int. J. Bifurcation Chaos*, Vol. 14, 2004, pp. 1195&ndash;1208. [[PDF Preprint](http://arxiv.org/abs/nlin/0304016)] [[http://dx.doi.org/10.1142/S0218127404009806](http://dx.doi.org/10.1142/S0218127404009806)]  
 2. A.D. Horchler, K.A. Daltorio, H.J. Chiel, and R.D Quinn &#8220;Designing Responsive Pattern Generators: Stable Heteroclinic Channel Cycles for Modeling and Control,&#8221; *Bioinspir. Biomim.*, Vol. 10, No. 2, 2015, pp. 1&ndash;16. [[http://dx.doi.org/10.1088/1748-3190/10/2/026001](http://dx.doi.org/10.1088/1748-3190/10/2/026001)]
 3. M.I. Rabinovich, R. Huerta, P. Varona, and V.S. Afraimovich, &#8220;Transient Cognitive Dynamics, Metastability, and Decision Making,&#8221; *PLoS Comp. Biol.*, Vol. 4, No. 5, 2008, e1000072. [[http://dx.doi.org/10.1371/journal.pcbi.1000072](http://dx.doi.org/10.1371/journal.pcbi.1000072)]  
 4. J.W. Reyn, &#8220;A Stability Criterion for Separatrix Polygons in the Phase Plane,&#8221; *Nieuw Archief Voor Wiskunde (3)*, Vol. 27, 1979, pp. 238&ndash;254. [[PDF](https://github.com/horchler/SHCTools/blob/master/references/reyn_criterion_1979.pdf?raw=true)]  
 5. K.M. Shaw, Y.-M. Park, H.J. Chiel, and P.J. Thomas, &#8220;Phase Resetting in an Asymptotically Phaseless System: On the Phase Response of Limit Cycles Verging on a Heteroclinic Orbit,&#8221; *SIAM J. Appl. Dyn. Syst.*, Vol. 11, No. 1, 2012, pp. 350&ndash;391. [[http://dx.doi.org/10.1137/110828976](http://dx.doi.org/10.1137/110828976)]  
 6. E. Stone and P. Holmes, &#8220;Random Perturbations of Heteroclinic Attractors,&#8221; *SIAM J. Appl. Math.*, Vol. 50, No. 3, 1990, pp. 726&ndash;743. [[http://jstor.org/stable/2101884](http://jstor.org/stable/2101884)]  
&nbsp;  

--------

Andrew D. Horchler, *adh9 @ case . edu*, [biorobots.case.edu](http://biorobots.case.edu/)  
Created: 1-4-12, Version: 1.0, 4-5-15  

This version tested with Matlab 8.5.0.173394 (R2015a)  
Mac OS X 10.10.2 (Build: 14C1514), Java 1.7.0_60-b19 
Compatibility maintained back through Matlab 7.4 (R2007a)  
&nbsp;  

--------

Acknowledgment of support: This material is based upon work supported by the [National Science Foundation](http://www.nsf.gov/) under [Grant No.&nbsp;1065489](http://www.nsf.gov/awardsearch/showAward.do?AwardNumber=1065489). Disclaimer: Any opinions, findings, and conclusions or recommendations expressed in this material are those of the author(s) and do not necessarily reflect the views of the National Science Foundation.  
&nbsp;  

Copyright &copy; 2012&ndash;2015, Andrew D. Horchler  
All rights reserved.  

Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
 * Neither the name of Case Western Reserve University nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL ANDREW D. HORCHLER BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.