## Description

This is a package for detecting the convexity of noisy functions. More specifically, given a sample function that can only be evaluated with noise (e.g. output of a simulation model), the MAIN script is able to return estimates of the posterior probability that the function is convex after each iteration of sampling. 

## To use

In sampleFcn.m there are already some test functions we used in the paper. You can either uncomment out your interested test functions, or put your own noisy function in there. When you do that, please follow the format in the script to make sure sampleFcn.m returns all its required outputs. Then run the MAIN.m script and get your estimates on how convex the function is!

The MATLAB console will prompt you with some configurations of the algorithm. We have provided the default configurations in paranthesis. For details, please read our paper below in the "Paper reference" section.

The final outputs of MAIN.m will be put under the "output" folder (please make sure this folder exists). All outputs will start with the month_day_year_hour_minute string of when the script stopped running. The outputs include a snapshot of the workspace at the end of the run, a text log of all the console outputs, and three plots for the estimator (estimated value, efficiency, and time vs. iteration).

All .m files are commented with heading and in-line explanations. The heading provides the user with information including inputs, outputs, and which functions are called by this one.

## Installation

You can run the package with any version of Matlab (R2013 or later recommended). There are two versions of the package -- one requires the additional installation of cvx and a compatible solver, e.g. Gurobi, and the other one that only needs Matlab. We recommend cvx + Gurobi because of the robustness of its linear program solver, as compared to linprog in Matlab. If you choose to stick with the basic version that only requires Matlab, you can skip the rest of this section.

Both cvx and Gurobi offer free academic licenses if you have a ".edu" mail address. To install cvx, visit http://cvxr.com/cvx/, which has detailed instruction of how to use Gurobi with cvx under http://web.cvxr.com/cvx/doc/gurobi.html#gurobi. The Gurobi academic license can be obtained at https://user.gurobi.com/download/licenses/free-academic.

After you have installed both cvx and Gurobi, move the files under "cvx_Gurobi" to the same directory as MAIN.m and other files and overwrite the files with the same name. Then you can run MAIN.m and try it out! =)

## Paper reference

The internal mechanism of this package is explained in our paper "Estimating the Probability that a Function Observed with Noise is Convex", which has been submitted to INFORMS Journal of Computing.

## Contributor

By Nanjing Jian (nj227@cornell.edu) at Cornell University
