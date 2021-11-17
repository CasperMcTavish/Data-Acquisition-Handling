
kdata.bin
===============

Data for analysis within the Project.


kdata-small.bin
===============

Data for analysis within the Project. Used in 'Q1.py', 'Q2-7.py', 'Q8B+B-.py' & 'Q9 Lifetime.py'


functions.py
==============

Holds key functions relating to the side-band subtraction method for finding events. Used in both 'Q2-7.py' and 'Q8 B+B-.py'


Q1.py
==============

Contains all code required to complete the tasks in Question 1.
Outputs:
-Saves 5 graphs of the histograms of each different data component;
	-Invariant mass
	-Transverse momentum
	-Pseudorapidity
	-Charge
	-Lifetime
 of the B+- Meson data.


Q2-7.py
==============

Contains all code required to complete the tasks in Question 2-7.

Outputs:
-Graphs the mass distribution of the data as a histogram.

-Prints the Total Number of Events.

-Prints the estimate number of events within the signal region with background removed.

-Graphs the mass distribution as a line graph, with the peak, sidebands and signal regions marked.

-Graphs the esimated fit of the normalised PDF using curve_fit's optimal parameters.

-Prints the result of the minimisation of scipy.optimise.minimize with the parameters given, using the Nelder-Mead method.

-Graphs the estimated fit of the normalised PDF using our minimisation engine's optimal parameters.

-(OPTIONAL)Prints the covariance for the mean mass.

-Graphs the residuals of the estimated fit of the normalised PDF using the nelder-mead method, with the initial histogram.

-Graphs the residuals with the initial nelder-mead and histogram graph.

Inputs:
-Asks for whether or not the covariance of the mean should be calculated (Y/N) WARNING: this may take some time on older computers.

Q8 B+B-.py
==============

Contains all code required to complete the tasks in Question 8. It is identical in almost every way to 'Q2-7.py', but it loops once to allow for the separated B+ and B- data sets to have their masses and yields calculated. The separation of the B+- to B+ and B- data sets is done within this code.

Q9 Lifetime.py
==============

Contains all code required to complete the task in Question 9: Finding the mean lifetime of B+- mesons.

Outputs:
-Graphs the histogram for the lifetime distribution of the B+- Mesons

-Graphs the esimated fit of the normalised PDF using curve_fit's optimal parameters.

-Prints the result of the minimisation of scipy.optimise.minimize with the parameters given, using the Nelder-Mead method.

-Graphs the estimated fit of the normalised PDF using our minimisation engine's optimal parameters.

-(OPTIONAL)Prints the covariance for the rate parameter.

Inputs:
-Asks for whether or not the covariance of the rate parameter should be calculated (Y/N) WARNING: this may take some time on older computers.

NOTE: The code in its current form implements the PDF of equation 7 as described in section 6.2 of the report. If you want to alter this to equation 6, simply comment out line 101 and uncomment 102, this will allow you to calculate the values found in the report for equation 6.