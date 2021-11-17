import numpy as np
import matplotlib.pyplot as plt
import math as m
import scipy.optimize as opt
from scipy.signal import find_peaks
from scipy import stats
import functions as fnc
import matplotlib.patches as mpatches



# import data
# xtime = np.loadtxt(sys.argv[1])
f = open("kdata-small.bin","r")
datalist = np.fromfile(f,dtype=np.float32)

# number of events
nevent = len(datalist)/5
xdata = np.split(datalist,nevent)


# make lists for all events
xtime = []

# For finding N values within signal region (Question 3)
sigvals = 0


# Giving each component of xfull a separate number of bins for better data evaluation, decided via eyeballing "suitable bins"
bns = 150
# Same as above, but with labels, titles and savefiles of graphs
label = 'Lifetime (ps)'
title = 'Lifetimes of $B^{\pm} Mesons$ against Count'
svfile = 'lifetimeplotb.pdf'



# Collect lifetime data
for i in range(0,int(nevent)):
    xtime.append(xdata[i][4])
#    This is for seeing all values of xtime
#    if i < int(nevent):
#        print(xtime[i])
print("No. of array components: " + str(len(xdata[0])))

#xtime = np.array(xtime)/1000

# Plot histograms
# Find binsize for each component
binsize = (max(xtime) - min(xtime))/bns
# Show histograms and save each one with unique title and labels
plt.hist(xtime, bins = bns)
plt.ylabel('Count')
plt.xlabel(label)
plt.title(title)
plt.show()



# Finding mass of B+ by determining the bins with the highest number of entries in the peak regions.

# Histogram of xtime, removed from plot to allow for more modification
xtimehist = np.histogram(xtime,bins = bns)
# 0th Component - Histogram values
# 1st Component - Bin edges

# Bin edges, first value sliced off so that graph works
binedges = np.around(xtimehist[1]).astype(int)[1:]
# Bin centers, no value sliced off. But keeping binedges as backup
bin_centers = 0.5*(xtimehist[1][1:]+xtimehist[1][:-1])


#Finding Lifetime and Variance via derivatives of NLL

tau = (1/len(xtime))*np.sum(xtime)

print("LIFETIME (PS)")
print(tau)
print("VARIANCE (PS)")
print(np.power(1/tau,2))


#

# Guesses at suitable values for each variable in PDF
# Composite Fraction#
lmbda = 0.6
c = 0

# new normalised histogram for PDF
xtimehistnorm = np.histogram(xtime,bins = bns, density = True)

# Defining x and y points to compare to
x = bin_centers
y = xtimehistnorm[0]

# Probability distribution function for exponential
def PDF(xp, lmb):
    # Commented out exponential is the fullly normalised derivative, but doesn't give an accurate or useful result.
    # The currently used exponential function gives the value same as tau above, but the fit looks wrong. Will be discussed in report
    expon = lmb*np.exp(-lmb*xp)
    #expon = np.exp(-lmb*xp)*(lmb/(np.exp(-lmb*min(xp))-np.exp(-lmb*max(xp))))
    return expon

# Negative Log Likelihood
def NLL(params, args):
    xi = args
    lmbl = params

    # METHOD 1
    pdfresult = PDF(xi, lmbl)
    #print("PDFresult")
    #print(pdfresult)
    L = -np.sum(np.log(pdfresult))
    return L

# constraint on f to keep it between 0 and 1
def conf(params):
    return 1-params[0]

def conf2(params):
    return params[0]-1


# TO DEMONSTRATE THAT OUR PDF IS CORRECT ###############

# Least squares fitting method to make PDF match data best
popt, pcov = opt.curve_fit(PDF, bin_centers, xtimehistnorm[0])
# 0th component - optimal parameters for PDF
# 1st component - variance of said parameters
print("Optimal Parameters: ")
print(popt)

print("Estimated covariant of the plot")
print(pcov)

# To demonstrate that the PDF correctly shows our function
ys = PDF(x, popt)
plt.plot(x,ys)
plt.xlabel(label)
plt.ylabel('Normalised Arbitrary Units')
plt.hist(xtime, bins = bns, density = True)
plt.title('PDF with curve_fit over normalised mass to demonstrate functionality')
plt.show()






# Optimise for NLL, these parameters
# lik_model = opt.minimize(NLL2, np.array([mu, F, lmbda, sigma]), method='Nelder-Mead')
#print(lik_model)

# [F, mu, lmbda, sigma] doesnt work
# [popt] doesn't work, both dont work in the exact same manner
con1 = {'type':'ineq','fun': conf}
con2 = {'type':'ineq','fun': conf2}
concol = [con1,con2]

lik_model = opt.minimize(NLL, np.array([popt]),args = (np.array(xtime)),method = 'Nelder-Mead')
print(lik_model)


# Collect optimized parameters, and apply them to initial PDF to get plot
# Using bin_centers to allow for correct parameters to apply to usable data.
lmbdamle = lik_model.x
ymle = PDF(bin_centers, lmbdamle)
plt.plot(bin_centers, ymle)
plt.xlabel(label)
plt.ylabel('Normalised Arbitrary Units')
plt.hist(np.array(xtime), bins = bns, density = True)
plt.title('Normalised Maximum Likelihood fit of Lifetime')
plt.show()



# FINDING VARIANCE ''

# Negative Log Likelihood - Rewritten again to hold mean mass as a constant value
def NLLparam(args):
    xi = args[0]
    lmbl = args[1]
    # METHOD 1
    pdfresult = PDF(xi, lmbl)
    L = -np.sum(np.log(pdfresult))
    return L



# Finding parameter uncertainties
# Take minimum values, but then change the values
# set initial value for NLL value

initval = lik_model.fun

# set initial value for mu

lambd = lik_model.x[0]
print("lambda value")
print(lambd)


# set up altered NLL for incrementing over
alteredNLL = NLLparam((np.array(xtime),lambd))

# As this takes some time, allow for this to be skipped via yes/no option

join = input("Calculate the covariance for mean of the mass? WARNING: This may take some time to complete. [Y/N]")
if join.lower() == 'y':
    while True:
        # As there are no values being minimised, can just alter NLL directly and find difference
        alteredNLL=NLLparam((np.array(xtime),lambd))
        if (alteredNLL-lik_model.fun) > 0.5 or (alteredNLL-lik_model.fun) < -0.5:
            break
        # change this value (lamb) to change speed, smaller number -> longer time.
        lambd += 0.000001
        print("Covariance not found. Incrementing to " + str(lambd))

    print("lambda value off normal")
    print(lambd)
    # Show variance in mu
    print("Variance in lambda")
    print(lambd-lik_model.x[0])
elif join.lower() == 'n':
    print("Lambda covariance will not be calculated.")
