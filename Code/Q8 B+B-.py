import numpy as np
import matplotlib.pyplot as plt
import math as m
import scipy.optimize as opt
from scipy.signal import find_peaks
from scipy import stats
import functions as fnc
import matplotlib.patches as mpatches




# import data
# xmass = np.loadtxt(sys.argv[1])
f = open("kdata-small.bin","r")
datalist = np.fromfile(f,dtype=np.float32)

# number of events
nevent = len(datalist)/5
xdata = np.split(datalist,nevent)
print(xdata[0])

# make lists for all events
xmass = []

# For finding N values within signal region (Question 3)
sigvals = 0


# Giving each component of xfull a separate number of bins for better data evaluation, decided via eyeballing "suitable bins"
bns = 200
# Same as above, but with labels, titles and savefiles of graphs
label = 'Invariant mass (GeV/$c^2$)'
title = 'Invariant masses of $J/ \psi K^+$ against Count'
svfile = ['peakmassplotB+.pdf','peakmassplotB-.pdf']







## B+ AND B- SECTION##

# to deal with

#Create a list of two arrays for B+ and B-, and separate data
xposneg=([],[])
# 0th component, positive data
# 1st component, negative data


# Create sorting loop
for i in range(0,int(nevent)):
    # If charge is positive
    if xdata[i][3] > 0:
        xposneg[0].append(xdata[i][0]/1000)
    else:
        xposneg[1].append(xdata[i][0]/1000)

# To allow for customising of axes in loop
labelb = ["B+ meson","B- meson"]

# create 2 times loop to apply for positive and negative
for b in range(2):
    #xmass = np.array(xmass)/1000

    # Plot histograms
    # Find binsize for each component
    binsize = (max(xposneg[b]) - min(xposneg[b]))/bns
    # Show histograms and save each one with unique title and labels
    plt.hist(xposneg[b], bins = bns)
    plt.ylabel('Count')
    plt.xlabel(label)
    plt.title(title + ' for ' + labelb[b])
    plt.show()


    # Finding mass of B+ by determining the bins with the highest number of entries in the peak regions.

    # Histogram of xposneg[b], removed from plot to allow for more modification
    xposneghist = np.histogram(xposneg[b],bins = bns)
    # 0th Component - Histogram values
    # 1st Component - Bin edges


    # Set up peaks function

    # Only includes peaks of the histogram where the prominence = 50, so no small scale minima/maxima.
    # 2nd graph includes when heigh = 700, to isolate local peaks at general maxima, but included an unwanted peak.
    peaks, _ = find_peaks(xposneghist[0], prominence = 50)

    # Bin edges, first value sliced off so that graph works
    binedges = np.around(xposneghist[1]).astype(int)[1:]
    # Bin centers, no value sliced off. But keeping binedges as backup
    bin_centers = 0.5*(xposneghist[1][1:]+xposneghist[1][:-1])



    # Print Maximum Peak value bin edges, the peak is within this bin
    print("Peak Value: ")
    print(xposneghist[1][peaks],xposneghist[1][peaks+1])
    print("Total Number of Events for " + labelb[b] + ":")
    print(len(xposneg[b]))

    # Find peak value
    peakval = (xposneghist[1][peaks] + xposneghist[1][peaks+1])/2
    # NOTE there is uncertainty here that isn't considered, due to bin size. Remember this!


    # find signal region
    sband = fnc.bandwidth(peakval,-25/1000,25/1000)

    # Find how many events within signal region #
    sigvals = fnc.events(xposneg[b],sband[0],sband[1])



    # Find upper and lower sidebands distributions #

    # lower band is -50,-62.5 from lower side of signal region
    lband = fnc.bandwidth(peakval, -75/1000, -50/1000)


    # upper band is +50,+62.5 from peakval
    uband = fnc.bandwidth(peakval, 50/1000, 75/1000)


    # Find how many events within upper and lower sidebands

    # lower band
    lbandevents = fnc.events(xposneg[b],lband[0],lband[1])

    # upper band
    ubandevents = fnc.events(xposneg[b],uband[0],uband[1])



    # Find how many events within signal region with background removed #

    truevals = sigvals - (ubandevents + lbandevents)
    print("The estimate number of events within the signal region with background removed: " + str(truevals))

    # Plot bin edges (x) against Count with Peak labelled
    plt.plot(bin_centers, xposneghist[0])
    plt.plot(bin_centers[peaks], xposneghist[0][peaks], "x")
    # Plot sidebands
    plt.axvspan(lband[0],lband[1],facecolor='b',alpha = 0.2)
    plt.axvspan(uband[0],uband[1],facecolor='b',alpha = 0.2)
    plt.axvspan(sband[0],sband[1],facecolor='g',alpha = 0.2)
    sidebands = mpatches.Patch(color='b', label = 'Sidebands')
    sigbands = mpatches.Patch(color='g', label = 'Signal Band')
    plt.legend(handles=[sidebands,sigbands])
    plt.ylabel('Count')
    plt.xlabel(label)
    plt.title('Invariant masses of ' + labelb[b] + ' against Count with marked Peak')
    plt.show()


    # Guesses at suitable values for each variable in PDF
    # Composite Fraction#
    F = 0.7
    # Lambda for exponential scaling
    lmbda = 0.001
    # Sigma for stdev
    sigma = 0.00794
    # mu for Mean
    mu = peakval[0]


    # new normalised histogram for PDF
    xmasshistnorm = np.histogram(xposneg[b],bins = bns, density = True)

    # Defining x and y points to compare to
    x = bin_centers

    y = xmasshistnorm[0]

    # Probability distribution function
    def PDF(xp, f, mew, lmb, sig):
        #print("Variables: ")
        #print(f, mew, lmb, sig)
        gauss = (1/(np.power(2*np.pi,0.5)*sig)) * np.exp(-0.5*np.power((xp-mew)/sig,2))
        expon = np.exp(-1*lmb*xp)*abs(lmb/(np.exp(-lmb*min(xp))-np.exp(-lmb*max(xp))))
        y = (1-f)*expon + (f)*gauss
        return y

    # Negative Log Likelihood
    def NLL(params, args):
        xi = args
        fl, mewl, lmbl, sigl = params

        # METHOD 1
        pdfresult = PDF(xi, fl, mewl, lmbl, sigl)
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
    popt, pcov = opt.curve_fit(PDF, bin_centers, xmasshistnorm[0], (F, mu, lmbda, sigma))
    # 0th component - optimal parameters for PDF
    # 1st component - variance of said parameters
    print("Optimal Parameters: ")
    print(popt)

    print("Estimated covariant of the plot")
    print(pcov)

    # To demonstrate that the PDF correctly shows our function
    ys = PDF(x, F, mu, lmbda, sigma)
    plt.plot(x,ys)
    plt.xlabel(label)
    plt.ylabel('Normalised Arbitrary Units')
    plt.hist(xposneg[b], bins = bns, density = True)
    plt.title('curve_fit PDF over normalised mass to demonstrate functionality')
    plt.show()






    # Optimise for NLL, these parameters
    # lik_model = opt.minimize(NLL2, np.array([mu, F, lmbda, sigma]), method='Nelder-Mead')
    #print(lik_model)

    # [F, mu, lmbda, sigma] doesnt work
    # [popt] doesn't work, both dont work in the exact same manner
    con1 = {'type':'ineq','fun': conf}
    con2 = {'type':'ineq','fun': conf2}
    concol = [con1,con2]

    lik_model = opt.minimize(NLL, np.array([popt]),args = (np.array(xposneg[b])),method = 'Nelder-Mead')
    print(lik_model)


    # Collect optimized parameters, and apply them to initial PDF to get plot
    # Using bin_centers to allow for correct parameters to apply to usable data.
    Fmle, mumle, lmbdamle, sigmamle = lik_model.x
    ymle = PDF(bin_centers, Fmle, mumle, lmbdamle, sigmamle)
    plt.plot(bin_centers, ymle)
    plt.xlabel(label)
    plt.ylabel('Normalised Arbitrary Units')
    plt.hist(np.array(xposneg[b]), bins = bns, density = True)
    plt.title('Maximum Likelihood fit via use of Minimize - Nelder-Mead of ' + labelb[b])
    plt.show()



    # FINDING VARIANCE ''

    # Negative Log Likelihood - Rewritten again to hold mean mass as a constant value
    def NLLparam(params, args):
        xi = args[0]
        mewl = args[1]
        fl, lmbl, mewfake, sigl = params

        # METHOD 1
        pdfresult = PDF(xi, fl, mewl, lmbl, sigl)
        L = -np.sum(np.log(pdfresult))
        return L



    # Finding parameter uncertainties
    # Take minimum values, but then change the values
    # set initial value for NLL value

    initval = lik_model.fun

    # set initial value for mu

    moo = lik_model.x[1]
    print("mu value")
    print(moo)


    # As this takes some time, allow for this to be skipped via yes/no option

    join = input("Calculate the covariance for mean of the mass? WARNING: This may take some time to complete. [Y/N]")
    if join.lower() == 'y':
        while True:
            lik_modelparam = opt.minimize(NLLparam, lik_model.x ,args = [(np.array(xposneg[b])),moo] ,method = 'Nelder-Mead')
            # if the difference in likelihoods is greater or less than 0.5, take that value of mu as uncertainty difference
            if (lik_modelparam.fun-lik_model.fun) > 0.5 or (lik_modelparam.fun-lik_model.fun) < -0.5:
                break
            # change this value (moo) to change speed, smaller number -> longer time.
            moo += 0.000001
            print("Covariance not found. Incrementing to " + str(moo))

        print("mu value off normal")
        print(moo)
        # Show variance in mu
        print("Variance in mu")
        print(moo-lik_model.x[1])
    elif join.lower() == 'n':
        print("Mean mass covariance will not be calculated, continuing to residuals...")




    # Finding Residuals - ymle already set for bin_centers

    yresiduals = (xmasshistnorm[0]-ymle)

    plt.bar(bin_centers, yresiduals, 0.001)
    plt.title('Nelder Mead NLL and Histogram Residuals of ' + labelb[b])
    plt.show()


    # Full plot with residuals underneath
    fig, axs = plt.subplots(2,gridspec_kw={
                               'height_ratios': [3, 1]})
    fig.suptitle('Normalised Nelder-Mead NLL fit with Data Histogram of ' + labelb[b])
    plt.xlabel(label)
    # done to force label to center, as it is (for some reason) attached to bottom subplot
    plt.ylabel('Residuals                                                                             Normalised Arbitrary Units')
    resid = mpatches.Patch(color='dodgerblue', label = 'Residuals')
    plt.legend(handles=[resid])
    axs[0].plot(bin_centers, ymle)
    axs[0].hist(np.array(xposneg[b]), bins = bns, density = True)
    axs[1].bar(bin_centers, yresiduals, 0.001)

    plt.show()
