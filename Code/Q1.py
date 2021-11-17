import numpy as np
import matplotlib.pyplot as plt

# Q1 is saved as a separate python file, due to its requirement for many plots to be displayed, which significantly slow down the speed at which


# import data
# xmass = np.loadtxt(sys.argv[1])
f = open("kdata-small.bin","r")
datalist = np.fromfile(f,dtype=np.float32)

# number of events
nevent = len(datalist)/5
xdata = np.split(datalist,nevent)


# make lists for all events
xmass = []
xmom = []
xpseu = []
xcharge = []
xlifetime = []

# list of lists method
xfull = []
xfull.append(xmass)
xfull.append(xmom)
xfull.append(xpseu)
xfull.append(xcharge)
xfull.append(xlifetime)

# xfull[i]
# i =;
# 0 - Mass
# 1 - Momentum
# 2 - Pseudorapidity
# 3 - Charge
# 4 - Lifetime

# Giving each component of xfull a separate number of bins for better data evaluation, decided via eyeballing "suitable bins"
bns = [200, 100, 75, 10, 150]
# Same as above, but with labels, titles and savefiles of graphs
label = ['Invariant mass (MeV/$c^2$)','Transverse momentum (MeV/c)','Pseudorapidity (MeV/c)','Charge (arbitrary units)','Lifetime (ps)']
title = ['Invariant masses of $J/ \psi K^+$ against Count','Transverse momenta of $J/ \psi K^+$ against Count','Pseudorapidities of $J/ \psi K^+$ against Count','Charges of $J/ \psi K^+$ against Count','Lifetimes of $J/ \psi K^+$ against Count']
svfile = ['massplot.pdf','momplot.pdf','Pseuplot.pdf','Chargeplot.pdf','Lifetimeplot.pdf']



# Allows for one full loop to apply all data to one list of components (0 -> 4)
# Easier for creating loops to apply the same tricks (histograms)
for l in range(0,len(xdata[0])):
    for i in range(0,int(nevent)):
        xfull[l].append(xdata[i][l])
    #    This is for seeing all values of xmass
    #    if i < int(nevent):
    #        print(xmass[i])
print("No. of array components: " + str(len(xdata[0])))


# Plot histograms
for l in range(0,len(xdata[0])):
    # Find binsize for each component
    binsize = (max(xfull[l]) - min(xfull[l]))/bns[l]
    # Show histograms and save each one with unique title and labels
    plt.hist(xfull[l], bins = bns[l])
    plt.ylabel('Count')
    plt.xlabel(label[l])
    plt.title(title[l])
    plt.savefig(svfile[l])
    plt.show()

    print(binsize)
