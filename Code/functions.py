import numpy as np
import matplotlib.pyplot as plt
import math as m


# Band Size
def bandwidth(peak, lower, upper):
    """
    Method to find the bandwidth from the peak value and two values that
    separate it from said peaks.

    :param peak: Peak value
    :param lower and upper: Lower and Upper values wrt the peak of our band
    :return: list of lower and upper values of our band
    """
    return([peak+lower,peak+upper])

# No. of events in region
def events(data, lower, upper):
    """
    Method to find how many events occur within a certain region of our data.

    :param data: list of data that is sorted through
    :param lower and upper: Lower and Upper values of the section we want to be counted
    :return: Total number of data points within selected region
    """
    # Counting no. of events
    val = 0
    # Loop over all data within region, add to count if said value is within region
    # /1000 due to conversion to TeV
    for i in range(0,int(len(data))):
        if (((data[i]) < (upper)) and ((data[i]) > (lower))):
            val += 1
    # return count
    return val

