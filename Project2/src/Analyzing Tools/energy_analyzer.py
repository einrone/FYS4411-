import numpy as np
from time import time

"""
Uses blocking to get mean and error of the energy for one simulation.

"""

def block(x):
    """
    The blocking function made by Marius. It is made so to give an error if
    the number if MC cycles times the number of processes are not on the form
    2^n.
    """
    # preliminaries
    n = len(x);
    if abs(np.log2(n) - int(np.log2(n))) > 1e-5:
        print("Number of MC_cycles are not on the form 2^n. Change that " \
                            +"and come back an other time")

        print(np.log2(n))
        exit()
    d = int(np.log2(n)); s, gamma = np.zeros(d), np.zeros(d);
    mu = np.mean(x); t0 = time()

    # estimate the auto-covariance and variances
    # for each blocking transformation
    for i in np.arange(0,d):
        n = len(x)
        # estimate autocovariance of x
        gamma[i] = (n)**(-1)*np.sum( (x[0:(n-1)]-mu)*(x[1:n]-mu) )
        # estimate variance of x
        s[i] = np.var(x)
        # perform blocking transformation
        x = 0.5*(x[0::2] + x[1::2])

    # generate the test observator M_k from the theorem
    M = (np.cumsum( ((gamma/s)**2*2**np.arange(1,d+1)[::-1])[::-1] )  )[::-1]

    # we need a list of magic numbers
    q =np.array([6.634897,9.210340, 11.344867, 13.276704, 15.086272, 16.811894, 18.475307, 20.090235, 21.665994, 23.209251, 24.724970, 26.216967, 27.688250, 29.141238, 30.577914, 31.999927, 33.408664, 34.805306, 36.190869, 37.566235, 38.932173, 40.289360, 41.638398, 42.979820, 44.314105, 45.641683, 46.962942, 48.278236, 49.587884, 50.892181])

    # use magic to determine when we should have stopped blocking
    for k in np.arange(0,d):
        if(M[k] < q[k]):
            break
    if (k >= d-1):
        print("Warning: Use more data")
    ans = s[k]/2**(d-k)
    # print("Runtime: %g sec" % (time()-t0)); print("Blocking Statistics :")
    # print("average            iterations      std. error")
    # print("%8g %20g %15g" % (mu, k, ans**.5))
    return ans


if __name__ == '__main__':
    data =  np.fromfile("../output/data.bin",sep=" ")
    print("Mean Energy: ",np.mean(data))
    print("Error: ",block(data))
