
# coding: utf-8

# In[1]:

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits as pyfits
import scipy.stats as sc

# In[2]:

def pause():
    '''
    This function pause the code until enter is pressed
    '''
#
# pauses until a keyboard entry (e.g. carrage return)
#
    print '\r'
    dummy = raw_input('Pause')
    print '\r'


# In[3]:

def npRan(seed1, seed2, seed3, it):
    '''
    Generates three sets of random numbers using np.random.uniform with a length of it
    -----------------------------------------------------------------------------
    inputs
    seed1, seed2, seed3: integers that start the random number generation
    it: integer that determines that length of the list of random numbers
    
    '''
    # it is the number of iterations the random number generator goes through
    np.random.seed(seed1)
    x = np.random.uniform(0,1,it)
    #
    np.random.seed(seed2)
    y = np.random.uniform(0,1,it)
    #
    np.random.seed(seed3)
    z = np.random.uniform(0,1,it)
    return x, y, z


# In[4]:

def johnRan(seed1,seed2,seed3,it):
    '''
    Generates three sets of random numbers using the John von Neumann method with a length     of it
    --------------------------------------------------------------------------------------
    inputs
    seed1, seed2, seed3: a 6 didget integer that start the random number generation
    it: integer that determines that length of the list of random numbers
    '''
    #makes empty arrays for the random numbers
    ran = np.array([])
    ran1 = np.array([])
    ran2 = np.array([])
    for n in range(it):
        g = str(seed1**2) #squares the seed and turns it to a string
        g = int(g[3:9]) #extracts the middle 6 numbers and makes it into a int
        if n > 0 and g < 1e4: #prevents zero's from existing
            g = (ran[n-1] * 1e6 + 123546)
        seed1 = g
        thing = g/1e6
        ran = np.append(ran, thing)
        if n > 0 and ran[n] in ran[:n]: #checks to make sure that it doesn't make a pattern
            seed1 = seed1 + 10
    for n in range(it):
        g = str(seed2**2) #squares the seed and turns it to a string
        g = int(g[3:9]) #extracts the middle 6 numbers and makes it into a int
        if n > 0 and g < 1e4: #prevents zero's from existing
            g = (ran1[n-1] * 1e6 + 123546)
        seed2 = g
        thing = g/1e6
        ran1 = np.append(ran1, thing)
        if n > 0 and ran1[n] in ran1[:n]: #checks to make sure that it doesn't make a pattern
            seed2 = seed2 + 10
    for n in range(it):
        g = str(seed3**2) #squares the seed and turns it to a string
        g = int(g[3:9]) #extracts the middle 6 numbers and makes it into a int
        if n > 0 and g < 1e4: #prevents zero's from existing
            g = (ran2[n-1] * 1e6 + 123546)
        seed3 = g
        thing = g/1e6
        ran2 = np.append(ran2, thing)
        if n > 0 and ran2[n] in ran2[:n]: #checks to make sure that it doesn't make a pattern
            seed3 = seed3 + 10
    return ran, ran1, ran2


# In[5]:
def chi(bins,obs,it,ex = 0): #finds the Chi^2 of the data
    '''
    calculates the chi^2 of a histogram
    -------------------------------------------------------------------------
    inputs:
    bins: integer, the amount of bins in the histogram
    obs: integer, the sizes of each bin in the histogram
    it: integer, how large is the set of data you are analyzing
    ex: list or integer defaults as 0, what is the expected value of the data
    '''
    T = 0
    Ohno = 0
    if type(ex) is int or type(ex) is float:
        if ex == 0:
            ex = it/bins #the expected value
            #the chi^2 equation
            t = (obs - ex)**2 / ex
            T = np.sum(t)
        else:
            #the chi^2 equation
            t = (obs - ex)**2 / ex
            T = np.sum(t)
    else:
        for n in range(bins-1): #the chi^2 equation
            if ex[n] < 5 or obs[n] < 5:
                Ohno +=1
            else:
                t = (obs[n] - ex[n])**2 / ex[n]
                T = T + t
    return T

# In[6]
def serial(ran,it): #does the serial test
    '''
    does the serial test on the set of data
    --------------------------------------------------------------
    inputs
    ran: the set of data that is being analyzed
    it: the sizes of the data set
    '''
    q = 0
    Q = 0
    for n in range(it): #sees how many times two identical numbers are next to eachother
        if n > 0 and ran[n-1]==ran[n]:
            q+=1
    for n in range(it): #sees how many times a number followed by the next number are next to eachother
        if n > 0 and ran[n-1]==ran[n] - 1e-11:
            Q+=1
    return q, Q
# In[7]
def csvArray(File): # makes an x, y, and z position array from a csv file
    '''
    takes apart a csv file and makes it into 3 arrays
    ---------------------------------------------------------------
    inputs
    File: string, the files location
    '''
    a = np.genfromtxt(File ,dtype=float,delimiter=',',skip_header=1) #generates an array from csv file
    x = np.array([])
    y = np.array([])
    z = np.array([])
    for n in range(len(a)):
        x = np.append(x, a[n][0]) #pulls out x values
        y = np.append(y, a[n][1]) #pulls out y values
        z = np.append(y, a[n][2]) #pulls out pixel intensity or z position
    return x , y, z

# In[8]
def figName(plotname, xaxis, yaxis):
    '''
    names the plot and the axes
    ---------------------------
    inputs
    plotname: string, name of the plot
    xaxis: string, name of the x axis
    yaxis: string, name of the y axis
    '''
    plt.title(plotname, fontsize=18)
    plt.xlabel(xaxis,fontsize=16)
    plt.ylabel(yaxis,fontsize=16)
    return ''

# In[9]
def nearestDistance(x,y):
    '''
    calculates the nearest and second nearestdistance between stars
    ---------------------------------------------------------------------
    inputs
    x: array, the x position of the star
    y: array, the y position of the star
    '''
    T = np.array([])
    O = np.array([])
    for n in range(len(x)): 
        xdis = x[n] - x #all the x distances for a single object, y bellow
        ydis = y[n] - y
        t = np.sqrt(xdis**2 + ydis**2) #Distance to other objects
        t= t[t != 0]
        T = np.append(T, min(t))
        t = t[t != min(t)] #part 1 in finding the second nearest value
        O = np.append(O, min(t)) #finds the second nearst value
    return T, O

# In[10]
def figNameSub(plotname, xaxis, yaxis, t, T):
    t.set_title(plotname, fontsize = 18)
    plt.xlabel(xaxis, fontsize = 16)
    t.set_ylabel(yaxis, fontsize = 16)
    T.set_ylabel(yaxis, fontsize = 16)
    return ''

# In[11]
def fitsHeader(filename):
    '''
    pulls the Juldate, AVGWIDTH, and RMSA values from an input fits file
    '''
    ki = pyfits.open(filename)
    zhdate = ki[0].header['JULDATE'] #pulls out the date of the observation
    zhavg = ki[0].header['AVGWIDTH'] #pulls out the average width
    zhrmsa = ki[0].header['RMSA'] #pull out the RMSA
    return zhdate, zhavg, zhrmsa

# In[12]
def leastSquare(x,y,sigma):
    '''
    calculates the slope b and the y intercept a from a list of data and the y error of that data
    --------------------------------------------------------------------------------
    inputs:
    x: array, the first set of data and the data in the x axis
    y: array, the second set of the data and the data in the y axis
    sigma: the error in the y axis
    '''
    s1 = np.sum(1/sigma**2)
    sx = np.sum(x/sigma**2)
    sy = np.sum(y/sigma**2)
    sxy = np.sum((x * y)/sigma**2)
    sxx = np.sum((x**2)/sigma**2)
    a = (sy * sxx - sx * sxy)/(s1 * sxx - sx**2) #the y intercept
    b = (s1 * sxy - sx * sy)/(s1 * sxx - sx**2) #the slope
    sigmaA = np.sqrt(sxx / (s1 * sxx - sx**2)) #the error in a
    sigmaB = np.sqrt(s1 / (s1 * sxx - sx**2)) # the error in b
    return a, b, sigmaA, sigmaB

# In[13]
def chiLLS(x,y,a,b,sigma= 5):
    '''
    calculates the chi squared value for the least square equation
    ------------------------------------------------------------------
    inputs:
    x: array, the first set of data and the data in the x axis
    y: array, the second set of the data and the data in the y axis
    a: float, the y intercept of the line of best fit
    b: float, the slope of the line of best fit
    sigma: float or array, the error in the y axis
    '''
    mu = b * x + a #mean of the data
    y = y - mu #ajusts the values to be around the line of best fit
    if sigma == 5:
        vernc = sc.moment(y,2) #the standard devation of the data
        chi = (y)**2 / vernc #the chi squared of the data
    else:
        chi = (y)**2 / sigma
    return np.sum(chi)

# In[14]
def red2Velocity(b):
    '''
    calculates the redshift and the effective velocity from the velocity of the object
    ----------------------------------------------------------------------------------
    inputs:
    b: array, the velocity data
    '''
    c = 3e5 #speed of light 
    ratio = b / c
    top = np.sqrt(1 + ratio)
    bottom = np.sqrt(1 - ratio)
    z = top / bottom - 1 #redshift
    veff = c * z #effective velocity
    return z, veff

# In[15]
def MLinearLeastSquare(distance,alpha,beta,velocity,sigma):
    '''
    calculates the linear least square throught the matrix method
    --------------------------------------------------------------------------------------
    inputs:
    distance: array, the distance of the objects
    alpha: float, the RA of the object in radians
    beta: float, the DEC of the object in radians
    velocity: array, the velocity of the objects
    sigma: array, the error in velocity
    '''
    R = distance # first part of equation
    X = np.cos(alpha)*np.cos(beta) # second part of equation
    Y = np.sin(alpha)*np.cos(beta) # third part of equation
    Z = np.sin(beta) # fourth part of equation
    A = np.array([
            [np.sum(R*R/sigma**2),np.sum(R*X/sigma**2),np.sum(R*Y/sigma**2),np.sum(R*Z/sigma**2)],
            [np.sum(X*R/sigma**2),np.sum(X*X/sigma**2),np.sum(X*Y/sigma**2),np.sum(X*Z/sigma**2)],
            [np.sum(Y*R/sigma**2),np.sum(Y*X/sigma**2),np.sum(Y*Y/sigma**2),np.sum(Y*Z/sigma**2)],
            [np.sum(Z*R/sigma**2),np.sum(Z*X/sigma**2),np.sum(Z*Y/sigma**2),np.sum(Z*Z/sigma**2)]
        ])
    b = np.sum(velocity * np.array([R,X,Y,Z]), axis = 1)
    a = np.linalg.inv(A).dot(b) #the constants
    return a