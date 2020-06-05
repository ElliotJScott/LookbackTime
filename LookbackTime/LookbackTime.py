import numpy
import scipy
import scipy.integrate as ig
import matplotlib.pyplot as pt

def IntegrationFunction(x):
    Ez = numpy.sqrt(OM*(1+x)**3 + Ok*(1+x)**2 + OL)
    func = 1 / ((1+x) * Ez)
    return func
def SecToGyr(s):
    return s / (3.1556926 * (10**16))

h = 0.72 #h from H0 = 100h km s^-1 Mpc^-1
tH = 3.0856 * (10**17) / h
OMVals = [1.0, 2.0, 0.04, 0.27]
OLVals = [0.0, 0.0, 0.0, 0.73]
OM = OMVals[3] #use decimal place if not present to force omega values to be floats
OL = OLVals[3]
Ok = 1.0 - (OM + OL)
zMaxValues = [0.1, 1.0, 2.0, 4.0, 6.7, numpy.inf] #numpy.inf added to get age of universe
for z in zMaxValues:
    integral = ig.quad(IntegrationFunction, 0, z)
    lookbackTime = tH * integral[0]
    lbtInGyr = SecToGyr(lookbackTime)
    print(str(z) + " " + str(lbtInGyr) + " " + "{:e}".format(lookbackTime))

plots = []
fig, ax = pt.subplots()
ax.set(xlabel='Redshift z', ylabel='Lookback Time $t_L$ / Gyr',
           title='A plot of Lookback Time against Redshift for various $\Omega_m$ and $\Omega_\Lambda$')
for i in range(4):
    OM = OMVals[i];
    OL = OLVals[i];
    Ok = 1.0 - (OM + OL)
    tLVals = []
    zVals = numpy.arange(0.0, 6.7, 0.01)
    for z in zVals:
        integral = ig.quad(IntegrationFunction, 0, z)
        lbTime = SecToGyr(tH * integral[0])
        tLVals = numpy.append(tLVals, lbTime)
    
    
    p = ax.plot(zVals, tLVals)
    plots = numpy.append(plots, p)
ax.legend((plots[0], plots[1], plots[2], plots[3]), ('$\Omega_m$ = 1.00, $\Omega_\Lambda$ = 0.00', '$\Omega_m$ = 2.00, $\Omega_\Lambda$ = 0.00', '$\Omega_m$ = 0.04, $\Omega_\Lambda$ = 0.00', '$\Omega_m$ = 0.27, $\Omega_\Lambda$ = 0.73'), loc='lower right', shadow=False)
pt.show()

