import numpy
from matplotlib.pylab import *

#reading file
data = []
f = open('project1datab1000.txt')
for line in f:
	data.append(line)
f.close()


#analytical solution
xanalytical = numpy.linspace(0,1,10000)
u = 1 - (1 - exp(-10)) * xanalytical - exp(-10 * xanalytical)

#plotting
xdata = numpy.linspace(0,1,len(data))
plot(xanalytical,u,'r', xdata,data,'b')
title(r"Analytical and numerical solution to $-u'' = f(x)$")
legend(["analytical", "numerical"])
xlabel('x(length)')
ylabel('u(electric potential)')
show()
