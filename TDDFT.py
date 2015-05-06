import os 
os.chdir('../')
print os.getcwd()

timeStep = 1
rydToSec = 4.8377687e-17
with open('./CO2.tddft-out', 'r') as f:
    data = f.readlines()
    x = []
    y = []
    z = []
    for line in data:
        words = line.split()
        if (len(words)>5 and words[0] == 'DIP'):
            x.append(float(words[3]))
            y.append(float(words[4]))
            z.append(float(words[5]))
        if (len(words) > 1 and words[0] == "Time"):
            timeStep = float(words[3])*rydToSec
            print "time step: ", timeStep, "sec"

N = len(x)
tpoints = []
fpoints = []
for n in range(N):
    tpoints.append(n*timeStep)
    fpoints.append(n/timeStep)

from pylab import plot, show, xlabel, ylabel, legend, subplot
import numpy as np
import numpy.fft as f
subplot(2,1,1)
plot(tpoints,x)
plot(tpoints,y)
plot(tpoints,z)
ylabel('dipole moment (arb units)')
legend(['x','y','z'])
xlabel('time (seconds)')
subplot(2,1,2)
print len(fpoints)/2, len(np.abs((f.rfft(x))))
plot(fpoints[0:len(fpoints)/2+1],np.abs((f.rfft(x))))
ylabel('fourier transform along E field direction')
legend(['fft(x)'])
xlabel('frequency (1/sec)')
show()
