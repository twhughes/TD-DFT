import os 
os.chdir('/Users/twh/Documents/Reed/Code')
os.chdir('../')
print os.getcwd()
from pylab import plot, show, xlabel, ylabel, legend, subplot, scatter
import numpy as np
import numpy.fft as fft

files =  ['./Graphene2.tddft-new2-1','./Graphene2.tddft-new2-2','./Graphene2.tddft-new2-3', \
         './Graphene2.tddft-new2-4','./Graphene2.tddft-new2-5','./Graphene2.tddft-new2-6', \
         './Graphene2.tddft-new2-7','./Graphene2.tddft-new2-8' \
         ]
         
timeStep = 1
auTime = 4.8377687e-17
auCharge = 1.6e-19       
auLength = 5e-11
auPotential = 20
auEfield = 5.14e11
d = 5*auLength
epoints = []
S = []
SF = []
SF2 =[]
for fname in files:
    print fname
    with open(fname, 'r') as f:
        data = f.readlines()
        x = []
        y = []
        z = []
        eChecked = 0
        
        for line in data:
            words = line.split()
            if (len(words)>5 and words[0] == 'DIP'):
                x.append(float(words[3]))
                y.append(float(words[4]))
                z.append(float(words[5])*auCharge*auLength)
            if (len(words) > 1 and words[0] == "Time"):
                timeStep = float(words[3])*auTime
                print "time step: ", timeStep, "sec"
            if (len(words) > 1 and words[0] == "E" and words[1] == "field" and words[2] == "amplitude" and eChecked == 0):
                print words[-1]
                eChecked = 1
                epoints.append(float(words[-1])*auEfield)
    N = len(x)
    tpoints = []
    fpoints = []
    for n in range(N):
        tpoints.append(n*timeStep)
        fpoints.append(n/timeStep)

    zt = []
    ztt = []
    #first time derivative
    for i in range(len(z)-1):
        zt.append((z[i+1] - z[i])/timeStep)

    #second time derivative
    for i in range(len(z)-2):
        ztt.append((zt[i+1]-zt[i])/timeStep)

    #time average
    timeAvg = 0
    for i in range(len(ztt)):
        timeAvg += zt[i]*timeStep
    timeAvg /= len(zt)*timeStep
    S.append(timeAvg)
    F_z = fft.rfft(z)
    F_dz = fft.rfft(zt)
    F_dzz = fft.rfft(ztt)
    s = 0
    for n in range(10):
        s += F_dzz[n]
    FNEW = F_dz
    FNEW[1:] = 0
    SF.append(s)
    dzNEW = fft.ifft(FNEW)
    #time average dznew
    ta = 0
    for i in range(len(dzNEW)):
        ta += dzNEW[i]
    ta /= len(zt)
    SF2.append(ta)

print auCharge*auLength/auTime

xlabel('time (seconds)')
#ylabel('dP/dt (current)')
#plot(tpoints[0:-1],zt)
#show()
"""  
plot(A)
show()
plot(B)
show()
plot(G)
xlabel('w (THz)')
ylabel('|F(P''(t))|')
show()
"""
epoints = [ 1e-4,2e-4,3e-4,4e-4,5e-4,6e-4,7e-4,8e-4]#,7,8,9,0.0011,0.00120,0.0013,0.0014,0.0015,0.0016,0.0017,0.0018,0.0019]  
for i in range(len(epoints)):
    epoints[i] *= -auPotential/auLength
print epoints
line = np.polyfit(epoints, SF2,1)
pls = []
lrange = np.arange(min(epoints), max(epoints), (max(epoints)-min(epoints))/100.0)
for l in lrange:
    pls.append(line[1]+l*line[0])
print 'slope = ', line[0]
print ''
print 'R = ', d*d/line[0], 'Ohms'
plot(lrange, pls)
xlabel('E field (Volts/m)')
ylabel('I(V) (Amps)')
scatter(epoints, SF2)
show()

