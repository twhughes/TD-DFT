import os 

#change to the directory where this file lives on your system
os.chdir('/Users/twh/Documents/Reed/Code') 

#change dir to where the tddft-out files live
os.chdir('../') 

#print working dir
print os.getcwd()

#import modules
from pylab import plot, show, xlabel, ylabel, legend, subplot, scatter
import numpy as np
import numpy.fft as fft

#list of tddft-out filenames (change to your needs)
files =  ['./Graphene2.tddft-new2-1','./Graphene2.tddft-new2-2','./Graphene2.tddft-new2-3', \
         './Graphene2.tddft-new2-4','./Graphene2.tddft-new2-5','./Graphene2.tddft-new2-6', \
         './Graphene2.tddft-new2-7','./Graphene2.tddft-new2-8' \
         ]

#initialize time step variable
timeStep = 1

#unit conversion values (au -> SI) [multiply an au value by these numbers to get SI]
auTime = 4.8377687e-17
auCharge = 1.6e-19       
auLength = 5e-11
auPotential = 20
auEfield = 5.14e11

#spacing between bilayers (specified in PW calculation, change for your case)
d = 5*auLength
print d, ' spacing'
#lists to store I-V values
epoints = []   #electric field (bias)
S = []         #straight time average list  (relevant)
SF = []        #fourier list                (doesn't work well) 
SF2 =[]        #fourier time average list   (relevant)

#loop through files
for fname in files:
    print fname
    with open(fname, 'r') as f:
        data = f.readlines()

        #lists to store dipole moment vs. time in x,y,z directions
        x = []
        y = []
        z = []
        
        #bool to specify whether we have stripped the bias from the output file already              #(only want to do it once)
        eChecked = 0

        #loop through lines in output file
        for line in data:
            words = line.split()

            #strip the dipole moment (and convert to SI)
            if (len(words)>5 and words[0] == 'DIP'):
                x.append(float(words[3])*auCharge*auLength)
                y.append(float(words[4])*auCharge*auLength)
                z.append(float(words[5])*auCharge*auLength)

            #strip the time step
            if (len(words) > 1 and words[0] == "Time"):
                timeStep = float(words[3])*auTime
                print "time step: ", timeStep, "sec"

            #strip the bias (if hasn't been done before)
            if (len(words) > 1 and words[0] == "E" and words[1] == "field" and words[2] == "amplitude" and eChecked == 0):
                print words[-1]
                eChecked = 1
                epoints.append(float(words[-1])*auEfield)

    #for plotting purposes, generates x axis (time and frequency) lists for plotting
    #stored in lists (tpoints, fpoints)
    N = len(x)
    tpoints = []
    fpoints = []
    for n in range(N):
        tpoints.append(n*timeStep)
        fpoints.append(n/timeStep)

    #new lists to hold time derivatives
    zt = []
    ztt = []
    
    #first time derivative
    for i in range(len(z)-1):
        zt.append((z[i+1] - z[i])/timeStep)

    #second time derivative
    for i in range(len(z)-2):
        ztt.append((zt[i+1]-zt[i])/timeStep)

    #time average of zt
    timeAvg = 0
    for i in range(len(zt)):
        timeAvg += zt[i]*timeStep
    timeAvg /= len(zt)*timeStep

    #add to S list
    S.append(timeAvg)

    #take fourier transforms
    F_z = fft.rfft(z)
    F_dz = fft.rfft(zt)
    F_dzz = fft.rfft(ztt)
    
    #for Fourier method (not relevant, but included) sum low freq. components of FFT p''(t)
    s = 0
    for n in range(10):
        s += F_dzz[n]
    #add to SF list
    SF.append(s)
    
    #for working Fourier smoothing method:
    #copy P'(t) fourier transform
    FNEW = F_dz
    
    #set all high frequency terms to 0 for smoothing
    highFreq = 1          #change to higher number if you want higher freq terms
    FNEW[highFreq:] = 0   

    #Inverse FT back to get current time series again
    dzNEW = fft.ifft(FNEW)
    
    #time average new time series
    ta = 0
    for i in range(len(dzNEW)):
        ta += dzNEW[i]
    ta /= len(zt)
    SF2.append(ta)

#list of electric field points (I manually entered but you can strip and append if you prefer)
epoints = [ 1e-4,2e-4,3e-4,4e-4,5e-4,6e-4,7e-4,8e-4]

#convert each epoint to SI (V/m)
for i in range(len(epoints)):
    epoints[i] *= -auPotential/auLength

#fit line to these points (line[0] = y intercept, line[1] = slope)
line = np.polyfit(epoints, SF2,1) #!!!!!Change SF2 to S or SF depending on what you want!!!!

#list to store line for plot
pls = []

#x axis for line plot
lrange = np.arange(min(epoints), max(epoints), (max(epoints)-min(epoints))/100.0)

#create y axis for line plot
for l in lrange:
    pls.append(line[1]+l*line[0])

#print slope
print 'slope = ', line[0]
print ''

#print resistance (all in SI, but need to divide by d^2 (in m^2, calc. above) to make units I(V)R = V)
print 'R = ', d*d/line[0], 'Ohms'

#plot the line
plot(lrange, pls)
xlabel('E field (Volts/m)')
ylabel('I(V) (Amps)')

#scatter plot the measured I-V points
scatter(epoints, SF2) #!!!!!Change SF2 to S or SF depending on what you want!!!!

#show the plot
show()