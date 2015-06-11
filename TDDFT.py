import os 
from pylab import plot, show, xlabel, ylabel, legend, subplot
import numpy as np
import numpy.fft as fft

print os.getcwd()

#clean and make tmp directory and navigate to it
os.chdir('/Users/twh/Documents/Reed/tmp/')
os.system('rm -rf tmp')
os.system('mkdir tmp')
os.chdir('/Users/twh/Documents/Reed/tmp/tmp')
print os.getcwd()

#transfer TDDFT output files from server to tmp directory
os.system('scp twhughes@centurion.stanford.edu:~/espresso/tddft-edit/Example/Graphene2/tmp/Graphene2.tddft-in*.0 ./')

#legend list
l = []

numfiles = 0
for fn in os.listdir('/Users/twh/Documents/Reed/tmp/tmp'):
    numfiles += 1


#loop through files in tmp directory and load data into x, y, z, & l lists
iter = 1
zs = []
for fn in os.listdir('/Users/twh/Documents/Reed/tmp/tmp'):

    print fn

    timeStep = 1                    #initialize time step variable
    rydToSec = 4.8377687e-17        #seconds / rydberg time units (convert rydberg t units by multiplying by this factor)
    dist = 1
    with open(fn, 'r') as f:        #open file
        dist = fn[18:]              #get the bilayer distance parameter from the last digits of the filename
        print 'distance: ', dist    #NOTE: should get correct units here (need to look up)
        data = f.readlines()        #load all file into 'data' list of lines
        x = []                      #dipole x / time
        y = []                      #dipole y / time
        z = []                      #dipole z / time
        for line in data:           #loop through lines
            words = line.split()                            #split lines by word / number
            if (len(words)>5 and words[0] == 'DIP'):        #if the line is specifying a dipole quantity
                x.append(float(words[3]))                   #load x into x
                y.append(float(words[4]))                   #load y into y
                z.append(float(words[5]))                   #load z into z
            if (len(words) > 1 and words[0] == "Time"):     #if the line has the time units used
                timeStep = float(words[3])*rydToSec         #strip that value into timeStep variable
                print "time step: ", timeStep, "sec"    

    N = len(x)                          #get length of time series
    tpoints = []                        #time points list
    fpoints = []                        #frequency points lists
    for n in range(N):                  #loop through values
        tpoints.append(n*timeStep)      #add correct time and frequency values into lists
        fpoints.append(n/timeStep)      # ... 


    #add this file's data into a plot
    print float(iter)/numfiles
    subplot(2,1,1)                                              #creat sublot 1
    #plot(tpoints,x,color=(0.5,0.5,float(iter)/numfiles) )
    #plot(tpoints,y,color=(float(dist)*0.05,0.5,0.5))   
    if (iter%3 == 0):
        plot(tpoints, z)
       #plot(tpoints,z,color=(1-float(iter)/numfiles,0,float(iter)/numfiles))                                             #plot dipole z vs t
        l.append(str(dist) + ' : z' )                               #append z specifier to list
    
    ylabel('dipole moment (arb units)')                         #plot labels (could be brought out of loop..)
    xlabel('time (seconds)')        
                                                                
    subplot(2,1,2)                                              #frequency sublot
    plot(fpoints[0:len(fpoints)/2+1],np.abs((fft.rfft(z))))     #plot fourier transform for first N/2 values
    ylabel('fourier transform along E field direction')         #labels
    xlabel('frequency (1/sec)')
    iter += 1
    zs.append(z)

''''   
import csv
resultFile = open("output.csv",'wb')
wr = csv.writer(resultFile, dialect='excel')
wr.writerows(zs)   
''' 
legend(l)               #append legend to plot
show()                  #show graph (comment out if you dont want it freezing)