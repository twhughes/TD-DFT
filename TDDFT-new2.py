import os 
os.chdir('../')
print os.getcwd()

#files = \
   #[#'./Graphene2.tddft-new2-1e-05','./Graphene2.tddft-new2-2e-05' ,'./Graphene2.tddft-new2-3e-05', \
         #'./Graphene2.tddft-new2-4e-05','./Graphene2.tddft-new2-5e-05','./Graphene2.tddft-new2-6e-05', \
         #'./Graphene2.tddft-new2-7e-05', './Graphene2.tddft-new2-8e-05'\
        #]
files =  ['./Graphene2.tddft-new2-1','./Graphene2.tddft-new2-2','./Graphene2.tddft-new2-3', \
         './Graphene2.tddft-new2-4','./Graphene2.tddft-new2-5','./Graphene2.tddft-new2-6', \
         './Graphene2.tddft-new2-7','./Graphene2.tddft-new2-8' \
         ]
timeStep = 1
rydToSec = 4.8377687e-17
S = []
for fname in files:
    print fname
    with open(fname, 'r') as f:
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

    from pylab import plot, show, xlabel, ylabel, legend, subplot, scatter
    import numpy as np
    import numpy.fft as f
    #subplot(2,1,1)
    #plot(tpoints,x)
    #plot(tpoints,y)
    #plot(tpoints,z)

    #take derivative in z
    dzz = []
    dz = []
    for i in range(len(z)-1):
        dzz.append(z[i+1] - z[i])
    for i in range(len(z)-2):
        dz.append(dzz[i+1] - dzz[i])
#    plot(tpoints,z)
#    plot(tpoints[:-1],dz)
   # show()

    ylabel('dipole moment (arb units)')
   # legend(['x','y','z'])
    xlabel('time (seconds)')
    F_dz = 0
    for n in range(len(dz)):
        F_dz += dz[n]
    F_dz /= len(dz)
#    F_dz = np.abs((f.rfft(dzz)))
    #plot(F)
    ylabel('derivative fourier')
   # show()
    #plot(np.abs((f.rfft(z))))
   # ylabel('fourier transform along E field direction')
    legend(['0.002','0.003','0.004','0.005','0.006', '0.007'])
#    xlabel('frequency (1/sec)')
   # show()
#    s = 0
#    for n in range(1):
#        s += F_dz[n]
#       # print s
    S.append(F_dz)
  
epoints = [ 1e-4,2e-4,3e-4,4e-4,5e-4,6e-4,7e-4,8e-4]#,7,8,9,0.0011,0.00120,0.0013,0.0014,0.0015,0.0016,0.0017,0.0018,0.0019]  
#show()
line = np.polyfit(epoints, S,1)
pls = []
lrange = np.arange(0, max(epoints), max(epoints)/100.0)
for l in lrange:
    pls.append(line[1]+l*line[0])
print 'slope = ', line[0]
plot(lrange, pls)
xlabel('E field (Atomic units)')
ylabel('Dipole moment (DC)')
scatter(epoints, S)
show()

