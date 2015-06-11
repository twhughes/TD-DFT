import os 
os.chdir('../')
print os.getcwd()

files = ['./Graphene2.tddft-in0.0001','./Graphene2.tddft-in0.0002' ,'./Graphene2.tddft-in0.0003', \
         './Graphene2.tddft-in0.0004','./Graphene2.tddft-in0.0005','./Graphene2.tddft-in0.0006', \
         './Graphene2.tddft-in0.0007','./Graphene2.tddft-in0.0008','./Graphene2.tddft-in0.0009', \
         './Graphene2.tddft-in0.0011','./Graphene2.tddft-in0.0012','./Graphene2.tddft-in0.0013', \
         './Graphene2.tddft-in0.0014','./Graphene2.tddft-in0.0015','./Graphene2.tddft-in0.0016', \
         './Graphene2.tddft-in0.0017','./Graphene2.tddft-in0.0018','./Graphene2.tddft-in0.0019', \
         ]
timeStep = 1
rydToSec = 4.8377687e-17
S = [0]
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
    dz = []
    for i in range(len(z)-1):
        dz.append(z[i+1] - z[i])
    for i in range(len(z)-1):
            dz.append(z[i+1] - z[i])
   # plot(tpoints,z)
#    plot(tpoints[:-1],dz)
   # show()

    ylabel('dipole moment (arb units)')
   # legend(['x','y','z'])
    xlabel('time (seconds)')
    F_dz = np.abs((f.rfft(dz)))
    #plot(F)
    ylabel('derivative fourier')
   # show()
    #plot(np.abs((f.rfft(z))))
   # ylabel('fourier transform along E field direction')
    legend(['0.002','0.003','0.004','0.005','0.006', '0.007'])
#    xlabel('frequency (1/sec)')
   # show()
    s = 0
    for n in range(1):
        s += F_dz[n]
       # print s
    S.append(s)
  
epoints = [0.0,0.0001,0.0002,0.0003,0.0004,0.0005,0.0006,0.0007,0.0008,0.0009,0.0011,0.00120,0.0013,0.0014,0.0015,0.0016,0.0017,0.0018,0.0019]  
#show()
line = np.polyfit(epoints, S,1)
pls = []
lrange = np.arange(0, max(epoints), max(epoints)/100.0)
for l in lrange:
    pls.append(line[1]+l*line[0])
print 'slope = ', line[0]
plot(lrange, pls)
xlabel('E field (AMU)')
ylabel('SUM (first 4) Abs(fourier(dz))')
scatter(epoints, S)
show()

