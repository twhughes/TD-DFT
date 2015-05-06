import os 
print os.getcwd()
with open('./GrapheneBilayer.tddft-out2', 'r') as f:
    data = f.readlines()
    x = []
    y = []
    z = []
    tot = []
    for line in data:
        words = line.split()
        if (len(words) > 4 and words[0] == 'DIP'):
            x.append(words[3])
            y.append(words[4])
            z.append(words[5])
            tot.append(float(words[3])+float(words[4])+float(words[5]))
            
from pylab import plot, show, xlabel, ylabel, legend

from numpy.fft import fft
plot(x)
plot(y)
plot(z)
#plot(tot)
ylabel('dipole moment (arb units)')
legend(['x','y','z'])
xlabel('time steps')
show()
'''
plot(fft(x))
plot(fft(y))
plot(fft(z))
ylabel('dipole moment (arb units)')
legend(['x','y','z'])
xlabel('time steps')
show()
'''
