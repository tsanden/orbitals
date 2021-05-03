import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri

# first Bohr radius
a0 = 5.29*10**(-11)

class Orbital:

    def __init__(self, r, theta):
        self.r = r
        self.theta = theta

    def prob_1s(self):
        return np.square(1/(np.sqrt(np.pi)*a0**(3/2))*np.exp(-self.r/a0))*(4*np.pi*self.r**2)

    def prob_2s(self):
        return (1/(8*a0)*(self.r/a0)**2*(2-self.r/a0)**2*np.exp(-self.r/a0))

    def prob_2p(self):
        return np.square((1/(4*np.sqrt(2*np.pi)*a0**(3/2)))*(self.r/a0)*np.exp(-self.r/(2*a0))*np.cos(self.theta))*(4*np.pi*self.r**2)
   
    def calc_prob(self, name: str, r, theta):
        calc = f"{name}"
        self.r = r
        self.theta = theta
        if hasattr(self, calc) and callable(func := getattr(self, calc)):
           return func()

# function to return probabilities given:
# n, l, m - 0 by default, quantum numbers
# N - 1 by default, sampling points
# radius - 10 by default, in number of Bohr radiuses
# linear - False by default, randomly distributed points
# r - [] by default, pass in radial values
# theta - [] by default, pass in theta values
def prob_orbital(n=0,l=0,m=0,N=1,radius=10, linear=False, r=[], theta=[]):

    orb_name = None

    if linear:
        if len(r) == 0:
            r = np.linspace(0,radius*a0,N)
        if len(theta) == 0:
            theta = np.linspace(0,2*np.pi,N)
    else:
        if len(r) == 0:
            r = radius*a0*np.sqrt(np.random.rand(N))
        if len(theta) == 0:
            theta = 2*np.pi*np.random.rand(N)
    
    orbital = Orbital(r, theta)

    if n > 2 or n < 1:
        print("Sorry, we cannot find those quantum numbers.")
    if n == 1:
        orb_name = 'prob_1s'
    elif (n == 2) & (l == 0):
        orb_name = 'prob_2s'
    elif (n == 2) & (l == 1):
        orb_name = 'prob_2p'

    probability = orbital.calc_prob(orb_name, r, theta)

    return probability, r, theta, orb_name

# function to plot orbitals:
# n, l, m - quantum numbers
# N - sampling points
# radius - in number of Bohr radiuses
# scatter - True by default, if false plot contour plot
# axis - turn axis on and off
def plot_orbital(n=0,l=0,m=0,sample=1,radius=10, sample_to_values=100, scatter=True, axis=True, cmap='cubehelix'):

    if scatter:
        
        probability, r, theta, orb_name = prob_orbital(n=n,l=l,m=m,N=sample_to_values*sample,radius=radius)
        
        # feature scaling
        probability = (probability-probability.min())/(probability.max()-probability.min())
        probability = probability/probability.sum()
        print(probability.max())
        
        elem = []
        elem = list(zip(r,theta))

        coord = np.random.choice([','.join(map(str,i)) for i in elem], size=sample, replace=False, p=probability)

        elem_mat = [i.split(',') for i in coord]

        r_c = [float(i[0]) for i in elem_mat] 
        theta_c = [float(i[1]) for i in elem_mat] 

        fig = plt.figure()
        ax = fig.add_subplot(projection='polar')
        ax.scatter(theta_c, r_c, alpha=1, s=2, c='#000000')

        ax.set_title("Hydrogen " + orb_name.split('_')[1] + " density")
        if not axis:
            plt.axis('off')
        plt.show()

    else:

        prob = np.zeros((sample,sample))
        
        r = np.linspace(0,radius*a0, sample)
        theta = np.linspace(0,2*np.pi, sample)
        
        for i, ir in enumerate(r):
            for j, jtheta in enumerate(theta):
                prob[i][j], r0, theta0, orb_name = prob_orbital(n=n,l=l,m=m, N=sample,radius=radius, linear=True, r=np.array([ir]), theta=np.array([jtheta]))

        # feature scaling
        prob = (prob-prob.min())/(prob.max()-prob.min())
        prob = prob/prob.sum()        

        fig, ax = plt.subplots(subplot_kw=dict(projection='polar'))
        
        ax.contourf(theta, r, prob, levels=255, cmap=cmap)
        if not axis:
            plt.axis('off')
        plt.show()

    return None