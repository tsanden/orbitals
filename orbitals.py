import numpy as np
import matplotlib.pyplot as plt

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
# n, l, m - quantum numbers
# N - sampling points
# radius - in number of Bohr radiuses
def prob_orbital(n=0,l=0,m=0,N=1,radius=10):

    orb_name = None

    # radial values
    # r = np.linspace(0,radius*a0,N)
    # uniform random values along
    r = radius*a0*np.sqrt(np.random.rand(N))

    # theta angles
    #theta = np.linspace(0,2*np.pi,N)
    theta = 2*np.pi*np.random.rand(N)
        
    orbital = Orbital(r, theta)

    if n > 2 or n < 1:
        print('Sorry n is too big or too small! Please try again.')
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
def plot_orbital(n=0,l=0,m=0,sample=1,radius=10, scatter=True):

    probability, r, theta, orb_name = prob_orbital(n=n,l=l,m=m,N=500000,radius=radius)
    probability = probability/probability.sum()

    elem = []
    elem = list(zip(r,theta))

    coord = np.random.choice([','.join(map(str,i)) for i in elem], size=sample, replace=False, p=probability)

    elem_mat = [i.split(',') for i in coord]

    r_c = [float(i[0]) for i in elem_mat] 
    theta_c = [float(i[1]) for i in elem_mat] 

    fig = plt.figure()

    if scatter:

        ax = fig.add_subplot(projection='polar')
        ax.scatter(theta_c, r_c, alpha=1, s=2)

        ax.set_title("Hydrogen " + orb_name.split('_')[1] + " density")
        plt.show()

    return None