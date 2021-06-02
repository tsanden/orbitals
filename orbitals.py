import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
import sys

# first Bohr radius
a0 = 5.29*10**(-11)

class Orbital:

    def __init__(self, r, theta):
        self.r = r
        self.theta = theta

    def prob_100(self):
        return np.square(1/(np.sqrt(np.pi)*a0**(3/2))*np.exp(-self.r/a0))*(4*np.pi*self.r**2)

    def prob_200(self):
        return np.square((1/(4*np.sqrt(2*np.pi)*a0**(3/2)))*(2-self.r/a0)*np.exp(-self.r/(2*a0)))*(4*np.pi*self.r**2)

    def prob_210(self):
        return np.square((1/(4*np.sqrt(2*np.pi)*a0**(3/2)))*(self.r/a0)*np.exp(-self.r/(2*a0))*np.cos(self.theta))*(4*np.pi*self.r**2)

    def prob_211(self):
        return np.square(1/(24*a0**3)*self.r**4*np.exp(-self.r/a0))        

    def prob_300(self):
        return np.square((1/(81*np.sqrt(3*np.pi)*a0**(3/2)))*(27-18*self.r/a0+2*(self.r/a0)**2)*np.exp(-self.r/(3*a0)))*(4*np.pi*self.r**2)

    def prob_310(self):
        return np.square((np.sqrt(2)/(81*np.sqrt(3*np.pi)*a0**(3/2)))*(6-self.r/a0)*self.r/a0*np.exp(-self.r/(3*a0))*np.cos(self.theta))*(4*np.pi*self.r**2)

    def prob_320(self):
        return np.square((1/(81*np.sqrt(6*np.pi)*a0**(3/2)))*(self.r/a0)**2*np.exp(-self.r/(3*a0))*(3*(np.cos(self.theta))**2-1))*(4*np.pi*self.r**2)

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
def prob_orbital(n=0, l=0, m=0, N=1, radius=10, linear=False, r=[], theta=[]):
    
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

    if n > 3 or n < 1:
        print("Sorry, we cannot find those quantum numbers.")
        sys.exit()
    if n == 1:
        orb_name = 'prob_100'
    elif (n == 2) & (l == 0) & (m == 0):
        orb_name = 'prob_200'
    elif (n == 2) & (l == 1) & (m == 0):
        orb_name = 'prob_210'
    elif (n == 2) & (l == 1) & (m == 1):
        orb_name = 'prob_211'      
    elif (n == 3) & (l == 0) & (m == 0):
        orb_name = 'prob_300'        
    elif (n == 3) & (l == 1) & (m == 0):
        orb_name = 'prob_310'   
    elif (n == 3) & (l == 2) & (m == 0):
        orb_name = 'prob_320'              

    probability = orbital.calc_prob(orb_name, r, theta)

    return probability, r, theta, orb_name

# function to plot orbitals:
# n, l, m - quantum numbers
# N - sampling points
# radius - in number of Bohr radiuses
# scatter - True by default, if false plot contour plot
# axis - turn axis on and off
def plot_orbital(n=[0], l=[0], m=[0], sample=1, radius=10, sample_to_values=10000, scatter=True, axis=True, cmap='cubehelix', rows=1, cols=1):

    if (len(n) != len(l)) or (len(n) != len(m)):
         print('Make sure that n EQUALS l.')
         sys.exit()

    if scatter:
        
        probability = np.zeros((len(n), sample_to_values*sample))
        r = np.zeros((len(n), sample_to_values*sample))
        theta = np.zeros((len(n), sample_to_values*sample))
        prev_prob = np.zeros((len(n), sample_to_values*sample))
        for k, q in enumerate(n):
            probability[k], r[k], theta[k], orb_name = prob_orbital(n=n[k], l=l[k], m=m[k], N=sample_to_values*sample, radius=radius)

        # feature scaling
        for k, g in enumerate(n):
            probability[k] = probability[k]/probability[k].sum()
        
        elem = np.zeros((len(n),sample_to_values*sample,2))
        r_c = np.zeros((len(n),sample))
        theta_c = np.zeros((len(n),sample))

        for k, g in enumerate(n):            
            elem[k] = list(zip(r[k],theta[k]))

        for k, g in enumerate(n):
            coord = np.random.choice([','.join(map(str,i)) for i in elem[k]], size=sample, replace=False, p=probability[k])
            elem_mat = [i.split(',') for i in coord]
            r_c[k] = [float(i[0]) for i in elem_mat]
            theta_c[k] = [float(i[1]) for i in elem_mat]

        figs, axs = plt.subplots(rows, cols, subplot_kw=dict(projection='polar'), figsize=(10,10))

        count = 0
        for i in range(rows):
            for j in range(cols):
                if isinstance(axs, np.ndarray):
                    if axs.ndim == 1:
                        if count <= len(n)-1:
                            axs[j].scatter(theta_c[count], r_c[count], alpha=1, s=2, c='#000000')
                        if not axis:
                            axs[j].axis('off')                  
                    elif axs.ndim == 2:
                        if count <= len(n)-1:
                            axs[i,j].scatter(theta_c[count], r_c[count], alpha=1, s=2, c='#000000')
                        if not axis:
                            axs[i, j].axis('off')                  
                else:
                    axs.scatter(theta_c[count], r_c[count], alpha=1, s=2, c='#000000')
                    if not axis:
                        axs.axis('off')
                count += 1

        plt.show()

    else:

        prob = np.zeros((len(n),sample,sample))

        r = np.linspace(0,radius*a0, sample)
        theta = np.linspace(0,2*np.pi, sample)
        
        for i, ir in enumerate(r):
            for j, jtheta in enumerate(theta):
                for k, q in enumerate(n):
                    prob[k][i][j], r0, theta0, orb_name = prob_orbital(n=n[k], l=l[k], m=m[k], N=sample,radius=radius, linear=True, r=np.array([ir]), theta=np.array([jtheta]))

        # feature scaling
        for k, g in enumerate(n):
            prob[k] = prob[k]/prob[k].sum()

        figs, axs = plt.subplots(rows, cols, subplot_kw=dict(projection='polar'), figsize=(10,10))
       
        count = 0
        for i in range(rows):
            for j in range(cols):
                if isinstance(axs, np.ndarray):
                    if axs.ndim == 1:
                        if count <= len(n)-1:
                            axs[j].contourf(theta, r, prob[count], levels=1000, cmap=cmap)  
                        if not axis:
                            axs[j].axis('off')                  
                    elif axs.ndim == 2:
                        if count <= len(n)-1:
                            axs[i,j].contourf(theta, r, prob[count], levels=1000, cmap=cmap)
                        if not axis:
                            axs[i, j].axis('off')                  
                else:
                    axs.contourf(theta, r, prob[count], levels=1000, cmap=cmap)
                    if not axis:
                        axs.axis('off')
                count += 1

        plt.show()

    return None