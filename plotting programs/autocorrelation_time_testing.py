import numpy as np
import matplotlib.pyplot as plt
import matplotlib.pyplot as plt
import csv
from astropy.io.votable import parse,from_table,writeto
from astropy.table import Table,vstack,join
import emcee
import matplotlib.patches as mpatches
import subprocess
import scipy.stats as st
from scipy.stats import kde
from os.path import exists
from numba import jit
import celerite
from celerite import terms
from scipy.optimize import minimize


sobject_id=131217003901075

reader=emcee.backends.HDFBackend('NGC_2682_reduction_fixed_photometric/_all_elements_long_run__prior_160106004101321_main_loop.h5')
auto_corr=reader.get_autocorr_time(tol=0,discard=0)


#calculate the autocorrelation time every 100 iterations
tau_steps=[]
for iteration in range(200,700000,100000):
    print(iteration)
    data_to_append=emcee.autocorr.integrated_time(reader.get_chain(discard=200,thin=1)[:iteration,:,:],tol=0)
    tau_steps.append(data_to_append)
tau_steps=np.array(tau_steps)



def autocorr_ml(y, thin=1, c=5.0):
    # Compute the initial estimate of tau using the standard method
    init = emcee.autocorr.integrated_time(y.T,tol=0)
    z = y[:, ::thin]
    N = z.shape[1]

    # Build the GP model
    tau = max(1.0, init / thin)
    kernel = terms.RealTerm(
        np.log(0.9 * np.var(z)),
        -np.log(tau),
        bounds=[(-10.0, 10.0), (-np.log(N), 0.0)],
    )
    kernel += terms.RealTerm(
        np.log(0.1 * np.var(z)),
        -np.log(0.5 * tau),
        bounds=[(-10.0, 10.0), (-np.log(N), 0.0)],
    )

    gp = celerite.GP(kernel, mean=np.mean(z))
    gp.compute(np.arange(z.shape[1]))

    # Define the objective
    def nll(p):
        # Update the GP model
        gp.set_parameter_vector(p)

        # Loop over the chains and compute likelihoods
        v, g = zip(*(gp.grad_log_likelihood(z0, quiet=True) for z0 in z))

        # Combine the datasets
        return -np.sum(v), -np.sum(g, axis=0)

    # Optimize the model
    p0 = gp.get_parameter_vector()
    bounds = gp.get_parameter_bounds()
    soln = minimize(nll, p0, jac=True, bounds=bounds)
    gp.set_parameter_vector(soln.x)

    # Compute the maximum likelihood tau
    a, c = kernel.coefficients[:2]
    tau = thin * 2 * np.sum(a / c) / np.sum(a)
    return tau

#calculate the autocorrelation time every 1000 iterations

# tau_steps_ml=[]
# for iteration in range(100,20000,500):
#     print(iteration)
#     #only input one chain
#     data_to_input=reader.get_chain(discard=100,thin=1)[:iteration,:,0]
#     tau_steps_ml.append(autocorr_ml(data_to_input.T,thin=1))

plt.figure()
plt.plot(np.arange(200,70000,10000),np.mean(tau_steps,axis=1),label='emcee')
# plt.plot(np.arange(100,3000,500),tau_steps_ml,label='auto_regression_model')
plt.xlabel('iteration')
plt.ylabel('autocorrelation time')
plt.legend()
#plt.savefig('autocorrelation_time_comparison.png')
