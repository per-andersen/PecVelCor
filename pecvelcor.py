import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt


# ToDo:
# - Check if lowz is defined by sample or redshift cut
# -

def get_sigma_redshift_pecvel():
	sigma_lineartheory = 150.
	sigma_shotnoise = 125.
	sigma_missingdata = 100.

	sigma_redshift_pecvel = np.sqrt(sigma_lineartheory**2 + sigma_shotnoise**2 + sigma_missingdata**2) / 2.99792e5

	return sigma_redshift_pecvel

def get_sigma_mu_pecvel(redshifts):

	sigma_redshift_pecvel = get_sigma_redshift_pecvel()

	sigma_mu_pecvel = sigma_redshift_pecvel * 5./np.log(10.)*((1.+redshifts)/redshifts/(1.+.5*redshifts))

	return sigma_mu_pecvel

if __name__ == '__main__':
	redshifts = np.linspace(0.001,0.1,100)

	sigma_mu_pecvel = get_sigma_mu_pecvel(redshifts)

	plt.figure()
	plt.xlim((0.,.1))
	plt.ylim((0,0.2))
	plt.xlabel('Cosmological redshift',size='large')
	plt.ylabel(r'$\sigma_{\mu ,pec}$',size='x-large')
	plt.plot(redshifts,sigma_mu_pecvel)
	plt.show()