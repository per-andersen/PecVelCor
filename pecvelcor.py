import numpy as np
from astropy.table import Table
from astropy.io import fits
import matplotlib.pyplot as plt

'''
Functions to derive the uncertainty associated with correcting observed redshifts in the CMB frame
for peculiar motion. Which function(s) is/are needed depends on the approach of the cosmological fit.

The simplest solution is to propagate the redshift uncertainty from get_sigma_redshift_pecvel() as
the uncertainty is constant in redshift. If that is not possible the uncertainty can be shifted
to the distance modulus (or whatever magnitude you have in your code, possibly m_B). For magnitudes
however the uncertainty is not constant, but decreases with redshift (corresponding to the slope of
the mu vs redshift plot). The function to call is then get_sigma_mu_pecvel() which takes as input a
numpy like array that contains the cosmological redshifts of the SNe (z_Hd or zcor) and returns a
numpy array with the magnitude uncertainties for each redshift to propagate. If you are working
with the full covariance matrix this array is added to the diagonal of that covariance matrix.

IMPORTANT: If you are already propagating an uncertainty from the linear theory parameter beta,
which scales linearly with the peculiar velocites, set sigma_lineartheory to zero.
'''

def get_sigma_redshift_pecvel():
	sigma_lineartheory = 150.
	sigma_shotnoise = 125.
	sigma_missingdata = 100.

	sigma_redshift_pecvel = np.sqrt(sigma_lineartheory**2 + sigma_shotnoise**2 + sigma_missingdata**2) / 2.99792e5

	return sigma_redshift_pecvel

def get_sigma_mu_pecvel(redshifts):

	sigma_redshift_pecvel = get_sigma_redshift_pecvel()

	sigma_mu_pecvel = sigma_redshift_pecvel * 5./np.log(10.)*(1.+redshifts)**2 / (redshifts*(1+0.5*redshifts))

	return sigma_mu_pecvel

if __name__ == '__main__':

	# Some basic tests

	redshifts = np.linspace(0.001,0.1,100)

	sigma_mu_pecvel = get_sigma_mu_pecvel(redshifts)

	plt.figure()
	plt.xlim((0.,.1))
	plt.ylim((0,0.2))
	plt.xlabel('Cosmological redshift',size='large')
	plt.ylabel(r'$\sigma_{\mu ,pec}$',size='x-large')
	plt.plot(redshifts,sigma_mu_pecvel)
	plt.show()