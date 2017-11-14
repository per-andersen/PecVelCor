
def read_fits_file():

	data = Table.read(root_dir + 'data/DES+nearby/DES+nearby_20171106.fits')
	#print data.keys()
	zcor = data['zcor']
	zcmb = data['zcmb']
	name = data['name']
	ra = data['ra']
	dec = data['dec']

	return zcor, zcmb, name, ra, dec

def read_fitres_file():
	fit_results = []
	fitres_file = '/Users/perandersen/Downloads/FITOPT000.FITRES'
	#fitres_file = '/Users/perandersen/Downloads/LOWZ.FITRES.TEXT'
	# scan header for variable names and number of lines to skip
	header = open(fitres_file, 'r').readlines()[:30]
	for i, line in enumerate(header):
		if line.startswith('NVAR'):
			nvar = int(line.split()[1])
			varnames = header[i+1].split()
		elif line.startswith('SN'):
			nskip = i
			break
	# read file
	

	fitres = np.genfromtxt(fitres_file, skip_header=nskip, usecols = tuple(range(1,nvar+1)), names = varnames[1:],dtype=None)
	fit_results.append(fitres)

	#print varnames
	data_fitres = fitres[:500]
	#print fitres
	zcmb_fitres = data_fitres['zCMB']
	cid_fitres = data_fitres['CID']
	vpec_fitres = data_fitres['VPEC']
	zhd_fitres = data_fitres['zHD']

	return zcmb_fitres, cid_fitres, vpec_fitres, zhd_fitres

def combine_data():
	zcor, zcmb, name, ra, dec = read_fits_file()
	zcmb_fitres, cid_fitres, vpec_fitres, zhd_fitres = read_fitres_file()
	bad_ones = 0
	
	combined_fitres = np.empty([0,3])
	combined_fits = np.empty([0,3])
	names_fitres = np.array([],dtype=np.string_)
	names_fits = np.array([],dtype=np.string_)

	for i, candidateid in enumerate(cid_fitres):
		for j, snname in enumerate(name):
			if candidateid in str(snname):

				combined_fitres = np.append(combined_fitres,np.array([[zcmb_fitres[i],zhd_fitres[i],vpec_fitres[i]]]),axis=0)
				vpec_temp = 2.99792e5 * ( (1.+zcor[j])/(1.+zcmb[j])-1. )
				combined_fits = np.append(combined_fits,np.array([[zcmb[j],zcor[j],vpec_temp]]),axis=0)

				names_fitres = np.append(names_fitres, candidateid)
				names_fits = np.append(names_fits, snname)

				#print name[j], zcmb[j]
				#print cid_fitres[i], zcmb_fitres[i]
				if np.abs(1. - zcmb[j] / zcmb_fitres[i]) > 0.001:
					print "ARGH ARGH"

				#print cid_fitres[i], zcmb_fitres[i], vpec_fitres[i], zcmb_fitres[i]+vpec_fitres[i]/2.99792e5, zhd_fitres[i]
				#print name[j], zcor[j]
				#print zcor[j] / (zcmb_fitres[i]+vpec_fitres[i]/2.99792e5), zcor[j] / zhd_fitres[i]

				if np.abs(1. - zcor[j] / zhd_fitres[i]) > 0.01:
					bad_ones += 1
					#print snname, np.round(ra[j],3), np.round(dec[j],3), np.round(np.abs(1. - zcor[j] / zhd_fitres[i]) * 100.,3)
					#print snname, zcmb[j], zcor[j], zcmb_fitres[i], zhd_fitres[i], vpec_fitres[i]
				
				#print ""

	return combined_fitres, combined_fits, names_fitres, names_fits

def plot_percentage_difference_zhd():
	plt.figure()
	plt.plot(combined_fits[:,1],1.-combined_fits[:,1]/combined_fitres[:,1],'bx')
	plt.grid(True)
	plt.xlabel('zcosmo')
	plt.ylabel(r'Relative difference (1 - $z_1 / z_2$)')
	plt.ylim((-0.1,0.1))

def plot_delta_zhd_vpec():
	f, axarr = plt.subplots(2, sharex=True)	
	axarr[0].set_ylabel(r'$\Delta v_{pec}$  [km/s]',size='x-large')
	axarr[0].set_ylim((-500,500))
	axarr[0].plot(combined_fits[:,0],combined_fits[:,2]-combined_fitres[:,2],'bx')

	axarr[1].set_ylabel(r'$\Delta z_{HD}$',size='x-large')
	axarr[1].set_ylim((-0.0016,0.0016))
	axarr[1].plot(combined_fits[:,0],combined_fits[:,1]-combined_fitres[:,1],'bx')
	plt.subplots_adjust(left=0.16,bottom=0.1,right=0.95,top=0.96,hspace=0.)
	plt.xlim((0,0.1))
	plt.xlabel('zcmb',size='large')
	#plt.savefig('compareplot.pdf',format='pdf')


def problematic_data():

	delta_vpec = combined_fits[:,2]-combined_fitres[:,2]
	problems_fitres = combined_fitres[:,2][np.abs(delta_vpec)>15.]
	problems_fits = combined_fits[:,2][np.abs(delta_vpec)>15.]
	problems = names_fits[np.abs(delta_vpec)>15.]
	problems_vel = combined_fitres[:,2][np.abs(delta_vpec)>15.]
	problems_delta = delta_vpec[np.abs(delta_vpec)>15.]

	for ii in range(len(problems[problems_vel!=0.])):
		print problems[problems_vel!=0.][ii], problems_delta[problems_vel!=0.][ii], problems_fits[problems_vel!=0.][ii],problems_fitres[problems_vel!=0.][ii]
	

if __name__ == '__main__':
	root_dir='/Users/perandersen/Github/SpectroSN-master/'
	combined_fitres, combined_fits, names_fitres, names_fits = combine_data()
	
	#problematic_data()
	#plot_percentage_difference_zhd()
	#plot_delta_zhd_vpec()
	#plt.show()