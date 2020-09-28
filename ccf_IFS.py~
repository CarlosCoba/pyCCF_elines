#!/usr/bin/python
from matplotlib import pyplot as plt
import numpy as np
import pyfits
from astropy.io import fits
from PyAstronomy import pyasl
from lmfit import Minimizer, Parameters, report_fit
import sys
#import scipy.interpolate as sci
from scipy.signal import square, sawtooth, correlate
import scipy
from r_ccf import r_ccf
import scipy.interpolate as sci

#import this scripts
#from set_params import set_parameters as params
#from fit_gaussians import fitting as fit







nargs=len(sys.argv)

def main_program():
	if (nargs==7):
		cube=sys.argv[1]
		z=float(sys.argv[2])
		fwhm=float(sys.argv[3])
		name=sys.argv[4]
		plot=int(sys.argv[5])
		dir_save=sys.argv[6]
		return cube,z,fwhm,name,plot,dir_save
		
	else:
		print 'USE: CCF_IFS.py GAS.cube.fits z fwhm OUTPUT plot dir_save'
		exit()




speed_light=299792.485#km/s


def ccf(cube,z,fwhm,name,plot,dir_save):
	#import this scripts
	from set_params import set_parameters as params
	from fit_gaussians import fitting as fit

 
	data, head = fits.getdata(cube, header=True)#,memmap=True)
	wave = head['CRVAL3'] + head['CDELT3'] * np.arange(head['NAXIS3'])
	#wave = head['CRVAL3'] + head['CD3_3'] * np.arange(head['NAXIS3'])


	[nz,ny,nx]=data.shape
	#elines to fit
	eline_params=params(z)
	elines=eline_params[0]
	#z=eline_params[1]
	n_elines=len(elines)
	n_elines=np.count_nonzero(~np.isnan(elines))


	#Remove NANs from cube:
	data[np.isnan(data)]=0

	#Prepare arrays to save data
	MEDIAN=np.ones((ny,nx))*np.nan
	MEAN=np.ones((ny,nx))*np.nan
	BISECTOR=np.ones((9,ny,nx))*np.nan
	RAD_VEL=np.ones((ny,nx))*np.nan
	DISPERSION=np.ones((ny,nx))*np.nan
	DISPERSION_ELINE=np.ones((ny,nx))*np.nan

	for jj in range(ny):
		for ii in range(nx):

			try:
				flux=data[:,jj,ii]
				print jj,ii


				l0,l1,l2=elines
				lamb_min,lamb_max=np.nanmin(elines),np.nanmax(elines)
				l0_z,l1_z,l2_z=l0*(1+z),l1*(1+z),l2*(1+z)



				waves=[l0,l1,l2]
				min_w=np.nanmin(waves)
				max_w=np.nanmax(waves)

				index=np.where((wave>(lamb_min*(1+z-1500./speed_light))) & (wave<(lamb_max*(1+z+1500./speed_light))))
				#wave range to estimate CCF
				wave_range=wave[index]
				flux_range=flux[index]
				max_flux=np.nanmax(flux_range)
				#Normalize the flux to the maximum value
				flux_range=flux_range/max_flux
				n_data=len(wave_range)
				min_w,max_w=np.min(wave_range),np.max(wave_range)


				#A guess of SN, is this necesary ?
				std_noise=np.std(flux_range[0:20])
				noise=std_noise*np.random.standard_normal(len(flux_range))*0

				SN=(max_flux)/std_noise
				mean_flux=np.nanmean(flux_range)

				#plt.plot(wave_range,flux_range)
				#plt.show()

				if np.isfinite(SN)==True:
					#Perform the fit to estimate the best values for A,sigma
					best_fit=fit(wave_range,flux_range,z,0)

					if n_elines==3:


						sigma=best_fit[0]
						A0,A1,A2=best_fit[1][0],best_fit[1][1],best_fit[1][2]
						z00=best_fit[2]
						l_0,l_1,l_2=l0*(1+z00),l1*(1+z00),l2*(1+z00)
						l_0,l_1,l_2=l0*(1+z),l1*(1+z),l2*(1+z)
						#model_spectra=best_fit[3]


						# Create the template
						xt=wave_range
						xt=np.linspace(lamb_min*(1+z-10000./speed_light),lamb_max*(1+z+10000./speed_light),100*n_data)
						def ccf_model(x):
							#A0,A1,A2=0.8/3.,1,0.8
							sigma_model=fwhm/2.345
							f0 =A0*np.exp(-(x-l_0)**2/(2.*sigma_model**2))
							f1 = A1*np.exp(-(x-l_1)**2/(2.*sigma_model**2))
							f2 = A2*np.exp(-(x-l_2)**2/(2.*sigma_model**2))
							return f0+f1+f2



	#########
					if n_elines==2:


						sigma=best_fit[0]
						A0,A1=best_fit[1][0],best_fit[1][1]
						z00=best_fit[2]
						model_spectra=best_fit[3]
						l_0,l_1=l0*(1+z00),l1*(1+z00)

						xt=wave_range
						xt=np.linspace(lamb_min*(1+z-10000./speed_light),lamb_max*(1+z+10000./speed_light),100*n_data)

						def ccf_model(x):
							sigma_model=fwhm/2.345
							f0 =A0*np.exp(-(x-l0_z)**2/(2.*sigma_model**2))
							f1 = A1*np.exp(-(x-l1_z)**2/(2.*sigma_model**2))
							return f0+f1



					# Carry out the cross-correlation.
					# The rv-range is -700 - +700 km/s in steps of 5 km/s.

					rv=np.arange(-700,700,5)
					ccf_vel=[]
					for j,delta_v in enumerate(rv):
						fi = sci.interp1d(xt*(1.0 + delta_v/speed_light), ccf_model(xt))
						# The values to interpolate must be in the range covered by xt*(1.0 + delta_v/speed_light)
						lags,r_=r_ccf(fi(wave_range),flux_range,False,1)
						ccf_vel.append(r_[0])

					# The cross correlation in velocity space:
					ccf_vel=np.asarray(ccf_vel)
					# Normalize the ccf:
					ccf_vel=ccf_vel/np.nanmax(ccf_vel)

					# Estimate the peak
					max_ccf_index=ccf_vel.argmax()
					max_ccf=ccf_vel[max_ccf_index]
					mean_vel=rv[max_ccf_index]
					v_peak=mean_vel

					# PERFORM A QUICK FIT TO ESTIMATE VELOCITY AND SIGMA 
					# AROUND THE PEAK
					from fit_gaussians import G1
					best_fit_ccf=G1(rv[rv>v_peak-500],ccf_vel[rv>v_peak-500],max_ccf,v_peak,100,0)
					A_ccf,v_ccf,sigma_ccf=best_fit_ccf[0],best_fit_ccf[1],abs(best_fit_ccf[2])

					v_centroid=v_ccf
					#shift velocity to the centroid
					rv=rv-v_centroid


					###########################################################
					#This is the real cross correlation function with 20 lags:#
					###########################################################
					lags,r_=r_ccf(ccf_model(wave_range),flux_range,True,20)			
					


					# shift cross correlation to the region of inrerenst -450km/s<rv<450 km/s
					index_vel=np.where((rv>-450) & (rv<450))
					rv=rv[index_vel]
					ccf_vel=ccf_vel[index_vel]

					#Calculate the bisectors
					x_bis=[]
					y_bis=[]
					delta_v=[]
					step=0.1
					for i in np.arange(step,1,step):
						li=(1-i)
						ls=(1+step-i)
						index=np.where((ccf_vel>li) & (ccf_vel<=ls))
						vels=np.mean(rv[index])
						x_bis.append(vels)
						y_bis.append(li)
						delta_v.append(vels)


					#SAVE CCF and BISECTORS
					J,I= int(jj),int(ii)
					MEAN[J][I]= np.nanmean(delta_v)
					MEDIAN[J][I]=np.nanmedian(delta_v)
					BISECTOR[:,J,I]=delta_v





					#SAVE VELOCITY PARAMETERS FROM FIT
					J,I= int(jj),int(ii)
					RAD_VEL[J][I]= v_ccf+0
					DISPERSION[J][I]=sigma_ccf+0
					DISPERSION_ELINE[J][I]=sigma+0


					if plot == 1:	
							plt.figure(figsize=(8,4))

							plt.subplot(131)
							plt.plot(xt,ccf_model(xt),label="template")
							plt.plot(wave_range,flux_range,'k-',label="spectra")
							plt.xlabel('wavelength')
							plt.legend(loc="upper left")

							plt.subplot(132)
							plt.plot(lags,r_,"r-",label='r_ccf')
							plt.axhline(y=0, color='r', linestyle='--')
							plt.xlabel('wavelength')
							plt.ylabel('r_ccf')
							plt.ylim(-1,1)
							plt.legend(loc="upper left")

							plt.subplot(133)
							plt.plot(rv,ccf_vel,'k+',label='norm. ccf_vel')
							G1(rv,ccf_vel,max_ccf,mean_vel,100,1)
							plt.xlabel('velocity')
							plt.ylabel('r_ccf')
							plt.axhline(y=0, color='r', linestyle='--')
							plt.ylim(-1,1.5)
							plt.plot(x_bis,y_bis,'g.')
							plt.legend(loc="upper left")
							plt.tight_layout()

							plt.show()


			except (ValueError,ZeroDivisionError,IOError,TypeError):
			#except (1):
				pass




	###
	# save all
	###

	hdu_1 = fits.PrimaryHDU()
	hdu_1.data=MEDIAN
	hdu_1.header=head
	hdu_1.writeto("%s%s.ccf.median.fits"%(dir_save,name))


	hdu_2 = fits.PrimaryHDU()
	hdu_2.data=MEAN
	hdu_2.header=head
	hdu_2.writeto("%s%s.ccf.mean.fits"%(dir_save,name))


	hdu_3 = fits.PrimaryHDU()
	hdu_3.data=BISECTOR
	hdu_3.header=head
	hdu_3.writeto("%s%s.ccf.bisector.cube.fits"%(dir_save,name))

	hdu_4 = fits.PrimaryHDU()
	hdu_4.data=RAD_VEL
	hdu_4.header=head
	hdu_4.writeto("%s%s.ccf.vel.fits"%(dir_save,name))

	hdu_5 = fits.PrimaryHDU()
	hdu_5.data=DISPERSION
	hdu_5.header=head
	hdu_5.writeto("%s%s.ccf.sigma_ccf.fits"%(dir_save,name))

	hdu_6 = fits.PrimaryHDU()
	hdu_6.data=DISPERSION_ELINE
	hdu_6.header=head
	hdu_6.writeto("%s%s.ccf.sigma_eline.fits"%(dir_save,name))


	print "done"

if __name__ == "__main__":
	main=main_program()
	ccf(main[0],main[1],main[2],main[3],main[4],main[5])
#else:
#	ccf()	







