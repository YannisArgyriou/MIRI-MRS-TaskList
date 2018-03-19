"""Function file miricap401"""

# import python modules
import os
import pickle
import numpy as np
import scipy.special as sp
from scipy.optimize import curve_fit
import scipy.interpolate as scp_interpolate
import matplotlib.pyplot as plt

# Definition
#-import cdps
def get_cdps(band,cdpDir,output='img'):
    """Returns fringe, photom, psf, and resolution CDP fits files """
    subchan_names = {'A':'SHORT','B':'MEDIUM','C':'LONG'}
    if int(band[0]) < 3: detnick = 'MIRIFUSHORT'
    else: detnick = 'MIRIFULONG'


    # fringe flat cdp, compiled by Fred Lahuis
    if int(band[0]) < 3: fringe_subchan_name = '12'
    else: fringe_subchan_name = '34'
    fringe_subchan_name += subchan_names[band[1]]
    fringe_file = os.path.join(cdpDir,\
      'MIRI_FM_%s_%s_FRINGE_06.02.00.fits' % (detnick,fringe_subchan_name))

    # photom cdp, compiled by Bart Vandenbussche
    photom_file = os.path.join(cdpDir,\
      'MIRI_FM_%s_%s_PHOTOM_06.03.02.fits' % (detnick,subchan_names[band[1]]))

    # psf cdp, compiled by Adrian Glauser
    psf_subchan_name = band[0]+subchan_names[band[1]]
    psf_file = os.path.join(cdpDir,\
      'MIRI_FM_%s_%s_PSF_05.02.00.fits' % (detnick,psf_subchan_name))

    # spatial and spectral resolution cdp, compiled by Alvaro Labiano
    if int(band[0]) < 3: resol_subchan_name = '12'
    else: resol_subchan_name = '34'
    resol_file = os.path.join(cdpDir,\
      'MIRI_FM_%s_%s_RESOL_7B.00.00.fits' % (detnick,resol_subchan_name))

    if output == 'filepath':
        return fringe_file,photom_file,psf_file,resol_file
    elif output == 'img':
        from astropy.io import fits
        fringe_img     = fits.open(fringe_file)[1].data        # [unitless]
        photom_img     = fits.open(photom_file)[1].data        # [DN/s * pixel/mJy]
        pixsiz_img     = fits.open(photom_file)[5].data        # [arcsec^2/pix]
        psffits        = fits.open(psf_file)                   # [unitless]
        specres_table  = fits.open(resol_file)[1].data         # [unitless]
        return fringe_img,photom_img,pixsiz_img,psffits,specres_table

# import MIRIM PSFs
def mirimpsfs(workDir=None):
    import glob
    # CV2 PSF measurements
    dataDir = workDir+'CV2_data/LVL2/'

    files = [os.path.basename(i) for i in glob.glob(dataDir+'*')]
    subchannels = ['SHORT','MED','LONG']
    pointings = ['P'+str(i) for i in range(17)]

    MIRIMPSF_dictionary = {}
    for pointing in pointings:
        mylist = []
        for subchannel in subchannels:
            sub = 'MIRM0363-{}-{}'.format(pointing,subchannel)
            mylist.extend([s for s in files if sub in s])
        MIRIMPSF_dictionary['CV2_'+pointing] = list(mylist)

    # CV3 PSF measurements
    dataDir = workDir+'CV3_data/LVL2/'

    files = [os.path.basename(i) for i in glob.glob(dataDir+'*')]

    pointings = ['Q'+str(i) for i in range(17)]
    for pointing in pointings:
        sub = 'MIRM103-{}-SHORT'.format(pointing)
        MIRIMPSF_dictionary['CV3_'+pointing] = [s for s in files if sub in s]

    # Dictionary keys are equivalent to PSF measurements in CV2 and CV3 tests (with different pointings).
    # Dictionary indeces within keys, for CV2 obs. are equivalent to:
    # [0,1,2,3] : SHORT_494,SHORT_495,SHORTB_494,SHORTB_495
    # [4,5,6,7] : MED_494,MED_495,MEDB_494,MEDB_495
    # [8,9,10,11] : LONG_494,LONG_495,LONGB_494,LONGB_495

    # Dictionary indeces within keys, for CV3 obs. are equivalent to:
    # [0,1,2,3] : SHORT_494,SHORT_495,SHORTB_494,SHORTB_495
    # There are only SHORT CV3 PSF measurements

    return MIRIMPSF_dictionary

#-corrections
def OddEvenRowSignalCorrection(sci_img,nRows=1024):
    copy_img = sci_img.copy()
    for nRow in range(nRows-2):
        copy_img[nRow+1,:] = (((sci_img[nRow,:]+sci_img[nRow+2,:])/2.)+sci_img[nRow+1,:])/2.
    return copy_img

# straylight correction
def Shepard2DKernel(R, k):
    """
    Calculates the kernel matrix of Shepard's modified algorithm
    R : Radius of influence
    k : exponent
    """
    xk, yk = np.meshgrid(np.arange(-R/2, R/2+1),np.arange(-R/2, R/2+1))
    d = np.sqrt(xk**2+yk**2)
    w = (np.maximum(0,R-d)/(R*d))**k
    w[d==0]=0
    return w

def straylightCorrection(sci_img,sliceMap,R=50, k=1):
    from astropy.convolution import convolve, Box2DKernel
    """
    Applies a modified version of the Shepard algorithm to remove straylight from the MIRI MRS detector
    img: Input image
    R: Radius of influence
    k: exponent of Shepard kernel
    sliceMap: Matrix indicating slice (band*100+slice_nr) and gap (0) pixels
    """
    w = Shepard2DKernel(R,k) # calculate kernel
    #mask where gap pixels are 1 and slice pixels are 0
    mask = np.zeros_like(sliceMap)
    mask[sliceMap == 0] = 1
    #apply mask to science img
    img_gap = sci_img*mask

    img_gap[img_gap>0.02*np.max(sci_img[sliceMap>0])] = 0 # avoid cosmic rays contaminating result
    img_gap[img_gap<0] = 0 # set pixels less than zero to 0
    img_gap = convolve(img_gap, Box2DKernel(3)) # smooth gap pixels with Boxkernel
    img_gap*=mask # reset sci pixels to 0
    # convolve gap pixel img with weight kernel
    astropy_conv = convolve(img_gap, w)
    # normalize straylight flux by weights
    norm_conv = convolve(mask, w)
    astropy_conv /= norm_conv
    # reinstate gap pixels to previous values
    #astropy_conv[sliceMap==0] = img_gap[sliceMap==0]
    return sci_img-astropy_conv

#-compute
def getSpecR(lamb0=None,band=None,specres_table=None):
    """Return spectral resolution (a.k.a. resolving power)"""
    res_select = 'res_avg'
    resDic = {'res_low':[2,3,4],'res_high':[5,6,7],'res_avg':[8,9,10]}
    subchan_names = {'A':'SHORT','B':'MEDIUM','C':'LONG'}
    resol_subchan_name = band[0]+subchan_names[band[1]]
    band_list = {'1SHORT':0,'1MEDIUM':1,'1LONG':2,\
                 '2SHORT':3,'2MEDIUM':4,'2LONG':5,\
                 '3SHORT':6,'3MEDIUM':7,'3LONG':8,\
                 '4SHORT':9,'4MEDIUM':10,'4LONG':11}

    specR = specres_table[band_list[resol_subchan_name]][resDic[res_select][0]] + \
            specres_table[band_list[resol_subchan_name]][resDic[res_select][1]]*(lamb0-specres_table[band_list[resol_subchan_name]][1]) + \
            specres_table[band_list[resol_subchan_name]][resDic[res_select][2]]*(lamb0-specres_table[band_list[resol_subchan_name]][1])**2

    return specR

def getSpecR_linearR(lamb0=None,band=None):
    """Return spectral resolution (a.k.a. resolving power) assuming a linear relation to wavelength"""
    import mrs_aux as maux
    bandlims = maux.MRS_bands[band]
    Rlims = maux.MRS_R[band]
    specR = (Rlims[1]-Rlims[0])/(bandlims[1]-bandlims[0]) * (lamb0-bandlims[0]) + Rlims[0]
    return specR

def spectral_gridding(band,d2cMaps,specres_table=None,oversampling = 1.):
    # Construct spectral (wavelength) grid
    lambdaMap = d2cMaps['lambdaMap']
    bandlims  = [lambdaMap[np.nonzero(lambdaMap)].min(),lambdaMap[np.nonzero(lambdaMap)].max()]
    #> initialize variables
    lambcens = []
    lambfwhms = []

    #> loop over wavelength bins (bin width defined based on MRS spectral resolution)
    lamb0   = bandlims[0]
    maxlamb = bandlims[1]

    # first iteration
    R       = getSpecR(lamb0=lamb0,band=band,specres_table=specres_table)
    fwhm    = lamb0 / R
    lambcen = lamb0 + (fwhm/2.)/oversampling

    lambcens.append(lambcen)
    lambfwhms.append(fwhm)

    # iterate over spectral range
    done = False
    while not done:
        R = getSpecR(lamb0=lamb0,band=band,specres_table=specres_table)
        fwhm = lamb0 / R
        lambcen = lamb0 + (fwhm/2.)/oversampling
        if (lambcen > maxlamb-(fwhm/2.)/oversampling):
            done = True
        else:
            lamb0 = lambcen + (fwhm/2.)/oversampling

        lambcens.append(lambcen)
        lambfwhms.append(fwhm)

    return np.array(lambcens),np.array(lambfwhms)

def spectral_gridding_linearR(band=None,d2cMaps=None,oversampling = 1.):
    # Construct spectral (wavelength) grid
    lambdaMap = d2cMaps['lambdaMap']
    bandlims  = [lambdaMap[np.nonzero(lambdaMap)].min(),lambdaMap[np.nonzero(lambdaMap)].max()]
    #> initialize variables
    lambcens = []
    lambfwhms = []

    #> loop over wavelength bins (bin width defined based on MRS spectral resolution)
    lamb0   = bandlims[0]
    maxlamb = bandlims[1]

    # first iteration
    R       = getSpecR_linearR(lamb0=lamb0,band=band)
    fwhm    = lamb0 / R
    lambcen = lamb0 + (fwhm/2.)/oversampling

    lambcens.append(lambcen)
    lambfwhms.append(fwhm)

    # iterate over spectral range
    done = False
    while not done:
        R = getSpecR_linearR(lamb0=lamb0,band=band)
        fwhm = lamb0 / R
        lambcen = lamb0 + (fwhm/2.)/oversampling
        if (lambcen > maxlamb-(fwhm/2.)/oversampling):
            done = True
        else:
            lamb0 = lambcen + (fwhm/2.)/oversampling

        lambcens.append(lambcen)
        lambfwhms.append(fwhm)

    return np.array(lambcens),np.array(lambfwhms)

def point_source_centroiding(sci_img,band,d2cMaps,spec_grid=None,fit='2D'):
    import mrs_aux as maux
    # distortion maps
    sliceMap  = d2cMaps['sliceMap']
    lambdaMap = d2cMaps['lambdaMap']
    alphaMap  = d2cMaps['alphaMap']
    betaMap   = d2cMaps['betaMap']
    nslices   = d2cMaps['nslices']
    mrs_fwhm  = maux.MRS_FWHM[band[0]]
    lambcens,lambfwhms = spec_grid[0],spec_grid[1]
    unique_betas = np.sort(np.unique(betaMap[(sliceMap>100*int(band[0])) & (sliceMap<100*(int(band[0])+1))]))
    fov_lims  = [alphaMap[np.nonzero(lambdaMap)].min(),alphaMap[np.nonzero(lambdaMap)].max()]

    print 'STEP 1: Rough centroiding'
    # premise> center of point source is located in slice with largest signal
    # across-slice center:
    sum_signals = np.zeros(nslices)
    for islice in xrange(1+nslices):
        sum_signals[islice-1] = sci_img[(sliceMap == 100*int(band[0])+islice) & (~np.isnan(sci_img))].sum()
    source_center_slice = np.argmax(sum_signals)+1

    print 'Slice {} has the largest summed flux'.format(source_center_slice)

    # along-slice center:
    det_dims = (1024,1032)
    img = np.full(det_dims,0.)
    sel = (sliceMap == 100*int(band[0])+source_center_slice)
    img[sel]  = sci_img[sel]

    first_nonzero_row = 0
    while all(img[first_nonzero_row,:][~np.isnan(img[first_nonzero_row,:])] == 0.): first_nonzero_row+=1
    source_center_alpha = alphaMap[first_nonzero_row,img[first_nonzero_row,:].argmax()]

    # summary:
    print 'Source position: beta = {}arcsec, alpha = {}arcsec \n'.format(round(unique_betas[source_center_slice-1],2),round(source_center_alpha,2))

    print 'STEP 2: 1D Gaussian fit'

    # Fit Gaussian distribution to along-slice signal profile
    sign_amp,alpha_centers,alpha_fwhms = [np.full((len(lambcens)),np.nan) for j in range(3)]
    failed_fits = []
    for ibin in xrange(len(lambcens)):
        coords = np.where((sliceMap == 100*int(band[0])+source_center_slice) & (np.abs(lambdaMap-lambcens[ibin])<=lambfwhms[ibin]/2.) & (~np.isnan(sci_img)))
        if len(coords[0]) == 0: failed_fits.append(ibin); continue
        try:popt,pcov = curve_fit(gauss1d_wBaseline, alphaMap[coords], sci_img[coords], p0=[sci_img[coords].max(),source_center_alpha,mrs_fwhm/2.355,0],method='lm')
        except: failed_fits.append(ibin); continue
        sign_amp[ibin]      = popt[0]+popt[3]
        alpha_centers[ibin] = popt[1]
        alpha_fwhms[ibin]   = 2.355*np.abs(popt[2])

    # omit outliers
    for i in xrange(len(np.diff(sign_amp))):
        if np.abs(np.diff(alpha_centers)[i]) > 0.05:
            sign_amp[i],sign_amp[i+1],alpha_centers[i],alpha_centers[i+1],alpha_fwhms[i],alpha_fwhms[i+1] = [np.nan for j in range(6)]

    print '[Along-slice fit] The following bins failed to converge:'
    print failed_fits

    # Fit Gaussian distribution to across-slice signal profile (signal brute-summed in each slice)
    summed_signal,beta_centers,beta_fwhms = [np.full((len(lambcens)),np.nan) for j in range(3)]
    failed_fits = []
    for ibin in range(len(lambcens)):
        if np.isnan(alpha_centers[ibin]): failed_fits.append(ibin);continue
        sel = (np.abs(lambdaMap-lambcens[ibin])<=lambfwhms[ibin]/2.) & (~np.isnan(sci_img))
        try:signals = np.array([sci_img[(sliceMap == 100*int(band[0])+islice) & sel][np.abs(alphaMap[(sliceMap == 100*int(band[0])+islice) & sel]-alpha_centers[ibin]).argmin()] for islice in range(1,1+nslices)])
        except ValueError: failed_fits.append(ibin); continue
        try:popt,pcov = curve_fit(gauss1d_wBaseline, unique_betas, signals, p0=[signals.max(),unique_betas[source_center_slice-1],mrs_fwhm/2.355,0],method='lm')
        except: failed_fits.append(ibin); continue
        summed_signal[ibin] = popt[0]+popt[3]
        beta_centers[ibin]  = popt[1]
        beta_fwhms[ibin]    = 2.355*np.abs(popt[2])

    # # omit outliers
    # for i in range(len(np.diff(summed_signal))):
    #     if np.abs(np.diff(beta_centers)[i]) > 0.05:
    #         summed_signal[i],summed_signal[i+1],beta_centers[i],beta_centers[i+1],beta_fwhms[i],beta_fwhms[i+1] = [np.nan for j in range(6)]

    print '[Across-slice fit] The following bins failed to converge:'
    print failed_fits
    print ''

    if fit == '1D':
        sigma_alpha, sigma_beta = alpha_fwhms/2.355, beta_fwhms/2.355
        return sign_amp,alpha_centers,beta_centers,sigma_alpha, sigma_beta

    elif fit == '2D':
        print 'STEP 3: 2D Gaussian fit'
        sign_amp2D,alpha_centers2D,beta_centers2D,sigma_alpha2D,sigma_beta2D = [np.full((len(lambcens)),np.nan) for j in range(5)]
        failed_fits = []

        for ibin in range(len(lambcens)):
            # initial guess for fitting, informed by previous centroiding steps
            amp,alpha0,beta0  = sign_amp[ibin],alpha_centers[ibin],beta_centers[ibin]
            sigma_alpha, sigma_beta = alpha_fwhms[ibin]/2.355, beta_fwhms[ibin]/2.355
            base = 0
            guess = [amp, alpha0, beta0, sigma_alpha, sigma_beta, base]
            bounds = ([0,-np.inf,-np.inf,0,0,-np.inf],[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf])

            # data to fit
            coords = (np.abs(lambdaMap-lambcens[ibin])<lambfwhms[ibin]/2.)
            alphas, betas, zobs   = alphaMap[coords],betaMap[coords],sci_img[coords]
            alphabetas = np.array([alphas,betas])

            # perform fitting
            try:popt,pcov = curve_fit(gauss2d, alphabetas, zobs, p0=guess,bounds=bounds)
            except: failed_fits.append(ibin); continue

            sign_amp2D[ibin]      = popt[0]
            alpha_centers2D[ibin] = popt[1]
            beta_centers2D[ibin]  = popt[2]
            sigma_alpha2D[ibin]   = popt[3]
            sigma_beta2D[ibin]    = popt[4]

        print 'The following bins failed to converge:'
        print failed_fits

        return sign_amp2D,alpha_centers2D,beta_centers2D,sigma_alpha2D,sigma_beta2D

def point_source_along_slice_centroiding(sci_img,band,d2cMaps,spec_grid=None,offset_slice=0,campaign=None,verbose=False):
    # same as "point_source_centroiding" function, however only performs 1D Gaussian fitting, in along-slice (alpha) direction
    # param. "offset slice" allows to perform the centroiding analysis in a neighboring slice
    import mrs_aux as maux
    # distortion maps
    sliceMap  = d2cMaps['sliceMap']
    lambdaMap = d2cMaps['lambdaMap']
    alphaMap  = d2cMaps['alphaMap']
    betaMap   = d2cMaps['betaMap']
    nslices   = maux.MRS_nslices[band[0]]
    mrs_fwhm  = maux.MRS_FWHM[band[0]]
    lambcens,lambfwhms = spec_grid[0],spec_grid[1]
    unique_betas = np.sort(np.unique(betaMap[(sliceMap>100*int(band[0])) & (sliceMap<100*(int(band[0])+1))]))
    fov_lims  = [alphaMap[np.nonzero(lambdaMap)].min(),alphaMap[np.nonzero(lambdaMap)].max()]

    # premise> center of point source is located in slice with largest signal
    # across-slice center:
    sum_signals = np.zeros(nslices)
    for islice in xrange(1+nslices):
        sum_signals[islice-1] = sci_img[(sliceMap == 100*int(band[0])+islice) & (~np.isnan(sci_img))].sum()
    source_center_slice = np.argmax(sum_signals)+1
    source_center_slice+=offset_slice

    if verbose==True:
        print 'Slice {} has the largest summed flux'.format(source_center_slice)

    # along-slice center:
    det_dims = (1024,1032)
    img = np.full(det_dims,0.)
    sel = (sliceMap == 100*int(band[0])+source_center_slice)
    img[sel]  = sci_img[sel]

    first_nonzero_row = 0
    while all(img[first_nonzero_row,:][~np.isnan(img[first_nonzero_row,:])] == 0.): first_nonzero_row+=1
    source_center_alpha = alphaMap[first_nonzero_row,img[first_nonzero_row,:].argmax()]
    if campaign=='CV1RR':
        source_center_alpha = alphaMap[np.where(img[~np.isnan(img)].max() == img)]

    # Fit Gaussian distribution to along-slice signal profile
    sign_amps,alpha_centers,alpha_sigmas,bkg_amps = [np.full((len(lambcens)),np.nan) for j in range(4)]
    failed_fits = []
    for ibin in xrange(len(lambcens)):
        coords = np.where((sliceMap == 100*int(band[0])+source_center_slice) & (np.abs(lambdaMap-lambcens[ibin])<=lambfwhms[ibin]/2.) & (~np.isnan(sci_img)))
        if len(coords[0]) == 0: failed_fits.append(ibin); continue
        try:popt,pcov = curve_fit(gauss1d_wBaseline, alphaMap[coords], sci_img[coords], p0=[sci_img[coords].max(),source_center_alpha,mrs_fwhm/2.355,0],method='lm')
        except: failed_fits.append(ibin); continue
        sign_amps[ibin]     = popt[0]
        alpha_centers[ibin] = popt[1]
        alpha_sigmas[ibin]  = np.abs(popt[2])
        bkg_amps[ibin]      = popt[3]

    # omit outliers
    for i in xrange(len(np.diff(sign_amps))):
        if np.abs(np.diff(alpha_centers)[i]) > 0.05:
            sign_amps[i],sign_amps[i+1],alpha_centers[i],alpha_centers[i+1],alpha_sigmas[i],alpha_sigmas[i+1],bkg_amps[i],bkg_amps[i+1] = [np.nan for j in range(8)]
    if campaign == 'CV1RR':
        for i in xrange(len(np.diff(sign_amps))):
            if np.abs(np.diff(sign_amps)[i]) > 10.:
                sign_amps[i],sign_amps[i+1],alpha_centers[i],alpha_centers[i+1],alpha_sigmas[i],alpha_sigmas[i+1],bkg_amps[i],bkg_amps[i+1] = [np.nan for j in range(8)]

    return source_center_slice,sign_amps,alpha_centers,alpha_sigmas,bkg_amps

def get_pixel_spatial_area(band=None,d2cMaps=None):
    # Calculate size map
    # The spatial area of a pixel (assumed quadrilateral) is calculated as the sum of two triangles
    # The two trangles have side lengths A,B,C, and side C is shared (i.e. equal in both triangles)

    #get dimensions
    alphaULMap = d2cMaps['alphaULMap']
    alphaLLMap = d2cMaps['alphaLLMap']
    alphaURMap = d2cMaps['alphaURMap']
    alphaLRMap = d2cMaps['alphaULMap']

    betaULMap = d2cMaps['betaULMap']
    betaLLMap = d2cMaps['betaLLMap']
    betaURMap = d2cMaps['betaURMap']
    betaLRMap = d2cMaps['betaULMap']

    A1 = np.sqrt( (alphaULMap-alphaLLMap)**2 + (betaULMap-betaLLMap)**2 )
    B1 = np.sqrt( (alphaURMap-alphaULMap)**2 + (betaURMap-betaULMap)**2 )
    C1 = np.sqrt( (alphaURMap-alphaLLMap)**2 + (betaURMap-betaLLMap)**2 )
    A2 = np.sqrt( (alphaURMap-alphaLRMap)**2 + (betaURMap-betaLRMap)**2 )
    B2 = np.sqrt( (alphaLRMap-alphaLLMap)**2 + (betaLRMap-betaLLMap)**2 )
    C2 = C1.copy()

    # The area of a triangle can be calculated from the length of its sides using Heron's formula
    s1 = (A1+B1+C1)/2. # half of triangle's perimeter
    s2 = (A2+B2+C2)/2. # " " "

    Area1 = np.sqrt(s1*(s1-A1)*(s1-B1)*(s1-C1))
    Area2 = np.sqrt(s2*(s2-A2)*(s2-B2)*(s2-C2))

    spaxelsizeMap = Area1 + Area2

    return spaxelsizeMap

def get_pixel_area_in_alphalambda(band=None,d2cMaps=None):
    # Calculate size map
    # The spatial area of a pixel (assumed quadrilateral) is calculated as the sum of two triangles
    # The two trangles have side lengths A,B,C, and side C is shared (i.e. equal in both triangles)

    #get dimensions
    alphaULMap = d2cMaps['alphaULMap']
    alphaLLMap = d2cMaps['alphaLLMap']
    alphaURMap = d2cMaps['alphaURMap']
    alphaLRMap = d2cMaps['alphaULMap']

    lambdaULMap = d2cMaps['lambdaULMap']
    lambdaLLMap = d2cMaps['lambdaLLMap']
    lambdaURMap = d2cMaps['lambdaURMap']
    lambdaLRMap = d2cMaps['lambdaULMap']

    A1 = np.sqrt( (alphaULMap-alphaLLMap)**2 + (lambdaULMap-lambdaLLMap)**2 )
    B1 = np.sqrt( (alphaURMap-alphaULMap)**2 + (lambdaURMap-lambdaULMap)**2 )
    C1 = np.sqrt( (alphaURMap-alphaLLMap)**2 + (lambdaURMap-lambdaLLMap)**2 )
    A2 = np.sqrt( (alphaURMap-alphaLRMap)**2 + (lambdaURMap-lambdaLRMap)**2 )
    B2 = np.sqrt( (alphaLRMap-alphaLLMap)**2 + (lambdaLRMap-lambdaLLMap)**2 )
    C2 = C1.copy()

    # The area of a triangle can be calculated from the length of its sides using Heron's formula
    s1 = (A1+B1+C1)/2. # half of triangle's perimeter
    s2 = (A2+B2+C2)/2. # " " "

    Area1 = np.sqrt(s1*(s1-A1)*(s1-B1)*(s1-C1))
    Area2 = np.sqrt(s2*(s2-A2)*(s2-B2)*(s2-C2))

    pixsiz_alphalambda = Area1 + Area2

    return pixsiz_alphalambda


def aperture_photometry(band=None,sci_img=None,apertureMask=None,aperture_area=None,d2cMaps=None,spec_grid=None,img_type='sci'):
    lambdaMap = d2cMaps['lambdaMap']
    lambcens,lambfwhms = spec_grid[0],spec_grid[1]

    pixsiz_img = get_pixel_spatial_area(band=band,d2cMaps=d2cMaps) # [arcsec*arcsec]
    pixelsiz_alphalambdaMap = get_pixel_area_in_alphalambda(band=band,d2cMaps=d2cMaps) # [arcsec*micron]
    pixvol_img = pixelsiz_alphalambdaMap*d2cMaps['bdel']

    signals_aper = np.zeros(len(lambcens))
    for ibin in range(len(lambcens)):
        # map containing only pixels within one spectral bin, within the defined aperture, omitting NaNs
        pixelsInBinInApertureNoNaN = np.where((np.abs(lambdaMap-lambcens[ibin])<lambfwhms[ibin]/2.) & (apertureMask!=0.) & (np.isnan(sci_img)==False) )
        nPixels = len(pixelsInBinInApertureNoNaN[0])

        # map containing only pixels within one spectral bin, within the defined aperture
        sci_img_masked    = sci_img[pixelsInBinInApertureNoNaN]
        pixsiz_img_masked = pixsiz_img[pixelsInBinInApertureNoNaN]
        pixvol_img_masked = pixvol_img[pixelsInBinInApertureNoNaN]

        # perform aperture photometry
        if img_type=='sci':
            signals_aper[ibin] = sci_img_masked.sum()/nPixels
            #signals_aper[ibin] = ((sci_img_masked*pixsiz_img_masked).sum()/pixvol_img_masked.sum()) * aperture_area * lambfwhms[ibin]
        elif img_type=='err':
            var_img = (sci_img_masked*pixsiz_img_masked)**2.
            signals_aper[ibin] = var_img.sum()/nPixels
            # signals_aper[ibin] = (var_img.sum()/pixvol_img_masked.sum()) * aperture_area * lambfwhms[ibin]
        elif img_type=='psf':
            psf = sci_img

            # enforce normalization of psf in every wavelength bin
            psf[np.abs(lambdaMap-lambcens[ibin])<lambfwhms[ibin]/2.] = psf[np.abs(lambdaMap-lambcens[ibin])<lambfwhms[ibin]/2.]/psf[np.abs(lambdaMap-lambcens[ibin])<lambfwhms[ibin]/2.].sum()

            psf_masked = psf[pixelsInBinInApertureNoNaN]
            signals_aper[ibin] = (psf_masked.sum()/pixvol_img_masked.sum()) * aperture_area * lambfwhms[ibin]

    return signals_aper

def optimal_extraction(band=None,sci_img=None,err_img=None,psf=None,spec_grid=None,d2cMaps=None):
    lambdaMap = d2cMaps['lambdaMap']
    lambcens,lambfwhms = spec_grid[0],spec_grid[1]
    pixsiz_img = get_pixel_spatial_area(band=band,d2cMaps=d2cMaps)
    var_img    = (err_img*pixsiz_img)**2.

    signals_opex,signals_error_opex = np.zeros(len(lambcens)),np.zeros(len(lambcens))
    for ibin in range(len(lambcens)):
        # enforce normalization of psf in every wavelength bin
        psf[np.abs(lambdaMap-lambcens[ibin])<lambfwhms[ibin]/2.] = psf[np.abs(lambdaMap-lambcens[ibin])<lambfwhms[ibin]/2.]/psf[np.abs(lambdaMap-lambcens[ibin])<lambfwhms[ibin]/2.].sum()

        # map containing only pixels within one spectral bin
        pixelsInBinNoNaN = np.where((np.abs(lambdaMap-lambcens[ibin])<lambfwhms[ibin]/2.) & (np.isnan(sci_img)==False) & (np.isnan(var_img)==False) )
        weights = (psf[pixelsInBinNoNaN] / var_img[pixelsInBinNoNaN]) / (psf[pixelsInBinNoNaN]**2./var_img[pixelsInBinNoNaN]).sum()

        signals_opex[ibin] = (weights*sci_img[pixelsInBinNoNaN]*pixsiz_img[pixelsInBinNoNaN]).sum()
        signals_error_opex[ibin] = np.sqrt( ( weights**2. * var_img[pixelsInBinNoNaN] ).sum() )

    return signals_opex,signals_error_opex

#-model
def elliptical_aperture(center=[0,0],r=1.,q=1.,pa=0,d2cMaps=None):
    """
    Elliptical aperture, written by Ruyman Azzollini (DIAS, ruyman.azzollini@gmail.com), edited by Ioannis Argyriou (KUL, ioannis.argyriou@kuleuven.be)

    centre: (alpha,beta)
    r : semi-major axis
    q : axis ratio. q=1 for circular apertures.
    pa : position angle, counter-clockwise from 'y' axis.

    """
    assert q <= 1.

    if q==1:
        spatial_extent = (d2cMaps['alphaMap']-center[0])**2. + (d2cMaps['betaMap']-center[1])**2.
    else:
        radeg = 180. / np.pi
        ang = pa / radeg

        cosang = np.cos(ang)
        sinang = np.sin(ang)

        xtemp = (d2cMaps['alphaMap']-centre[0]) * cosang + (d2cMaps['betaMap']-centre[1]) * sinang
        ytemp = -(d2cMaps['alphaMap']-centre[0]) * sinang + (d2cMaps['betaMap']-centre[1]) * cosang

        spatial_extent = (xtemp/q)**2. + ytemp**2.

    elliptical_aperture_area = np.pi * r**2. * q

    return spatial_extent <= r**2.,elliptical_aperture_area

def rectangular_aperture(center=[0,0],width=1.,height=1.,d2cMaps=None):
    """
    Rectangular aperture, written by Ruyman Azzollini (DIAS, ruyman.azzollini@gmail.com), edited by Ioannis Argyriou (KUL, ioannis.argyriou@kuleuven.be)

    centre: (alpha,beta)
    width : width (alpha dimension)
    height : height (beta dimension)
    """

    pixels_inside_rectangle = (np.abs(d2cMaps['alphaMap']-center[0])<=(width/2.) ) & (np.abs(d2cMaps['betaMap']-center[1])<=(height/2.) )

    rectangular_aperture_area = width * height

    return pixels_inside_rectangle,rectangular_aperture_area

def evaluate_psf_cdp(psffits=None,source_center=[0,0],d2cMaps=None):
    # PSF CDP is provided as a spectral cube
    #>get values
    psf_values = psffits[1].data.transpose(2,1,0) # flip data from Z,Y,X to X,Y,Z

    #>get grid
    NAXIS1,NAXIS2,NAXIS3 = psf_values.shape

    alphastpix = psffits[1].header['CRPIX1'] # pixel nr
    alpha_step = psffits[1].header['CDELT1'] # arcsec/pix
    stalpha    = psffits[1].header['CRVAL1']-(alphastpix-1)*alpha_step # arcsec

    betastpix = psffits[1].header['CRPIX2'] # pixel nr
    beta_step = psffits[1].header['CDELT2'] # arcsec/pix
    stbeta    = psffits[1].header['CRVAL2']-(betastpix-1)*beta_step # arcsec

    stwavl = psffits[1].header['CRVAL3'] # microns
    wavl_step   = psffits[1].header['CDELT3'] # microns/pix

    alpha_slices = np.linspace(stalpha,stalpha+ (NAXIS1-1.5)*alpha_step,NAXIS1)
    beta_slices  = np.linspace(stbeta,stbeta+ (NAXIS2-1.5)*beta_step,NAXIS2)
    wvl_slices   = np.linspace(stwavl ,stwavl+NAXIS3*wavl_step,NAXIS3)

    #> center psf to source
    alpha_slices += source_center[0]
    beta_slices  += source_center[1]

    #> create interpolant based on regular grid
    interpolpsf = scp_interpolate.RegularGridInterpolator((alpha_slices,beta_slices,wvl_slices),psf_values)
    interpolpsf.fill_value=0.
    interpolpsf.bounds_error=False

    # evaluate psf at each pixel center and pixel corner
    alphaULMap = d2cMaps['alphaULMap']
    alphaURMap = d2cMaps['alphaURMap']
    alphaLLMap = d2cMaps['alphaLLMap']
    alphaLRMap = d2cMaps['alphaLRMap']
    alphaMap   = d2cMaps['alphaMap']

    betaULMap = d2cMaps['betaULMap']
    betaURMap = d2cMaps['betaURMap']
    betaLLMap = d2cMaps['betaLLMap']
    betaLRMap = d2cMaps['betaLRMap']
    betaMap   = d2cMaps['betaMap']

    lambdaULMap = d2cMaps['lambdaULMap']
    lambdaURMap = d2cMaps['lambdaURMap']
    lambdaLLMap = d2cMaps['lambdaLLMap']
    lambdaLRMap = d2cMaps['lambdaLRMap']
    lambdaMap = d2cMaps['lambdaMap']

    #> interpolate psf to science image pixel centers and corners
    #-- assume no significant change in wavelength over one pixel size
    psfUL  = interpolpsf((alphaULMap,betaULMap,lambdaULMap))
    psfUR  = interpolpsf((alphaURMap,betaURMap,lambdaURMap))
    psfLL  = interpolpsf((alphaLLMap,betaLLMap,lambdaLLMap))
    psfLR  = interpolpsf((alphaLRMap,betaLRMap,lambdaLRMap))
    psfCEN = interpolpsf((alphaMap,betaMap,lambdaMap))

    #> evaluate psf as a weighted average
    w = np.array([0.125,0.125,0.125,0.125,0.5]) # WARNING: ARBITRARY!
    sumweights = w.sum()

    psf = (w[0]*psfUL+w[1]*psfUR+w[2]*psfLL+w[3]*psfLR+w[4]*psfCEN)/sumweights

    return psf

#-fit
#--1d
def straight_line(x,a,b):
    return a*x+b

def order2polyfit(x,a,b,c):
    return a*x**2 + b*x +c

def order3polyfit(x,a,b,c,d):
    return a*x**3 + b*x**2 + c*x +d

def order4polyfit(x,a,b,c,d,e):
    return a*x**4 + b*x**3 + c*x**2 + d*x + e

def gauss1d_wBaseline(x, A, mu, sigma, baseline):
    """1D Gaussian distribution function"""
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))+baseline

def gauss1d_woBaseline(x, A, mu, sigma):
    """1D Gaussian distribution function"""
    return A*np.exp(-(x-mu)**2/(2.*sigma**2))

def skewnorm_func(x, A, mu, sigmag, alpha):
    #normal distribution
    normpdf = (1/(sigmag*np.sqrt(2*np.pi)))*np.exp(-(np.power((x-mu),2)/(2*np.power(sigmag,2))))
    normcdf = (0.5*(1+sp.erf((alpha*((x-mu)/sigmag))/(np.sqrt(2)))))
    return 2*A*normpdf*normcdf

def lorentzian_profile(x,A,mu,sigma):
    # Wrong amplitude?
    L_nu = (2*A/(np.pi*sigma))/(1+4*((x-mu)/sigma)**2)
    return L_nu

def voight_profile(x,A,mu,sigma,f):
    G_nu = (A/sigma) * np.sqrt(4*np.log(2)/np.pi) * np.exp(-4*np.log(2)*((x-mu)/sigma)**2)
    L_nu = (2*A/(np.pi*sigma))/(1+4*((x-mu)/sigma)**2)
    return f*L_nu + (1-f)*G_nu

def skewed_voight(x,A,mu,sigma0,f,a):
    """ According to Stancik and Brauns (2008)"""
    A /= ( ((1-f)/sigma0) * np.sqrt(4*np.log(2)/np.pi) + f*(2./(np.pi*sigma0)))
    sigma = 2*sigma0 / (1+np.exp(a*(x-mu)))
    G_nu = (A/sigma) * np.sqrt(4*np.log(2)/np.pi) * np.exp(-4*np.log(2)*((x-mu)/sigma)**2)
    L_nu = (2*A/(np.pi*sigma))/(1+4*((x-mu)/sigma)**2)
    return f*L_nu + (1-f)*G_nu

def FPfunc(wavenumber,R,D,phi,theta=0):
    # T = 1-R-A
    return (1 + (4*R/(1-R)**2) * np.sin(2*np.pi*D*wavenumber*np.cos(theta) - (phi-np.pi))**2 )**-1

def FPfunc_noPhaseShift(wavenumber,R,D,theta=0):
    # T = 1-R-A
    return (1 + (4*R/(1-R)**2) * np.sin(2*np.pi*D*wavenumber*np.cos(theta))**2 )**-1

def reflectivity_from_continuum(y):
    """
    Spectral continuum where: sin(delta/2) = 0 ==> 1 / (1 + F) = y
    ==> (1 + F) = 1/y ==> F = 1/y -1

    F = 4*R/(1-R)**2 ==> F * (1 -2*R + R**2) - 4*R = 0
    ==> F -2*F*R + F*R**2  - 4*R = 0 ==> R**2 -(2 + 4/F)*R + 1 = 0 (second order equation)

    a = 1; b = -(2 + 4/F) ; c = 1
    Discr = b**2 - 4*a*c = (2 + 4/F)**2 - 4 = 4 + 16/F + 16/F**2 - 4 = 16/F + 16/F**2
    R1 = -b +sqrt(Discr) / 2*a = ((2 + 4/F) + sqrt(16/F + 16/F**2)) / 2
    R2 = -b -sqrt(Discr) / 2*a = ((2 + 4/F) - sqrt(16/F + 16/F**2)) / 2
    """
    F = 1/y -1
    R1 = ((2. + 4./F) + np.sqrt(16./F + 16./F**2)) / 2.
    R2 = ((2. + 4./F) - np.sqrt(16./F + 16./F**2)) / 2.
    return F,R1,R2

#--2d
def gauss2d(xy, amp, x0, y0, sigma_x, sigma_y, base):
    x, y = xy
    a = 1/(2*sigma_x**2)
    b = 1/(2*sigma_y**2)
    inner = a * (x - x0)**2
    inner += b * (y - y0)**2
    return amp * np.exp(-inner) + base

def voight_profile2d(xy,amp,x0,y0,sigma_x,sigma_y,f):
    x, y = xy
    a = 1/(2*sigma_x**2)
    b = 1/(2*sigma_y**2)
    inner = a * (x - x0)**2
    inner += b * (y - y0)**2
    G_nu = amp * np.exp(-inner)
    L_nu = (2*amp/(np.pi*(sigma_x+sigma_y)))/(1+4*inner)
    return f*L_nu + (1-f)*G_nu

#--etalon lines
def fit_etalon_lines(x,y,peaks,fit_func='skewed_voight',sigma0=3.5,f0=0.5,a0=0.1):
    # Available fitting functions: 'gauss1d','skewnorm_func','voight_profile','skewed_voight'
    xdata = x.copy()
    xdata[np.isnan(xdata)] = 0.
    ydata = y.copy()
    ydata[np.isnan(ydata)] = 0.

    bounds_gauss = ([0,0,0],[np.inf,np.inf,np.inf])
    bounds_skewnorm = ([0,0,0,-np.inf],[np.inf,np.inf,np.inf,np.inf])
    bounds_lorentzian = ([0,0,0],[np.inf,np.inf,np.inf])
    bounds_voight = ([0,0,0,0],[np.inf,np.inf,np.inf,1])
    bounds_skewvoight = ([0,0,0,0,-np.inf],[np.inf,np.inf,np.inf,1,np.inf])

    fitparams = []
    fiterrors = []
    fitting_flag = []
    range_ini = np.full(len(peaks),np.nan) # fit range first position
    range_fin = np.full(len(peaks),np.nan) # fit range last position
    for i in range(len(peaks)):
        if i == len(peaks)-1: N = np.diff(peaks)[i-1]/2 - 1
        else: N = np.diff(peaks)[i]/2

        peak_idx = peaks[i]

        guess_gauss = [ydata[peak_idx],xdata[peak_idx],sigma0]
        guess_skewnorm = guess_gauss+[a0]
        guess_lorentzian = [ydata[peak_idx],xdata[peak_idx],sigma0]
        guess_voight = guess_gauss+[f0]
        guess_skewvoight = guess_gauss+[f0,a0]

        if peak_idx<N:
            range_ini[i] = xdata[0]
            range_fin[i] = xdata[peak_idx+N]
            if fit_func == 'gauss1d':
                popt,pcov = curve_fit(gauss1d_woBaseline,xdata[0:peak_idx+N],ydata[0:peak_idx+N],p0=guess_gauss,absolute_sigma=True,bounds=bounds_gauss)
                fitting_flag.append('gauss1d')
            elif fit_func == 'skewnorm_func':
                try:
                    popt,pcov = curve_fit(skewnorm_func,xdata[0:peak_idx+N],ydata[0:peak_idx+N],p0=guess_skewnorm,absolute_sigma=True,bounds=bounds_skewnorm)
                    fitting_flag.append('skewnorm_func')
                except RuntimeError:
                    popt,pcov = curve_fit(gauss1d_woBaseline,xdata[0:peak_idx+N],ydata[0:peak_idx+N],p0=guess_gauss,absolute_sigma=True,bounds=bounds_gauss)
                    fitting_flags.append('gauss1d')
            elif fit_func == 'lorentzian_profile':
                popt,pcov = curve_fit(lorentzian_profile,xdata[0:peak_idx+N],ydata[0:peak_idx+N],p0=guess_lorentzian,absolute_sigma=True,bounds=bounds_lorentzian)
                fitting_flag.append('lorentzian_profile')
            elif fit_func == 'voight_profile':
                popt,pcov = curve_fit(voight_profile,xdata[0:peak_idx+N],ydata[0:peak_idx+N],p0=guess_voight,absolute_sigma=True,bounds=bounds_voight)
                fitting_flag.append('voight_profile')
            elif fit_func == 'skewed_voight':
                try:
                    popt,pcov = curve_fit(skewed_voight,xdata[0:peak_idx+N],ydata[0:peak_idx+N],p0=guess_skewvoight,absolute_sigma=True,bounds=bounds_skewvoight)
                    fitting_flag.append('skewed_voight')
                except RuntimeError:
                    try:
                        popt,pcov = curve_fit(voight_profile,xdata[0:peak_idx+N],ydata[0:peak_idx+N],p0=guess_voight,absolute_sigma=True,bounds=bounds_voight)
                        fitting_flag.append('voight_profile')
                    except RuntimeError:
                        popt,pcov = curve_fit(gauss1d_woBaseline,xdata[0:peak_idx+N],ydata[0:peak_idx+N],p0=guess_gauss,absolute_sigma=True,bounds=bounds_gauss)
                        fitting_flag.append('gauss1d')
        elif len(ydata)-peak_idx<N:
            range_ini[i] = xdata[peak_idx-N]
            range_fin[i] = xdata[-1]
            if fit_func == 'gauss1d':
                popt,pcov = curve_fit(gauss1d_woBaseline,xdata[peak_idx-N:],ydata[peak_idx-N:],p0=guess_gauss,absolute_sigma=True,bounds=bounds_gauss)
                fitting_flag.append('gauss1d')
            elif fit_func == 'skewnorm_func':
                try:
                    popt,pcov = curve_fit(skewnorm_func,xdata[peak_idx-N:],ydata[peak_idx-N:],p0=guess_skewnorm,absolute_sigma=True,bounds=bounds_skewnorm)
                    fitting_flag.append('skewnorm_func')
                except RuntimeError:
                    popt,pcov = curve_fit(gauss1d_woBaseline,xdata[peak_idx-N:],ydata[peak_idx-N:],p0=guess_gauss,absolute_sigma=True,bounds=bounds_gauss)
                    fitting_flag.append('gauss1d')
            elif fit_func == 'lorentzian_profile':
                popt,pcov = curve_fit(lorentzian_profile,xdata[peak_idx-N:],ydata[peak_idx-N:],p0=guess_lorentzian,absolute_sigma=True,bounds=bounds_lorentzian)
                fitting_flag.append('lorentzian_profile')
            elif fit_func == 'voight_profile':
                popt,pcov = curve_fit(voight_profile,xdata[peak_idx-N:],ydata[peak_idx-N:],p0=guess_voight,absolute_sigma=True,bounds=bounds_voight)
                fitting_flag.append('voight_profile')
            elif fit_func == 'skewed_voight':
                try:
                    popt,pcov = curve_fit(skewed_voight,xdata[peak_idx-N:],ydata[peak_idx-N:],p0=guess_skewvoight,absolute_sigma=True,bounds=bounds_skewvoight)
                    fitting_flag.append('skewed_voight')
                except RuntimeError:
                    try:
                        popt,pcov = curve_fit(voight_profile,xdata[peak_idx-N:],ydata[peak_idx-N:],p0=guess_voight,absolute_sigma=True,bounds=bounds_voight)
                        fitting_flag.append('voight_profile')
                    except RuntimeError:
                        popt,pcov = curve_fit(gauss1d_woBaseline,xdata[0:peak_idx+N],ydata[0:peak_idx+N],p0=guess_gauss,absolute_sigma=True,bounds=bounds_gauss)
                        fitting_flag.append('gauss1d')
        else:
            range_ini[i] = xdata[peak_idx-N]
            range_fin[i] = xdata[peak_idx+N]
            if fit_func == 'gauss1d':
                popt,pcov = curve_fit(gauss1d_woBaseline,xdata[peak_idx-N:peak_idx+N],ydata[peak_idx-N:peak_idx+N],p0=guess_gauss,absolute_sigma=True,bounds=bounds_gauss)
                fitting_flag.append('gauss1d')
            elif fit_func == 'skewnorm_func':
                try:
                    popt,pcov = curve_fit(skewnorm_func,xdata[peak_idx-N:peak_idx+N],ydata[peak_idx-N:peak_idx+N],p0=guess_skewnorm,absolute_sigma=True,bounds=bounds_skewnorm)
                    fitting_flag.append('skewnorm_func')
                except RuntimeError:
                    popt,pcov = curve_fit(gauss1d_woBaseline,xdata[peak_idx-N:peak_idx+N],ydata[peak_idx-N:peak_idx+N],p0=guess_gauss,absolute_sigma=True,bounds=bounds_gauss)
                    fitting_flag.append('gauss1d')
            elif fit_func == 'lorentzian_profile':
                popt,pcov = curve_fit(lorentzian_profile,xdata[peak_idx-N:peak_idx+N],ydata[peak_idx-N:peak_idx+N],p0=guess_lorentzian,absolute_sigma=True,bounds=bounds_lorentzian)
                fitting_flag.append('gauss1d')
            elif fit_func == 'voight_profile':
                    popt,pcov = curve_fit(voight_profile,xdata[peak_idx-N:peak_idx+N],ydata[peak_idx-N:peak_idx+N],p0=guess_voight,absolute_sigma=True,bounds=bounds_voight)
                    fitting_flag.append('voight_profile')
            elif fit_func == 'skewed_voight':
                try:
                    popt,pcov = curve_fit(skewed_voight,xdata[peak_idx-N:peak_idx+N],ydata[peak_idx-N:peak_idx+N],p0=guess_skewvoight,absolute_sigma=True,bounds=bounds_skewvoight)
                    fitting_flag.append('skewed_voight')
                except RuntimeError:
                    popt,pcov = curve_fit(voight_profile,xdata[peak_idx-N:peak_idx+N],ydata[peak_idx-N:peak_idx+N],p0=guess_voight,absolute_sigma=True,bounds=bounds_voight)
                    fitting_flag.append('voight_profile')
        fitparams.append(popt)
        fiterrors.append(pcov)

    return fitparams,fiterrors,fitting_flag,range_ini,range_fin

def get_amplitude(fitparams,fitting_flag):
    amplitude = np.full(len(fitparams),np.nan)
    for i in range(len(fitparams)):
        if fitting_flag[i] == 'gauss1d':
            amplitude[i] = fitparams[i][0]
        elif fitting_flag[i] == 'skewnorm_func':
            plotx = np.linspace(fitparams[i][1]-1*fitparams[i][2],fitparams[i][1]+1*fitparams[i][2],500)
            ploty = skewnorm_func(plotx,*fitparams[i])
            amplitude[i] = np.max(ploty)
        elif fitting_flag[i] == 'lorentzian_profile':
            amplitude[i] = fitparams[i][0]
        elif fitting_flag[i] == 'voight_profile':
            amplitude[i] = fitparams[i][0]
        elif fitting_flag[i] == 'skewed_voight':
            plotx = np.linspace(fitparams[i][1]-1*fitparams[i][2],fitparams[i][1]+1*fitparams[i][2],500)
            ploty = skewed_voight(plotx,*fitparams[i])
            amplitude[i] = np.max(ploty)
    return amplitude

def get_linecenter(fitparams,fitting_flag):
    linecenter = np.full(len(fitparams),np.nan)
    for i in range(len(fitparams)):
        if fitting_flag[i] == 'gauss1d':
            linecenter[i] = fitparams[i][1]
        elif fitting_flag[i] == 'skewnorm_func':
            plotx = np.linspace(fitparams[i][1]-1*fitparams[i][2],fitparams[i][1]+1*fitparams[i][2],500)
            ploty = skewnorm_func(plotx,*fitparams[i])
            linecenter[i] = plotx[np.argmax(ploty)]
        elif fitting_flag[i] == 'lorentzian_profile':
            linecenter[i] = fitparams[i][1]
        elif fitting_flag[i] == 'voight_profile':
            linecenter[i] = fitparams[i][1]
        elif fitting_flag[i] == 'skewed_voight':
            plotx = np.linspace(fitparams[i][1]-1*fitparams[i][2],fitparams[i][1]+1*fitparams[i][2],500)
            ploty = skewed_voight(plotx,*fitparams[i])
            linecenter[i] = plotx[np.argmax(ploty)]
    return linecenter

def get_FWHM(fitparams,fitting_flag):
    fwhm = np.full(len(fitparams),np.nan)
    for i in range(len(fitparams)):
        if fitting_flag[i] == 'gauss1d':
            fwhm[i] = 2.355*fitparams[i][2]
        elif fitting_flag[i] == 'skewnorm_func':
            plotx = np.linspace(fitparams[i][1]-3*fitparams[i][2],fitparams[i][1]+3*fitparams[i][2],500)
            ploty = skewnorm_func(plotx,*fitparams[i])
            fwhm_idxs = np.abs(ploty-ploty.max()/2.).argsort()[:2]
            fwhm[i] = np.abs(plotx[fwhm_idxs[1]]-plotx[fwhm_idxs[0]])
        elif fitting_flag[i] == 'lorentzian_profile':
            plotx = np.linspace(fitparams[i][1]-3*fitparams[i][2],fitparams[i][1]+3*fitparams[i][2],500)
            ploty = lorentzian_profile(plotx,*fitparams[i])
            fwhm_idxs = np.abs(ploty-ploty.max()/2.).argsort()[:2]
            fwhm[i] = np.abs(plotx[fwhm_idxs[1]]-plotx[fwhm_idxs[0]])
        elif fitting_flag[i] == 'voight_profile':
            plotx = np.linspace(fitparams[i][1]-3*fitparams[i][2],fitparams[i][1]+3*fitparams[i][2],500)
            ploty = voight_profile(plotx,*fitparams[i])
            fwhm_idxs = np.abs(ploty-ploty.max()/2.).argsort()[:2]
            fwhm[i] = np.abs(plotx[fwhm_idxs[1]]-plotx[fwhm_idxs[0]])
        elif fitting_flag[i] == 'skewed_voight':
            plotx = np.linspace(fitparams[i][1]-3*fitparams[i][2],fitparams[i][1]+3*fitparams[i][2],500)
            ploty = skewed_voight(plotx,*fitparams[i])
            fwhm_idxs = np.abs(ploty-ploty.max()/2.).argsort()[:2]
            fwhm[i] = np.abs(plotx[fwhm_idxs[1]]-plotx[fwhm_idxs[0]])
    return fwhm

def get_skewness(fitparams,fitting_flag):
    skewparam = np.full(len(fitparams),np.nan)
    for i in range(len(fitparams)):
        if fitting_flag[i] == 'gauss1d':
            skewparam[i] = 0.
        elif fitting_flag[i] == 'skewnorm_func':
            skewparam[i] = fitparams[i][3]
        elif fitting_flag[i] == 'lorentzian_profile':
            skewparam[i] = 0.
        elif fitting_flag[i] == 'voight_profile':
            skewparam[i] = 0.
        elif fitting_flag[i] == 'skewed_voight':
            skewparam[i] = fitparams[i][4]
    return skewparam

def sum_etalon_lines(xdata,fitparams,fitting_flag):
    summed_signal = np.zeros(len(xdata))
    for i in range(len(fitparams)):
        if fitting_flag[i] == 'gauss1d':
            signal = gauss1d_woBaseline(xdata,*fitparams[i])
        elif fitting_flag[i] == 'skewnorm_func':
            signal = skewnorm_func(xdata,*fitparams[i])
        elif fitting_flag[i] == 'lorentzian_profile':
            signal = lorentzian_profile(xdata,*fitparams[i])
        elif fitting_flag[i] == 'voight_profile':
            signal = voight_profile(xdata,*fitparams[i])
        elif fitting_flag[i] == 'skewed_voight':
            signal = skewed_voight(xdata,*fitparams[i])
        summed_signal+= signal
    return summed_signal


def plot_etalon_fit(fitparams,fitting_flag):
    for i in range(len(fitparams)):
        if fitting_flag[i] == 'gauss1d':
            plotx = np.linspace(fitparams[i][1]-10*fitparams[i][2],fitparams[i][1]+10*fitparams[i][2],100)
            ploty = gauss1d_woBaseline(plotx,*fitparams[i])
        elif fitting_flag[i] == 'skewnorm_func':
            plotx = np.linspace(fitparams[i][1]-10*fitparams[i][2],fitparams[i][1]+10*fitparams[i][2],100)
            ploty = skewnorm_func(plotx,*fitparams[i])
        elif fitting_flag[i] == 'lorentzian_profile':
            plotx = np.linspace(fitparams[i][1]-10*fitparams[i][2],fitparams[i][1]+10*fitparams[i][2],500)
            ploty = lorentzian_profile(plotx,*fitparams[i])
        elif fitting_flag[i] == 'voight_profile':
            plotx = np.linspace(fitparams[i][1]-10*fitparams[i][2],fitparams[i][1]+10*fitparams[i][2],500)
            ploty = voight_profile(plotx,*fitparams[i])
        elif fitting_flag[i] == 'skewed_voight':
            plotx = np.linspace(fitparams[i][1]-10*fitparams[i][2],fitparams[i][1]+10*fitparams[i][2],500)
            ploty = skewed_voight(plotx,*fitparams[i])
        plt.plot(plotx,ploty,'r')

# normalize fringes
def norm_fringe(sci_data,thres=0,min_dist=2,k=3,ext=3):
    # determine peaks
    sci_data_noNaN = sci_data.copy()
    sci_data_noNaN[np.isnan(sci_data_noNaN)] = 0.
    peaks = find_peaks(sci_data_noNaN,thres=thres,min_dist=min_dist)
    if peaks[0] == np.nonzero(sci_data_noNaN)[0][0]:
        peaks = np.delete(peaks,0)
    # determine fringe continuum
    arr_interpolator = scp_interpolate.InterpolatedUnivariateSpline(peaks,sci_data_noNaN[peaks],k=k,ext=ext)
    sci_data_profile = arr_interpolator(range(len(sci_data_noNaN)))

    return sci_data_noNaN,peaks,sci_data_profile

# def fringe_scanning(wvnr,sci_data,xpos,thres=0):
#     # fringe peaks are found according to two parameters
#     # 1) threshold
#     # 2) minimum distance between subsequent peaks
#     detsubstrfringe_wvnr_spacing = (1/0.1710)/2. # assuming the AR coating to pixel aluminization fringe is always the dominant source of fringing
#     min_dist = 0
#     while wvnr[len(wvnr)/2 - min_dist]-wvnr[len(wvnr)/2]<detsubstrfringe_wvnr_spacing :
#         min_dist+=1
#     min_dist = int(round(min_dist/1.6))
#     print 'Min_dist between subsequent peaks: {} elements (pixels)'.format(min_dist)
#
#     # jumps/drops in signal due to displacement in pixel column
#     discont_idxs = np.unique(np.sort(len(sci_data)-1-1-np.where((np.abs(np.diff(xpos))>=1))[0]))+1
#     discont_idxs = np.concatenate([discont_idxs,np.array([len(sci_data)-1])])
#
#     # omit sections with too few elements
#     bad_sections = []
#     for i in range(len(discont_idxs)-1):
#         if discont_idxs[i+1]-discont_idxs[i] <= min_dist:
#             if discont_idxs[i+1] == -1:
#                 sci_data[discont_idxs[i]:] = np.nan
#                 bad_sections.extend([i,i+1])
#             elif discont_idxs[i] == 0:
#                 sci_data[:discont_idxs[i+1]] = np.nan
#                 bad_sections.extend([i,i+1])
#             else:
#                 sci_data[discont_idxs[i]-1:discont_idxs[i+1]+1] = np.nan
#                 bad_sections.extend([i,i+1])
#     discont_idxs = np.delete(discont_idxs,bad_sections)
#
#     omit_nan = ~np.isnan(sci_data).copy()
#     wvnr_noNaN = wvnr[omit_nan]
#     sci_data_noNaN = sci_data[omit_nan]
#
#     # second iteration
#     discont_idxs = np.unique(np.sort(len(sci_data_noNaN)-1-1-np.where((np.abs(np.diff(xpos[np.flipud(omit_nan)]))>=1))[0]))+1
#     if discont_idxs[0] != 0:
#         discont_idxs = np.concatenate((np.array([0]),discont_idxs))
#     if discont_idxs[-1] != len(sci_data_noNaN)-1:
#         discont_idxs = np.concatenate((discont_idxs,np.array([len(sci_data_noNaN)-1])))
#
#     outlier_condition = 5*np.std(np.abs(np.diff(sci_data_noNaN))) # five sigma
#     outliers = np.where(np.abs(np.diff(sci_data_noNaN)) > outlier_condition)[0]
#     invalid_outlier_idxs = [i for i in range(len(outliers)) if outliers[i] in discont_idxs-1]
#     outliers = np.delete(outliers,invalid_outlier_idxs)
#
#     test_case = find_peaks(np.abs(np.diff(sci_data_noNaN)),thres=outlier_condition/np.max(np.abs(np.diff(sci_data_noNaN))),min_dist=min_dist)
#     invalid_testcase_idxs = [i for i in range(len(test_case)) if test_case[i] in outliers]
#     test_case = np.delete(test_case,invalid_testcase_idxs)
#     if (len(test_case) != 0) & (set(test_case).issubset(discont_idxs[(discont_idxs!=0)]-1)) & (alpha_position < 0):
#         case = '1'
#     elif (len(test_case) != 0) & (set(test_case).issubset(discont_idxs[(discont_idxs!=0)]-1)) & (alpha_position > 0):
#         case = '2'
#     else:
#         case = '3'
#     print 'Case {}'.format(case)
#
#     if case in ['1','2']:
#         peaks_idxs = []
#         pseudo_wvnr = []
#         pseudo_contin = []
#         for i in range(len(discont_idxs)-1):
#             idx1,idx2 = discont_idxs[i],discont_idxs[i+1]
#             pseudo_continuum = ((sci_data_noNaN[idx2-1]-sci_data_noNaN[idx1])/(wvnr_noNaN[idx2-1]-wvnr_noNaN[idx1]))*(np.linspace(wvnr_noNaN[idx1],wvnr_noNaN[idx2],idx2-idx1)-wvnr_noNaN[idx1])+sci_data_noNaN[idx1]
#             # could try a 4th order polynomial for the pseudo_continuum, but too time consuming?..
#             peak_idxs = idx1+funcs.find_peaks(sci_data_noNaN[idx1:idx2]-pseudo_continuum,thres=thres,min_dist=min_dist)
#             invalid_peaks = [j for j in range(len(peak_idxs)) if peak_idxs[j] in outliers+1]
#             peak_idxs = np.delete(peak_idxs,invalid_peaks)
#             peaks_idxs.extend(peak_idxs)
#             pseudo_wvnr.extend(wvnr_noNaN[idx1:idx2])
#             pseudo_contin.extend(pseudo_continuum)
#     elif case == '3':
#         peaks_idxs = funcs.find_peaks(sci_data_noNaN,thres=thres,min_dist=min_dist)

def cleanRD(R,D):
    # take care of numerical instabilities
    cleanR = R.copy()
    cleanD = D.copy()
    numerics = [] # list of indexes were "cleaning" required

    # reflectivity
    #--have not found exceptions yet

    # optical thickness
    diffD = np.diff(cleanD)
    offset = np.mean(np.abs(diffD)[np.where(np.abs(diffD)*10000.>1)[0]])
    while len(np.where(np.abs(diffD)*10000.>1)[0] != 0):
        clean_idx_pos = np.where(diffD*10000.>1)[0]
        numerics.extend(clean_idx_pos+1)
        cleanD[clean_idx_pos+1] -= offset
        diffD = np.diff(cleanD)

    numerics = np.sort(np.unique(np.array(numerics)))
    return cleanR,cleanD,numerics

# find
def find_nearest(array,value):
    return np.abs(array-value).argmin()

def find_peaks(ydata, thres=0.3, min_dist=1):
    """Peak detection routine.

    Finds the numeric index of the peaks in *y* by taking its first order difference. By using
    *thres* and *min_dist* parameters, it is possible to reduce the number of
    detected peaks. *y* must be signed.

    Parameters
    ----------
    y : ndarray (signed)
        1D amplitude ydata to search for peaks.
    thres : float between [0., 1.]
        Normalized threshold. Only the peaks with amplitude higher than the
        threshold will be detected.
    min_dist : int
        Minimum distance between each detected peak. The peak with the highest
        amplitude is preferred to satisfy this constraint.

    Returns
    -------
    ndarray
        Array containing the numeric indexes of the peaks that were detected
    """
    if isinstance(ydata, np.ndarray) and np.issubdtype(ydata.dtype, np.unsignedinteger):
        raise ValueError("ydata must be signed")

    y = ydata.copy()
    y[np.isnan(y)] = 0

    thres = thres * (np.max(y) - np.min(y)) + np.min(y)
    min_dist = int(min_dist)

    # compute first order difference
    dy = np.diff(y)

    # propagate left and right values successively to fill all plateau pixels (0-value)
    zeros,=np.where(dy == 0)

    while len(zeros):
        # add pixels 2 by 2 to propagate left and right value onto the zero-value pixel
        zerosr = np.hstack([dy[1:], 0.])
        zerosl = np.hstack([0., dy[:-1]])

        # replace 0 with right value if non zero
        dy[zeros]=zerosr[zeros]
        zeros,=np.where(dy == 0)

        # replace 0 with left value if non zero
        dy[zeros]=zerosl[zeros]
        zeros,=np.where(dy == 0)

    # find the peaks by using the first order difference
    peaks = np.where((np.hstack([dy, 0.]) < 0.)
                     & (np.hstack([0., dy]) > 0.)
                     & (y > thres))[0]

    if peaks.size > 1 and min_dist > 1:
        highest = peaks[np.argsort(y[peaks])][::-1]
        rem = np.ones(y.size, dtype=bool)
        rem[peaks] = False

        for peak in highest:
            if not rem[peak]:
                sl = slice(max(0, peak - min_dist), peak + min_dist + 1)
                rem[sl] = True
                rem[peak] = False

        peaks = np.arange(y.size)[~rem]

    return peaks

def detpixel_trace(band,d2cMaps,sliceID=None,alpha_pos=None):
    # detector dimensions
    det_dims=(1024,1032)
    # initialize placeholders
    ypos,xpos = np.arange(det_dims[0]),np.zeros(det_dims[0])
    slice_img,alpha_img = [np.full(det_dims,0.) for j in range(2)]
    # create pixel masks
    sel_pix = (d2cMaps['sliceMap'] == 100*int(band[0])+sliceID) # select pixels with correct slice number
    slice_img[sel_pix] = d2cMaps['sliceMap'][sel_pix]           # image containing single slice
    alpha_img[sel_pix] = d2cMaps['alphaMap'][sel_pix]           # image containing alpha positions in single slice

    # find pixel trace
    for row in ypos:
        if band[0] in ['1','4']:
            try:xpos[row] = np.argmin(alpha_img[row,:])+find_nearest(alpha_img[row,:][(slice_img[row,:]!=0)],alpha_pos)
            except ValueError:
                """Band 1C has missing pixel values (zeros instead of valid values)"""
                continue
        elif band[0] in ['2','3']:
            xpos[row] = np.argmax(alpha_img[row,:])+find_nearest(alpha_img[row,:][(slice_img[row,:]!=0)],alpha_pos)
    xpos = xpos.astype(int)

    return ypos,xpos

def detpixel_trace_compactsource(sci_img,band,d2cMaps,offset_slice=0,verbose=False):
    # detector dimensions
    det_dims  = (1024,1032)
    nslices   = d2cMaps['nslices']
    sliceMap  = d2cMaps['sliceMap']
    ypos,xpos = np.arange(det_dims[0]),np.zeros(det_dims[0])

    sum_signals = np.zeros(nslices)
    for islice in xrange(1+nslices):
        sum_signals[islice-1] = sci_img[(sliceMap == 100*int(band[0])+islice) & (~np.isnan(sci_img))].sum()
    source_center_slice = np.argmax(sum_signals)+1
    if verbose==True:
        print 'Source center slice ID: {}'.format(source_center_slice)

    signal_img = np.full(det_dims,0.)
    sel_pix = (sliceMap == 100*int(band[0])+source_center_slice+offset_slice)
    signal_img[sel_pix] = sci_img[sel_pix]
    for row in ypos:
        xpos[row] = np.argmax(signal_img[row,:])
    xpos = xpos.astype(int)
    return ypos,xpos

def slice_alphapositions(band,d2cMaps,sliceID=None):
    import mrs_aux as maux
    # find how many alpha positions fill an entire slice
    det_dims = (1024,1032)
    mrs_fwhm  = maux.MRS_FWHM[band[0]]

    ypos = np.arange(det_dims[0])
    slice_img,alpha_img,alpha_img2 = np.full(det_dims,0.),np.full(det_dims,0.),np.full(det_dims,0.)
    sel_pix = (d2cMaps['sliceMap'] == 100*int(band[0])+sliceID) # select pixels with correct slice number
    slice_img[sel_pix] = d2cMaps['sliceMap'][sel_pix]           # image containing single slice
    alpha_img[sel_pix] = d2cMaps['alphaMap'][sel_pix]           # image containing alpha positions in single slice

    alpha_pos = alpha_img[(alpha_img!=0)].min() # arcsec
    step = mrs_fwhm/2.
    increment = mrs_fwhm/40.
    while (alpha_img2-alpha_img).any() != 0:
        xpos = np.zeros(len(ypos))
        for row in ypos:
            if band[0] in ['1','4']:
                xpos[row] = np.argmin(alpha_img[row,:])+find_nearest(alpha_img[row,:][(slice_img[row,:]!=0)],alpha_pos)
            elif band[0] in ['2','3']:
                xpos[row] = np.argmax(alpha_img[row,:])+find_nearest(alpha_img[row,:][(slice_img[row,:]!=0)],alpha_pos)
        xpos = xpos.astype(int)

        # Spectrum origination on detector
        alpha_img2[ypos,xpos] = d2cMaps['alphaMap'][ypos,xpos]
        alpha_pos += step

        if (alpha_pos > alpha_img[(alpha_img!=0)].max() + 2*mrs_fwhm):
            alpha_img2 = np.full(det_dims,0.)
            alpha_pos = alpha_img[(alpha_img!=0)].min()
            step -= increment

        if step <= 0:
            alpha_img2 = np.full(det_dims,0.)
            alpha_pos = alpha_img[(alpha_img!=0)].min()
            step = 0.2
            increment /= 2.

    alpha_positions = np.arange(alpha_img[(alpha_img!=0)].min(),alpha_pos,step)

    rmv_positions = []
    for j in range(len(alpha_positions)-1):
        xpos1 = np.zeros(len(ypos))
        for row in ypos:
            if band[0] in ['1','3']:
                xpos1[row] = np.argmin(alpha_img[row,:])+find_nearest(alpha_img[row,:][(slice_img[row,:]!=0)],alpha_positions[j])
            elif band[0] in ['2','4']:
                xpos1[row] = np.argmax(alpha_img[row,:])+find_nearest(alpha_img[row,:][(slice_img[row,:]!=0)],alpha_positions[j])

        xpos2 = np.zeros(len(ypos))
        for row in ypos:
            if band[0] in ['1','3']:
                xpos2[row] = np.argmin(alpha_img[row,:])+find_nearest(alpha_img[row,:][(slice_img[row,:]!=0)],alpha_positions[j+1])
            elif band[0] in ['2','4']:
                xpos2[row] = np.argmax(alpha_img[row,:])+find_nearest(alpha_img[row,:][(slice_img[row,:]!=0)],alpha_positions[j+1])

        if np.array_equal(xpos1, xpos2):
            # print j
            rmv_positions.append(j+1)
    new_alpha_positions = np.delete(alpha_positions,rmv_positions)

    return new_alpha_positions


# plot
def plot_point_source_centroiding(band=None,sci_img=None,d2cMaps=None,spec_grid=None,centroid=None,ibin=None,data=None):
    import mrs_aux as maux
    # distortion maps
    sliceMap  = d2cMaps['sliceMap']
    lambdaMap = d2cMaps['lambdaMap']
    alphaMap  = d2cMaps['alphaMap']
    betaMap   = d2cMaps['betaMap']
    nslices   = maux.MRS_nslices[band[0]]
    mrs_fwhm  = maux.MRS_FWHM[band[0]]
    unique_betas = np.sort(np.unique(betaMap[(sliceMap>100*int(band[0])) & (sliceMap<100*(int(band[0])+1))]))
    fov_lims  = [alphaMap[np.nonzero(lambdaMap)].min(),alphaMap[np.nonzero(lambdaMap)].max()]
    lambcens,lambfwhms = spec_grid[0],spec_grid[1]
    sign_amp,alpha_centers,beta_centers,sigma_alpha,sigma_beta = centroid[0],centroid[1],centroid[2],centroid[3],centroid[4]

    # across-slice center:
    sum_signals = np.zeros(nslices)
    for islice in xrange(1+nslices):
        sum_signals[islice-1] = sci_img[(sliceMap == 100*int(band[0])+islice) & (~np.isnan(sci_img))].sum()
    source_center_slice = np.argmax(sum_signals)+1

    # along-slice center:
    det_dims = (1024,1032)
    img = np.full(det_dims,0.)
    sel = (sliceMap == 100*int(band[0])+source_center_slice)
    img[sel]  = sci_img[sel]

    first_nonzero_row = 0
    while all(img[first_nonzero_row,:][~np.isnan(img[first_nonzero_row,:])] == 0.): first_nonzero_row+=1
    source_center_alpha = alphaMap[first_nonzero_row,img[first_nonzero_row,:].argmax()]

    # plot centroiding process in a single bin
    fig,axs = plt.subplots(2,1,figsize=(12,10))
    coords = np.where((sliceMap == 100*int(band[0])+source_center_slice) & (np.abs(lambdaMap-lambcens[ibin])<=lambfwhms[ibin]/2.) & (~np.isnan(sci_img)))
    popt,pcov = curve_fit(gauss1d_wBaseline, alphaMap[coords], sci_img[coords], p0=[sci_img[coords].max(),alpha_centers[ibin],mrs_fwhm/2.355,0.],method='lm')
    testx = np.linspace(alphaMap[coords].min(),alphaMap[coords].max(),1000)
    testy = gauss1d_wBaseline(testx,*popt)

    axs[0].plot(alphaMap[coords], sci_img[coords],'bo',label='along-slice data')
    axs[0].plot(testx, testy,'r',label='1D Gauss fit')
    axs[0].plot(testx,gauss1d_wBaseline(testx,popt[0],popt[1],0.31*(lambcens[ibin]/8.)/2.355,popt[3]),alpha=0.4,label='diffraction-limited PSF')
    axs[0].vlines([alpha_centers[ibin]-3*sigma_alpha[ibin].max(),alpha_centers[ibin]+3*sigma_alpha[ibin].max()],testy.min(),testy.max(),label=r'3$\sigma$ lines')
    axs[0].set_xlim(fov_lims[0],fov_lims[1])
    axs[0].tick_params(axis='both',labelsize=20)
    axs[0].set_xlabel(r'along-slice direction $\alpha$ [arcsec]',fontsize=20)
    if data == 'slope': axs[0].set_ylabel('Signal [DN/sec]',fontsize=20)
    elif data == 'surfbright': axs[0].set_ylabel('Signal [mJy/arcsec^2]',fontsize=20)
    axs[0].legend(loc='best',fontsize=14)

    # across-slice source centroiding
    sel = (np.abs(lambdaMap-lambcens[ibin])<=lambfwhms[ibin]/2.) & (~np.isnan(sci_img))
    signals = np.zeros(nslices)
    for islice in range(1,1+nslices):
        signals[islice-1] = sci_img[(sliceMap == 100*int(band[0])+islice) & sel][np.abs(alphaMap[(sliceMap == 100*int(band[0])+islice) & sel]-popt[1]).argmin()]
    popt,pcov = curve_fit(gauss1d_wBaseline, unique_betas, signals, p0=[signals.max(),beta_centers[ibin],mrs_fwhm/2.355,0],method='lm')
    testx = np.linspace(unique_betas.min(),unique_betas.max(),1000)
    testy = gauss1d_wBaseline(testx,*popt)

    axs[1].plot(unique_betas, signals,'bo',label='across-slice data')
    axs[1].plot(testx, testy,'r',label='1D Gauss fit')
    axs[1].plot(testx,gauss1d_wBaseline(testx,popt[0],popt[1],0.31*(lambcens[ibin]/8.)/2.355,popt[3]),alpha=0.4,label='diffraction-limited PSF')
    axs[1].vlines([beta_centers[ibin]-3*sigma_beta[ibin].max(),beta_centers[ibin]+3*sigma_beta[ibin].max()],testy.min(),testy.max(),label=r'3$\sigma$ lines')
    axs[1].tick_params(axis='both',labelsize=20)
    axs[1].set_xlabel(r'across-slice direction $\beta$ [arcsec]',fontsize=20)
    if data == 'slope': axs[1].set_ylabel('Signal [DN/sec]',fontsize=20)
    elif data == 'surfbright': axs[1].set_ylabel('Signal [mJy/arcsec^2]',fontsize=20)
    axs[1].legend(loc='best',fontsize=14)
    plt.suptitle('1D Gaussian Fitting',fontsize=20)
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])

    # initial guess for fitting, informed by previous centroiding steps
    amp,alpha0,beta0  = sign_amp[ibin],alpha_centers[ibin],beta_centers[ibin]
    sigma_alpha0, sigma_beta0 = sigma_alpha[ibin], sigma_beta[ibin]
    base = 0.
    guess = [amp, alpha0, beta0, sigma_alpha0, sigma_beta0, base]
    bounds = ([0,-np.inf,-np.inf,0,0,-np.inf],[np.inf,np.inf,np.inf,np.inf,np.inf,np.inf])

    # data to fit
    coords = (np.abs(lambdaMap-lambcens[ibin])<lambfwhms[ibin]/2.)
    alphas, betas, zobs   = alphaMap[coords],betaMap[coords],sci_img[coords]
    alphabetas = np.array([alphas,betas])

    # projected grid
    betai, alphai = np.mgrid[unique_betas.min():unique_betas.max():300j, fov_lims[0]:fov_lims[1]:300j]
    alphabetai = np.vstack([alphai.ravel(), betai.ravel()])

    zpred = gauss2d(alphabetai, *guess)
    zpred.shape = alphai.shape

    # plot result
    fig, ax = plt.subplots(figsize=(8,6))
    ax.scatter(alphas, betas, c=zobs, s=50)
    im = ax.imshow(zpred, extent=[alphai.min(), alphai.max(), betai.max(), betai.min()],
                   aspect='auto')
    cbar = fig.colorbar(im)
    cbar.ax.tick_params(labelsize=14)
    if data == 'slope': cbar.set_label(r'Signal [DN/sec]', labelpad=30,rotation=270,fontsize=16)
    elif data == 'surfbright': cbar.set_label(r'Signal [mJy/arcsec^2]', labelpad=30,rotation=270,fontsize=16)
    plt.xlabel(r'Along-slice direction (alpha) [arcsec]',fontsize=16)
    plt.ylabel(r'Across-slice direction (beta) [arcsec]',fontsize=16)
    plt.tick_params(axis='both',labelsize=20)
    ax.invert_yaxis()
    plt.suptitle('2D Gaussian Fitting',fontsize=20)
    plt.tight_layout(rect=[0, 0.03, 1, 0.96])

    # make wireframe plot
    from mpl_toolkits.mplot3d import Axes3D

    fig = plt.figure(figsize=(10,7))
    ax = fig.gca(projection='3d')
    ax.scatter(alphas, betas, zobs)
    ax.plot_wireframe(alphai,betai, zpred,color='r',alpha=0.15)
    ax.set_xlim(fov_lims[0],fov_lims[1])
    ax.set_ylim(unique_betas.min(),unique_betas.max())
    ax.set_zlim(0)
    ax.set_xlabel('Along-slice direction [arcsec]',fontsize=16)
    ax.set_ylabel('Across-slice direction [arcsec]',fontsize=16)
    if data == 'slope': ax.set_zlabel('Signal [DN/sec]',fontsize=16)
    if data == 'surfbright': ax.set_zlabel('Signal [mJy/arcsec^2]',fontsize=16)
    ax.text2D(0.15, 0.85, str(round(lambcens[ibin],2))+'um', transform=ax.transAxes,fontsize=20)
    ax.tick_params(axis='both',labelsize=10)
    plt.suptitle('Wireframe plot',fontsize=20)
    plt.tight_layout(rect=[0, 0.03, 1, 0.98])
    plt.show()

def tick_function(X):
    # for plotting wavelengths and wavenumbers on the same plot (two x-axes)
    V = 10000./X
    return ["%.2f" % z for z in V]

# optical coefficients
def indexOfRefractionZnS(wav):
    """ Index of refraction of Zinc Sulfide (AR coatings) according to M. R. Querry. "Optical constants of minerals and other materials from the millimeter to the ultraviolet"
    Data source: https://refractiveindex.info """
    wav_data,n_data,k_data = np.genfromtxt('/Users/ioannisa/Desktop/python/miri_devel/refractive_index_ZnS.txt',usecols=(0,1,2),delimiter=',',unpack=True)

    interp_n  = scp_interpolate.interp1d(wav_data,n_data)
    interp_k  = scp_interpolate.interp1d(wav_data,k_data)

    try:
        n,k = [],[]
        for wvl in wav:
            n.append(interp_n(wvl))
            k.append(interp_k(wvl))
        n,k = np.array(n),np.array(k)
    except TypeError:
        n,k = interp_n(wav),interp_k(wav)

    return n+k*1j

def indexOfRefractionSi(wav):
    # Salzberg and Villa 1957: n 1.36-11 microns
    wav2= wav**2
    C1 = 10.6684293
    C2 = (0.301516485)**2
    C3 = 0.003043475
    C4 = (1.13475115)**2
    C5 = 1.54133408
    C6 = (1104.0)**2
    n = np.sqrt( 1 + C1*wav2/(wav2-C2) + C3*wav2/(wav2-C4) + C5*wav2/(wav2-C6) )
    wav0 = (10.6)**2
    n_10 = np.sqrt( 1 + C1*wav0/(wav0-C2) + C3*wav0/(wav0-C4) + C5*wav0/(wav0-C6) )
    n = n * 3.38966 / n_10

    # Chandler-Horowitz and Amirtharaj 2005
    n = np.sqrt(11.67316 + (1/wav**2) + (0.004482633/(wav**2 - 1.108205**2)) )

    wav_data,k_data = np.genfromtxt('/Users/ioannisa/Desktop/python/miri_devel/extinction_coeff_silicon.txt',usecols=(0,1),delimiter=',',unpack=True)
    interp_k  = scp_interpolate.interp1d(wav_data,k_data)

    try:
        k = []
        for wvl in wav:
            if (wvl <wav_data[0]) or (wvl >wav_data[-1]):
                k.append(0)
            else:
                k.append(interp_k(wvl))
        k = np.array(k)
    except TypeError:
        if (wav <wav_data[0]) or (wav >wav_data[-1]):
            k = 0
        else:
            k = interp_k(wav)

    return n+k*1j

def indexOfRefractionSiAs(wav):
    # real component of index of refraction (assume refractive index of pure Silicon, due to lack of data)
    n = indexOfRefractionSi(wav)
    # imaginary component of index of refraction (data from "qe_report_new_rev.pdf", sent to me by George Rieke)
    absorption_coeff = 102.*(wav/7.)**2 # [cm-1]
    k = absorption_coeff*wav*1e-4/(4*np.pi)
    # or directly
    k = 5.69*10**-3 * (wav/7.)**3
    return n+k*1j

def indexOfRefractionAl(wav):
    wav_data,n_data,k_data = np.genfromtxt('/Users/ioannisa/Desktop/python/miri_devel/extinction_coeff_aluminium.txt',usecols=(0,1,2),delimiter=',',unpack=True)
    interp_n  = scp_interpolate.interp1d(wav_data,n_data)
    interp_k  = scp_interpolate.interp1d(wav_data,k_data)
    try:
        n,k = [],[]
        for wvl in wav:
            n.append(interp_n(wvl))
            k.append(interp_k(wvl))
        n,k = np.array(n),np.array(k)
    except TypeError:
        n,k = interp_n(wav),interp_k(wav)

    return n # +k*1j

def indexOfRefractionTransCont(wav):
    # Index of refraction of transparent contact
    n = indexOfRefractionSi(wav)-0.1
    k = 0.
    return n+k*1j

def indexOfRefractionCdTe(wav):
    """ Index of refraction of Cadmium Telluride (dichroic)
    Data source: refractiveindex.info"""
    # wav is wavelength in microns
    return np.sqrt(1 + (6.0599879*wav**2)/(wav**2 - 0.1004272) + (3.7564378*wav**2)/(wav**2 - 6138.789))

def indexOfRefractionAl2O3(wav):
    """ Index of refraction of Aluminium Oxide (front surface of BiB detector (after buried layer and front contact, as seen in Woods et al. 2011))
    Data source: https://refractiveindex.info """
    wav_data,n_data,k_data = np.genfromtxt('/Users/ioannisa/Desktop/python/miri_devel/opticalconstants_Al2O3.txt',usecols=(0,1,2),delimiter=',',unpack=True)

    interp_n  = scp_interpolate.interp1d(wav_data,n_data)
    interp_k  = scp_interpolate.interp1d(wav_data,k_data)

    try:
        n,k = [],[]
        for wvl in wav:
            n.append(interp_n(wvl))
            k.append(interp_k(wvl))
        n,k = np.array(n),np.array(k)
    except TypeError:
        n,k = interp_n(wav),interp_k(wav)

#     if k<0:
#         k=np.abs(k)
    return n+k*1j

def indexOfRefractionZnSe(wav):
    """ Index of refraction of Zinc Selenide (etalons used by INTA to produce FTS measurements)
    Data source: https://refractiveindex.info """
    wav_data,n_data,k_data = np.genfromtxt('/Users/ioannisa/Desktop/python/miri_devel/opticalconstants_ZnSe.txt',usecols=(0,1,2),delimiter=',',unpack=True)

    interp_n  = scp_interpolate.interp1d(wav_data,n_data)
    interp_k  = scp_interpolate.interp1d(wav_data,k_data)

    try:
        n,k = [],[]
        for wvl in wav:
            n.append(interp_n(wvl))
            k.append(interp_k(wvl))
        n,k = np.array(n),np.array(k)
    except TypeError:
        n,k = interp_n(wav),interp_k(wav)

    return n+k*1j

def indexOfRefractionSiO2(wav):
    """ Index of refraction of Silicon Oxide
    Data source: https://refractiveindex.info """
    wav_data,n_data,k_data = np.genfromtxt('/Users/ioannisa/Desktop/python/miri_devel/refractive_index_SiO2.txt',usecols=(0,1,2),delimiter=',',unpack=True)

    interp_n  = scp_interpolate.interp1d(wav_data,n_data)
    interp_k  = scp_interpolate.interp1d(wav_data,k_data)

    try:
        n,k = [],[]
        for wvl in wav:
            if (wvl <wav_data[0]):
                n.append(np.sqrt( 1 + (0.6961663*wvl**2 / (wvl**2 - 0.0684043**2)) + (0.4079426*wvl**2 / (wvl**2 - 0.1162414**2)) + (0.8974794*wvl**2 / (wvl**2 - 9.896161**2))) )
                k.append(0)
            else:
                n.append(interp_n(wvl))
                k.append(interp_k(wvl))
        n,k = np.array(n),np.array(k)
    except TypeError:
        if (wav <wav_data[0]):
            n = np.sqrt( 1 + (0.6961663*wav**2 / (wav**2 - 0.0684043**2)) + (0.4079426*wav**2 / (wav**2 - 0.1162414**2)) + (0.8974794*wav**2 / (wav**2 - 9.896161**2)))
            k = 0
        else:
            n,k = interp_n(wav),interp_k(wav)

    return n+k*1j

def ALrefract_index(dosage):
    # Data from "Refractive index and extinction coefficient of doped polycrystalline silicon films in infrared spectrum"
    # Note that phosphorus is used as the doping element (i.e. not Arsenic as in MIRI's detectors)
    Implanted_dosage = np.array([1.,6.,16.,51.])*10**14
    ALrefract_index = np.array([3.3,3.27,3.14,2.67])
    ALextinct_coeff = np.array([3.48e-3,1.09e-2,2.23e-2,2.04e-1])
    ALrefract_index_interpolator = scp_interpolate.InterpolatedUnivariateSpline(Implanted_dosage,ALrefract_index,k=2,ext=0)
    ALextinct_coeff_interpolator = scp_interpolate.InterpolatedUnivariateSpline(Implanted_dosage,ALextinct_coeff,k=2,ext=0)
    return ALrefract_index_interpolator(dosage) + ALextinct_coeff_interpolator(dosage)*1j

def buriedelectrode_transmission(workDir=None):
    wav_data,transmission = np.genfromtxt(workDir+'transp_contact_transm_5e14implant_poly.txt',skip_header=3,usecols=(0,1),delimiter='',unpack=True)
    transmission /= 100.
    # wav_data in micron
    # transmission normalized to 1
    return wav_data,transmission

def SW_ARcoat_reflectance(workDir=None):
    wav_data,reflectance = np.genfromtxt(workDir+'SW_ARcoat_reflectance.txt',skip_header=4,usecols=(0,1),delimiter=',',unpack=True)
    # wav_data in micron
    # reflectance normalized to 1
    return wav_data,reflectance

# transfer matrix method
def simple_tmm(n_list,d_list,th_0,lambda_vacuum):
    from scipy import arcsin
    #------------------
    num_layers = n_list.size
    # th_list is a list with, for each layer, the angle that the light travels
    # through the layer. Computed with Snell's law. Note that the "angles" may be
    # complex!
    th_list = arcsin(n_list[0]*np.sin(th_0) / n_list)
    # kz is the z-component of (complex) angular wavevector for forward-moving
    # wave. Positive imaginary part means decaying.
    kz_list = 2 * np.pi * n_list * np.cos(th_list) / lambda_vacuum
    # delta is the total phase accrued by traveling through a given layer.
    # Ignore warning about inf multiplication
    delta = kz_list * d_list
    #------------------
    # t_list[i,j] and r_list[i,j] are transmission and reflection amplitudes,
    # respectively, coming from i, going to j. Only need to calculate this when
    # j=i+1. (2D array is overkill but helps avoid confusion.)

    # s-polarization
    t_list_spol = np.zeros((num_layers, num_layers), dtype=complex)
    r_list_spol = np.zeros((num_layers, num_layers), dtype=complex)
    for i in range(num_layers-1):
        t_list_spol[i,i+1] = 2 * n_list[i] * np.cos(th_list[i]) / (n_list[i] * np.cos(th_list[i]) + n_list[i+1] * np.cos(th_list[i+1]))
        r_list_spol[i,i+1] = ((n_list[i] * np.cos(th_list[i]) - n_list[i+1] * np.cos(th_list[i+1])) / (n_list[i] * np.cos(th_list[i]) + n_list[i+1] * np.cos(th_list[i+1])))

    # p-polarization
    t_list_ppol = np.zeros((num_layers, num_layers), dtype=complex)
    r_list_ppol = np.zeros((num_layers, num_layers), dtype=complex)
    for i in range(num_layers-1):
        t_list_ppol[i,i+1] = 2 * n_list[i] * np.cos(th_list[i]) / (n_list[i+1] * np.cos(th_list[i]) + n_list[i] * np.cos(th_list[i+1]))
        r_list_ppol[i,i+1] = ((n_list[i+1] * np.cos(th_list[i]) - n_list[i] * np.cos(th_list[i+1])) / (n_list[i+1] * np.cos(th_list[i]) + n_list[i] * np.cos(th_list[i+1])))
    #------------------
    def make_2x2_array(a, b, c, d, dtype=float):
        my_array = np.empty((2,2), dtype=dtype)
        my_array[0,0] = a
        my_array[0,1] = b
        my_array[1,0] = c
        my_array[1,1] = d
        return my_array

    # At the interface between the (n-1)st and nth material, let v_n be the
    # amplitude of the wave on the nth side heading forwards (away from the
    # boundary), and let w_n be the amplitude on the nth side heading backwards
    # (towards the boundary). Then (v_n,w_n) = M_n (v_{n+1},w_{n+1}). M_n is
    # M_list[n]. M_0 and M_{num_layers-1} are not defined.
    # My M is a bit different than Sernelius's, but Mtilde is the same.
    M_list_spol = np.zeros((num_layers, 2, 2), dtype=complex)
    for i in range(1, num_layers-1):
        M_list_spol[i] = (1/t_list_spol[i,i+1]) * np.dot(
            make_2x2_array(np.exp(-1j*delta[i]), 0, 0, np.exp(1j*delta[i]),
                           dtype=complex),
            make_2x2_array(1, r_list_spol[i,i+1], r_list_spol[i,i+1], 1, dtype=complex))
    Mtilde_spol = make_2x2_array(1, 0, 0, 1, dtype=complex)
    for i in range(1, num_layers-1):
        Mtilde_spol = np.dot(Mtilde_spol, M_list_spol[i])
    Mtilde_spol = np.dot(make_2x2_array(1, r_list_spol[0,1], r_list_spol[0,1], 1,
                                   dtype=complex)/t_list_spol[0,1], Mtilde_spol)

    M_list_ppol = np.zeros((num_layers, 2, 2), dtype=complex)
    for i in range(1, num_layers-1):
        M_list_ppol[i] = (1/t_list_ppol[i,i+1]) * np.dot(
            make_2x2_array(np.exp(-1j*delta[i]), 0, 0, np.exp(1j*delta[i]),
                           dtype=complex),
            make_2x2_array(1, r_list_ppol[i,i+1], r_list_ppol[i,i+1], 1, dtype=complex))
    Mtilde_ppol = make_2x2_array(1, 0, 0, 1, dtype=complex)
    for i in range(1, num_layers-1):
        Mtilde_ppol = np.dot(Mtilde_ppol, M_list_ppol[i])
    Mtilde_ppol = np.dot(make_2x2_array(1, r_list_ppol[0,1], r_list_ppol[0,1], 1,
                                   dtype=complex)/t_list_ppol[0,1], Mtilde_ppol)
    #------------------
    # Net complex transmission and reflection amplitudes
    r_spol = Mtilde_spol[1,0]/Mtilde_spol[0,0]
    t_spol = 1/Mtilde_spol[0,0]

    r_ppol = Mtilde_ppol[1,0]/Mtilde_ppol[0,0]
    t_ppol = 1/Mtilde_ppol[0,0]
    #------------------
    # Net transmitted and reflected power, as a proportion of the incoming light
    # power.
    R_spol = abs(r_spol)**2
    T_spol = abs(t_spol**2) * (((n_list[-1]*np.cos(th_list[-1])).real) / (n_list[0]*np.cos(th_0)).real)
    power_entering_spol = ((n_list[0]*np.cos(th_0)*(1+np.conj(r_spol))*(1-r_spol)).real
                         / (n_list[0]*np.cos(th_0)).real)

    R_ppol = abs(r_ppol)**2
    T_ppol = abs(t_ppol**2) * (((n_list[-1]*np.conj(np.cos(th_list[-1]))).real) / (n_list[0]*np.conj(np.cos(th_0))).real)
    power_entering_ppol = ((n_list[0]*np.conj(np.cos(th_0))*(1+r_ppol)*(1-np.conj(r_ppol))).real
                          / (n_list[0]*np.conj(np.cos(th_0))).real)
    #------------------
    # Calculates reflected and transmitted power for unpolarized light.
    R = (R_spol + R_ppol) / 2.
    T = (T_spol + T_ppol) / 2.
    A = 1-R-T

    return R,T,A

# save and load objects
def save_obj(obj,name,path='' ):
    with open(path+name + '.pkl', 'wb') as f:
        pickle.dump(obj, f, pickle.HIGHEST_PROTOCOL)

def load_obj(name,path='' ):
    with open(path+name + '.pkl', 'rb') as f:
        return pickle.load(f)
