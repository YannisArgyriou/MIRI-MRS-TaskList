# -*- coding: utf-8 -*-
"""

Auxiliary data for mrs.

Probably this could be superseded with access to a similar module elsewhere
in "miri" package: TBD.


CAVEAT: Some values may need update.

:History:

Created on Fri May 23 18:19:48 2014

@author: Ruyman Azzollini (DIAS, ruyman.azzollini@gmail.com)
"""


allbands = ['1A','1B','1C','2A','2B','2C','3A','3B','3C','4A','4B','4C']
allchannels = ['1','2','3','4']
allsubchannels = ['A','B','C']

MRS_bands = {'1A':[4.83,5.82],
    '1B':[5.62,6.73],
    '1C':[6.46,7.76],
    '2A':[7.44,8.90],
    '2B':[8.61,10.28],
    '2C':[9.94,11.87],
    '3A':[11.47,13.67],
    '3B':[13.25,15.80],
    '3C':[15.30,18.24],
    '4A':[17.54,21.10],
    '4B':[20.44,24.72],
    '4C':[23.84,28.82]} # microns

MRS_R = {'1A':[3320.,3710.],
    '1B':[3190.,3750.],
    '1C':[3100.,3610.],
    '2A':[2990.,3110.],
    '2B':[2750.,3170.],
    '2C':[2860.,3300.],
    '3A':[2530.,2880.],
    '3B':[1790.,2640.],
    '3C':[1980.,2790.],
    '4A':[1460.,1930.],
    '4B':[1680.,1770.],
    '4C':[1630.,1330.]} # R = lambda / delta_lambda


MRS_alphapix = {'1':0.196,'2':0.196,'3':0.245,'4':0.273} # arcseconds

MRS_slice = {'1':0.176,'2':0.277,'3':0.387,'4':0.645} # arcseconds

MRS_nslices = {'1':21,'2':17,'3':16,'4':12} # Number of slices in channel

MRS_FOV = {'1':[3.70,3.70],'2':[4.51,4.71],'3':[6.13,6.19],
           '4':[7.74,7.74]} # arcseconds along and across slices

MRS_FWHM = {'1':2.16*MRS_alphapix['1'],'2':3.30*MRS_alphapix['2'],
            '3':4.04*MRS_alphapix['3'],'4':5.56*MRS_alphapix['4']} # MRS PSF
