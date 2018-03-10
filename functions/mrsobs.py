# -*- coding: utf-8 -*-
"""
File paths to mrs observations.
Author: Ioannis Argyriou (Institute of Astronomy, KU Leuven, ioannis.argyriou@kuleuven.be)
"""

# some trivia
allbands = ['1A','1B','1C','2A','2B','2C','3A','3B','3C','4A','4B','4C']
allchannels = ['1','2','3','4']
allsubchannels = ['A','B','C']

Etalons_names = ['ET_1A','ET_1B','ET_2A','ET_2B']

# filepaths
#-- FM campaign
def FM_MTS_BB_extended_source(lvl2path,band,bb_temp=None,output='img'):
    if bb_temp == '800K':
        sci_imgs = {"1A":lvl2path +'FM1T00011282/MIRFM1T00011282_1_495_SE_2011-05-31T02h15m32_LVL2.fits',
                    "1B":lvl2path +'FM1T00011283/MIRFM1T00011283_1_495_SE_2011-05-31T03h12m30_LVL2.fits',
                    "1C":lvl2path +'FM1T00011284/MIRFM1T00011284_1_495_SE_2011-05-31T04h09m25_LVL2.fits',
                    "2A":lvl2path +'FM1T00011282/MIRFM1T00011282_1_495_SE_2011-05-31T02h15m32_LVL2.fits',
                    "2B":lvl2path +'FM1T00011283/MIRFM1T00011283_1_495_SE_2011-05-31T03h12m30_LVL2.fits',
                    "2C":lvl2path +'FM1T00011284/MIRFM1T00011284_1_495_SE_2011-05-31T04h09m25_LVL2.fits',
                    "3A":lvl2path +'FM1T00011282/MIRFM1T00011282_1_494_SE_2011-05-31T02h15m02_LVL2.fits',
                    "3B":lvl2path +'FM1T00011283/MIRFM1T00011283_1_494_SE_2011-05-31T03h11m59_LVL2.fits',
                    "3C":lvl2path +'FM1T00011284/MIRFM1T00011284_1_494_SE_2011-05-31T04h08m55_LVL2.fits',
                    "4A":lvl2path +'FM1T00011282/MIRFM1T00011282_1_494_SE_2011-05-31T02h15m02_LVL2.fits',
                    "4B":lvl2path +'FM1T00011283/MIRFM1T00011283_1_494_SE_2011-05-31T03h11m59_LVL2.fits',
                    "4C":lvl2path +'FM1T00011284/MIRFM1T00011284_1_494_SE_2011-05-31T04h08m55_LVL2.fits'}

        bkg_imgs = {"1A":lvl2path +'FM1T00011285/MIRFM1T00011285_1_495_SE_2011-05-31T05h06m47_LVL2.fits',
                    "1B":lvl2path +'FM1T00011286/MIRFM1T00011286_1_495_SE_2011-05-31T06h03m43_LVL2.fits',
                    "1C":lvl2path +'FM1T00011287/MIRFM1T00011287_1_495_SE_2011-05-31T07h00m44_LVL2.fits',
                    "2A":lvl2path +'FM1T00011285/MIRFM1T00011285_1_495_SE_2011-05-31T05h06m47_LVL2.fits',
                    "2B":lvl2path +'FM1T00011286/MIRFM1T00011286_1_495_SE_2011-05-31T06h03m43_LVL2.fits',
                    "2C":lvl2path +'FM1T00011287/MIRFM1T00011287_1_495_SE_2011-05-31T07h00m44_LVL2.fits',
                    "3A":lvl2path +'FM1T00011285/MIRFM1T00011285_1_494_SE_2011-05-31T05h06m17_LVL2.fits',
                    "3B":lvl2path +'FM1T00011286/MIRFM1T00011286_1_494_SE_2011-05-31T06h03m14_LVL2.fits',
                    "3C":lvl2path +'FM1T00011287/MIRFM1T00011287_1_494_SE_2011-05-31T07h00m15_LVL2.fits',
                    "4A":lvl2path +'FM1T00011285/MIRFM1T00011285_1_494_SE_2011-05-31T05h06m17_LVL2.fits',
                    "4B":lvl2path +'FM1T00011286/MIRFM1T00011286_1_494_SE_2011-05-31T06h03m14_LVL2.fits',
                    "4C":lvl2path +'FM1T00011287/MIRFM1T00011287_1_494_SE_2011-05-31T07h00m15_LVL2.fits'}
    elif bb_temp == '600K':
        sci_imgs = {"1A":lvl2path +'FM1T00011982/MIRFM1T00011982_1_495_SE_2011-06-19T03h57m48_LVL2.fits',
                    "1B":lvl2path +'FM1T00011984/MIRFM1T00011984_1_495_SE_2011-06-19T05h55m41_LVL2.fits',
                    "1C":lvl2path +'FM1T00011986/MIRFM1T00011986_1_495_SE_2011-06-19T07h53m34_LVL2.fits',
                    "2A":lvl2path +'FM1T00011982/MIRFM1T00011982_1_495_SE_2011-06-19T03h57m48_LVL2.fits',
                    "2B":lvl2path +'FM1T00011984/MIRFM1T00011984_1_495_SE_2011-06-19T05h55m41_LVL2.fits',
                    "2C":lvl2path +'FM1T00011986/MIRFM1T00011986_1_495_SE_2011-06-19T07h53m34_LVL2.fits',
                    "3A":lvl2path +'FM1T00011982/MIRFM1T00011982_1_494_SE_2011-06-19T03h57m18_LVL2.fits',
                    "3B":lvl2path +'FM1T00011984/MIRFM1T00011984_1_494_SE_2011-06-19T05h55m11_LVL2.fits',
                    "3C":lvl2path +'FM1T00011986/MIRFM1T00011986_1_494_SE_2011-06-19T07h53m05_LVL2.fits',
                    "4A":lvl2path +'FM1T00011982/MIRFM1T00011982_1_494_SE_2011-06-19T03h57m18_LVL2.fits',
                    "4B":lvl2path +'FM1T00011984/MIRFM1T00011984_1_494_SE_2011-06-19T05h55m11_LVL2.fits',
                    "4C":lvl2path +'FM1T00011986/MIRFM1T00011986_1_494_SE_2011-06-19T07h53m05_LVL2.fits'}

        bkg_imgs = {"1A":lvl2path +'FM1T00011983/MIRFM1T00011983_1_495_SE_2011-06-19T04h57m29_LVL2.fits',
                    "1B":lvl2path +'FM1T00011985/MIRFM1T00011985_1_495_SE_2011-06-19T06h55m53_LVL2.fits',
                    "1C":lvl2path +'FM1T00011987/MIRFM1T00011987_1_495_SE_2011-06-19T08h50m05_LVL2.fits',
                    "2A":lvl2path +'FM1T00011983/MIRFM1T00011983_1_495_SE_2011-06-19T04h57m29_LVL2.fits',
                    "2B":lvl2path +'FM1T00011985/MIRFM1T00011985_1_495_SE_2011-06-19T06h55m53_LVL2.fits',
                    "2C":lvl2path +'FM1T00011987/MIRFM1T00011987_1_495_SE_2011-06-19T08h50m05_LVL2.fits',
                    "3A":lvl2path +'FM1T00011983/MIRFM1T00011983_1_494_SE_2011-06-19T04h57m00_LVL2.fits',
                    "3B":lvl2path +'FM1T00011985/MIRFM1T00011985_1_494_SE_2011-06-19T06h55m23_LVL2.fits',
                    "3C":lvl2path +'FM1T00011987/MIRFM1T00011987_1_494_SE_2011-06-19T08h49m36_LVL2.fits',
                    "4A":lvl2path +'FM1T00011983/MIRFM1T00011983_1_494_SE_2011-06-19T04h57m00_LVL2.fits',
                    "4B":lvl2path +'FM1T00011985/MIRFM1T00011985_1_494_SE_2011-06-19T06h55m23_LVL2.fits',
                    "4C":lvl2path +'FM1T00011987/MIRFM1T00011987_1_494_SE_2011-06-19T08h49m36_LVL2.fits'}
    elif bb_temp == '400K':
        sci_imgs = {"1A":lvl2path +'FM1T00011819/MIRFM1T00011819_1_495_SE_2011-06-17T14h32m35_LVL2.fits',
                    "1B":lvl2path +'FM1T00011822/MIRFM1T00011822_1_495_SE_2011-06-17T16h41m44_LVL2.fits',
                    "1C":lvl2path +'FM1T00011817/MIRFM1T00011817_1_495_SE_2011-06-17T12h27m28_LVL2.fits',
                    "2A":lvl2path +'FM1T00011819/MIRFM1T00011819_1_495_SE_2011-06-17T14h32m35_LVL2.fits',
                    "2B":lvl2path +'FM1T00011822/MIRFM1T00011822_1_495_SE_2011-06-17T16h41m44_LVL2.fits',
                    "2C":lvl2path +'FM1T00011817/MIRFM1T00011817_1_495_SE_2011-06-17T12h27m28_LVL2.fits',
                    "3A":lvl2path +'FM1T00011819/MIRFM1T00011819_1_494_SE_2011-06-17T14h32m05_LVL2.fits',
                    "3B":lvl2path +'FM1T00011822/MIRFM1T00011822_1_494_SE_2011-06-17T16h41m14_LVL2.fits',
                    "3C":lvl2path +'FM1T00011817/MIRFM1T00011817_1_494_SE_2011-06-17T12h26m58_LVL2.fits',
                    "4A":lvl2path +'FM1T00011819/MIRFM1T00011819_1_494_SE_2011-06-17T14h32m05_LVL2.fits',
                    "4B":lvl2path +'FM1T00011822/MIRFM1T00011822_1_494_SE_2011-06-17T16h41m14_LVL2.fits',
                    "4C":lvl2path +'FM1T00011817/MIRFM1T00011817_1_494_SE_2011-06-17T12h26m58_LVL2.fits'}

        bkg_imgs = {"1A":lvl2path +'FM1T00011820/MIRFM1T00011820_1_495_SE_2011-06-17T15h29m07_LVL2.fits',
                    "1B":lvl2path +'FM1T00011823/MIRFM1T00011823_1_495_SE_2011-06-17T17h39m40_LVL2.fits',
                    "1C":lvl2path +'FM1T00011818/MIRFM1T00011818_1_495_SE_2011-06-17T13h34m18_LVL2.fits',
                    "2A":lvl2path +'FM1T00011820/MIRFM1T00011820_1_495_SE_2011-06-17T15h29m07_LVL2.fits',
                    "2B":lvl2path +'FM1T00011823/MIRFM1T00011823_1_495_SE_2011-06-17T17h39m40_LVL2.fits',
                    "2C":lvl2path +'FM1T00011818/MIRFM1T00011818_1_495_SE_2011-06-17T13h34m18_LVL2.fits',
                    "3A":lvl2path +'FM1T00011820/MIRFM1T00011820_1_494_SE_2011-06-17T15h28m37_LVL2.fits',
                    "3B":lvl2path +'FM1T00011823/MIRFM1T00011823_1_494_SE_2011-06-17T17h39m10_LVL2.fits',
                    "3C":lvl2path +'FM1T00011818/MIRFM1T00011818_1_494_SE_2011-06-17T13h33m48_LVL2.fits',
                    "4A":lvl2path +'FM1T00011820/MIRFM1T00011820_1_494_SE_2011-06-17T15h28m37_LVL2.fits',
                    "4B":lvl2path +'FM1T00011823/MIRFM1T00011823_1_494_SE_2011-06-17T17h39m10_LVL2.fits',
                    "4C":lvl2path +'FM1T00011818/MIRFM1T00011818_1_494_SE_2011-06-17T13h33m48_LVL2.fits'}
    if output == 'filename':
        return sci_imgs[band],bkg_imgs[band]
    elif output == 'img':
        from astropy.io import fits
        return fits.open(sci_imgs[band])[0].data[0,:,:],fits.open(bkg_imgs[band])[0].data[0,:,:]

def MIRI_internal_calibration_source(lvl2path,band,campaign=None,output='img'):
    """
    * 800K BB source
    * observed in different test campaigns
    * CCC closed during internal calibration source observation (no background observations)
    """
    if campaign == 'FM':
        sci_imgs = {"1A":lvl2path +'FM1T00012668/MIRFM1T00012668_1_495_SE_2011-07-06T17h08m40_LVL2.fits',
                    "1B":lvl2path +'FM1T00012668/MIRFM1T00012668_5_495_SE_2011-07-06T17h32m15_LVL2.fits',
                    "1C":lvl2path +'FM1T00012668/MIRFM1T00012668_9_495_SE_2011-07-06T17h55m46_LVL2.fits',
                    "2A":lvl2path +'FM1T00012668/MIRFM1T00012668_1_495_SE_2011-07-06T17h08m40_LVL2.fits',
                    "2B":lvl2path +'FM1T00012668/MIRFM1T00012668_5_495_SE_2011-07-06T17h32m15_LVL2.fits',
                    "2C":lvl2path +'FM1T00012668/MIRFM1T00012668_9_495_SE_2011-07-06T17h55m46_LVL2.fits',
                    "3A":lvl2path +'FM1T00012668/MIRFM1T00012668_1_494_SE_2011-07-06T17h08m24_LVL2.fits',
                    "3B":lvl2path +'FM1T00012668/MIRFM1T00012668_5_494_SE_2011-07-06T17h31m59_LVL2.fits',
                    "3C":lvl2path +'FM1T00012668/MIRFM1T00012668_9_494_SE_2011-07-06T17h55m30_LVL2.fits',
                    "4A":lvl2path +'FM1T00012668/MIRFM1T00012668_1_494_SE_2011-07-06T17h08m24_LVL2.fits',
                    "4B":lvl2path +'FM1T00012668/MIRFM1T00012668_5_494_SE_2011-07-06T17h31m59_LVL2.fits',
                    "4C":lvl2path +'FM1T00012668/MIRFM1T00012668_9_494_SE_2011-07-06T17h55m30_LVL2.fits'}
    elif campaign == 'OTIS':
        sci_imgs = {"1A":lvl2path +'MIRV00331001001P0000000002103_1_495_SE_2017-08-25T19h09m24_LVL2.fits',
                    "1B":lvl2path +'',
                    "1C":lvl2path +'',
                    "2A":lvl2path +'MIRV00331001001P0000000002103_1_495_SE_2017-08-25T19h09m24_LVL2.fits',
                    "2B":lvl2path +'',
                    "2C":lvl2path +'',
                    "3A":lvl2path +'',
                    "3B":lvl2path +'',
                    "3C":lvl2path +'',
                    "4A":lvl2path +'',
                    "4B":lvl2path +'',
                    "4C":lvl2path +''}
    elif campaign == 'CV3':
        sci_imgs = {"1A":lvl2path +'MIRM33541-A-A-8MA-6019093539_1_495_SE_2016-01-19T09h59m18_LVL2.fits',
                    "1B":lvl2path +'',
                    "1C":lvl2path +'',
                    "2A":lvl2path +'MIRM33541-A-A-8MA-6019093539_1_495_SE_2016-01-19T09h59m18_LVL2.fits',
                    "2B":lvl2path +'',
                    "2C":lvl2path +'',
                    "3A":lvl2path +'',
                    "3B":lvl2path +'',
                    "3C":lvl2path +'',
                    "4A":lvl2path +'',
                    "4B":lvl2path +'',
                    "4C":lvl2path +''}
    if output == 'filename':
        return sci_imgs[band]
    elif output == 'img':
        from astropy.io import fits
        return fits.open(sci_imgs[band])[0].data[0,:,:]

def FM_MTS_800K_BB_extended_source_through_etalon(lvl2path,band,etalon=None,output='img'):
    if etalon == 'ET1A':
        sci_imgs = {"1A":lvl2path +'FM1T00010937/MIRFM1T00010937_1_495_SE_2011-05-19T22h24m06_LVL2.fits',
                    "1B":lvl2path +'FM1T00010939/MIRFM1T00010939_1_495_SE_2011-05-20T00h04m44_LVL2.fits',
                    "1C":lvl2path +'FM1T00010843/MIRFM1T00010843_1_495_SE_2011-05-18T00h47m06_LVL2.fits',
                    "2A":lvl2path +'FM1T00010937/MIRFM1T00010937_1_495_SE_2011-05-19T22h24m06_LVL2.fits',
                    "2B":lvl2path +'FM1T00010939/MIRFM1T00010939_1_495_SE_2011-05-20T00h04m44_LVL2.fits',
                    "2C":lvl2path +'FM1T00010843/MIRFM1T00010843_1_495_SE_2011-05-18T00h47m06_LVL2.fits',
                    "3A":lvl2path +'FM1T00010937/MIRFM1T00010937_1_494_SE_2011-05-19T22h23m43_LVL2.fits',
                    "3B":lvl2path +'FM1T00010939/MIRFM1T00010939_1_494_SE_2011-05-20T00h04m21_LVL2.fits',
                    "3C":lvl2path +'FM1T00010843/MIRFM1T00010843_1_494_SE_2011-05-18T00h46m35_LVL2.fits',
                    "4A":lvl2path +'FM1T00010937/MIRFM1T00010937_1_494_SE_2011-05-19T22h23m43_LVL2.fits',
                    "4B":lvl2path +'FM1T00010939/MIRFM1T00010939_1_494_SE_2011-05-20T00h04m21_LVL2.fits',
                    "4C":lvl2path +'FM1T00010843/MIRFM1T00010843_1_494_SE_2011-05-18T00h46m35_LVL2.fits'}

        bkg_imgs = {"1A":lvl2path +'FM1T00010938/MIRFM1T00010938_1_495_SE_2011-05-19T23h16m52_LVL2.fits',
                    "1B":lvl2path +'FM1T00010940/MIRFM1T00010940_1_495_SE_2011-05-20T00h53m40_LVL2.fits',
                    "1C":lvl2path +'FM1T00010844/MIRFM1T00010844_1_495_SE_2011-05-18T01h53m59_LVL2.fits',
                    "2A":lvl2path +'FM1T00010938/MIRFM1T00010938_1_495_SE_2011-05-19T23h16m52_LVL2.fits',
                    "2B":lvl2path +'FM1T00010940/MIRFM1T00010940_1_495_SE_2011-05-20T00h53m40_LVL2.fits',
                    "2C":lvl2path +'FM1T00010844/MIRFM1T00010844_1_495_SE_2011-05-18T01h53m59_LVL2.fits',
                    "3A":lvl2path +'FM1T00010938/MIRFM1T00010938_1_494_SE_2011-05-19T23h16m30_LVL2.fits',
                    "3B":lvl2path +'FM1T00010940/MIRFM1T00010940_1_494_SE_2011-05-20T00h53m17_LVL2.fits',
                    "3C":lvl2path +'FM1T00010844/MIRFM1T00010844_1_494_SE_2011-05-18T01h53m27_LVL2.fits',
                    "4A":lvl2path +'FM1T00010938/MIRFM1T00010938_1_494_SE_2011-05-19T23h16m30_LVL2.fits',
                    "4B":lvl2path +'FM1T00010940/MIRFM1T00010940_1_494_SE_2011-05-20T00h53m17_LVL2.fits',
                    "4C":lvl2path +'FM1T00010844/MIRFM1T00010844_1_494_SE_2011-05-18T01h53m27_LVL2.fits'}
    elif etalon == 'ET1B':
        sci_imgs = {"1A":lvl2path +'FM1T00010918/MIRFM1T00010918_1_495_SE_2011-05-19T15h17m35_LVL2.fits',
                    "1B":lvl2path +'FM1T00010920/MIRFM1T00010920_1_495_SE_2011-05-19T16h52m14_LVL2.fits',
                    "1C":lvl2path +'FM1T00010943/MIRFM1T00010943_1_495_SE_2011-05-20T03h23m04_LVL2.fits',
                    "2A":lvl2path +'FM1T00010918/MIRFM1T00010918_1_495_SE_2011-05-19T15h17m35_LVL2.fits',
                    "2B":lvl2path +'FM1T00010920/MIRFM1T00010920_1_495_SE_2011-05-19T16h52m14_LVL2.fits',
                    "2C":lvl2path +'FM1T00010943/MIRFM1T00010943_1_495_SE_2011-05-20T03h23m04_LVL2.fits',
                    "3A":lvl2path +'FM1T00010918/MIRFM1T00010918_1_494_SE_2011-05-19T15h17m12_LVL2.fits',
                    "3B":lvl2path +'FM1T00010920/MIRFM1T00010920_1_494_SE_2011-05-19T16h51m50_LVL2.fits',
                    "3C":lvl2path +'FM1T00010943/MIRFM1T00010943_1_494_SE_2011-05-20T03h22m41_LVL2.fits',
                    "4A":lvl2path +'FM1T00010918/MIRFM1T00010918_1_494_SE_2011-05-19T15h17m12_LVL2.fits',
                    "4B":lvl2path +'FM1T00010920/MIRFM1T00010920_1_494_SE_2011-05-19T16h51m50_LVL2.fits',
                    "4C":lvl2path +'FM1T00010943/MIRFM1T00010943_1_494_SE_2011-05-20T03h22m41_LVL2.fits'}

        bkg_imgs = {"1A":lvl2path +'FM1T00010919/MIRFM1T00010919_1_495_SE_2011-05-19T16h04m35_LVL2.fits',
                    "1B":lvl2path +'FM1T00010921/MIRFM1T00010921_1_495_SE_2011-05-19T17h38m33_LVL2.fits',
                    "1C":lvl2path +'FM1T00010942/MIRFM1T00010942_1_495_SE_2011-05-20T02h30m42_LVL2.fits',
                    "2A":lvl2path +'FM1T00010919/MIRFM1T00010919_1_495_SE_2011-05-19T16h04m35_LVL2.fits',
                    "2B":lvl2path +'FM1T00010921/MIRFM1T00010921_1_495_SE_2011-05-19T17h38m33_LVL2.fits',
                    "2C":lvl2path +'FM1T00010942/MIRFM1T00010942_1_495_SE_2011-05-20T02h30m42_LVL2.fits',
                    "3A":lvl2path +'FM1T00010919/MIRFM1T00010919_1_494_SE_2011-05-19T16h04m13_LVL2.fits',
                    "3B":lvl2path +'FM1T00010921/MIRFM1T00010921_1_494_SE_2011-05-19T17h38m10_LVL2.fits',
                    "3C":lvl2path +'FM1T00010942/MIRFM1T00010942_1_494_SE_2011-05-20T02h30m20_LVL2.fits',
                    "4A":lvl2path +'FM1T00010919/MIRFM1T00010919_1_494_SE_2011-05-19T16h04m13_LVL2.fits',
                    "4B":lvl2path +'FM1T00010921/MIRFM1T00010921_1_494_SE_2011-05-19T17h38m10_LVL2.fits',
                    "4C":lvl2path +'FM1T00010942/MIRFM1T00010942_1_494_SE_2011-05-20T02h30m20_LVL2.fits'}
    elif etalon == 'ET2A':
        sci_imgs = {"1A":lvl2path +'FM1T00010912/MIRFM1T00010912_1_495_SE_2011-05-19T09h52m51_LVL2.fits',
                    "1B":lvl2path +'FM1T00010914/MIRFM1T00010914_1_495_SE_2011-05-19T12h06m04_LVL2.fits',
                    "1C":lvl2path +'FM1T00010916/MIRFM1T00010916_1_495_SE_2011-05-19T13h40m32_LVL2.fits',
                    "2A":lvl2path +'FM1T00010912/MIRFM1T00010912_1_495_SE_2011-05-19T09h52m51_LVL2.fits',
                    "2B":lvl2path +'FM1T00010914/MIRFM1T00010914_1_495_SE_2011-05-19T12h06m04_LVL2.fits',
                    "2C":lvl2path +'FM1T00010916/MIRFM1T00010916_1_495_SE_2011-05-19T13h40m32_LVL2.fits',
                    "3A":lvl2path +'FM1T00010912/MIRFM1T00010912_1_494_SE_2011-05-19T09h52m28_LVL2.fits',
                    "3B":lvl2path +'FM1T00010914/MIRFM1T00010914_1_494_SE_2011-05-19T12h05m42_LVL2.fits',
                    "3C":lvl2path +'FM1T00010916/MIRFM1T00010916_1_494_SE_2011-05-19T13h40m09_LVL2.fits',
                    "4A":lvl2path +'FM1T00010912/MIRFM1T00010912_1_494_SE_2011-05-19T09h52m28_LVL2.fits',
                    "4B":lvl2path +'FM1T00010914/MIRFM1T00010914_1_494_SE_2011-05-19T12h05m42_LVL2.fits',
                    "4C":lvl2path +'FM1T00010916/MIRFM1T00010916_1_494_SE_2011-05-19T13h40m09_LVL2.fits'}

        bkg_imgs = {"1A":lvl2path +'FM1T00010913/MIRFM1T00010913_1_495_SE_2011-05-19T11h07m18_LVL2.fits',
                    "1B":lvl2path +'FM1T00010915/MIRFM1T00010915_1_495_SE_2011-05-19T12h52m36_LVL2.fits',
                    "1C":lvl2path +'FM1T00010917/MIRFM1T00010917_1_495_SE_2011-05-19T14h28m24_LVL2.fits',
                    "2A":lvl2path +'FM1T00010913/MIRFM1T00010913_1_495_SE_2011-05-19T11h07m18_LVL2.fits',
                    "2B":lvl2path +'FM1T00010915/MIRFM1T00010915_1_495_SE_2011-05-19T12h52m36_LVL2.fits',
                    "2C":lvl2path +'FM1T00010917/MIRFM1T00010917_1_495_SE_2011-05-19T14h28m24_LVL2.fits',
                    "3A":lvl2path +'FM1T00010913/MIRFM1T00010913_1_494_SE_2011-05-19T11h06m56_LVL2.fits',
                    "3B":lvl2path +'FM1T00010915/MIRFM1T00010915_1_494_SE_2011-05-19T12h52m13_LVL2.fits',
                    "3C":lvl2path +'FM1T00010917/MIRFM1T00010917_1_494_SE_2011-05-19T14h28m01_LVL2.fits',
                    "4A":lvl2path +'FM1T00010913/MIRFM1T00010913_1_494_SE_2011-05-19T11h06m56_LVL2.fits',
                    "4B":lvl2path +'FM1T00010915/MIRFM1T00010915_1_494_SE_2011-05-19T12h52m13_LVL2.fits',
                    "4C":lvl2path +'FM1T00010917/MIRFM1T00010917_1_494_SE_2011-05-19T14h28m01_LVL2.fits'}
    elif etalon == 'ET2B':
        sci_imgs = {"1A":lvl2path +'FM1T00010963/MIRFM1T00010963_1_495_SE_2011-05-20T23h07m01_LVL2.fits',
                    "1B":lvl2path +'FM1T00010952/MIRFM1T00010952_1_495_SE_2011-05-20T11h03m18_LVL2.fits',
                    "1C":lvl2path +'FM1T00010954/MIRFM1T00010954_1_495_SE_2011-05-20T12h38m33_LVL2.fits',
                    "2A":lvl2path +'FM1T00010963/MIRFM1T00010963_1_495_SE_2011-05-20T23h07m01_LVL2.fits',
                    "2B":lvl2path +'FM1T00010952/MIRFM1T00010952_1_495_SE_2011-05-20T11h03m18_LVL2.fits',
                    "2C":lvl2path +'FM1T00010954/MIRFM1T00010954_1_495_SE_2011-05-20T12h38m33_LVL2.fits',
                    "3A":lvl2path +'FM1T00010963/MIRFM1T00010963_1_494_SE_2011-05-20T23h06m39_LVL2.fits',
                    "3B":lvl2path +'FM1T00010952/MIRFM1T00010952_1_494_SE_2011-05-20T11h02m55_LVL2.fits',
                    "3C":lvl2path +'FM1T00010954/MIRFM1T00010954_1_494_SE_2011-05-20T12h38m09_LVL2.fits',
                    "4A":lvl2path +'FM1T00010963/MIRFM1T00010963_1_494_SE_2011-05-20T23h06m39_LVL2.fits',
                    "4B":lvl2path +'FM1T00010952/MIRFM1T00010952_1_494_SE_2011-05-20T11h02m55_LVL2.fits',
                    "4C":lvl2path +'FM1T00010954/MIRFM1T00010954_1_494_SE_2011-05-20T12h38m09_LVL2.fits'}

        bkg_imgs = {"1A":lvl2path +'FM1T00010964/MIRFM1T00010964_1_495_SE_2011-05-21T00h03m53_LVL2.fits',
                    "1B":lvl2path +'FM1T00010953/MIRFM1T00010953_1_495_SE_2011-05-20T11h50m14_LVL2.fits',
                    "1C":lvl2path +'FM1T00010955/MIRFM1T00010955_1_495_SE_2011-05-20T13h25m22_LVL2.fits',
                    "2A":lvl2path +'FM1T00010964/MIRFM1T00010964_1_495_SE_2011-05-21T00h03m53_LVL2.fits',
                    "2B":lvl2path +'FM1T00010953/MIRFM1T00010953_1_495_SE_2011-05-20T11h50m14_LVL2.fits',
                    "2C":lvl2path +'FM1T00010955/MIRFM1T00010955_1_495_SE_2011-05-20T13h25m22_LVL2.fits',
                    "3A":lvl2path +'FM1T00010964/MIRFM1T00010964_1_494_SE_2011-05-21T00h03m30_LVL2.fits',
                    "3B":lvl2path +'FM1T00010953/MIRFM1T00010953_1_494_SE_2011-05-20T11h49m51_LVL2.fits',
                    "3C":lvl2path +'FM1T00010955/MIRFM1T00010955_1_494_SE_2011-05-20T13h24m59_LVL2.fits',
                    "4A":lvl2path +'FM1T00010964/MIRFM1T00010964_1_494_SE_2011-05-21T00h03m30_LVL2.fits',
                    "4B":lvl2path +'FM1T00010953/MIRFM1T00010953_1_494_SE_2011-05-20T11h49m51_LVL2.fits',
                    "4C":lvl2path +'FM1T00010955/MIRFM1T00010955_1_494_SE_2011-05-20T13h24m59_LVL2.fits'}
    if output == 'filename':
        return sci_imgs[band],bkg_imgs[band]
    elif output == 'img':
        from astropy.io import fits
        return fits.open(sci_imgs[band])[0].data[0,:,:],fits.open(bkg_imgs[band])[0].data[0,:,:]

def FM_MTS_800K_BB_extended_source_through_etalon_through_pinhole(lvl2path,band,etalon=None,output='img'):
    if etalon == 'ET1A':
        sci_imgs = {"1A":lvl2path +'FM1T00012163/MIRFM1T00012163_1_495_SE_2011-06-25T08h53m56_LVL2.fits',
                    "1B":lvl2path +'',
                    "1C":lvl2path +'',
                    "2A":lvl2path +'',
                    "2B":lvl2path +'',
                    "2C":lvl2path +'',
                    "3A":lvl2path +'',
                    "3B":lvl2path +'',
                    "3C":lvl2path +'',
                    "4A":lvl2path +'',
                    "4B":lvl2path +'',
                    "4C":lvl2path +''}

        bkg_imgs = {"1A":lvl2path +'',
                    "1B":lvl2path +'',
                    "1C":lvl2path +'',
                    "2A":lvl2path +'',
                    "2B":lvl2path +'',
                    "2C":lvl2path +'',
                    "3A":lvl2path +'',
                    "3B":lvl2path +'',
                    "3C":lvl2path +'',
                    "4A":lvl2path +'',
                    "4B":lvl2path +'',
                    "4C":lvl2path +''}
    elif etalon == 'ET1B':
        sci_imgs = {"1A":lvl2path +'',
                    "1B":lvl2path +'',
                    "1C":lvl2path +'',
                    "2A":lvl2path +'',
                    "2B":lvl2path +'',
                    "2C":lvl2path +'',
                    "3A":lvl2path +'',
                    "3B":lvl2path +'',
                    "3C":lvl2path +'',
                    "4A":lvl2path +'',
                    "4B":lvl2path +'',
                    "4C":lvl2path +''}

        bkg_imgs = {"1A":lvl2path +'',
                    "1B":lvl2path +'',
                    "1C":lvl2path +'',
                    "2A":lvl2path +'',
                    "2B":lvl2path +'',
                    "2C":lvl2path +'',
                    "3A":lvl2path +'',
                    "3B":lvl2path +'',
                    "3C":lvl2path +'',
                    "4A":lvl2path +'',
                    "4B":lvl2path +'',
                    "4C":lvl2path +''}
    elif etalon == 'ET2A':
        sci_imgs = {"1A":lvl2path +'',
                    "1B":lvl2path +'',
                    "1C":lvl2path +'',
                    "2A":lvl2path +'',
                    "2B":lvl2path +'',
                    "2C":lvl2path +'',
                    "3A":lvl2path +'',
                    "3B":lvl2path +'',
                    "3C":lvl2path +'',
                    "4A":lvl2path +'',
                    "4B":lvl2path +'',
                    "4C":lvl2path +''}

        bkg_imgs = {"1A":lvl2path +'',
                    "1B":lvl2path +'',
                    "1C":lvl2path +'',
                    "2A":lvl2path +'',
                    "2B":lvl2path +'',
                    "2C":lvl2path +'',
                    "3A":lvl2path +'',
                    "3B":lvl2path +'',
                    "3C":lvl2path +'',
                    "4A":lvl2path +'',
                    "4B":lvl2path +'',
                    "4C":lvl2path +''}
    elif etalon == 'ET2B':
        sci_imgs = {"1A":lvl2path +'',
                    "1B":lvl2path +'',
                    "1C":lvl2path +'',
                    "2A":lvl2path +'',
                    "2B":lvl2path +'',
                    "2C":lvl2path +'',
                    "3A":lvl2path +'',
                    "3B":lvl2path +'',
                    "3C":lvl2path +'',
                    "4A":lvl2path +'',
                    "4B":lvl2path +'',
                    "4C":lvl2path +''}

        bkg_imgs = {"1A":lvl2path +'',
                    "1B":lvl2path +'',
                    "1C":lvl2path +'',
                    "2A":lvl2path +'',
                    "2B":lvl2path +'',
                    "2C":lvl2path +'',
                    "3A":lvl2path +'',
                    "3B":lvl2path +'',
                    "3C":lvl2path +'',
                    "4A":lvl2path +'',
                    "4B":lvl2path +'',
                    "4C":lvl2path +''}
    if output == 'filename':
        return sci_imgs[band],bkg_imgs[band]
    elif output == 'img':
        from astropy.io import fits
        return fits.open(sci_imgs[band])[0].data[0,:,:],0 #fits.open(bkg_imgs[band])[0].data[0,:,:]

def RAL_FTS_ET_observations(obsDir,etalon=None):
    import numpy as np
    if etalon == 'ET1A':
        etalon1A_file = obsDir+'MIRI_b_E1av_etalon1A_77K.txt'
        wvnrs,etalon_data = np.genfromtxt(etalon1A_file, skip_footer = 1, usecols=(0,1), delimiter = '',unpack='True')
    elif etalon == 'ET1B':
        etalon1B_file = obsDir+'MIRI_b_G3av_etalon1B_80K.txt'
        wvnrs,etalon_data = np.genfromtxt(etalon1B_file, skip_footer = 1, usecols=(0,1), delimiter = '',unpack='True')
    elif etalon == 'ET2A':
        etalon2A_file = obsDir+'MIRI_b_K1av_etalon2A_80K.txt'
        wvnrs,etalon_data = np.genfromtxt(etalon2A_file, skip_footer = 1, usecols=(0,1), delimiter = '',unpack='True')
    elif etalon == 'ET2B':
        import csv
        etalon2B_file = obsDir+'MIRI_etalon2B_80K.csv'
        wvnrs,etalon_data = [],[]
        with open(etalon2B_file, 'rb') as csvfile:
            spamreader = csv.reader(csvfile, delimiter=',')
            for row in spamreader:
                wvnrs.append(np.float(row[0]) )
                etalon_data.append(np.float(row[1]) )
        wvnrs,etalon_data = np.array(wvnrs),np.array(etalon)
    return wvnrs,etalon_data

def FM_MTS_800K_BB_point_source_raster(lvl2path,position=None,pointing='all',output='img'):
    """
    * all pointings are in DGA setup 1A/2A
    """
    if position == 'left':
        # Left side of FOV
        sci_imgs = {'P1':lvl2path+'FM1T00012204/MIRFM1T00012204_1_495_SE_2011-06-27T22h27m07_LVL2.fits',
                    'P2':lvl2path+'FM1T00012204/MIRFM1T00012204_2_495_SE_2011-06-27T22h42m57_LVL2.fits',
                    'P3':lvl2path+'FM1T00012204/MIRFM1T00012204_3_495_SE_2011-06-27T22h58m48_LVL2.fits',
                    'P4':lvl2path+'FM1T00012204/MIRFM1T00012204_4_495_SE_2011-06-27T23h14m33_LVL2.fits',
                    'P5':lvl2path+'FM1T00012204/MIRFM1T00012204_5_495_SE_2011-06-27T23h30m18_LVL2.fits',
                    'P6':lvl2path+'FM1T00012204/MIRFM1T00012204_6_495_SE_2011-06-27T23h46m09_LVL2.fits',
                    'P7':lvl2path+'FM1T00012204/MIRFM1T00012204_7_495_SE_2011-06-28T00h01m59_LVL2.fits',
                    'P8':lvl2path+'FM1T00012204/MIRFM1T00012204_8_495_SE_2011-06-28T00h17m50_LVL2.fits',
                    'P9':lvl2path+'FM1T00012204/MIRFM1T00012204_9_495_SE_2011-06-28T00h33m35_LVL2.fits',
                    'P10':lvl2path+'FM1T00012204/MIRFM1T00012204_10_495_SE_2011-06-28T00h49m31_LVL2.fits',
                    'P11':lvl2path+'FM1T00012204/MIRFM1T00012204_11_495_SE_2011-06-28T01h05m16_LVL2.fits',
                    'P12':lvl2path+'FM1T00012204/MIRFM1T00012204_12_495_SE_2011-06-28T01h21m02_LVL2.fits',
                    'P13':lvl2path+'FM1T00012204/MIRFM1T00012204_13_495_SE_2011-06-28T01h36m52_LVL2.fits',
                    'P14':lvl2path+'FM1T00012204/MIRFM1T00012204_14_495_SE_2011-06-28T01h52m42_LVL2.fits',
                    'P15':lvl2path+'FM1T00012204/MIRFM1T00012204_15_495_SE_2011-06-28T02h08m33_LVL2.fits',
                    'P16':lvl2path+'FM1T00012204/MIRFM1T00012204_16_495_SE_2011-06-28T02h24m23_LVL2.fits',
                    'P17':lvl2path+'FM1T00012204/MIRFM1T00012204_17_495_SE_2011-06-28T02h40m09_LVL2.fits',
                    'P18':lvl2path+'FM1T00012204/MIRFM1T00012204_18_495_SE_2011-06-28T02h55m54_LVL2.fits',
                    'P19':lvl2path+'FM1T00012204/MIRFM1T00012204_19_495_SE_2011-06-28T03h11m50_LVL2.fits',
                    'P20':lvl2path+'FM1T00012204/MIRFM1T00012204_20_495_SE_2011-06-28T03h27m35_LVL2.fits',
                    'P21':lvl2path+'FM1T00012204/MIRFM1T00012204_21_495_SE_2011-06-28T03h43m20_LVL2.fits',
                    'P22':lvl2path+'FM1T00012204/MIRFM1T00012204_22_495_SE_2011-06-28T03h59m11_LVL2.fits',
                    'P23':lvl2path+'FM1T00012204/MIRFM1T00012204_23_495_SE_2011-06-28T04h14m56_LVL2.fits'}
        bkg_file = lvl2path+'FM1T00012203/MIRFM1T00012203_1_495_SE_2011-06-27T22h03m11_LVL2.fits'

    elif position == 'middle':
        # Middle of FOV
        sci_imgs = {'P1':lvl2path+'FM1T00012308/MIRFM1T00012308_1_495_SE_2011-06-30T22h43m36_LVL2.fits',
                    'P2':lvl2path+'FM1T00012308/MIRFM1T00012308_2_495_SE_2011-06-30T22h59m21_LVL2.fits',
                    'P3':lvl2path+'FM1T00012308/MIRFM1T00012308_3_495_SE_2011-06-30T23h15m12_LVL2.fits',
                    'P4':lvl2path+'FM1T00012308/MIRFM1T00012308_4_495_SE_2011-06-30T23h31m02_LVL2.fits',
                    'P5':lvl2path+'FM1T00012308/MIRFM1T00012308_5_495_SE_2011-06-30T23h46m53_LVL2.fits',
                    'P6':lvl2path+'FM1T00012308/MIRFM1T00012308_6_495_SE_2011-07-01T00h02m38_LVL2.fits',
                    'P7':lvl2path+'FM1T00012308/MIRFM1T00012308_7_495_SE_2011-07-01T00h18m29_LVL2.fits',
                    'P8':lvl2path+'FM1T00012308/MIRFM1T00012308_8_495_SE_2011-07-01T00h34m19_LVL2.fits',
                    'P9':lvl2path+'FM1T00012308/MIRFM1T00012308_9_495_SE_2011-07-01T00h50m05_LVL2.fits',
                    'P10':lvl2path+'FM1T00012308/MIRFM1T00012308_10_495_SE_2011-07-01T01h06m05_LVL2.fits',
                    'P11':lvl2path+'FM1T00012308/MIRFM1T00012308_11_495_SE_2011-07-01T01h21m50_LVL2.fits',
                    'P12':lvl2path+'FM1T00012308/MIRFM1T00012308_12_495_SE_2011-07-01T01h37m36_LVL2.fits',
                    'P13':lvl2path+'FM1T00012308/MIRFM1T00012308_13_495_SE_2011-07-01T01h53m21_LVL2.fits',
                    'P14':lvl2path+'FM1T00012308/MIRFM1T00012308_14_495_SE_2011-07-01T02h09m17_LVL2.fits',
                    'P15':lvl2path+'FM1T00012308/MIRFM1T00012308_15_495_SE_2011-07-01T02h25m02_LVL2.fits',
                    'P16':lvl2path+'FM1T00012308/MIRFM1T00012308_16_495_SE_2011-07-01T02h40m52_LVL2.fits',
                    'P17':lvl2path+'FM1T00012308/MIRFM1T00012308_17_495_SE_2011-07-01T02h56m38_LVL2.fits',
                    'P18':lvl2path+'FM1T00012308/MIRFM1T00012308_18_495_SE_2011-07-01T03h12m23_LVL2.fits',
                    'P19':lvl2path+'FM1T00012308/MIRFM1T00012308_19_495_SE_2011-07-01T03h28m19_LVL2.fits',
                    'P20':lvl2path+'FM1T00012308/MIRFM1T00012308_20_495_SE_2011-07-01T03h44m04_LVL2.fits',
                    'P21':lvl2path+'FM1T00012308/MIRFM1T00012308_21_495_SE_2011-07-01T03h59m55_LVL2.fits',
                    'P22':lvl2path+'FM1T00012308/MIRFM1T00012308_22_495_SE_2011-07-01T04h15m40_LVL2.fits',
                    'P23':lvl2path+'FM1T00012308/MIRFM1T00012308_23_495_SE_2011-07-01T04h31m30_LVL2.fits',
                    'P24':lvl2path+'FM1T00012308/MIRFM1T00012308_24_495_SE_2011-07-01T04h47m16_LVL2.fits',
                    'P25':lvl2path+'FM1T00012308/MIRFM1T00012308_25_495_SE_2011-07-01T05h03m06_LVL2.fits',
                    'P26':lvl2path+'FM1T00012308/MIRFM1T00012308_26_495_SE_2011-07-01T05h18m52_LVL2.fits',
                    'P27':lvl2path+'FM1T00012308/MIRFM1T00012308_27_495_SE_2011-07-01T05h34m43_LVL2.fits'}
        bkg_file = lvl2path+'FM1T00012203/MIRFM1T00012203_1_495_SE_2011-06-27T22h03m11_LVL2.fits'

    elif position == 'right':
        # Right side of FOV
        sci_imgs = {'P1':lvl2path+'FM1T00012206/MIRFM1T00012206_1_495_SE_2011-06-28T04h51m37_LVL2.fits',
                    'P2':lvl2path+'FM1T00012206/MIRFM1T00012206_2_495_SE_2011-06-28T05h07m28_LVL2.fits',
                    'P3':lvl2path+'FM1T00012206/MIRFM1T00012206_3_495_SE_2011-06-28T05h23m13_LVL2.fits',
                    'P4':lvl2path+'FM1T00012206/MIRFM1T00012206_4_495_SE_2011-06-28T05h39m04_LVL2.fits',
                    'P5':lvl2path+'FM1T00012206/MIRFM1T00012206_5_495_SE_2011-06-28T05h54m49_LVL2.fits',
                    'P6':lvl2path+'FM1T00012206/MIRFM1T00012206_6_495_SE_2011-06-28T06h10m34_LVL2.fits',
                    'P7':lvl2path+'FM1T00012206/MIRFM1T00012206_7_495_SE_2011-06-28T06h26m30_LVL2.fits',
                    'P8':lvl2path+'FM1T00012206/MIRFM1T00012206_8_495_SE_2011-06-28T06h42m15_LVL2.fits',
                    'P9':lvl2path+'FM1T00012206/MIRFM1T00012206_9_495_SE_2011-06-28T06h58m11_LVL2.fits',
                    'P10':lvl2path+'FM1T00012206/MIRFM1T00012206_10_495_SE_2011-06-28T07h13m56_LVL2.fits',
                    'P11':lvl2path+'FM1T00012206/MIRFM1T00012206_11_495_SE_2011-06-28T07h29m41_LVL2.fits',
                    'P12':lvl2path+'FM1T00012206/MIRFM1T00012206_12_495_SE_2011-06-28T07h45m32_LVL2.fits',
                    'P13':lvl2path+'FM1T00012206/MIRFM1T00012206_13_495_SE_2011-06-28T08h01m17_LVL2.fits',
                    'P14':lvl2path+'FM1T00012206/MIRFM1T00012206_14_495_SE_2011-06-28T08h17m03_LVL2.fits',
                    'P15':lvl2path+'FM1T00012206/MIRFM1T00012206_15_495_SE_2011-06-28T08h32m58_LVL2.fits',
                    'P16':lvl2path+'FM1T00012206/MIRFM1T00012206_16_495_SE_2011-06-28T08h48m44_LVL2.fits',
                    'P17':lvl2path+'FM1T00012206/MIRFM1T00012206_17_495_SE_2011-06-28T09h04m34_LVL2.fits',
                    'P18':lvl2path+'FM1T00012206/MIRFM1T00012206_18_495_SE_2011-06-28T09h20m24_LVL2.fits',
                    'P19':lvl2path+'FM1T00012206/MIRFM1T00012206_19_495_SE_2011-06-28T09h36m10_LVL2.fits',
                    'P20':lvl2path+'FM1T00012206/MIRFM1T00012206_20_495_SE_2011-06-28T09h52m00_LVL2.fits',
                    'P21':lvl2path+'FM1T00012206/MIRFM1T00012206_21_495_SE_2011-06-28T10h07m46_LVL2.fits',
                    'P22':lvl2path+'FM1T00012206/MIRFM1T00012206_22_495_SE_2011-06-28T10h23m36_LVL2.fits',
                    'P23':lvl2path+'FM1T00012206/MIRFM1T00012206_23_495_SE_2011-06-28T10h39m27_LVL2.fits'}
        bkg_file = lvl2path+'FM1T00012203/MIRFM1T00012203_1_495_SE_2011-06-27T22h03m11_LVL2.fits'
    if pointing == 'all':
        return sci_imgs,bkg_file
    elif (pointing != 'all') & (output=='img'):
        from astropy.io import fits
        return fits.open(sci_imgs[pointing])[0].data[0,:,:],fits.open(bkg_file)[0].data[0,:,:]

#-- CV2 & CV3 data
def CV_800K_BB_point_source_raster(dataDir,campaign=None):
    """ load MIRIM PSFs"""
    import os
    import glob

    files = [os.path.basename(i) for i in glob.glob(dataDir+'*')]
    subchannels = ['SHORT','MED','LONG']
    MIRIMPSF_dictionary = {}
    if campaign == 'CV2':
        pointings = ['P'+str(i) for i in range(17)]
        for pointing in pointings:
            mylist = []
            for subchannel in subchannels:
                sub = 'MIRM0363-{}-{}'.format(pointing,subchannel)
                mylist.extend([s for s in files if sub in s])
            MIRIMPSF_dictionary[campaign+'_'+pointing] = list(mylist)

    if campaign == 'CV3':
        pointings = ['Q'+str(i) for i in range(17)]
        for pointing in pointings:
            sub = 'MIRM103-{}-SHORT'.format(pointing)
            MIRIMPSF_dictionary[campaign+'_'+pointing] = [s for s in files if sub in s]

    """
    Dictionary keys are equivalent to PSF measurements in CV2 and CV3 tests (with different pointings).
    Dictionary indeces within keys, for CV2 obs. are equivalent to:
    [0,1,2,3] : SHORT_494,SHORT_495,SHORTB_494,SHORTB_495
    [4,5,6,7] : MED_494,MED_495,MEDB_494,MEDB_495
    [8,9,10,11] : LONG_494,LONG_495,LONGB_494,LONGB_495

    Dictionary indeces within keys, for CV3 obs. are equivalent to:
    [0,1,2,3] : SHORT_494,SHORT_495,SHORTB_494,SHORTB_495
    There are only SHORT CV3 PSF measurements
    """

    return MIRIMPSF_dictionary

#-- OTIS campaign
def OTIS_ASPA_semiextended_source(lvl2path,band,pointing=None,output='img'):
    if pointing == 'v03':
        sci_imgs = {"1A":lvl2path +'MIRM32313-SS-V03-7249024654_1_495_SE_2017-09-06T02h54m21_LVL2.fits',
                    "1B":lvl2path +'MIRM32313-MM-V03-7249025513_1_495_SE_2017-09-06T03h03m01_LVL2.fits',
                    "1C":lvl2path +'MIRM32313-LL-V03-7249030323_1_495_SE_2017-09-06T03h08m31_LVL2.fits',
                    "2A":lvl2path +'MIRM32313-SS-V03-7249024654_1_495_SE_2017-09-06T02h54m21_LVL2.fits',
                    "2B":lvl2path +'MIRM32313-MM-V03-7249025513_1_495_SE_2017-09-06T03h03m01_LVL2.fits',
                    "2C":lvl2path +'MIRM32313-LL-V03-7249030323_1_495_SE_2017-09-06T03h08m31_LVL2.fits',
                    "3A":lvl2path +'MIRM32313-SS-V03-7249024654_1_494_SE_2017-09-06T02h54m21_LVL2.fits',
                    "3B":lvl2path +'MIRM32313-MM-V03-7249025513_1_494_SE_2017-09-06T03h03m01_LVL2.fits',
                    "3C":lvl2path +'MIRM32313-LL-V03-7249030323_1_494_SE_2017-09-06T03h08m31_LVL2.fits',
                    "4A":lvl2path +'MIRM32313-SS-V03-7249024654_1_494_SE_2017-09-06T02h54m21_LVL2.fits',
                    "4B":lvl2path +'MIRM32313-MM-V03-7249025513_1_494_SE_2017-09-06T03h03m01_LVL2.fits',
                    "4C":lvl2path +'MIRM32313-LL-V03-7249030323_1_494_SE_2017-09-06T03h08m31_LVL2.fits'}

        bkg_imgs = {"1A":lvl2path +'MIRM32313-SS-V03B-7249024337_1_495_SE_2017-09-06T02h49m41_LVL2.fits',
                    "1B":lvl2path +'MIRM32313-MM-V03B-7249025018_1_495_SE_2017-09-06T02h57m31_LVL2.fits',
                    "1C":lvl2path +'MIRM32313-LL-V03B-7249025837_1_495_SE_2017-09-06T03h05m21_LVL2.fits',
                    "2A":lvl2path +'MIRM32313-SS-V03B-7249024337_1_495_SE_2017-09-06T02h49m41_LVL2.fits',
                    "2B":lvl2path +'MIRM32313-MM-V03B-7249025018_1_495_SE_2017-09-06T02h57m31_LVL2.fits',
                    "2C":lvl2path +'MIRM32313-LL-V03B-7249025837_1_495_SE_2017-09-06T03h05m21_LVL2.fits',
                    "3A":lvl2path +'MIRM32313-SS-V03B-7249024337_1_494_SE_2017-09-06T02h49m41_LVL2.fits',
                    "3B":lvl2path +'MIRM32313-MM-V03B-7249025018_1_494_SE_2017-09-06T02h57m31_LVL2.fits',
                    "3C":lvl2path +'MIRM32313-LL-V03B-7249025837_1_494_SE_2017-09-06T03h05m21_LVL2.fits',
                    "4A":lvl2path +'MIRM32313-SS-V03B-7249024337_1_494_SE_2017-09-06T02h49m41_LVL2.fits',
                    "4B":lvl2path +'MIRM32313-MM-V03B-7249025018_1_494_SE_2017-09-06T02h57m31_LVL2.fits',
                    "4C":lvl2path +'MIRM32313-LL-V03B-7249025837_1_494_SE_2017-09-06T03h05m21_LVL2.fits'}
    elif pointing == 'v05':
        sci_imgs = {"1A":lvl2path +'MIRM32313-SS-V05-7249031538_1_495_SE_2017-09-06T03h23m41_LVL2.fits',
                    "1B":lvl2path +'MIRM32313-MM-V05-7249032353_1_495_SE_2017-09-06T03h31m21_LVL2.fits',
                    "1C":lvl2path +'MIRM32313-LL-V05-7249033205_1_495_SE_2017-09-06T03h37m41_LVL2.fits',
                    "2A":lvl2path +'MIRM32313-SS-V05-7249031538_1_495_SE_2017-09-06T03h23m41_LVL2.fits',
                    "2B":lvl2path +'MIRM32313-MM-V05-7249032353_1_495_SE_2017-09-06T03h31m21_LVL2.fits',
                    "2C":lvl2path +'MIRM32313-LL-V05-7249033205_1_495_SE_2017-09-06T03h37m41_LVL2.fits',
                    "3A":lvl2path +'MIRM32313-SS-V05-7249031538_1_494_SE_2017-09-06T03h23m41_LVL2.fits',
                    "3B":lvl2path +'MIRM32313-MM-V05-7249032353_1_494_SE_2017-09-06T03h31m21_LVL2.fits',
                    "3C":lvl2path +'MIRM32313-LL-V05-7249033205_1_494_SE_2017-09-06T03h37m41_LVL2.fits',
                    "4A":lvl2path +'MIRM32313-SS-V05-7249031538_1_494_SE_2017-09-06T03h23m41_LVL2.fits',
                    "4B":lvl2path +'MIRM32313-MM-V05-7249032353_1_494_SE_2017-09-06T03h31m21_LVL2.fits',
                    "4C":lvl2path +'MIRM32313-LL-V05-7249033205_1_494_SE_2017-09-06T03h37m41_LVL2.fits'}

        bkg_imgs = {"1A":lvl2path +'MIRM32313-SS-V05B-7249031039_1_495_SE_2017-09-06T03h18m01_LVL2.fits',
                    "1B":lvl2path +'MIRM32313-MM-V05B-7249031906_1_495_SE_2017-09-06T03h26m01_LVL2.fits',
                    "1C":lvl2path +'MIRM32313-LL-V05B-7249032716_1_495_SE_2017-09-06T03h34m51_LVL2.fits',
                    "2A":lvl2path +'MIRM32313-SS-V05B-7249031039_1_495_SE_2017-09-06T03h18m01_LVL2.fits',
                    "2B":lvl2path +'MIRM32313-MM-V05B-7249031906_1_495_SE_2017-09-06T03h26m01_LVL2.fits',
                    "2C":lvl2path +'MIRM32313-LL-V05B-7249032716_1_495_SE_2017-09-06T03h34m51_LVL2.fits',
                    "3A":lvl2path +'MIRM32313-SS-V05B-7249031039_1_494_SE_2017-09-06T03h18m01_LVL2.fits',
                    "3B":lvl2path +'MIRM32313-MM-V05B-7249031906_1_494_SE_2017-09-06T03h26m01_LVL2.fits',
                    "3C":lvl2path +'MIRM32313-LL-V05B-7249032716_1_494_SE_2017-09-06T03h34m51_LVL2.fits',
                    "4A":lvl2path +'MIRM32313-SS-V05B-7249031039_1_494_SE_2017-09-06T03h18m01_LVL2.fits',
                    "4B":lvl2path +'MIRM32313-MM-V05B-7249031906_1_494_SE_2017-09-06T03h26m01_LVL2.fits',
                    "4C":lvl2path +'MIRM32313-LL-V05B-7249032716_1_494_SE_2017-09-06T03h34m51_LVL2.fits'}
    if output == 'filename':
        return sci_imgs[band],bkg_imgs[band]
    elif output == 'img':
        from astropy.io import fits
        return fits.open(sci_imgs[band])[0].data[0,:,:],fits.open(bkg_imgs[band])[0].data[0,:,:]
