from astropy.io import fits
from astropy.table import Table
import numpy as np
import glob
import os

# set memmap=True for large files
fits_files = sorted(glob.glob('/storage/home/nxt5197/work/530_stellar_atmospheres/project/NLTT43564/*.fits'))
for f in fits_files:
    # filename = raw_input(f)
    # output = raw_input(os.path.splitext(os.path.basename(f))[0] + '.csv')

    # # Open the given fits file
    # hdulist = fits.open(filename)
    # scidata = hdulist[0].data

    # # save your new file as a csv
    # numpy.savetxt(output, scidata, fmt='%d', delimiter=',')
 
    with fits.open(f, memmap=True):

        # select the HDU you want
        hdu_list = fits.open(f)
        hdul1=hdu_list[1].header['EXTNAME']
        hdul2=hdu_list[2].header['EXTNAME']
        hdul3=hdu_list[3].header['EXTNAME']
        hdul4=hdu_list[4].header['EXTNAME']
        hdul5=hdu_list[5].header['EXTNAME']
        hdul6=hdu_list[6].header['EXTNAME']
        hdul7=hdu_list[7].header['EXTNAME']
        hdul8=hdu_list[8].header['EXTNAME']
        hdul9=hdu_list[9].header['EXTNAME']
        # hdu_flux = hdul[1].data[26]
        # hdu_wavelength = hdul[7].data[26]
        # hdu_var = hdul[4].data[26]

        #flux values in science, sky, and calibrator data
        #each one contains 28 lists of flux values corresponding to the 28 spectral orders
        sci=hdu_list[1].data
        sky=hdu_list[2].data
        cal=hdu_list[3].data
        scivar=hdu_list[4].data
        skyvar=hdu_list[5].data
        calvar=hdu_list[6].data

        #x values for science and sky fiber data
        scix=hdu_list[7].data
        skyx=hdu_list[8].data
        calx=hdu_list[9].data

        #convert NaN pixels to 0
        nans=np.isnan(sky)
        sky[nans]=0

        nans=np.isnan(sci)
        sci[nans]=0

        dt=Table()
        dt[hdul1]=sci
        dt[hdul2]=sky
        dt[hdul3]=cal
        dt[hdul4]=scivar
        dt[hdul5]=skyvar
        dt[hdul6]=calvar
        dt[hdul7]=scix
        dt[hdul8]=skyx
        dt[hdul9]=calx

        # write to a CSV file
        filname = os.path.splitext(os.path.basename(f))[0] + '.csv'
        dt.write(filname, delimiter='\t', format='ascii', overwrite=True)

#         import numpy
# from astropy.io import fits

# # ask for an input file 
# filename = raw_input('Enter a FITS file to CSV-ize: ')

# # ask for an output file name
# output = raw_input('What would you like to call your new file?: ')

# # Open the given fits file
# hdulist = fits.open(filename)
# scidata = hdulist[0].data

# # save your new file as a csv
# numpy.savetxt(output, scidata, fmt='%d', delimiter=',')