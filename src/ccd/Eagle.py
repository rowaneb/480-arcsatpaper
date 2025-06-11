
from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np
from astropy.coordinates import get_body
from astropy.visualization import ImageNormalize, LinearStretch, ZScaleInterval
import glob
from astropy.stats import sigma_clip
import os
import pathlib
from astropy.table import Table
import matplotlib.pyplot as plt
from pathlib import Path

#########

#Creating a Bias File:

#outlining directory path
input_dir = Path("data")
output_file = "eagle_bias.fits"

#identifying bias images
bias_data = []

for file in sorted(input_dir.glob("Bias*.fits")):
    data = fits.getdata(file).astype('f4')

    bias_data.append(data)

#converting bias into an array 
bias_3d = np.array(bias_data)
clipping = sigma_clip(bias_3d, cenfunc='median', sigma=3, axis=0)

#taking the median
median_bias = np.ma.median(clipping, axis=0)

#Plot Formatting and Creation
plot_data = median_bias.filled(np.nan)
vmin = np.nanpercentile(data_for_plot, 5)
vmax = np.nanpercentile(data_for_plot, 95)

plt.imshow(plot_data, origin='lower', cmap = plt.cm.plasma, vmin=vmin, vmax=vmax)
plt.savefig(Path("Eagle_Bias.png"), dpi=150)
plt.show()
    
#Saving the final image
output_path = Path("Eagle Images") / output_file
primary = fits.PrimaryHDU(data=data_for_plot, header=fits.Header())
hdul = fits.HDUList([primary])
hdul.writeto(output_path, overwrite=True)

#########

#Creating a Dark File:

#Referencing bias frame
bias_filename= "work/480-arcsatpaper/eagle_bias.fits"
bias_frame = fits.getdata("eagle_bias.fits").astype('f4')

#identifying dark images
dark_data = []
    
Calling the data from dark_list
for file in sorted(input_dir.glob("Dark*.fits")):
        data = fits.getdata(file).astype('f4')
        
        #Setting exposure time to what was used during observation
        exptime = 300
        
        #Subtracting bias from dark frame
        dark_sub = data - bias_frame
        
        #divinding the bias by exp time, becomes 2d array and is appended into dark_bias_data list
        dark_data.append(dark_sub / exptime)

#turning dark back into array, sigma clipping
dark_3d = np.array(dark_data)
clipping = sigma_clip(dark_3d, cenfunc='median', sigma=3, axis=0)

#taking the median
median_dark = np.ma.median(clipping, axis=0)

#Plot Formatting and Creation
data_for_plot = median_dark.filled(np.nan)
vmin = np.nanpercentile(data_for_plot, 5)
vmax = np.nanpercentile(data_for_plot, 95)

#making the actual plot into a png
plt.imshow(data_for_plot, origin='lower', cmap = plt.cm.plasma, vmin=vmin, vmax=vmax)
plt.savefig(Path("Eagle_Dark.png"), dpi=150)
plt.show()

#Saving the final image
output_file= "eagle_dark.fits"
output_path = Path("Eagle Images") / output_file
primary = fits.PrimaryHDU(data=data_for_plot, header=fits.Header())
hdul = fits.HDUList([primary])
hdul.writeto(output_path, overwrite=True)

#########

#Creating a Flat File:

#Referencing bias frame
bias_frame = fits.getdata("eagle_bias.fits").astype('f4')

#identifying flat images
flat_data = []

for file in sorted(input_dir.glob("domeflat*.fits")):
    data = fits.getdata(file).astype('f4')

    #Subtracting bias from flat frame
    sub_bias = data - bias_frame

    flat_data.append(sub_bias)

#turning back into array, sigma clipping
array = np.array(flat_data)
clipping = sigma_clip(array, cenfunc='median', sigma=3, axis=0)

#taking median
me_flat = np.ma.median(clipping, axis=0)

#normalizing
median_flat = me_flat / np.ma.median(me_flat)

#Plot Formatting and Creation
data_for_plot = median_flat.filled(np.nan)
vmin = np.nanpercentile(data_for_plot, 5)
vmax = np.nanpercentile(data_for_plot, 95)

#making the actual plot
plt.imshow(data_for_plot, origin='lower', cmap = plt.cm.plasma, vmin=vmin, vmax=vmax)
plt.savefig(Path("Eagle_Flat.png"), dpi=150)
plt.show()

#Saving the final image
output_file= "eagle_flat.fits"
output_path = Path("Eagle Images") / output_file
primary = fits.PrimaryHDU(data=data_for_plot, header=fits.Header())
hdul = fits.HDUList([primary])
hdul.writeto(output_path, overwrite=True)

#########

#Creating a Science File:

#referencing the previous reduction files
bias_frame = fits.getdata("eagle_bias.fits").astype('f4')
flat_frame = fits.getdata("eagle_flat.fits").astype('f4')
dark_frame = fits.getdata("eagle_dark.fits").astype('f4')
#referencing the image I will be reducing
sci_data = fits.getdata("EAGLE_H-Alpha_20250603_102159.fits").astype('f4')

#setting the exposure time to that used during observation
exptime = 600

#creating an array, sigma clipping
array = np.array(sci_data)
sci_clip = sigma_clip(array, cenfunc='median', sigma=3, axis=0)
  
#subtracting bias from science frame
sub_bias = sci_clip - bias_frame

#subtracting dark frame, accounting for the exposure time
dark_corrected = sub_bias - dark_frame * exptime

#normalizing (flat)
flat_norm = flat_frame / np.mean(flat_frame)

#normalzing (science)
final_sci = dark_corrected / flat_norm
    
#Plot Formatting and Creation
data_for_plot = final_sci.filled(np.nan)
vmin = np.nanpercentile(data_for_plot, 5)
vmax = np.nanpercentile(data_for_plot, 95)

#making the actual plot
plt.imshow(data_for_plot, origin='lower', cmap = plt.cm.plasma, vmin=vmin, vmax=vmax)
plt.savefig(Path("Eagle_Final.png"), dpi=150)
plt.show()

#Saving the final final final image
output_file= "eagle_science.fits"
output_path = Path("Eagle Images") / output_file
primary = fits.PrimaryHDU(data=data_for_plot, header=fits.Header())
hdul = fits.HDUList([primary])
hdul.writeto(output_path, overwrite=True)