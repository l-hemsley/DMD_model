import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from functions import *
import time

## this model is used firstly to calculate the diffraction pattern from a DMD, and secondly to calculate the transmission efficiency with wavelength when collected with a lens.
#The modelled configuration assumes point source at infinity -> focussing lens -> DMD -> collection lens -> view diffraction pattern/calculate collection efficiency

#useful constants
mm=10**-3
um=10**-6
nm=10**-9

#initialize system using the paramters from the prototype as an example

#DMD_parameters(pitch, fill_factor, tilt_angle, mirror axis tilt direction)
DMD=DMD_parameters(10.8*um,0.98,np.radians(12),'diagonal')
a=8.73

#input_parameters(wavelength, axis_angle_x, axis_angle_y,lens,lens_diameter,focal_length):
#use lens='ideal' as default to calculate MTF function for unknown lenses
#axis_angles_x/y describe the direction of the lens axis to the DMD.
#assume the DMD is at the focal plane of the lens

input=input_parameters(600*nm,np.radians(a),np.radians(-a),'TL200',0.05*2*150*mm,150*mm)

#output_parameters(lens_NA, angle_x_centre, angle_y_centre, datapoints, DMD,lens,lens_diameter,focal_length, input)
#datapoints gives number of data points used to caluclate diffraction pattern - more points is slower but more accurate
output=output_parameters(np.radians(a),np.radians(-a),200,DMD,'TL200',0.05*2*150*mm,150*mm,input)

 #ANTOINES SYSTEM - out of date
# DMD=DMD_parameters(5.4*um,0.98,np.radians(17),'vertical')
# input=input_parameters(600*nm,np.radians(0),np.radians(0),'TL165')
# output=output_parameters(0.05,np.radians(34),np.radians(0),400,DMD,input,'TL165')

#wavelength range of interest
wavelengths=np.arange(420*nm,700*nm,20*nm)

transmission_collected=np.zeros((np.size(wavelengths)))
image_collected=np.zeros([output.angle_x_array_meshed.shape[0],output.angle_x_array_meshed.shape[1]])
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)

#The model is run over the wavelength range, producing a diffraction pattern for each wavelength
for i in np.arange(np.size(wavelengths)):
     ax1.cla()
     ax2.cla()
     ax3.cla()
     ax4.cla()

     t = time.time()
     input.wavelength=wavelengths[i]

     #calculation diffraction image
     [diffraction_image,total_power_collected,E2_grating,E2_envelope,image_collected]=calculate_diffraction_pattern_image(input, output, DMD)
     #image_collected=image_collected+image_collected2*wavelengths[i]
     elapsed = time.time() - t
     print('calculation time elapsed = '+str(elapsed))
     t = time.time()

     # TODO - how to correctly normalize the tranmission?
     transmission_collected[i] = total_power_collected / np.sum(np.sum(diffraction_image))
     print('Wavelength =' + str(np.round(input.wavelength/nm,0))+'nm')

     #plot results

     im1=ax1.contourf(np.degrees(output.angle_x_array_meshed),np.degrees(output.angle_y_array_meshed),E2_grating,50)
     im2=ax2.contourf(np.degrees(output.angle_x_array_meshed),np.degrees(output.angle_y_array_meshed),E2_envelope,50)
     im3=ax3.contourf(np.degrees(output.angle_x_array_meshed),np.degrees(output.angle_y_array_meshed),diffraction_image,50)
     im4=ax4.contourf(np.degrees(output.angle_x_array_meshed),np.degrees(output.angle_y_array_meshed),image_collected,50)
     ax1.set_title('Grating Orders')
     ax2.set_title('Mirror Envelope')
     ax3.set_title('Combined Diffraction Pattern')
     ax4.set_title('Collected Image (w/ MTF)')
     ax1.set_xlabel('output angle x degrees')
     ax1.set_ylabel('output angle y degrees')

     fig.tight_layout()

     fig.suptitle('Wavelength =' +str(np.round(input.wavelength/nm,0))+'nm', fontsize=16)
     fig.savefig('Figures/wavelength' +str(np.round(input.wavelength/nm,0))+'nm' + '.png')

     plt.pause(0.000001)
     elapsed = time.time() - t
     print('plotting time elapsed = '+str(elapsed))

#plt.show()
#import some experimental data for comparison

experimental_data=pd.read_excel(r'experimental.xlsx')
wavelengths_experimental=experimental_data.loc[:,'Wavelength'].to_numpy()
transmission_experimental=experimental_data.loc[:,'Transmission'].to_numpy()

experimental_data = pd.read_excel(r'TL200-transmission.xlsx')
wavelengths_lens = experimental_data.loc[:, 'Wavelength'].to_numpy()
transmission_lens= experimental_data.loc[:, 'Transmission'].to_numpy()/100
transmission_lens=np.interp(wavelengths/nm,wavelengths_lens,transmission_lens)

fig=plt.figure(2)
plt.plot(wavelengths/nm,transmission_collected*transmission_lens**2)
plt.plot(wavelengths/nm,transmission_collected)
plt.plot(wavelengths_experimental,transmission_experimental)
plt.grid()
plt.ylim((0,1.2))
plt.xlim((420,700))
plt.xlabel('Wavelengths (nm)')
plt.ylabel('Transmission')
fig.savefig('Figures/Transmission.png')
plt.show()
#
