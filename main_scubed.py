import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from functions import *

import time

#useful constants
mm=10**-3
um=10**-6
nm=10**-9

#initialize system using the paramters from the prototype

#DMD_parameters(pitch, fill_factor, tilt_angle, mirror axis tilt direction)
DMD=DMD_parameters(7.6*um,0.98,np.radians(12),'diagonal')
#input_parameters(wavelength, angle_x_centre, angle_y_centre)
a=8.54
input=input_parameters(700*nm,0,np.radians(-a),'ideal',22*mm,50*mm)
# output_parameters(lens_NA, angle_x_centre, angle_y_centre, datapoints, DMD, input)
output=output_parameters(np.radians(2*a),np.radians(-a),200,DMD,'ideal',44*mm,100*mm,input)


#wavelength range of interest
wavelengths=np.arange(700*nm,850*nm,10*nm)
#wavelengths=[732*nm]
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
     fig.savefig('Figures/scubed_wavelength' +str(np.round(input.wavelength/nm,0))+'nm' + '.png')

     plt.pause(0.000001)
     elapsed = time.time() - t
     print('plotting time elapsed = '+str(elapsed))

#plt.show()
#import some experimental data for comparison


fig=plt.figure(2)
plt.plot(wavelengths/nm,transmission_collected)
plt.grid()
plt.ylim((0,1.2))
plt.xlim((700,850))
plt.xlabel('Wavelengths (nm)')
plt.ylabel('Transmission')
fig.savefig('Figures/scubed_Transmission.png')
plt.show()
#
