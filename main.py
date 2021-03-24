import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from functions import *
<<<<<<< HEAD
=======
import time
>>>>>>> 2f6e068ce994767c094452733a178d539cedcccd

#useful constants
mm=10**-3
um=10**-6
nm=10**-9

#initialize system
dmd=DMD_parameters('Standard',10.8*um,0.98,np.radians(12))
input=input_parameters(600*nm,0.05,np.radians(8.54),np.radians(-8.54),150*mm)
output=output_parameters(0.05,np.radians(8.54),np.radians(-8.54),100)
wavelengths=np.arange(400*nm,750*nm,5*nm)
transmission_collected=np.zeros((np.size(wavelengths)))
<<<<<<< HEAD
transmission_collected_triangle=np.zeros((np.size(wavelengths)))
=======

#simulation run over wavelengths

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
>>>>>>> 2f6e068ce994767c094452733a178d539cedcccd

for i in np.arange(np.size(wavelengths)):

     t = time.time()
     input.wavelength=wavelengths[i]
<<<<<<< HEAD
     [diffraction_image,total_power_collected,total_power_collected_triangle,E2_grating,E2_envelope]=calculate_diffraction_pattern_image(input, output, dmd)
     transmission_collected[i] = total_power_collected#*wavelengths[i]
     transmission_collected_triangle[i]=total_power_collected_triangle#*wavelengths[i]
     print('Wavelength =' +str(input.wavelength/nm)+'nm')
=======
     [diffraction_image,total_power_collected,E2_grating,E2_envelope,image_collected]=calculate_diffraction_pattern_image(input, output, dmd)
     transmission_collected[i] = total_power_collected#
     print('Wavelength =' + str(np.round(input.wavelength/nm,0))+'nm')

     #plot results

     im1=ax1.contourf(np.degrees(output.angle_x_array_meshed),np.degrees(output.angle_y_array_meshed),E2_grating,100)
     im2=ax2.contourf(np.degrees(output.angle_x_array_meshed),np.degrees(output.angle_y_array_meshed),E2_envelope,100)
     im3=ax3.contourf(np.degrees(output.angle_x_array_meshed),np.degrees(output.angle_y_array_meshed),diffraction_image,100)
     im4=ax4.contourf(np.degrees(output.angle_x_array_meshed),np.degrees(output.angle_y_array_meshed),image_collected,100)

     fig.suptitle('Wavelength =' +str(np.round(input.wavelength/nm,0))+'nm', fontsize=16)
     fig.savefig('Figures/wavelength' +str(np.round(input.wavelength/nm,0))+'nm' + '.png')
     ax1.set_title('Grating Orders')
     ax2.set_title('Mirror Envelope')
     ax3.set_title('Combined Diffraction Pattern')
     ax4.set_title('Collected Image')
     if i==0:
          cbar1=plt.colorbar(im1, ax=ax1)
          cbar2 = plt.colorbar(im2, ax=ax2)
          cbar3 = plt.colorbar(im3, ax=ax3)
          cbar4 = plt.colorbar(im4, ax=ax4)
     else:
           cbar1.update_normal(im1)
           cbar2.update_normal(im2)
           cbar3.update_normal(im3)
           cbar4.update_normal(im4)

     fig.tight_layout()
     plt.pause(0.0001)
     elapsed = time.time() - t
     print('elapsed = '+str(elapsed))
     ax1.cla()
     ax2.cla()
     ax3.cla()
     ax4.cla()
>>>>>>> 2f6e068ce994767c094452733a178d539cedcccd

#import experimental data

experimental_data=pd.read_excel(r'experimental.xlsx')
wavelengths_experimental=experimental_data.loc[:,'Wavelength'].to_numpy()
transmission_experimental=experimental_data.loc[:,'Transmission'].to_numpy()

<<<<<<< HEAD
#plot results
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
im1=ax1.contourf(np.degrees(output.angle_x_array_meshed),np.degrees(output.angle_y_array_meshed),E2_grating,500)
im2=ax2.contourf(np.degrees(output.angle_x_array_meshed),np.degrees(output.angle_y_array_meshed),E2_envelope,500)
im3=ax3.contourf(np.degrees(output.angle_x_array_meshed),np.degrees(output.angle_y_array_meshed),diffraction_image,500)

ax4.plot(wavelengths/nm,transmission_collected/max(transmission_collected))
ax4.plot(wavelengths/nm,transmission_collected_triangle/max(transmission_collected_triangle))
ax4.plot(wavelengths_experimental,transmission_experimental)
ax4.set_xlim((400,750))
ax4.set_ylim((0,1.2))
ax4.legend(['T','T (triangle)','exp.'])

ax1.set_title('Grating Orders')
ax2.set_title('Mirror Envelope')
ax3.set_title('Combined Diffraction Pattern')
fig.savefig('Figure_single.png')
plt.show()

=======
plt.clf()
plt.plot(wavelengths/nm,transmission_collected/max(transmission_collected))
plt.plot(wavelengths_experimental,transmission_experimental)
plt.ylim((0,1.2))
plt.xlim((400,700))

fig.savefig('Figures/Transmission.png')
plt.show()

>>>>>>> 2f6e068ce994767c094452733a178d539cedcccd
