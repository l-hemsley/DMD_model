
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from functions import *


#useful constants
mm=10**-3
um=10**-6
nm=10**-9

#initialize system
dmd=DMD_parameters('Standard',10.8*um,0.99,np.radians(12))
input=input_parameters(600*nm,0.05,np.radians(8.54),np.radians(-8.54),150*mm)
output=output_parameters(0.05,np.radians(8.54),np.radians(-8.54),50)
wavelengths=np.arange(400*nm,750*nm,5*nm)
#wavelengths=np.arange(400*nm,720*nm,20*nm)
transmission_collected_integrated=np.zeros((np.size(wavelengths)))
transmission_collected=np.zeros((np.size(wavelengths)))


#simulation run over wavelengths
# - compare integrated over input beam VS single beam model

for i in np.arange(np.size(wavelengths)):
     input.wavelength=wavelengths[i]
     [diffraction_image,total_power_collected,E2_grating,E2_envelope]=calculate_diffraction_pattern_image(input, output, dmd)
     transmission_collected[i] = total_power_collected#*wavelengths[i]

     #[diffraction_image_int,total_power_collected_integrated,E2_grating_int,E2_envelope_int]=diff_image_integrated_input_NA(input, output, dmd, 20) # 20-50 ok
     #transmission_collected_integrated[i] = total_power_collected_integrated
     print('Wavelength =' +str(input.wavelength/nm)+'nm')

#import experimental data
experimental_data=pd.read_excel(r'experimental.xlsx')
wavelengths_experimental=experimental_data.loc[:,'Wavelength'].to_numpy()
transmission_experimental=experimental_data.loc[:,'Transmission'].to_numpy()

#plot results
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
im1=ax1.contourf(np.degrees(output.angle_x_array_meshed),np.degrees(output.angle_y_array_meshed),E2_grating,500)
im2=ax2.contourf(np.degrees(output.angle_x_array_meshed),np.degrees(output.angle_y_array_meshed),E2_envelope,500)
im3=ax3.contourf(np.degrees(output.angle_x_array_meshed),np.degrees(output.angle_y_array_meshed),diffraction_image,500)

ax4.plot(wavelengths/nm,transmission_collected/max(transmission_collected))
ax4.plot(wavelengths_experimental,transmission_experimental)
ax4.set_xlim((400,750))
ax4.set_ylim((0,1.2))

fig.colorbar(im1, ax=ax1)
fig.colorbar(im2, ax=ax2)
fig.colorbar(im3, ax=ax3)

ax1.set_title('Grating Orders')
ax2.set_title('Mirror Envelope')
ax3.set_title('Combined Diffraction Pattern')
#fig.savefig('Figure_single.png')
plt.show()

# fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
# im1=ax1.contourf(np.degrees(output.angle_x_array_meshed),np.degrees(output.angle_y_array_meshed),E2_grating_int,100)
# im2=ax2.contourf(np.degrees(output.angle_x_array_meshed),np.degrees(output.angle_y_array_meshed),E2_envelope_int,100)
# im3=ax3.contourf(np.degrees(output.angle_x_array_meshed),np.degrees(output.angle_y_array_meshed),diffraction_image_int,100)
# ax4.plot(wavelengths/nm,transmission_collected_integrated/max(transmission_collected_integrated))
# ax4.plot(wavelengths_experimental,transmission_experimental)
# ax4.set_xlim((420,700))
# ax4.set_ylim((0,1.2))
#
# fig.colorbar(im1, ax=ax1)
# fig.colorbar(im2, ax=ax2)
# fig.colorbar(im3, ax=ax3)
#
# ax1.set_title('Grating Orders')
# ax2.set_title('Mirror Envelope')
# ax3.set_title('Combined Diffraction Pattern')
# fig.savefig('Figure_integrated.png')
#
# fig=plt.figure()
# plt.plot(wavelengths/nm,transmission_collected/max(transmission_collected))
# plt.plot(wavelengths/nm,transmission_collected_integrated/max(transmission_collected_integrated))
# plt.plot(wavelengths_experimental,transmission_experimental)
# plt.xlim((400,750))
# plt.ylim((0,1.2))
# plt.title('Transmission Function')
# fig.savefig('Transmission_Plot.png')
# plt.show()

