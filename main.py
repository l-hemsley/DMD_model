
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from functions import *

#useful constants
um=10**-6
nm=10**-9

#initialize system
dmd=DMD_parameters('Standard',10.8*um,0.99,np.radians(12))
input=input_parameters(600*nm,0.05,np.radians(8.54),np.radians(-8.54),10*um)
output=output_parameters(0.045,np.radians(8.54),np.radians(-8.54),200)
wavelengths=np.arange(400*nm,750*nm,10*nm)
transmission_collected=np.zeros((np.size(wavelengths)))


#simulation run over wavelengths

for i in np.arange(np.size(wavelengths)):
     input.wavelength=wavelengths[i]
     #[diffraction_image,total_power_collected,E2_grating,E2_envelope]=calculate_diffraction_pattern_image(input, output, dmd)
     [diffraction_image,total_power_collected,E2_grating,E2_envelope]=diff_image_integrated_input_NA(input, output, dmd, 10)
     transmission_collected[i] = total_power_collected
     print('Wavelength =' +str(input.wavelength)+'nm')

#plot results
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
im1=ax1.imshow(E2_grating,extent=(output.angle_x_array_deg[0],output.angle_x_array_deg[-1],output.angle_y_array_deg[0],output.angle_y_array_deg[-1]))
im2=ax2.imshow(E2_envelope,extent=(output.angle_x_array_deg[0],output.angle_x_array_deg[-1],output.angle_y_array_deg[0],output.angle_y_array_deg[-1]))
im3=ax3.imshow(diffraction_image,extent=(output.angle_x_array_deg[0],output.angle_x_array_deg[-1],output.angle_y_array_deg[0],output.angle_y_array_deg[-1]))

fig.colorbar(im1, ax=ax1)
fig.colorbar(im2, ax=ax2)
fig.colorbar(im3, ax=ax3)

ax1.set_title('Grating Orders')
ax2.set_title('Mirror Envelope')
ax3.set_title('Combined Diffraction Pattern')

#import experimental data
experimental_data=pd.read_excel(r'experimental.xlsx')
wavelengths_experimental=experimental_data.loc[:,'Wavelength'].to_numpy()
transmission_experimental=experimental_data.loc[:,'Transmission'].to_numpy()

im4=ax4.plot(wavelengths/nm,transmission_collected/np.max(transmission_collected))
im4=ax4.plot(wavelengths_experimental,transmission_experimental)
ax4.set_xlim((400,750))
ax4.set_ylim((0,1.2))
ax4.set_title('Transmission Function')

plt.show()
