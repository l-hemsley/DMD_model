
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from functions import *

#useful constants
um=10**-6
nm=10**-9

dmd=DMD_parameters('Standard',10.8*um,0.95,np.radians(-12))
input=input_parameters(600*nm,0.05,np.radians(8.54),np.radians(-8.54),20*um)
output=output_parameters(0.05,np.radians(8.54),np.radians(-8.54),500)
envelope_function(input, output, dmd)

# find orders for DMD and input
wavelengths=np.arange(400*nm,700*nm,10*nm)
#wavelengths=np.array([400*nm])
transmission=np.zeros((np.size(wavelengths)))
transmission_collected=np.zeros((np.size(wavelengths)))

fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
#im2=ax1.imshow(np.zeros((1,1)))
#fig.colorbar(im2, ax=ax2)
# fig.colorbar(im2, ax=ax2)
# fig.colorbar(im3, ax=ax3)

for i in np.arange(np.size(wavelengths)):
     input.wavelength=wavelengths[i]
     E_grating=grating_function(input, output, dmd)**2
     E_envelope=envelope_function(input, output, dmd)**2
     E_total=E_grating*E_envelope
     transmission[i]=np.sum(np.sum(E_total))
     collected_vectors = abs(output.effective_NA_of_vector) < output.half_angle

     E_total_collected=E_total*collected_vectors
     transmission_collected[i] = np.sum(np.sum(E_total_collected))

im1=ax1.imshow(E_grating,extent=(output.angle_x_array_deg[0],output.angle_x_array_deg[-1],output.angle_y_array_deg[0],output.angle_y_array_deg[-1]))

im2=ax2.imshow(E_envelope,extent=(output.angle_x_array_deg[0],output.angle_x_array_deg[-1],output.angle_y_array_deg[0],output.angle_y_array_deg[-1]))


im3=ax3.imshow(E_total,extent=(output.angle_x_array_deg[0],output.angle_x_array_deg[-1],output.angle_y_array_deg[0],output.angle_y_array_deg[-1]))

ax4.clear()

#im4=ax4.plot(wavelengths/nm,transmission/np.max(transmission))
im4=ax4.plot(wavelengths/nm,transmission_collected/np.max(transmission_collected))

experimental_data=pd.read_excel(r'experimental.xlsx')
wavelengths_experimental=experimental_data.loc[:,'Wavelength'].to_numpy()
transmission_experimental=experimental_data.loc[:,'Transmission'].to_numpy()

im4=ax4.plot(wavelengths_experimental,transmission_experimental)
ax4.set_xlim((400,700))
ax4.set_ylim((0,1.2))
ax1.text(output.angle_x_array_deg[0],output.angle_y_array_deg[0],'wavelength='+ str(int(input.wavelength/nm))+ 'nm, beam size='+str(input.effective_beam_size/um)+'um',color='w')

#
ax1.set_title('Diffraction Orders, input ('+ str(np.degrees(input.angle_x_centre)) +','+ str(np.degrees(input.angle_y_centre))+')')
ax2.set_title('Envelope, input ('+ str(np.degrees(input.angle_x_centre)) +','+ str(np.degrees(input.angle_y_centre))+')')
ax3.set_title('Combined, input ('+ str(np.degrees(input.angle_x_centre)) +','+ str(np.degrees(input.angle_y_centre))+')')
#plt.pause(0.01)

print('done')
plt.show()
