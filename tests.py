import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from functions import *


#useful constants
um=10**-6
nm=10**-9

#initialize system
dmd=DMD_parameters('Standard',10.8*um,0.99,np.radians(12))
input=input_parameters(600*nm,0.05,np.radians(8.54),np.radians(-8.54),15*um)
wavelengths=np.arange(420*nm,700*nm,5*nm)
transmission_collected_integrated=np.zeros((np.size(wavelengths)))
transmission_collected=np.zeros((np.size(wavelengths)))

#import experimental data
experimental_data=pd.read_excel(r'experimental.xlsx')
wavelengths_experimental=experimental_data.loc[:,'Wavelength'].to_numpy()
transmission_experimental=experimental_data.loc[:,'Transmission'].to_numpy()

fig=plt.figure()
output=output_parameters(0.04,np.radians(8.54),np.radians(-8.54),100)
for i in np.arange(np.size(wavelengths)):
     input.wavelength=wavelengths[i]
     [diffraction_image,total_power_collected,E2_grating,E2_envelope]=calculate_diffraction_pattern_image(input, output, dmd)
     transmission_collected[i] = total_power_collected*wavelengths[i]

plt.plot(wavelengths/nm,transmission_collected/max(transmission_collected))
plt.plot(wavelengths_experimental,transmission_experimental-0.1)
plt.xlim((420,720))
plt.ylim((0,1.2))
plt.title('Transmission Function - multiplied by wavelength')
plt.show()
fig.savefig('NA.png')

fig=plt.figure()

# #### run over NA
# NAs=np.arange(0.03,0.05,0.005)
#
# for x in np.arange(np.size(NAs)):
#
#     output=output_parameters(NAs[x],np.radians(8.54),np.radians(-8.54),100)
#     for i in np.arange(np.size(wavelengths)):
#          input.wavelength=wavelengths[i]
#          [diffraction_image,total_power_collected,E2_grating,E2_envelope]=calculate_diffraction_pattern_image(input, output, dmd)
#          transmission_collected[i] = total_power_collected
#
#     plt.plot(wavelengths/nm,transmission_collected/max(transmission_collected))
#
#
# plt.legend(np.around(NAs,3))
# plt.plot(wavelengths_experimental,transmission_experimental)
# plt.xlim((420,720))
# plt.ylim((0,1.2))
# plt.title('Transmission Function - NA')
# plt.show()
# fig.savefig('NA.png')
#
# #### run overbeam width
# effective_beam_width=np.arange(5,30,5)
#
# output=output_parameters(0.04,np.radians(8.54),np.radians(-8.54),100)
# for x in np.arange(np.size(effective_beam_width)):
#
#     input = input_parameters(600 * nm, 0.05, np.radians(8.54), np.radians(-8.54), effective_beam_width[x] * um)
#     for i in np.arange(np.size(wavelengths)):
#          input.wavelength=wavelengths[i]
#          [diffraction_image,total_power_collected,E2_grating,E2_envelope]=calculate_diffraction_pattern_image(input, output, dmd)
#          transmission_collected[i] = total_power_collected
#
#     plt.plot(wavelengths/nm,transmission_collected/max(transmission_collected))
#
#
# plt.legend(np.around(effective_beam_width,3))
# plt.plot(wavelengths_experimental,transmission_experimental)
# plt.xlim((420,720))
# plt.ylim((0,1.2))
# plt.title('Transmission Function - effective_beam_width')
# plt.show()
# fig.savefig('beam_width.png')
#
# #### run over beam angle
# beam_angle_x=np.arange(8.54-5,8.54+5,1)
# output=output_parameters(0.04,np.radians(8.54),np.radians(-8.54),100)
#
# for x in np.arange(np.size(beam_angle_x)):
#
#     input = input_parameters(600 * nm, 0.05, np.radians(beam_angle_x[x]), np.radians(-8.54), 10 * um)
#     for i in np.arange(np.size(wavelengths)):
#          input.wavelength=wavelengths[i]
#          [diffraction_image,total_power_collected,E2_grating,E2_envelope]=calculate_diffraction_pattern_image(input, output, dmd)
#          transmission_collected[i] = total_power_collected
#
#     plt.plot(wavelengths/nm,transmission_collected)
#
#
# plt.legend(np.around(beam_angle_x,3))
# #plt.plot(wavelengths_experimental,transmission_experimental)
# plt.xlim((420,720))
# #plt.ylim((0,1.2))
# plt.title('Transmission Function - effective_beam_width y')
# plt.show()
# fig.savefig('beam_angle_x.png')
#
# #### run over beam angle y
# beam_angle_y=np.arange(8.54-5,8.54+5,1)
# output=output_parameters(0.04,np.radians(8.54),np.radians(-8.54),100)
#
# for x in np.arange(np.size(beam_angle_y)):
#
#     input = input_parameters(600 * nm, 0.05, np.radians(8.54), np.radians(beam_angle_y[x]), 10 * um)
#     for i in np.arange(np.size(wavelengths)):
#          input.wavelength=wavelengths[i]
#          [diffraction_image,total_power_collected,E2_grating,E2_envelope]=calculate_diffraction_pattern_image(input, output, dmd)
#          transmission_collected[i] = total_power_collected
#
#     plt.plot(wavelengths/nm,transmission_collected)
#
#
# plt.legend(np.around(beam_angle_x,3))
# #plt.plot(wavelengths_experimental,transmission_experimental)
# plt.xlim((420,720))
# #plt.ylim((0,1.2))
# plt.title('Transmission Function - effective_beam_width x')
# plt.show()
# fig.savefig('beam_angle_y.png')
#
