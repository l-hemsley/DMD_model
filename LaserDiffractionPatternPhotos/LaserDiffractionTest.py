from functions import *
from matplotlib import colors, ticker, cm
from matplotlib.colors import LogNorm

#useful constants
um=10**-6
mm=10**-3
nm=10**-9
cm=10**-2

wavelength=633*nm
input=input_parameters(wavelength,0.05,np.radians(8.54),np.radians(-8.54),20*um)
output=output_parameters(0.2,np.radians(8.54),np.radians(-8.54),200)

screen_distance=7
x_distance=np.degrees(output.angle_x_array_meshed)
y_distance=np.degrees(output.angle_y_array_meshed)
x_distance=screen_distance*np.tan(output.angle_x_array_meshed)
y_distance=screen_distance*np.tan(output.angle_y_array_meshed)

#Red Laser On Direction

dmd=DMD_parameters('Standard',10.8*um,0.98,np.radians(12))
[diffraction_image_on, _, E2_grating_on, E2_envelope_on] = calculate_diffraction_pattern_image(input, output,
                                                                                                    dmd)
#plot results
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
im1=ax1.contourf(x_distance,y_distance,E2_grating_on,500)
im2=ax2.contourf(x_distance,y_distance,E2_envelope_on,500)
im3=ax3.contourf(x_distance,y_distance,diffraction_image_on,500,vmax=5)
#ax2.grid(color='w', linestyle='-', linewidth=2)
fig.colorbar(im1, ax=ax1)
fig.colorbar(im2, ax=ax2)
fig.colorbar(im3, ax=ax3)

ax1.set_title('Grating Orders')
ax2.set_title('Mirror Envelope')
ax3.set_title('Combined Diffraction Pattern')
plt.show()
fig.savefig('Red_on.png')

 #Red Laser Off Direction
dmd=DMD_parameters('Standard',10.8*um,0.98,np.radians(-12))
output=output_parameters(0.2,np.radians(-22),np.radians(22),200)
screen_distance=7
x_distance=np.degrees(output.angle_x_array_meshed)
y_distance=np.degrees(output.angle_y_array_meshed)
# x_distance=screen_distance*np.tan(output.angle_x_array_meshed)
# y_distance=screen_distance*np.tan(output.angle_y_array_meshed)
[diffraction_image_off, _, E2_grating_off, E2_envelope_off] = calculate_diffraction_pattern_image(input, output,
                                                                                                     dmd)

#plot results
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
im1=ax1.contourf(x_distance,y_distance,E2_grating_off,500)
im2=ax2.contourf(x_distance,y_distance,E2_envelope_off,500)
im3=ax3.contourf(x_distance,y_distance,diffraction_image_off,500,vmax=5)
#ax2.grid(color='w', linestyle='-', linewidth=2)
fig.colorbar(im1, ax=ax1)
fig.colorbar(im2, ax=ax2)
fig.colorbar(im3, ax=ax3)

ax1.set_title('Grating Orders')
ax2.set_title('Mirror Envelope')
ax3.set_title('Combined Diffraction Pattern')
plt.show()
fig.savefig('Red_off.png')

####################Green

wavelength=533*nm
input=input_parameters(wavelength,0.05,np.radians(8.54),np.radians(-8.54),20*um)
output=output_parameters(0.2,np.radians(8.54),np.radians(-8.54),200)

screen_distance=7
x_distance=np.degrees(output.angle_x_array_meshed)
y_distance=np.degrees(output.angle_y_array_meshed)
x_distance=screen_distance*np.tan(output.angle_x_array_meshed)
y_distance=screen_distance*np.tan(output.angle_y_array_meshed)

#Green Laser On Direction

dmd=DMD_parameters('Standard',10.8*um,0.98,np.radians(12))
[diffraction_image_on, _, E2_grating_on, E2_envelope_on] = calculate_diffraction_pattern_image(input, output,
                                                                                                    dmd)
#plot results
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
im1=ax1.contourf(x_distance,y_distance,E2_grating_on,500)
im2=ax2.contourf(x_distance,y_distance,E2_envelope_on,500)
im3=ax3.contourf(x_distance,y_distance,diffraction_image_on,500,vmax=5)
#ax2.grid(color='w', linestyle='-', linewidth=2)
fig.colorbar(im1, ax=ax1)
fig.colorbar(im2, ax=ax2)
fig.colorbar(im3, ax=ax3)

ax1.set_title('Grating Orders')
ax2.set_title('Mirror Envelope')
ax3.set_title('Combined Diffraction Pattern')
plt.show()
fig.savefig('Green_on.png')

 #Green Laser Off Direction
dmd=DMD_parameters('Standard',10.8*um,0.98,np.radians(-12))
output=output_parameters(0.2,np.radians(-22),np.radians(22),200)
screen_distance=7
x_distance=np.degrees(output.angle_x_array_meshed)
y_distance=np.degrees(output.angle_y_array_meshed)
# x_distance=screen_distance*np.tan(output.angle_x_array_meshed)
# y_distance=screen_distance*np.tan(output.angle_y_array_meshed)
[diffraction_image_off, _, E2_grating_off, E2_envelope_off] = calculate_diffraction_pattern_image(input, output,
                                                                                                     dmd)

#plot results
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
im1=ax1.contourf(x_distance,y_distance,E2_grating_off,500)
im2=ax2.contourf(x_distance,y_distance,E2_envelope_off,500)
im3=ax3.contourf(x_distance,y_distance,diffraction_image_off,500,vmax=5)
#ax2.grid(color='w', linestyle='-', linewidth=2)
fig.colorbar(im1, ax=ax1)
fig.colorbar(im2, ax=ax2)
fig.colorbar(im3, ax=ax3)

ax1.set_title('Grating Orders')
ax2.set_title('Mirror Envelope')
ax3.set_title('Combined Diffraction Pattern')
plt.show()
fig.savefig('Green_off.png')