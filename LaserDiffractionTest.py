from functions import *
from matplotlib import colors, ticker, cm
from matplotlib.colors import LogNorm

#useful constants
um=10**-6
mm=10**-3
nm=10**-9
cm=10**-2

#Red Laser Off Direction
wavelength=633*nm
dmd=DMD_parameters('Standard',10.8*um,0.98,np.radians(12))
input=input_parameters(wavelength,0.05,np.radians(8.54),np.radians(-8.54),20*um)
output=output_parameters(0.2,np.radians(8.54),np.radians(-8.54),200)
#output=output_parameters(0.2,np.radians(-3*8.54),np.radians(3*8.54),200)

[diffraction_image, total_power_collected, E2_grating, E2_envelope] = calculate_diffraction_pattern_image(input, output,
                                                                                                    dmd)

screen_distance=7
x_distance=screen_distance*np.tan(output.angle_x_array)
y_distance=screen_distance*np.tan(output.angle_y_array)
x_distance=output.angle_x_array_deg
y_distance=output.angle_y_array_deg

#plot results
fig, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
im1=ax1.imshow(E2_grating,extent=(x_distance[0],x_distance[-1],y_distance[0],y_distance[-1]))
im2=ax2.imshow(E2_envelope,extent=(x_distance[0],x_distance[-1],y_distance[0],y_distance[-1]))
im3=ax3.imshow(diffraction_image,extent=(x_distance[0],x_distance[-1],y_distance[0],y_distance[-1]), vmin=0, vmax=10)
#ax2.grid(color='w', linestyle='-', linewidth=2)
fig.colorbar(im1, ax=ax1)
fig.colorbar(im2, ax=ax2)
fig.colorbar(im3, ax=ax3)

ax1.set_title('Grating Orders')
ax2.set_title('Mirror Envelope')
ax3.set_title('Combined Diffraction Pattern')
plt.show()
#fig.savefig('Figure_single.png')