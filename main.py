
import matplotlib.pyplot as plt
import numpy as np

#useful constants
um=10**-6
nm=10**-9

#### Python 3.7.9

class DMD_parameters:

    def __init__(self, name, pitch, fill_factor,  tilt_angle):
        self.name=name
        self.pitch = pitch
        self.fill_factor=fill_factor
        self.mirror_width=pitch*fill_factor
        self.gap=pitch-self.mirror_width
        self.tilt_angle=tilt_angle

class input_parameters:

    def __init__(self, wavelength,lens_NA,angle_x_centre,angle_y_centre, effective_beam_size):
        self.wavelength=wavelength
        self.lens_NA=lens_NA
        self.half_angle=np.arcsin(lens_NA)
        self.angle_x_centre=angle_x_centre
        self.angle_y_centre=angle_y_centre
        self.effective_beam_size=effective_beam_size
        self.vector_z=1/np.sqrt(np.tan(angle_x_centre**2+np.tan(angle_y_centre)**2+1))
        self.vector_x=self.vector_z*np.tan(angle_x_centre)
        self.vector_y= self.vector_z * np.tan(angle_y_centre)
        self.vector = np.array([self.vector_x, self.vector_y, self.vector_z])

class output_parameters:

    def __init__(self, lens_NA, angle_x_centre, angle_y_centre, datapoints):
        self.lens_NA = lens_NA
        self.half_angle=np.arcsin(lens_NA)
        self.angle_x_centre = angle_x_centre
        self.angle_y_centre = angle_y_centre
        self.datapoint=datapoints
        self.angle_x_array=np.linspace(angle_x_centre-self.half_angle,angle_x_centre+self.half_angle,datapoints)
        self.angle_y_array = np.linspace(angle_y_centre - self.half_angle, angle_y_centre + self.half_angle, datapoints)
        self.angle_x_array_deg=np.degrees(self.angle_x_array)
        self.angle_y_array_deg = np.degrees(self.angle_y_array)
        [self.angle_x_array_meshed,self.angle_y_array_meshed]=np.meshgrid(self.angle_x_array,self.angle_y_array)
        self.vector_z = -1 / np.sqrt(np.tan(self.angle_x_array_meshed ** 2 + np.tan(self.angle_y_array_meshed) ** 2 + 1))
        self.vector_x = self.vector_z * np.tan(self.angle_x_array_meshed)
        self.vector_y = self.vector_z * np.tan(self.angle_y_array_meshed)
        self.vector = np.array([self.vector_x, self.vector_y, self.vector_z])
        self.effective_NA_of_vector=np.sqrt((self.angle_x_array_meshed-angle_x_centre)**2+(self.angle_y_array_meshed-self.angle_y_centre)**2)
        self.collected_vectors=self.effective_NA_of_vector<self.half_angle

def envelope_function(input, output, dmd):
    w=dmd.mirror_width
    c=np.cos(dmd.tilt_angle)
    s = np.sin(dmd.tilt_angle)
    a=input.vector
    b=output.vector
    a=np.reshape(a,(3,1,1))
    #print(c,s)
    diff=a-b
    f1=diff[0]*0.5*(1+c)+diff[1]*(0.5*(1-c))+diff[2]*(s/np.sqrt(2))
    f2= diff[0] *( 0.5 * (1 - c))+ diff[1] * (0.5 * (1 + c)) + diff[2] * ( - s / np.sqrt(2))
    A=2*np.pi*f1/input.wavelength
    B = 2 * np.pi * f2/ input.wavelength
    data=(2/A)*(2/B)*np.sin(A*w/2)*np.sin(B*w/2)
    return data

def calculate_orders(input, output, DMD):
    #find range of orders which fit in the output lens

    alpha_x=input.angle_x_centre
    beta_x=output.angle_x_centre
    half_angle=output.half_angle
    order_x_max=np.ceil(DMD.pitch*(np.sin(alpha_x)+np.sin(beta_x+half_angle))/input.wavelength)
    order_x_min=np.floor(DMD.pitch*(np.sin(alpha_x)+np.sin(beta_x-half_angle))/input.wavelength)
    order_array_x=np.arange(order_x_min-1,order_x_max+1,1)
    order_angles_x=np.arcsin(order_array_x*input.wavelength/DMD.pitch-np.sin(alpha_x))


    alpha_y = input.angle_y_centre
    beta_y = output.angle_y_centre
    order_y_max = np.ceil(DMD.pitch * (np.sin(alpha_y) + np.sin(beta_y + half_angle)) / input.wavelength)
    order_y_min = np.floor(DMD.pitch * (np.sin(alpha_y) + np.sin(beta_y - half_angle)) / input.wavelength)
    order_array_y = np.arange(order_y_min-1, order_y_max+1, 1)
    order_angles_y = np.arcsin(order_array_y * input.wavelength / DMD.pitch - np.sin(alpha_y))

    return order_angles_x, order_angles_y

def gaussian2D_normalized(x,x0,y,y0,w):
    value=np.exp(-0.5*((x-x0)**2+(y-y0)**2)/w**2)/(np.pi*w*np.sqrt(2))
    return value

def grating_function(input, output, dmd):
    [order_angles_x, order_angles_y]=calculate_orders(input, output, dmd)
    data=np.zeros((output.datapoint,output.datapoint))

    sigma=input.wavelength/(2*input.effective_beam_size*np.pi)

    for order_x in order_angles_x:
        for order_y in order_angles_y:
            data=data+gaussian2D_normalized(output.angle_x_array_meshed,order_x,output.angle_y_array_meshed,order_y,sigma)
    return data

########################## main

dmd=DMD_parameters('Standard',10.8*um,0.95,np.radians(12))
input=input_parameters(600*nm,0.05,np.radians(-8.54),np.radians(8.54),10*um)
output=output_parameters(0.05,np.radians(-8.54),np.radians(8.54),200)
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
    im4=ax4.plot(wavelengths/nm,transmission)
    im4=ax4.plot(wavelengths/nm,transmission_collected)

    ax1.text(output.angle_x_array_deg[0],output.angle_y_array_deg[0],'wavelength='+ str(int(input.wavelength/nm))+ 'nm, beam size='+str(input.effective_beam_size/um)+'um',color='w')

    #
    ax1.set_title('Diffraction Orders, input ('+ str(np.degrees(input.angle_x_centre)) +','+ str(np.degrees(input.angle_y_centre))+')')
    ax2.set_title('Envelope, input ('+ str(np.degrees(input.angle_x_centre)) +','+ str(np.degrees(input.angle_y_centre))+')')
    ax3.set_title('Combined, input ('+ str(np.degrees(input.angle_x_centre)) +','+ str(np.degrees(input.angle_y_centre))+')')
    plt.pause(0.01)

plt.show()

# input=input_parameters(np.arange(400*nm,700*nm,5*nm),0.05,8.54,-8.54,10*um)#
# print(input.wavelength)