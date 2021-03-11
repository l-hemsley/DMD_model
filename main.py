
import matplotlib.pyplot as plt
import numpy as np

#useful constants
um=10**-6
nm=10**-9

####

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
        self.half_angle=np.degrees(np.arcsin(lens_NA))
        self.angle_x_centre=angle_x_centre
        self.angle_y_centre=angle_y_centre
        self.effective_beam_size=effective_beam_size

class output_parameters:

    def __init__(self, lens_NA, angle_x_centre, angle_y_centre, datapoints):
        self.lens_NA = lens_NA
        self.half_angle=np.degrees(np.arcsin(lens_NA))
        self.angle_x_centre = angle_x_centre
        self.angle_y_centre = angle_y_centre
        self.datapoint=datapoints
        self.angle_x_array=np.linspace(angle_x_centre-self.half_angle,angle_x_centre+self.half_angle,datapoints)
        self.angle_y_array = np.linspace(angle_y_centre - self.half_angle, angle_y_centre + self.half_angle, datapoints)
        [self.angle_x_array_meshed,self.angle_y_array_meshed]=np.meshgrid(self.angle_x_array,self.angle_y_array)

def calculate_orders(input, output, DMD):
    #find range of orders which fit in the output lens

    alpha_x=np.radians(input.angle_x_centre)
    beta_x=np.radians(output.angle_x_centre)
    half_angle=np.radians(output.half_angle)
    order_x_max=np.ceil(DMD.pitch*(np.sin(alpha_x)+np.sin(beta_x+half_angle))/input.wavelength)
    order_x_min=np.floor(DMD.pitch*(np.sin(alpha_x)+np.sin(beta_x-half_angle))/input.wavelength)
    order_array_x=np.arange(order_x_min,order_x_max,1)
    order_angles_x=np.arcsin(order_array_x*input.wavelength/DMD.pitch-np.sin(alpha_x))


    alpha_y = np.radians(input.angle_y_centre)
    beta_y = np.radians(output.angle_y_centre)
    order_y_max = np.ceil(DMD.pitch * (np.sin(alpha_y) + np.sin(beta_y + half_angle)) / input.wavelength)
    order_y_min = np.floor(DMD.pitch * (np.sin(alpha_y) + np.sin(beta_y - half_angle)) / input.wavelength)
    order_array_y = np.arange(order_y_min, order_y_max, 1)
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
            data=data+gaussian2D_normalized(np.radians(output.angle_x_array_meshed),order_x,np.radians(output.angle_y_array_meshed),order_y,sigma)
    return data



########################## main

dmd=DMD_parameters('Standard',10.8*um,0.95,12)
input=input_parameters(600*nm,0.05,8.54,-8.54,10*um)
output=output_parameters(0.05,8.54,-8.54,100)

# find orders for DMD and input
wavelengths=np.arange(400*nm,700*nm,5*nm)
for w in wavelengths:
    input.wavelength=w;
    data=grating_function(input, output, dmd)

    # fig, ax = plt.subplots()
    # im=ax.imshow(data,extent=(output.angle_x_array[0],output.angle_x_array[-1],output.angle_y_array[0],output.angle_y_array[-1]))
    # ax.text(output.angle_x_array[0],output.angle_y_array[0],'wavelength='+ str(int(input.wavelength/nm))+ 'nm, beam size='+str(input.effective_beam_size/um)+'um',color='w')
    # plt.xlabel('Angle Out x')
    # plt.ylabel('Angle Out y')
    # plt.title('Diffraction Orders, input ('+ str(input.angle_x_centre) +','+ str(input.angle_y_centre)+')')
    # fig.colorbar(im)
    # plt.show()

input=input_parameters(np.arange(400*nm,700*nm,5*nm),0.05,8.54,-8.54,10*um)#
print(input.wavelength)