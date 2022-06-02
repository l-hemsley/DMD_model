import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from functions import *
import time

# Python 3.7.9

#TODO - how to consider spatial offset, e.g. pixels that are off centre/corner of image
#TODO - MTF for different lenses

#useful constants
mm=10**-3
um=10**-6
nm=10**-9

class DMD_parameters:

# class for parameters related to the DMD.
# TODO - Could include DMD mask pattern info here?

    def __init__(self, pitch, fill_factor,  tilt_angle, mirror_axis_angle):
        self.pitch = pitch
        self.fill_factor = fill_factor
        self.mirror_width = pitch*fill_factor
        self.gap = pitch-self.mirror_width
        self.tilt_angle = tilt_angle
        self.mirror_axis_angle=mirror_axis_angle

class vector:

# class makes normalized 3D unit vector from x and y angles

    def __init__(self, angle_x, angle_y, direction):
        self.z = direction * 1 / \
                 np.sqrt(1+np.tan(angle_x) ** 2 + np.tan(angle_y) ** 2)
        self.x = self.z * np.tan(angle_x)
        self.y = self.z * np.tan(angle_y)

class input_parameters:
    #class related to the input parameters on the DMD
    def __init__(self, wavelength, axis_angle_x, axis_angle_y,lens,lens_diameter,focal_length):
        self.wavelength = wavelength
        self.axis_angle_x = axis_angle_x #angle of optical axis to DMD perpendicular direction, =0 if normal incidence
        self.axis_angle_y= axis_angle_y
        self.beam_vector = vector(self.axis_angle_x, self.axis_angle_y, 1)
        self.lens=lens
        self.lens_diameter=lens_diameter
        self.focal_length=focal_length

class output_parameters:

# class related to the 'output parameters' of the system after the DMD, for example the collection lens.
# gives the 2D array of angles used to build the diffraction pattern image

    def __init__(self, axis_angle_x, axis_angle_y, datapoints,DMD,input):

        self.half_angle = input.wavelength/DMD.mirror_width # double approx the first minimum of the diffraction envelope due to the mirror
        self.axis_angle_x = axis_angle_x  # angle of optical axis to DMD perpendicular direction, =0 if normal incidence
        self.axis_angle_y = axis_angle_y
        self.datapoint = datapoints #number of datapoints to use in diffraction image
        self.angle_x_array = np.linspace(
            axis_angle_x -  self.half_angle, axis_angle_x +  self.half_angle, datapoints)
        self.angle_y_array = np.linspace(
            axis_angle_y -  self.half_angle, axis_angle_y + self.half_angle, datapoints)
        [self.angle_x_array_meshed, self.angle_y_array_meshed] = np.meshgrid(
            self.angle_x_array, self.angle_y_array)
        self.beam_vector = vector(self.angle_x_array_meshed,
                             self.angle_y_array_meshed, -1)
        self.effective_angle_of_vector = np.sqrt(
            (self.angle_x_array_meshed - axis_angle_x) ** 2 + (self.angle_y_array_meshed - axis_angle_y ) ** 2) #angle compared to optical axis

def envelope_function(input, output, DMD):

    # this calculated the diffraction pattern for the tilted mirror, which is used as the envelope function for the entire diffraction pattern
    #Calculated using equation given in; 'Simulating digital micromirror devices for patterning coherent excitation light in structured illumination microscopy'
    #link = https://www.biorxiv.org/content/10.1101/2020.10.02.323527v1.supplementary-material

    w = DMD.mirror_width
    c = np.cos(DMD.tilt_angle)
    s = np.sin(DMD.tilt_angle)
    a = input.beam_vector
    b = output.beam_vector
    diff = [a.x-b.x, a.y-b.y, a.z-b.z]

    n=R_matrix(DMD)
    nx=n[0]
    ny=n[1]

    f1 = diff[0]*(nx**2*(1-c)+c)+diff[1]*nx*ny*(1-c)-diff[2]*ny*s
    f2 = diff[0] *nx*ny* (1 - c) + diff[1] *(ny**2 * (1 - c) +c)+ diff[2] *nx* s

    #f1 = diff[0]*0.5*(1+c)+diff[1]*0.5*(1-c)+diff[2]*(-s/np.sqrt(2))
    #f2 = diff[0] * 0.5 * (1 - c) + diff[1] * 0.5 * (1 + c) + diff[2] * ( s / np.sqrt(2))

    A = np.pi*f1*w/input.wavelength
    B = np.pi * f2*w/ input.wavelength
    data =w**2*np.sin(A)*np.sin(B)/(A*B)

    return data

def R_matrix(DMD):

    if DMD.mirror_axis_angle =='diagonal':
        n=[1,1]/np.sqrt(2)
    elif DMD.mirror_axis_angle =='vertical':
        n = [0, 1]
    elif DMD.mirror_axis_angle == 'horizontal' \
                                  '':
        n = [1, 0]

    return n


def calculate_orders(input, output, DMD):

    #this function calculates the angular location of diffraction orders in x and y
    #keeps only the orders which fall within the 'half angle' defined by the collection lens

    # for x
    alpha_x = input.axis_angle_x
    beta_x = output.axis_angle_x
    half_angle = output.half_angle
    order_mx_max = np.ceil(DMD.pitch*(np.sin(alpha_x)+np.sin(beta_x+half_angle))/input.wavelength)
    order_mx_min = np.floor(DMD.pitch*(np.sin(alpha_x)+np.sin(beta_x-half_angle))/input.wavelength)
    order_array_mx = np.arange(order_mx_min, order_mx_max, 1)
    order_angles_mx = np.arcsin(
        order_array_mx*input.wavelength/DMD.pitch-np.sin(alpha_x))

    #for y
    alpha_y = input.axis_angle_y
    beta_y = output.axis_angle_y
    order_my_max = np.ceil(DMD.pitch * (np.sin(alpha_y) + np.sin(beta_y + half_angle)) / input.wavelength)
    order_my_min = np.floor(DMD.pitch * (np.sin(alpha_y) + np.sin(beta_y - half_angle)) / input.wavelength)
    order_array_my = np.arange(order_my_min, order_my_max, 1)
    order_angles_my = np.arcsin(
        order_array_my * input.wavelength / DMD.pitch - np.sin(alpha_y))

    #print(order_array_mx, order_array_my)

    return order_angles_mx, order_angles_my

def MTF_function(spatial_frequency,input):

    lens=input.lens
    f=spatial_frequency
    MTF=[]
    if lens == 'ideal':
        f0=(input.lens_diameter/(input.wavelength*input.focal_length))
        MTF=((2/np.pi)*(np.arccos(f/f0)-(f/f0)*np.sqrt(1-(f/f0)**2)))
        zerovals=f>f0
        MTF[zerovals]=0
    else:
        file=lens+'-MTF.xlsx'
        experimental_data = pd.read_excel(file)
        spatial_frequency_data= experimental_data.loc[:, 'Spatial Frequency'].to_numpy()
        MTF_data = experimental_data.loc[:, 'MTF'].to_numpy()
        spatial_frequency_data=spatial_frequency_data/mm
        MTF=np.interp(f,spatial_frequency_data,MTF_data)

    return MTF

def grating_function(input, output, DMD):

    #gives diffraction order convoluted by FT of aperture function (or PSF on DMD)

    [order_angles_mx, order_angles_my] = calculate_orders(input, output, DMD)
    data = np.zeros((output.datapoint, output.datapoint))

    for angle_mx in order_angles_mx:
        for angle_my in order_angles_my:
            #The fourier transform of the MTF of the lens gives the point spread function (PSF) on the DMD (fopr ideal lens?)
            #The diffraction orders are convoluted with the FT of the 'aperture'  to give the diffraction pattern
            #Therefore the diffraction due to the spot on the DMD is the fourier transform of the PSF = the MTF, where k=2*pi*sin(angle)/wavelength
            k=2*np.pi*np.sin(np.sqrt((output.angle_x_array_meshed-angle_mx)**2+(output.angle_y_array_meshed-angle_my)**2))/input.wavelength
            data=data+MTF_function(k,input)**2
            #not sure about whether MTF needs to be squared? As there are two lenses in the system before the DMD
            #the influence of the DMD extent is negligible - but how to include the masks pattern?
    return data

def calculate_diffraction_pattern_image(input, output, DMD):

    #gives the diffraction pattern in intensity

    E2_grating = grating_function(input, output, DMD) ** 2
    E2_envelope = envelope_function(input, output, DMD) **2
    diffraction_image = E2_grating  * E2_envelope

    #calculate total power collected by lens
    spatial_frequencies = np.sin(output.effective_angle_of_vector)/input.wavelength #check this
    image_collected_MTF = diffraction_image * MTF_function(spatial_frequencies,input)
    total_power_collected_MTF= np.sum(np.sum(image_collected_MTF))

    return [diffraction_image,total_power_collected_MTF,E2_grating,E2_envelope,image_collected_MTF]

