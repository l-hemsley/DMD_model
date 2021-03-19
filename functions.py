import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
import pandas as pd
from functions import *
from scipy.special import j1
import sys
from tqdm import tqdm
import time

# Python 3.7.9

###TO DO
# sort normalization of E so that can calculate total power, and transmission as a percentage.
# speed up timing
# camera pixel location

class DMD_parameters:

    def __init__(self, name, pitch, fill_factor,  tilt_angle):
        self.name = name
        self.pitch = pitch
        self.fill_factor = fill_factor
        self.mirror_width = pitch*fill_factor
        self.gap = pitch-self.mirror_width
        self.tilt_angle = tilt_angle


class input_parameters:

    def __init__(self, wavelength, lens_NA, angle_x_centre, angle_y_centre, focal_length):
        self.wavelength = wavelength
        self.lens_NA = lens_NA
        self.focal_length=focal_length
        self.half_angle = np.arcsin(lens_NA)
        self.angle_x_centre = angle_x_centre
        self.angle_y_centre = angle_y_centre
        self.beam_vector = vector(self.angle_x_centre, self.angle_y_centre, 1)

class output_parameters:

    def __init__(self, lens_NA, angle_x_centre, angle_y_centre, datapoints):
        self.lens_NA = lens_NA
        self.half_angle = 2*np.arcsin(lens_NA)
        self.angle_x_centre = angle_x_centre
        self.angle_y_centre = angle_y_centre
        self.datapoint = datapoints
        self.angle_x_array = np.linspace(
            angle_x_centre-self.half_angle, angle_x_centre+self.half_angle, datapoints)
        self.angle_y_array = np.linspace(
            angle_y_centre - self.half_angle, angle_y_centre + self.half_angle, datapoints)
        self.angle_x_array_deg = np.degrees(self.angle_x_array)
        self.angle_y_array_deg = np.degrees(self.angle_y_array)
        [self.angle_x_array_meshed, self.angle_y_array_meshed] = np.meshgrid(
            self.angle_x_array, self.angle_y_array)
        self.beam_vector = vector(self.angle_x_array_meshed,
                             self.angle_y_array_meshed, -1)
        self.effective_angle_of_vector = np.sqrt(
            (self.angle_x_array_meshed-angle_x_centre)**2+(self.angle_y_array_meshed-self.angle_y_centre)**2)
        self.collected_vectors = self.effective_angle_of_vector < self.half_angle
        self.collected_vectors_triangle=abs(self.half_angle-self.effective_angle_of_vector)


class vector:
    def __init__(self, angle_x, angle_y, direction):
        self.z = direction*1 / \
            np.sqrt(1-np.tan(angle_x)** 2-np.tan(angle_y) ** 2)
        self.x = self.z * np.tan(angle_x)
        self.y = self.z * np.tan(angle_y)


def envelope_function(input, output, dmd):
    #https://www.biorxiv.org/content/10.1101/2020.10.02.323527v1.supplementary-material
    #https://www.biorxiv.org/content/10.1101/2020.07.27.223941v3.supplementary-material

    w = dmd.mirror_width
    c = np.cos(dmd.tilt_angle)
    s = np.sin(dmd.tilt_angle)
    a = input.beam_vector
    b = output.beam_vector
    diff = [a.x-b.x, a.y-b.y, a.z-b.z]

    f1 = diff[0]*0.5*(1+c)+diff[1]*0.5*(1-c)+diff[2]*(-s/np.sqrt(2))
    f2 = diff[0] * 0.5 * (1 - c) + diff[1] * 0.5 * (1 + c) + diff[2] * ( s / np.sqrt(2))

    A = np.pi*f1*w/input.wavelength
    B = np.pi * f2*w/ input.wavelength
    data =np.sin(A)*np.sin(B)/(A*B)
    return data


def calculate_orders(input, output, DMD):
    # find range of orders which fit in the output lens

    # for x angles
    alpha_x = input.angle_x_centre
    beta_x = output.angle_x_centre
    half_angle = output.half_angle
    order_x_max = np.ceil(DMD.pitch*(np.sin(alpha_x)+np.sin(beta_x+half_angle))/input.wavelength)
    order_x_min = np.floor(DMD.pitch*(np.sin(alpha_x)+np.sin(beta_x-half_angle))/input.wavelength)
    order_array_x = np.arange(order_x_min-1, order_x_max+1, 1)
    order_angles_x = np.arcsin(
        order_array_x*input.wavelength/DMD.pitch-np.sin(alpha_x))

    #for y angles
    alpha_y = input.angle_y_centre
    beta_y = output.angle_y_centre
    order_y_max = np.ceil(DMD.pitch * (np.sin(alpha_y) + np.sin(beta_y + half_angle)) / input.wavelength)
    order_y_min = np.floor(DMD.pitch * (np.sin(alpha_y) + np.sin(beta_y - half_angle)) / input.wavelength)
    order_array_y = np.arange(order_y_min-1, order_y_max+1, 1)
    order_angles_y = np.arcsin(
        order_array_y * input.wavelength / DMD.pitch - np.sin(alpha_y))

    return order_angles_x, order_angles_y


def gaussian2D_normalized(x, x0, y, y0, w):
    value = np.exp(-0.5*((x-x0)**2+(y-y0)**2)/w**2)/(np.pi*w*np.sqrt(2))
    return value


def jinc(x):
    if x.all() == 0.0:
        return 0.5
    return j1(x) / x

def jinc_functions(x, x0, y, y0, a,wavelength,f):
    print(a)
    value =(2*np.pi*a**2/(4*wavelength*f))*jinc((np.pi*a/wavelength)*np.sqrt(np.tan(x-x0)**2+np.tan(y-y0)**2))
    return value


def grating_function(input, output, dmd):
    [order_angles_x, order_angles_y] = calculate_orders(input, output, dmd)
    data = np.zeros((output.datapoint, output.datapoint))

    effective_beam_size=4*input.wavelength/(np.pi*input.lens_NA)
    sigma = input.wavelength/(2*effective_beam_size*np.pi)

    for order_x in order_angles_x:
        for order_y in order_angles_y:
            data = data+gaussian2D_normalized(output.angle_x_array_meshed,
                                              order_x, output.angle_y_array_meshed, order_y, sigma)
            #data = data+jinc_functions(output.angle_x_array_meshed,
             #                                 order_x, output.angle_y_array_meshed, order_y, 2*(input.focal_length)*input.lens_NA, input.wavelength, input.focal_length)
    return data

def calculate_diffraction_pattern_image(input, output, dmd):
    E2_grating = grating_function(input, output, dmd) ** 2
    E2_envelope = envelope_function(input, output, dmd) **2
    diffraction_image = E2_grating * E2_envelope
    #image_collected = diffraction_image * output.collected_vectors
    image_collected = diffraction_image * output.collected_vectors_triangle
    total_power_collected=np.sum(np.sum(image_collected))


    return [diffraction_image,total_power_collected,E2_grating,E2_envelope]



def diff_image_integrated_input_NA(input, output, dmd,integration_points):

    input_angle_x_array = np.linspace(
        input.angle_x_centre - input.half_angle, input.angle_x_centre + input.half_angle, integration_points)
    input_angle_y_array = np.linspace(
        input.angle_y_centre - input.half_angle, input.angle_y_centre + input.half_angle, integration_points)
    [input_angle_x_array_meshed, input_angle_y_array_meshed] = np.meshgrid(
        input_angle_x_array, input_angle_y_array)

    effective_angle_of_vector = np.sqrt(
        (input_angle_x_array_meshed - input.angle_x_centre) ** 2 + (input_angle_y_array_meshed - input.angle_y_centre) ** 2)
    collected_input_vectors = abs(effective_angle_of_vector) < input.half_angle

    input_angles_x=input_angle_x_array_meshed*collected_input_vectors
    input_angles_x=input_angles_x.flatten()
    input_angles_x=input_angles_x[input_angles_x!=0]
    input_angles_y = input_angle_y_array_meshed * collected_input_vectors
    input_angles_y = input_angles_y.flatten()
    input_angles_y= input_angles_y[input_angles_y != 0]

    total_power_collected=0
    diffraction_image=E2_grating=E2_envelope=np.zeros((np.shape(output.angle_x_array_meshed)))

    # plt.plot(input_angle_x_array_meshed.flatten(),input_angle_y_array_meshed.flatten(),'*')
    # plt.plot(input_angles_x, input_angles_y, '*')
    # plt.show()
    # toolbar_width=np.size(input_angles_x)
    #
    # sys.stdout.write("[%s]" % (" " * toolbar_width))
    # sys.stdout.flush()
    # sys.stdout.write("\b" * (toolbar_width + 1))  # return to start of line, after '['

    for x in tqdm(input_angles_x):

        for y in input_angles_y:
            input_adjusted = input_parameters(input.wavelength, input.lens_NA, x, y, input.effective_beam_size)
            [diffraction_image_temp, total_power_collected_temp, E2_grating_temp,E2_envelope_temp] = calculate_diffraction_pattern_image(input_adjusted,
                                                                                                              output,
                                                                                                              dmd)
            diffraction_image=diffraction_image+diffraction_image_temp
            total_power_collected=total_power_collected+total_power_collected_temp
            E2_grating=E2_grating+E2_grating_temp
            E2_envelope=E2_envelope+E2_envelope_temp

    #     sys.stdout.write("-")
    #     sys.stdout.flush()
    #
    # sys.stdout.write("]\n")  # this ends the progress bar

    return [diffraction_image,total_power_collected,E2_grating,E2_envelope]