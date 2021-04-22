import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from functions import *
import time

# Python 3.7.9

###TO DO
# sort normalization of E so that can calculate total power, and transmission as a percentage.
# camera pixel location

class DMD_parameters:

# class for parameters related to the DMD.
# TODO - Could include DMD mask pattern info here?
# TODO - Also should include mask x/y dimensions but so far this is neglected in the diffraction pattern as the effect is small.

    def __init__(self, pitch, fill_factor,  tilt_angle):
        self.pitch = pitch
        self.fill_factor = fill_factor
        self.mirror_width = pitch*fill_factor
        self.gap = pitch-self.mirror_width
        self.tilt_angle = tilt_angle


class vector:

# class makes normalized 3D unit vector from x and y angles

    def __init__(self, angle_x, angle_y, direction):
        self.z = direction * 1 / \
                 np.sqrt(1 - np.tan(angle_x) ** 2 - np.tan(angle_y) ** 2)
        self.x = self.z * np.tan(angle_x)
        self.y = self.z * np.tan(angle_y)

class input_parameters:

# class related to so-called 'input parameters' of the system before the DMD, for example, the focusing lens parameters and wavelength.
# TODO - What if there is no input lens?

    def __init__(self, wavelength, lens_NA, angle_x_centre, angle_y_centre, focal_length):
        self.wavelength = wavelength
        self.lens_NA = lens_NA
        self.focal_length=focal_length
        self.half_angle = np.arcsin(lens_NA)
        self.angle_x_centre = angle_x_centre
        self.angle_y_centre = angle_y_centre
        self.beam_vector = vector(self.angle_x_centre, self.angle_y_centre, 1)


class output_parameters:

# class related to the 'output parameters' of the system after the DMD, for example the collection lens.
# This is useful for working out which bit of the diffraction pattern is relevant
# gives the 2D array of discretized angles used to build the diffraction pattern image


    def __init__(self, lens_NA, angle_x_centre, angle_y_centre, datapoints):
        self.lens_NA = lens_NA
        self.half_angle = np.arcsin(lens_NA)
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

        # params below for finding which angles are within the NA of the lens
        self.effective_angle_of_vector = np.sqrt(
            (self.angle_x_array_meshed-angle_x_centre)**2+(self.angle_y_array_meshed-self.angle_y_centre)**2)
        # collected vectors are the angles which fall within the 'half angle' defined by the NA of the output lens. i.e. which vectors can be collected.
        self.collected_vectors = self.effective_angle_of_vector < self.half_angle
        # this weights the vectors using a triangle function (max in the centre, to zero at 2* the half angle) to account for the lens MTF
        # TODO - check this is correct method to account for lens MTF
        self.collected_vectors_triangle=abs(2*self.half_angle-self.effective_angle_of_vector)*self.collected_vectors

def envelope_function(input, output, dmd):

    # this calculated the diffraction pattern for the tilted mirror, which is used as the envelope function for the entire diffraction pattern
    #Calculated using equation given in; 'Simulating digital micromirror devices for patterning coherent excitation light in structured illumination microscopy'
    #link = https://www.biorxiv.org/content/10.1101/2020.10.02.323527v1.supplementary-material

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
    data =w**2*np.sin(A)*np.sin(B)/(A*B)

    return data


def calculate_orders(input, output, DMD):

    #this function calculates the angular location of diffraction orders due to the DMD periodicity in x and y
    #keeps only the orders which fall within the 'half angle' defined by the collection lens

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


def grating_function(input, output, dmd):

    # this function gives the diffraction pattern due to the DMD periodicity
    # this function puts a guassian function (FT of a the assumed gaussian spot on the DMD) wherever there is an order
    # equivalent to FT of beam convoluted with the 2D array of delta functions
    # TODO - could include mask here? better way to do convolution?

    [order_angles_x, order_angles_y] = calculate_orders(input, output, dmd)
    data = np.zeros((output.datapoint, output.datapoint))

    # we assume that the effective beam size depends on the input lens NA - given by minimum beam waist of focused gaussian beam.
    # TODO - another way to define this?
    effective_beam_size=4*input.wavelength/(np.pi*input.lens_NA)
    sigma = input.wavelength/(2*effective_beam_size*np.pi)

    for order_x in order_angles_x:
        for order_y in order_angles_y:
            data = data+gaussian2D_normalized(output.angle_x_array_meshed,
                                              order_x, output.angle_y_array_meshed, order_y, sigma)

    return data

def calculate_diffraction_pattern_image(input, output, dmd):

    #gives the diffraction pattern

    E2_grating = grating_function(input, output, dmd) ** 2
    E2_envelope = envelope_function(input, output, dmd) **2
    diffraction_image = E2_grating * E2_envelope

    #calculate total power collected by lens

    #image_collected = diffraction_image * output.collected_vectors
    image_collected_triangle = diffraction_image * output.collected_vectors_triangle
    #total_power_collected=np.sum(np.sum(image_collected))
    total_power_collected_triangle= np.sum(np.sum(image_collected_triangle))
    total_power_collected=total_power_collected_triangle

    return [diffraction_image,total_power_collected,E2_grating,E2_envelope,image_collected ]

