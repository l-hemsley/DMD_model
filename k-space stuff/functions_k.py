import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from functions import *
import time

# Python 3.7.9

class DMD_parameters:
    def __init__(self, name, pitch, fill_factor,  tilt_angle,mirrors_x,mirrors_y):
        self.name = name
        self.pitch = pitch
        self.fill_factor = fill_factor
        self.mirror_width = pitch*fill_factor
        self.gap = pitch-self.mirror_width
        self.tilt_angle = tilt_angle
        self.grating_period=pitch
        self.K=2*np.pi/self.grating_period
        self.mirrors_x=mirrors_x
        self.mirrors_y=mirrors_y
        self.W=mirrors_x*pitch
        self.H = mirrors_y * pitch

def grating order location(k_i,DMD,m,n):
    k_x=k_i.k_x+m*DMD.K
    k_y=k_i.k_y+n*DMD.K
    return [k_x,k_y]

def diffraction_pattern_DMD_extent(DMD,k_x,k_y):
    I=(DMD.H*DMD.W)*np.sinc(np.pi*DMD.W*k_x/2)*np.sinc(-np.pi*DMD.H*k_y/2)
    return I


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

#
class input_parameters:
    def __init__(self, angle_x_centre, angle_y_centre):
        self.angle_x_centre = angle_x_centre
        self.angle_y_centre = angle_y_centre
        self.k_vector = vector(self.angle_x_centre, self.angle_y_centre, 1,wavelength)

#
class output_parameters:
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
        self.k_vector = vector(self.angle_x_centre, self.angle_y_centre, 1, wavelength)
        self.effective_angle_of_vector = np.sqrt(
            (self.angle_x_array_meshed-angle_x_centre)**2+(self.angle_y_array_meshed-self.angle_y_centre)**2)
        self.collected_vectors = self.effective_angle_of_vector < self.half_angle
        self.collected_vectors_triangle=abs(2*self.half_angle-self.effective_angle_of_vector)


class k_vector:
    def __init__(self, angle_x, angle_y, direction,wavelength):
        self.k=(2*np.pi/wavelength)
        self.k_z = direction*self.k/ np.sqrt(1-np.tan(angle_x)** 2-np.tan(angle_y) ** 2)
        self.k_x = self.k_z * np.tan(angle_x)
        self.k_y = self.k_z * np.tan(angle_y)
#
#
def envelope_function_k(input, output, dmd):
     w = dmd.mirror_width
     c = np.cos(dmd.tilt_angle)
     s = np.sin(dmd.tilt_angle)
     a = input.beam_vector
     b = output.beam_vector
     diff = [a.x-b.x, a.y-b.y, a.z-b.z]

     f1 = diff[0]*0.5*(1+c)+diff[1]*0.5*(1-c)+diff[2]*(-s/np.sqrt(2))
     f2 = diff[0] * 0.5 * (1 - c) + diff[1] * 0.5 * (1 + c) + diff[2] * ( s / np.sqrt(2))

    A = f1*w/2
    B =  f2*w/2
    data =w**2*np.sin(A)*np.sin(B)/(A*B)

    return data
#

def gaussian2D_normalized(x, x0, y, y0, w):
    value = np.exp(-0.5*((x-x0)**2+(y-y0)**2)/w**2)/(np.pi*w*np.sqrt(2))
    return value


def grating_function(input, output, dmd):
    [order_angles_x, order_angles_y] = calculate_orders(input, output, dmd)
    data = np.zeros((output.datapoint, output.datapoint))
    #effective beam size depends on the lens NA - given by minimum beam waist of focused gaussian beam
    m=1
    effective_beam_size=4*input.wavelength/(np.pi*input.lens_NA)/m
    sigma = input.wavelength/(2*effective_beam_size*np.pi)

    for order_x in order_angles_x:
        for order_y in order_angles_y:
            data = data+gaussian2D_normalized(output.angle_x_array_meshed,
                                              order_x, output.angle_y_array_meshed, order_y, sigma)

    return data

def calculate_diffraction_pattern_image(input, output, dmd):
    E2_grating = grating_function(input, output, dmd) ** 2
    E2_envelope = envelope_function(input, output, dmd) **2
    diffraction_image = E2_grating * E2_envelope
    image_collected = diffraction_image * output.collected_vectors
    image_collected_triangle = diffraction_image * output.collected_vectors_triangle
    total_power_collected=np.sum(np.sum(image_collected))
    total_power_collected_triangle= np.sum(np.sum(image_collected_triangle))

    return [diffraction_image,total_power_collected,E2_grating,E2_envelope,image_collected ]

