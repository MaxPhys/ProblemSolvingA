# -*- coding: utf-8 -*-
"""
Created on Fri Feb 04

@author: Eduard Vives
febrer 2022
"""

import pandas as pnd
import numpy as np
import matplotlib.pyplot as plt

# Read the column from the csv file
# (note: you'll have to manually correct the error on
# the energy exponent directly on the file)

dades = pnd.read_csv('charcoal-2023-02-21-01-HIT.TXT', sep='\s+', skiprows=7)

# Select the column which contains the energies, all ranges. Since it is a dataframe,
# it must be done by writing the title of the column as it is shown below.
# Same thing goes for the duration of the pulseswhich will be used to filtrate the energies
# We also select the times when the hits occured.

energies = dades['ABS-ENERGY']
durades = dades['DURATION']
temps = dades['SSSSSSSS.mmmuuun']

energia_minima = np.min(energies)
energia_maxima = np.max(energies)

# Since it is easier to work with numpy arrays, a conversion must be done.
# The energies must be filtered: all values smaller than 0,5 attoJ are considered noise,
# not relevant enough to be identified as signals, and they must be removed.
# Some of the durations are 0, and they also must be erased

# We filter energies above a desired threshold (in aJ) that will be changed .
# We propose 7 values for the threshold

th_energia = [0.5, 1.0, 5.0, 10.0, 50.0, 100.0, 500.0]

energia = energies.to_numpy()
durada = durades.to_numpy()
temp = temps.to_numpy()
dimensio = np.size(energia)

# The bins of the histogram must be logarithmic, and the number of bins dim_bins may be
# modified
# We histogram the scaled waiting times between 10^-8 and 10^4

dim_bins = 49
bins_log = np.logspace(-8, 4, num=dim_bins, base=10.0)

# The number of intervals is dim_bins-1. The width of the bins is defined below,
# and it is calculated by substracting the value of two consecutive bins

amplada_bins = np.zeros(dim_bins)
for ii in range(0, dim_bins - 1):
    amplada_bins[ii] = bins_log[ii + 1] - bins_log[ii]

# We prepare the vector for the histogram heights.
# They will be computed for the different threshold values

alcada = np.zeros((dim_bins, 7))

# This index th, scans the values of the thresholds

for th in range(0, 6):
    energia_filtrada = []
    durada_filtrada = []
    temp_filtrat = []

    for ii in range(0, dimensio):
        if durada[ii] > 0 and energia[ii] > th_energia[th]:
            energia_filtrada.append(energia[ii])
            durada_filtrada.append(energia[ii])
            temp_filtrat.append(temp[ii])

    # Computing waiting times, and the mean waiting time

    dim = np.size(temp_filtrat) - 1

    meanwaiting = 0
    count = 0
    waiting = []
    for ii in range(0, dim - 1):
        waiting.append(temp_filtrat[ii + 1] - temp_filtrat[ii])
        meanwaiting = meanwaiting + temp_filtrat[ii + 1] - temp_filtrat[ii]
        count = count + 1

    meanwaiting = meanwaiting / count

    dim2 = np.size(waiting)

    # Computing main waiting time as total time divided by number of events
    temps_minim = np.min(temp_filtrat)
    temps_maxim = np.max(temp_filtrat)
    meanwaiting2 = (temps_maxim - temps_minim) / count

    # In theory it should be dim2=dim-1

    print(count)
    print(dim2, dim - 1)
    print(meanwaiting2, meanwaiting)

    waiting_minim = np.min(waiting)
    waiting_maxim = np.max(waiting)
    print(waiting_minim, waiting_maxim)

    # Compute the scaled values of the waiting times
    for ii in range(0, dim - 1):
        waiting[ii] = waiting[ii] / meanwaiting

    # The next loop is the one that actually makes the histogram. All energies must
    # be evaluated, and if they are between the highest and the lowest value of the bin,
    # 1 is added to the counter. The counter will dictate the height of the bin.

    comptador = np.zeros(dim_bins)
    comptador_total = 0
    for ii in range(0, dim_bins - 1):

        for jj in range(0, dim2):
            if waiting[jj] >= bins_log[ii] and waiting[jj] < bins_log[ii + 1]:
                comptador[ii] = comptador[ii] + 1
                comptador_total = comptador_total + 1

    # The height of the histogram is normalized by dividing the value of the counter
    # by the product of the total number of counts and the width of the bin,
    # in order to estimate the probability density.

    for ii in range(0, dim_bins - 1):
        alcada[ii, th] = comptador[ii] / (amplada_bins[ii] * comptador_total)

# This computes the normalized double power-law probability density
a = 0.91
b = 2.29
x0 = (2.0 - a) * (b - 2.0) / ((1.0 - a) * (b - 1.0))
print(x0)
fac1 = (1.0 - a) * (b - 1.0) / ((b - a) * (x0 ** (1.0 - a)))
fac2 = fac1 * (x0 ** (b - a))
exponencial = np.zeros(dim_bins)
doublepower = np.zeros(dim_bins)
for ii in range(0, dim_bins):
    exponencial[ii] = np.exp(-1.0 * bins_log[ii])
    if bins_log[ii] <= x0:
        doublepower[ii] = fac1 * bins_log[ii] ** (-a)
    else:
        doublepower[ii] = fac2 * bins_log[ii] ** (-b)

# Plotting

plt.xscale('log')
plt.yscale('log')
plt.xlabel(' Waiting time/average waiting time')
plt.ylabel('Prob. density ')
plt.title('Histogram, bins = 0.25 decades')

axes = plt.gca()
axes.set_ylim([10 ** -10, 10 ** 6])
axes.set_xlim([10 ** -8, 10 ** 4])

# plot bars
# plt.bar(bins_log, alcada, amplada_bins, color='none',edgecolor='black', align='edge')

# plot steps
plt.step(bins_log, alcada[:, 0], where='post', color='green')
plt.step(bins_log, alcada[:, 1], where='post', color='blue')
plt.step(bins_log, alcada[:, 2], where='post', color='black')
plt.step(bins_log, alcada[:, 3], where='post', color='red')
plt.step(bins_log, alcada[:, 4], where='post', color='orange')
plt.step(bins_log, alcada[:, 5], where='post', color='magenta')
plt.step(bins_log, alcada[:, 6], where='post', color='cyan')
plt.plot(bins_log, exponencial, color='black')
plt.plot(bins_log, doublepower, color='black')

plt.show()
