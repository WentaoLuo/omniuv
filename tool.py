#!/usr/bin/env python

import numpy as np
import sys
from matplotlib import pyplot as plt, rc
from mpl_toolkits.mplot3d import Axes3D

# Input:
# uvw:      in units of wave number, 2-D array
# cellsize: in mas
# nc:       grid number along one axis
# name:     saved figure file name
def plot_uv(uvw, cellsize, nc, name):

    uv  =   uvw[:, :2]

    urange  =   1. / mas2rad(cellsize)
    umin    =   -urange * 0.5  

    sc  =   1E6
    plt.clf()
    fig =   plt.figure()
    ax  =   fig.add_subplot()
    ax.set_aspect('equal')
    plt.plot(uv[:, 0]/sc, uv[:, 1]/sc, marker = '.', c = 'steelblue', \
            ls = 'none', ms = 2)
    ax.set_xlabel('$u [10^6\lambda]$')
    ax.set_ylabel('$v [10^6\lambda]$')
    ax.set_xlim(umin/sc, -umin/sc) 
    ax.set_ylim(umin/sc, -umin/sc) 
    plt.savefig(name)

# Input:
# beam:     2-D array
# cellsize: in mas
# nc:       grid number along one axis
# name:     saved figure file name
def plot_beam(beam, cellsize, nc, name):

    cellsize    /=  60E3

    hcs     =   cellsize * 0.5
    s       =   cellsize * nc * 0.5
    vmin    =   np.min(beam)
    vmax    =   np.max(beam)

    plt.clf()
    fig =   plt.figure()
    ax  =   fig.add_subplot()
    im  =   ax.imshow(beam, vmin = vmin, vmax = vmax, \
                origin = 'lower', cmap = plt.get_cmap('rainbow'), \
                extent = (-s, s, -s, s))
    ax.set_xlabel('X [arcmin]')
    ax.set_ylabel('Y [arcmin]')
#    cb  =   plt.colorbar(im, ax = ax)
#    cb.ax.set_ylabel('Flux [Jy]', rotation=90, va='bottom')
    cb  =   plt.colorbar(im, orientation='vertical')
    cb.ax.set_ylabel('Strength')

    plt.savefig(name)

# Input:
# image:    2-D array
# cellsize: in mas
# nc:       grid number along one axis
# name:     saved figure file name
def plot_image(image, cellsize, nc, name):

    cellsize    /=  60E3

    hcs     =   cellsize * 0.5
    s       =   cellsize * nc * 0.5
    vmin    =   np.min(image)
    vmax    =   np.max(image)

    plt.clf()
    fig =   plt.figure()
    ax  =   fig.add_subplot()
    im  =   ax.imshow(image, vmin = vmin, vmax = vmax, \
                origin = 'lower', cmap = plt.get_cmap('rainbow'), \
                extent = (-s, s, -s, s))
    ax.set_xlabel('X [arcmin]')
    ax.set_ylabel('Y [arcmin]')
#    cb  =   plt.colorbar(im, ax = ax)
#    cb.ax.set_ylabel('Flux [Jy]', rotation=90, va='bottom')
    cb  =   plt.colorbar(im, orientation='vertical')
    cb.ax.set_ylabel('Flux [Jy]')

    plt.savefig(name)

def plot_crs(task):

    fig =   plt.figure()
    ax  =   fig.add_subplot(projection = '3d')
    for stn in task.stns:
        
#        if stn.cb != 'Earth':
#            print('Skip non Earth crs...')
#            continue

# Retrive the CRS coordinates of each station
        x, y, z =   zip(*stn.p_crs.tolist())
        ax.plot(x, y, z)
       
    plt.show()

def plot_lcs(task):

    fig =   plt.figure()
    ax  =   fig.add_subplot(projection = '3d')
    for stn in task.stns:
        if stn.cb != 'Moon':
            print('plot_lcs(): skip none lunar station %s' % (stn.name))
            continue
# The LCS coordinates are only available for Moon related stations 
# (orbit, surface)
        x, y, z =   zip(*stn.p_lcs.tolist())
        ax.plot(x, y, z)

    Rm  =   1737.4E3 # Moon radius, in m
    lim =   (-Rm, Rm)

    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')

    plt.show()

def mas2rad(mas):
    return mas / 1E3 / 3600 / 180. * np.pi

# Input:
# hl, hh:   orbital height in m
# R:        celestial object radius in m
# Output:
# a:        semi major axis in m
# e:        eccentricity
def h2ae(hl, hh, R):
    
    a   =   (hl + hh + 2 * R) * 0.5
    c   =   a - hl - R
    e   =   c / a
    return a, e

# Calculate the beam shape using TPJ's algorithm in DIFMAP
# Input:
# uv:   2-D array, uv in wave number
# Output:
# Major and minor size in rad, angle in rad
def calc_beam_param(uv):

    uu  =   uv[:, 0]
    vv  =   uv[:, 1]

# u, v should have been used! Otherwise multiply with lam
    muu =   np.average(uu**2)
    mvv =   np.average(vv**2)
    muv =   np.average(uu*vv)
    
    fudge   =   0.7
    ftmp    =   np.sqrt((muu-mvv)**2 + 4.*muv**2)
    e_bpa   =   -0.5 * np.arctan2(2.*muv, muu-mvv)
    e_bmin  =   fudge/np.sqrt(2.*(muu+mvv) + 2.*ftmp)
    e_bmaj  =   fudge/np.sqrt(2.*(muu+mvv) - 2.*ftmp)

    if e_bmin > e_bmaj:
        e_bmaj, e_bmin  =   e_bmin, e_bmaj

    return e_bmaj, e_bmin, e_bpa

def rad2mas(rad):
    return rad / np.pi * 180. * 3600E3

def rad2min(rad):
    return rad / np.pi * 180. * 60.

def rad2deg(rad):
    return rad / np.pi * 180.

def mas2rad(mas):
    return mas / 1E3 / 3600 / 180. * np.pi


