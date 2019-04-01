# -*- coding: utf-8 -*-
"""
Created on Sat Mar 30 16:28:35 2019

@author: pnola

This script shows how to compute the attraction and repulsion rate fields (s_1 & s_2) and their associated
eigenvectors for a 3 dimensional flow.

This script Assumes that the flow is presented in a MESHGRID format, i.e. arrays are structued as [i,j,...]
where i is the index for rows or y-dimension and j is the index for columns or x-dimension.

Note! while this script is written for a 2-dimensional flow, it can be easily extended to n-dimensions.
"""
import numpy as np

#Load your velocity field
#A random sample is provided for demonstration
u = np.array([[1,3],[1,2]])
v = np.array([[1,1],[2,2]])

#This line isn't necessarary, BUT it allows the script to be easily reused for
#different sized data sets
[ydim,xdim] = u.shape

#set dx and dy, this can be done dynamically with 
#dx = x[1]-x[0] and dy = y[1]-y[0]
dx=1
dy=1

#Calculate the gradients of the velocity field
dudy,dudx = np.gradient(u,dy,dx)
dvdy,dvdx = np.gradient(v,dy,dx)

#Initialize arrays for the attraction rate, repullsion rate, and eigenvectors
#Using masked arrays can be very useful when dealing with geophysical data and
#data with gaps in it.
s1 = np.ma.empty([ydim,xdim])
Xi1 = np.ma.empty([ydim,xdim,2])
s2 = np.ma.empty([ydim,xdim])
Xi2 = np.ma.empty([ydim,xdim,2])

# For each point in the domain (i,j)
for i in range(ydim):
    for j in range(xdim):
        #Make sure the data is not masked, masked gridpoints do not work with
        #Python's linalg module
        if (dudx[i,j] and dudy[i,j] and dvdx[i,j] and dvdy[i,j]) is not np.ma.masked:
            #If the data is not masked, compute s_1, s_2 and eigenvectors
            Grad = np.array([[dudx[i,j], dudy[i,j]], [dvdx[i,j], dvdy[i,j]]])
            S = 0.5*(Grad + np.transpose(Grad))
            eigenValues, eigenVectors = np.linalg.eig(S)
            idx = eigenValues.argsort()
            s1[i,j] = eigenValues[idx[0]]
            s2[i,j] = eigenValues[idx[-1]] #The use of -1 here allows this to be more easily extended to n-dimensions.
            Xi1[i,j,:] = eigenVectors[:,idx[0]]
            Xi2[i,j,:] = eigenVectors[:,idx[-1]]
            
        else:
            #If the data is masked, then mask the grid point in the output.
            s1[i,j] = np.ma.masked
            s2[i,j] = np.ma.masked
            Xi1[i,j,1] = np.ma.masked
            Xi1[i,j,2] = np.ma.masked
            Xi2[i,j,1] = np.ma.masked
            Xi2[i,j,2] = np.ma.masked

#Done
