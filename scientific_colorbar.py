# -*- coding: utf-8 -*-
"""
Created on Fri Jul 27 12:10:58 2018

A function to put a colorbar into scientific notation.

@author: pnola
"""

# This function is for putting the colorbar into scientific notation
def fmt(x, pos):
    a, b = '{:.2e}'.format(x).split('e')
    b = int(b)
    return r'${} \cdot 10^{{{}}}$'.format(a, b)
