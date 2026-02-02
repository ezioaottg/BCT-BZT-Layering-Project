#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 13:30:39 2023

@author: reidmarkland
"""

import os
import PIL
path = r'/Users/natanii/Desktop/Latticevector/'
file = r'IFFT of 1065D_3_bw_new.png'



image = PIL.Image.open(os.path.join(path, file)+'.tif').convert('L')
image.save(os.path.join(path,file)+'_bw.png')