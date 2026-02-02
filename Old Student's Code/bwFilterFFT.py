#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 13:30:39 2023

@author: reidmarkland
"""

import os
import PIL
import os
import hyperspy.api as hs
import matplotlib.pyplot as plt
import temul.api as tml
import atomap.api as am
import numpy as np

# path = r'/Users/natanii/Desktop/Latticevector'
#
# file = r'IFFT of 1064A_200kV_15MX_JEOL  ADF_0076-2 cropped 2_bw-1.png'
#
# image = hs.load(os.path.join(path,file))
#
# # Resize the image to desired dimensions
# image_resized = image.rebin((934, 934))
# print(image_resized)
#
# import temul.api as tml
# tml.visualise_dg_filter(image_resized)
# filtered_image = tml.double_gaussian_fft_filter(image_resized, 40, 200)
# image_resized.plot()
# filtered_image.plot()
# plt.show()

import hyperspy.api as hs
import os
import numpy as np
import temul.api as tml
import matplotlib.pyplot as plt

import hyperspy.api as hs
import os
import numpy as np
import temul.api as tml
import matplotlib.pyplot as plt
import hyperspy.api as hs
import numpy as np


# Load the image
path = r'/Users/natanii/Desktop/Latticevector'
file = r'2.38V_009_IFFT of Untitled 1.png'
image = hs.load(os.path.join(path, file))


# # Convert RGB to grayscale
# # Using weighted average to account for human perception (standard formula)
# grayscale_data = 0.2989 * image.data['R'] + 0.5870 * image.data['G'] + 0.1140 * image.data['B']
#
# # Convert to float32
# grayscale_data = grayscale_data.astype(np.float32)
#
# # Replace the image data with the grayscale data
# image.data = grayscale_data

# Now you can proceed with your processing


# Resize the image to desired dimensions
image_resized = image.rebin((770, 770))

# Visualize and filter the image
tml.visualise_dg_filter(image_resized)
filtered_image = tml.double_gaussian_fft_filter(image_resized, 40, 200)

# Shift all values to be positive
min_val = filtered_image.data.min()
if min_val < 0:
    filtered_image.data -= min_val

# Plot the images
image_resized.plot()
filtered_image.plot()
plt.show()
