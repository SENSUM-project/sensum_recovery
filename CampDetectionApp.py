'''
-----------------------------------------------------------------------------
This program is used to extract the camp from image. The ideas is to first segment
the image and then classify the image.
The dependend libraries as follows:
1. python 2.7
2. numpy 1.6.1
3. scipy 0.9.0
4. gdal 1.10.1
5. pymorph library (http://luispedro.org/software/pymorph)
6. wavelet library (http://www.pybytes.com/pywavelets/)
author: Dr. Shifeng Wang
create: 22 Oct 2013
modified: 22 Oct 2013
-----------------------------------------------------------------------------
'''

from ImageProcLib import *

# tent detection
strInputFile = 'I:\\Shifeng_Wang\\program\\camp\\camp_data_from_thailand\\BNK_Camp.tif'
strOutputFile = 'I:\\Shifeng_Wang\\program\\camp\\camp_data_from_thailand\\BNK_Camp_Detection.tif'
tentDetection_MM(strInputFile, 1, 60, strOutputFile)
