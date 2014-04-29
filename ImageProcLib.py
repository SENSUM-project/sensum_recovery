'''
------------------------------------------------------------
description: the library is for image processing and can be used to camp detection.
The idea behind is to first use pansharpened, segmentation and classification.
author: Dr. Shifeng Wang, wangsf1013@gmail.com
time:18 Oct 2013

The dependend libraries as follows:
1. python 2.7
2. numpy 1.6.1
3. scipy 0.9.0
4. gdal 1.10.1
5. pymorph library (http://luispedro.org/software/pymorph)
6. wavelet library (http://www.pybytes.com/pywavelets/)
------------------------------------------------------------
'''

import os, sys
import string
import osgeo.gdal, gdal
from osgeo.gdalconst import *
from gdalconst import *
import scipy as sp
import scipy.stats
import numpy as np
from numpy import unravel_index
import osgeo.osr
import osgeo.ogr


'''
---------------------------------------------------------------
Description: The function generate the slash for different system
author: Dr. Shifeng Wang
create: 18 Oct 2013
modified:
---------------------------------------------------------------
'''
def sys_sep():
    if os.name == 'posix':
       separator = '/'
    else:
       separator = '\\'
    return separator

'''
---------------------------------------------------------------
Description: Uses a gdal geomatrix  to calculate  the geospatial
             coordinates of top-left and down-right pixel
----------------------------------------------------------------
    '''
def Pixel2world(gt, cols, rows):
        
    minx = gt[0]
    miny = gt[3] + cols * gt[4] + rows * gt[5] 
    maxx = gt[0] + cols * gt[1] + rows * gt[2]
    maxy = gt[3] 
    
    
    return (maxx, miny)

'''
--------------------------------------------------------------------------------------------------
Description: The function is used to write out data to image. It is modified from EU center
Input: strOutFilePath: the full path of the output file
       iBands: the number of bands to be saved
       array_list: list containing all the data to be written; each element of the list should be a matrix
       icols: number of columns, in case set to 0 the number of columns is taken from the reference image
       irows: number of rows, in case set to 0 the number of rows is taken from the reference image
       itype: type of data to be written into the output file, if 0 the default is GDT_FLoat32
       strRefProjFilePath: the full path of the reference image used to get the projection
Output:
author: Dr. Shifeng Wang
Create:
Modified:
--------------------------------------------------------------------------------------------------
'''
def WriteOutputImage(strOutFilePath, iBands, array_list, iCols, iRows, iType, strRefProjFilePath):
    
    # create the output image using a reference image for the projection
    # type is the type of data
    # array_list is a list containing all the data matrixes; a list is used because could be more than one matrix (more than one band)
    # if cols and rows are not provided, the algorithm uses values from the reference image
    # nbands contains the number of bands in the output image
    # print ('len(array_list[0]',len(array_list[0]))
    
    if iType == 0:
        iType = GDT_Float32
        
    inb = osgeo.gdal.Open(strRefProjFilePath, GA_ReadOnly)
    
    driver = inb.GetDriver()
    
    if iRows == 0 or iCols == 0:
        iRows = inb.RasterYSize
        iCols = inb.RasterXSize
   
    outDs = driver.Create(strOutFilePath, iCols, iRows, iBands, iType)
    
    if outDs is None:
        print 'Could not create file'
        sys.exit(1)
        
    for i in range(iBands):
        outBand = outDs.GetRasterBand(i + 1)
        outmatrix = array_list[i].reshape(iRows, iCols)
        outBand.WriteArray(outmatrix, 0, 0)
        
    # georeference the image and set the projection
    outDs.SetGeoTransform(inb.GetGeoTransform())
    outDs.SetProjection(inb.GetProjection())


''''
-------------------------------------------------------------------------------------------
Description: Detect the tent from high resolution remote sensing image,
             using pymorph library (http://luispedro.org/software/pymorph)
Inputs: strInputFile: the full path of the high resolution image which contains tents
        resolution: the resolution of image. It should be converted into meter
        maxTentArea: the maximual area of tent included in image. It should be m^2
        strOutputFile: the full path of the output image
outputs:
created: 21 Nov 2013
modified: 2 Dec 2013
Author: Dr. Shifeng Wang
-------------------------------------------------------------------------------------------
'''              
def tentDetection_MM(strInputFile, resolution, maxTentArea, strOutputFile, strShape='box', iThresh_coeff=0):
    
    # five step to do this
    # 1. opening-determine the square structure element (6-60 m2/resolution)
    # 2. opening by reconstruction
    # 3. top-hat by reconstruction
    # 4. lower threshold
    # 5. double threshold
    import pymorph as pymm
        
    objImg = osgeo.gdal.Open(strInputFile, GA_ReadOnly)
    nRasterCount = objImg.RasterCount
    poDataset = objImg.ReadAsArray().astype(np.float) 
    # NoDataValue = objImg.GetRasterBand(1).GetNoDataValue()
    
    # gray scale image
    if (nRasterCount == 1):  
        objnImg = pymm.to_int32(poDataset)
    # RGB image   
    elif(nRasterCount == 3):
        objnImg = pymm.to_gray(poDataset) 
    else:
        print 'it only supports gray-scale or RGB image'
        sys.exit(1)
        
    # determine the structure element
    iNum = int(np.sqrt(maxTentArea) / resolution) + 1
    if (strShape == 'box'):
        objStructureElement = pymm.sebox(iNum)
    elif (strShape == 'cross'):
        objStructureElement = pymm.secross(iNum)
    else:
        objStructureElement = pymm.sedisk(iNum)
          
    # opening
    objOpen = pymm.open(objnImg, objStructureElement)
                   
    # opening by reconstruction
    objOpenRec = pymm.openrec(objOpen, objStructureElement, objStructureElement)
        
    objtophat = pymm.openrecth(objnImg, objStructureElement, objStructureElement)
    # objtophat = pymm.subm(objnImg, objOpenRec)
              
    # objTent = pymm.threshad(objtophat, 0.25 * objnImg, 0.40 * objnImg)
    # y = mean + k*std
    (minValue, maxValue, meanValue, stdValue) = objImg.GetRasterBand(1).GetStatistics(0, 1)
    
    if (nRasterCount == 3):
       (minValue2, maxValue2, meanValue2, stdValue2) = objImg.GetRasterBand(2).GetStatistics(0, 1)
       (minValue3, maxValue3, meanValue3, stdValue3) = objImg.GetRasterBand(3).GetStatistics(0, 1)
       meanValue = 0.2989 * meanValue + 0.5870 * meanValue2 + 0.1140 * meanValue3
       maxValue = 0.2989 * maxValue + 0.5870 * maxValue2 + 0.1140 * maxValue3
       
    # meanValue = 438
    # maxValue = 2047
    threshad = meanValue + iThresh_coeff * stdValue
    
    objTent = pymm.threshad(objtophat, threshad, maxValue)
            
    data_list = []
    data_list.append(objTent)
   
    WriteOutputImage(strOutputFile, 1, data_list, 0, 0, 0, strInputFile)

    '''
-------------------------------------------------------------------------------------------
Description: Detect the tent from high resolution remote sensing image,
             using wavelet+ pymorph library (http://www.pybytes.com/pywavelets/ and 
             http://luispedro.org/software/pymorph)
Inputs: strInputFile: the full path of the high resolution image which contains tents
        resolution: the resolution of image. It should be converted into meter
        maxTentArea: the maximual area of tent included in image. It should be m^2
        strOutputFile: the full path of the output image
outputs:
created: 29 Nov 2013
modified: 2 Dec 2013
Author: Dr. Shifeng Wang
-------------------------------------------------------------------------------------------
'''                
def tentDetection_wt_mm(strInputFile, resolution, maxTentArea, strOutputFile, strShape='box', iThresh_coeff=0):
    import pywt
    import pymorph as pymm
    
    objImg = osgeo.gdal.Open(strInputFile, GA_ReadOnly)
    nRasterCount = objImg.RasterCount
    poDataset = objImg.ReadAsArray().astype(np.float) 
    # NoDataValue = objImg.GetRasterBand(1).GetNoDataValue()
    
    # gray scale image
    if (nRasterCount == 1):  
        objnImg = pymm.to_int32(poDataset)
    # RGB image   
    elif(nRasterCount == 3):
        objnImg = pymm.to_gray(poDataset) 
    else:
        print 'it only supports gray-scale or RGB image'
        sys.exit(1)
        
    # determine the structure element
    iNum = int(np.sqrt(maxTentArea) / resolution) + 1
    if (strShape == 'box'):
        objStructureElement = pymm.sebox(iNum)
    elif (strShape == 'cross'):
        objStructureElement = pymm.secross(iNum)
    else:
        objStructureElement = pymm.sedisk(iNum)
          
    # decomposition until 1 level
    wp = pywt.WaveletPacket2D(data=objnImg, wavelet='db4', mode='sym', maxlevel=1)
    # iMaxLevel = wp.maxlevel()
    # top-hat
    wp['h'].data = pymm.openrecth(pymm.to_int32(wp['h'].data), objStructureElement, objStructureElement)
    wp['v'].data = pymm.openrecth(pymm.to_int32(wp['v'].data), objStructureElement, objStructureElement) 
    wp['d'].data = pymm.openrecth(pymm.to_int32(wp['d'].data), objStructureElement, objStructureElement)
    wp['a'].data = 0.5 * wp['a'].data
    # reconstruction 
    wp.reconstruct(update=True)
    
    # top-hat for reconstructed image
    objtophat = pymm.openrecth(pymm.to_int32(wp.data), objStructureElement, objStructureElement)
    
    # y = mean + k*std
    (minValue, maxValue, meanValue, stdValue) = objImg.GetRasterBand(1).GetStatistics(0, 1)
    
    if (nRasterCount == 3):
       (minValue2, maxValue2, meanValue2, stdValue2) = objImg.GetRasterBand(2).GetStatistics(0, 1)
       (minValue3, maxValue3, meanValue3, stdValue3) = objImg.GetRasterBand(3).GetStatistics(0, 1)
       meanValue = 0.2989 * meanValue + 0.5870 * meanValue2 + 0.1140 * meanValue3
       maxValue = 0.2989 * maxValue + 0.5870 * maxValue2 + 0.1140 * maxValue3
    
    # meanValue = 438
    # maxValue = 2047

    threshad = meanValue + iThresh_coeff * stdValue
    
    objTent = pymm.threshad(objtophat, stdValue, maxValue)
            
    data_list = []
    data_list.append(objTent)
   
    WriteOutputImage(strOutputFile, 1, data_list, 0, 0, 0, strInputFile)
    
'''
-------------------------------------------------------------------------------------------
Description: Detect the possible camp from the tent detection function, then this function
             to calcuate the average density of density for each window             
Inputs: strInputFile: the full path of the high resolution image which contains tents
        resolution: the resolution of image. It should be converted into meter
        maxTentArea: the maximual area of tent included in image. It should be m^2
        strOutputFile: the full path of the output image
outputs:
created: 23 Apr 2014
modified: 
Author: Dr. Shifeng Wang
-------------------------------------------------------------------------------------------
'''     
def densityWindow(objImgData, iWinH, iWinW, fThreshold):
    
    iImgRows = len(objImgData)  # rows of image
    iImgCols = len(objImgData[0])  # cols of image
    
    iMaxPosRow = np.int(iImgRows / iWinH)
    iMaxPosCol = np.int(iImgCols / iWinW)
    
    for i in range(0, iMaxPosRow - 1):
        for j in range(0, iMaxPosCol - 1):
            tempImgData = objImgData[i * iWinH:(i + 1) * iWinH - 1, j * iWinW:(j + 1) * iWinW - 1 ]  # extract the data
            density = np.average(tempImgData)
            if (density >= fThreshold):
                tempImgData[:, :] = 1
            else:
                tempImgData[:, :] = 0
            objImgData[i * iWinH:(i + 1) * iWinH - 1, j * iWinW:(j + 1) * iWinW - 1 ] = tempImgData[:, :]
    
    # deal with that the window is not cover the whole image
    bRow = False
    bCol = False
    if (iMaxPosRow * iWinH < iImgRows):
       bRow = True
    if(iMaxPosCol * iWinW < iImgCols):
        bCol = True
    
    if(bRow) and (bCol):
        for i in range(0, iMaxPosRow - 1):
            tempImgData = objImgData[i * iWinH:(i + 1) * iWinH - 1, iMaxPosCol * iWinW:iImgCols - 1 ]  # extract the data
            density = np.average(tempImgData)
            if (density >= fThreshold):
               tempImgData[:, :] = 1
            else:
               tempImgData[:, :] = 0
            objImgData[i * iWinH:(i + 1) * iWinH - 1, iMaxPosCol * iWinW:iImgCols - 1] = tempImgData[:, :]
        for j in range(0, iMaxPosCol - 1):
            tempImgData = objImgData[iMaxPosRow * iWinH:iImgRows - 1, j * iWinW:(j + 1) * iWinW - 1]
            density = np.average(tempImgData)
            if (density >= fThreshold):
               tempImgData[:, :] = 1
            else:
               tempImgData[:, :] = 0
            objImgData[iMaxPosRow * iWinH:iImgRows - 1, j * iWinW:(j + 1) * iWinW - 1] = tempImgData[:, :]
        # the last row and column
        tempImgData = objImgData[iMaxPosRow * iWinH:iImgRows - 1, iMaxPosCol * iWinW:iImgCols - 1] 
        density = np.average(tempImgData)
        if (density >= fThreshold):
            tempImgData[:, :] = 1
        else:
            tempImgData[:, :] = 0
        objImgData[iMaxPosRow * iWinH:iImgRows - 1, iMaxPosCol * iWinW:iImgCols - 1] = tempImgData[:, :]    
    elif(bRow):
        for j in range(0, iMaxPosCol - 1):
            tempImgData = objImgData[iMaxPosRow * iWinH:iImgRows - 1, j * iWinW:(j + 1) * iWinW - 1]
            density = np.average(tempImgData)
            if (density >= fThreshold):
               tempImgData[:, :] = 1
            else:
               tempImgData[:, :] = 0
            objImgData[iMaxPosRow * iWinH:iImgRows - 1, j * iWinW:(j + 1) * iWinW - 1] = tempImgData[:, :]
    elif(bCol):
         for i in range(0, iMaxPosRow - 1):
            tempImgData = objImgData[i * iWinH:(i + 1) * iWinH - 1, iMaxPosCol * iWinW:iImgCols - 1 ]  # extract the data
            density = np.average(tempImgData)
            if (density >= fThreshold):
               tempImgData[:, :] = 1
            else:
               tempImgData[:, :] = 0
            objImgData[i * iWinH:(i + 1) * iWinH - 1, iMaxPosCol * iWinW:iImgCols - 1] = tempImgData[:, :]
    return objImgData
'''
###################################################################################################################
Conversion from ESRI shapefile to raster
Input:
- input_shape: path of the input shapefile
- output_image: path and name of the output raster file
- rows: rows of the output raster
- cols: columns of the output raster
- field_name: name of the field from the shapefile used to differenciate segments (for example DN)
Output:
Nothing is returned. Output image is automatically saved.
###################################################################################################################
'''
def Shp2Rast(input_shape, output_image, rows, cols, field_name, px_W, px_H, x_min, x_max, y_min, y_max):
           
    driver_shape = osgeo.ogr.GetDriverByName('ESRI Shapefile')
    data_source = driver_shape.Open(input_shape)
    source_layer = data_source.GetLayer()
    source_srs = source_layer.GetSpatialRef()
    if x_min == 0 or x_max == 0 or y_min == 0 or y_max == 0:
        x_min, x_max, y_min, y_max = source_layer.GetExtent()
    
    if rows != 0 and cols != 0 and px_W != 0 and px_H != 0 and x_min != 0 and y_max != 0:
        pixel_size_x = px_W
        pixel_size_y = abs(px_H)
        
    else:
        if rows != 0 and cols != 0:
            pixel_size_x = float((x_max - x_min)) / float(cols)
            pixel_size_y = float((y_max - y_min)) / float(rows)
        else:
            pixel_size_x = px_W
            pixel_size_y = abs(px_H)
            cols = int(float((x_max - x_min)) / float(pixel_size_x))
            rows = int(float((y_max - y_min)) / float(pixel_size_y))
    if rows != 0 and cols != 0:
        target_ds = osgeo.gdal.GetDriverByName('GTiff').Create(output_image, cols, rows, 1, GDT_Float32)
        target_ds.SetGeoTransform((x_min, pixel_size_x, 0, y_max, 0, -pixel_size_y))
        if source_srs:
            # Make the target raster have the same projection as the source
            target_ds.SetProjection(source_srs.ExportToWkt())
        else:
            # Source has no projection (needs GDAL >= 1.7.0 to work)
            target_ds.SetProjection('LOCAL_CS["arbitrary"]')
        
        # Rasterize
        err = osgeo.gdal.RasterizeLayer(target_ds, [1], source_layer, burn_values=[0], options=["ATTRIBUTE=" + field_name])
        if err != 0:
            raise Exception("error rasterizing layer: %s" % err)
        
    return x_min, x_max, y_min, y_max
'''
###################################################################################################################
Conversion from raster to ESRI shapefile
Input:
- input_image: path of the input raster
- output_shape: path and name of the output shapefile
Output:
Nothing is returned. Output shapefile is automatically saved.
###################################################################################################################
'''    
def Rast2Shp(input_image, output_shape):
       
    src_image = osgeo.gdal.Open(input_image)
    src_band = src_image.GetRasterBand(1)
    projection = src_image.GetProjection()
    # mask = np.equal(src_band,1)
    
    driver_shape = osgeo.ogr.GetDriverByName('ESRI Shapefile')
    outfile = driver_shape.CreateDataSource(output_shape)
    outlayer = outfile.CreateLayer('Conversion', geom_type=osgeo.ogr.wkbPolygon)
    dn = osgeo.ogr.FieldDefn('DN', osgeo.ogr.OFTInteger)
    outlayer.CreateField(dn)
    
    # Polygonize
    osgeo.gdal.Polygonize(src_band, src_band.GetMaskBand(), outlayer, 0)
    
    outprj = osgeo.osr.SpatialReference(projection)
    outprj.MorphToESRI()
    file_prj = open(output_shape[:-4] + '.prj', 'w')
    file_prj.write(outprj.ExportToWkt())
    file_prj.close()
    
'''
-------------------------------------------------------------------------------------------
Description: vector statistics            
Inputs: strInputFile: the full path of the high resolution image which contains tents
        outputs:
created: 28 Apr 2014
modified: 
Author: Dr. Shifeng Wang
-------------------------------------------------------------------------------------------
'''     
def vecStat(strInputFile, strOutputFile):
    driver = osgeo.ogr.GetDriverByName('ESRI Shapefile')
    dataSource = driver.Open(strInputFile, 0)
    if dataSource is None:
        print 'Could not open %s' % (strInputFile)
        sys.exit(1)
        
    text_file = open(strOutputFile, 'w')
    if text_file is None:
        print 'could not create %' % (strOutputFile)
        sys.exit(1)
        
    layer = dataSource.GetLayer()
    featureCount = layer.GetFeatureCount()
    text_file.write('Total number of feature:    ' + str(featureCount) + '\n')
    
    text_file.write('ID     length      area\n')
    for feature in layer:
        text_file.write(str(feature.GetFID()) + '     ')
        geom = feature.GetometryRef()
        length = geom.Length()
        text_file.write(str(length) + '     ')
        area = geom.GetArea()
        text_file.write(str(area) + '     \n')
    text_file.close()
