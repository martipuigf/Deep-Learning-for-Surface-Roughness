#!/usr/bin/env python3

# Import packages
import os
import numpy as np
import matplotlib.pyplot as plt
from tifffile import imwrite # Write TIFF files


from tifffile import imwrite
from tifffile import imread
import ipywidgets as widgets
from cil.io import TIFFStackReader, TIFFWriter
from cil.framework import AcquisitionGeometry
from cil.framework import ImageGeometry
from cil.recon import FDK
from gvxrPython3.JSON2gVXRDataReader import *
from cil.processors import TransmissionAbsorptionConverter
from cil.utilities.display import show_geometry, show2D
from cil.recon import FBP, FDK
from cil.plugins.astra.processors.FDK_Flexible import FDK_Flexible

has_cil = True
has_tigre = False
has_rtk = False
current_path = os.getcwd()

def reconstructFDKWithCIL(data, ig, verbose):
    if verbose > 0: print("Cone beam detected")

    if has_tigre:
        if verbose > 0: print("Backend: Tigre")
        reconstruction:ImageData | None = FDK(data, ig).run()
    else:
        if verbose > 0: print("Backend: Astra-Toolbox")
        fbk = FDK_Flexible(ig, data.geometry)
        fbk.set_input(data)
        reconstruction:ImageData | None = fbk.get_output()
    
    return reconstruction

def reconstruct(JSON_fname, verbose=0):
    """
    

    Parameters
    ----------
    JSON_fname : TYPE
        DESCRIPTION.
    verbose : TYPE, optional
        DESCRIPTION. The default is 0.

    Raises
    ------
    ValueError
        DESCRIPTION.

    Returns
    -------
    data : TYPE
        DESCRIPTION.
    reconstruction : TYPE
        DESCRIPTION.

    """
    
    data = None
    reconstruction = None
    
    source_shape = json2gvxr.params["Source"]["Shape"]

    if verbose > 0:
        print("Source shape:", source_shape)

    # Use CIL
    if has_cil:
    
        if verbose > 0: print("Use CIL")

        reader = JSON2gVXRDataReader(file_name=JSON_fname)
        data = reader.read()

        print("data.geometry", data.geometry)
        
        if has_tigre:
            data.reorder(order='tigre')
            data.geometry.set_angles(-data.geometry.angles)
        else:
            data.reorder("astra")

        ig = data.geometry.get_ImageGeometry()

        data_corr = TransmissionAbsorptionConverter(white_level=data.max(), min_intensity=0.000001)(data)

        if type(source_shape) == str:

            if source_shape.upper() == "PARALLELBEAM" or source_shape.upper() == "PARALLEL":
                reconstruction:ImageData | None = reconstructFBPWithCIL(data_corr, ig, verbose)

            elif source_shape.upper() == "POINTSOURCE" or source_shape.upper() == "POINT" or source_shape.upper() == "CONE" or source_shape.upper() == "CONEBEAM":
                reconstruction:ImageData | None = reconstructFDKWithCIL(data_corr, ig, verbose)

            else:
                raise ValueError("Unknown source shape:" + source_shape)

        elif type(source_shape) == type([]):
            if source_shape[0].upper() == "FOCALSPOT":
                reconstruction:ImageData | None = reconstructFDKWithCIL(data_corr, ig, verbose)

            else:
                raise ValueError("Unknown source shape:" + source_shape)

        else:
            raise ValueError("Unknown source shape:" + source_shape)    

    # Use ITK-RTK
    elif has_rtk:
        if verbose > 0: print("Use RTK")

        if type(source_shape) == str:

            if source_shape.upper() == "PARALLELBEAM" or source_shape.upper() == "PARALLEL":
                reconstruction = reconstructFBPWithRTK(verbose)
                
            elif source_shape.upper() == "POINTSOURCE" or source_shape.upper() == "POINT" or source_shape.upper() == "CONE" or source_shape.upper() == "CONEBEAM":
                reconstruction = reconstructFDKWithRTK(verbose)
                
            else:
                raise ValueError("Unknown source shape:" + source_shape)

        elif type(source_shape) == type([]):
            if source_shape[0].upper() == "FOCALSPOT":
                reconstruction = reconstructFDKWithRTK(verbose)

            else:
                raise ValueError("Unknown source shape:" + source_shape)

        else:
            raise ValueError("Unknown source shape:" + source_shape)    

    else:
        raise ValueError("CIL and RTK are not installed")    

    return data, reconstruction


"""Parameters of the cubes to be read in"""

A_min = 3
A_max = 5
A_int = 1

WL_min = 4
WL_max = 8
WL_int = 4

T_min = 10
T_max = 20
T_int = 10

repeats = 1

file_list = [] # creates a variable containing the name of the files to be scanned by gVXR
               # make sure to adapt this to the be consistent with the naming convention employed
for A in np.arange(A_min, A_max+A_int, A_int):
    for wavelength in np.arange(WL_min, WL_max+WL_int, WL_int):
        for thickness in np.arange(T_min, T_max+T_int, T_int):
            for repeat in np.arange(1, repeats+1, 1):
                file_title = ('A='+str(A)+', I='+str(wavelength)+', T='
                              +str(thickness)+'-'+str(repeat)+'(Cleaned)_Small_Array(Inverted)')
                file_list.append(file_title)

"""From here onwards, the code iterates through all the files conained in
the file_list variable, creates projections trough gVXR and reconstructs them 
using CIL"""

for entry in file_list:
    
    from gvxrPython3 import gvxr # Simulate X-ray images
    from gvxrPython3 import json2gvxr
    
    fname = entry+'.stl' 
    
    json2gvxr.initGVXR('notebook.json', "OPENGL") # make sure to use the correct file name for the JSON file
    
    gvxr.loadMeshFile("cube", fname, "mm") # reads in the STL file
    gvxr.setCompound("cube", "C3H4O2") # defines the object's composition
    gvxr.setDensity("cube", 1.24, "g/cm3") # defines the object's density
    
    json2gvxr.initSourceGeometry() # reads in geometry settings from JSON file
    
    spectrum, unit, k, f= json2gvxr.initSpectrum(verbose=0) # creates source spectrum

    json2gvxr.initDetector("notebook.json") 
    gvxr.moveToCentre("cube") # the object is placed at the middle
    
    raw_projection_output_dir = json2gvxr.getFilePath(json2gvxr.params["Scan"]["OutFolder"]) #output directory for the projections
    print("The raw projections were saved in", raw_projection_output_dir)
    
    angles = json2gvxr.initScan()
    
    angles = np.array(json2gvxr.doCTScan())
    gvxr.removePolygonMeshesFromSceneGraph() # removes object from simulation to add following iteration
    print('gVXR simulation completed!')
    
    data, reconstruction = reconstruct(current_path+"/notebook.json") # CIL reconstruction
    writer = TIFFWriter(data=reconstruction, file_name=(current_path+'/reconstructions'+'/'+entry+'/'+entry+'_recon')) #saves reconstruction as TIFF image sequence
    print('The reconstruction will be saved on', (current_path+'/reconstructions'+'/'+entry+'/'+entry+'_recon'))
    writer.write()
    
    print('CIL reconstruction completed!')
    
    
    
