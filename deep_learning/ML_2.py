# -*- coding: utf-8 -*-
"""U_Net_4_superresolution_in_EM.ipynb

Automatically generated by Colaboratory.

Original file is located at
    https://colab.research.google.com/drive/1bTbbAyUQHMTCQ7nR3nqjWMnRAXa3K5-j

# Deep Learning example: U-Net for super-resolution

---
## Introduction
This is a notebook that shows how to design and train a U-Net-like network for super-resolution on Electron Miscroscopy (EM) images. The aim is to train the network using low resolution versions of the images as input, and the high resolution versions as output.

<img src="https://drive.google.com/uc?id=1KUCwas63FD6AfiKOetiGjLRR6oldNkKC" width="450">




## Data
The image data used in the notebook was produced by [Lichtman Lab at Harvard University](https://lichtmanlab.fas.harvard.edu/) (Daniel R. Berger, Richard Schalek, Narayanan "Bobby" Kasthuri, Juan-Carlos Tapia, Kenneth Hayworth, Jeff W. Lichtman). Their corresponding biological findings were published in [Cell (2015)](https://www.ncbi.nlm.nih.gov/pubmed/26232230).
The training and test data sets are both 3D stacks of 100 sections from a serial section Scanning Electron Microscopy (ssSEM) data set of mouse cortex. The microcube measures 6 x 6 x 3 microns approx., with a resolution of 6 x 6 x 30 nm/voxel. For simplicity, in this notebook we will only use 10 sections of the test set.

## Getting started
First, we make sure we are using Tensorflow version compatible with DeepImageJ (<= 1.13).
"""

# Commented out IPython magic to ensure Python compatibility.
# Use Tensorflow and Keras versions compatible with DeepImageJ
#%pip install tensorflow-gpu==2.5.0
#%pip install keras==2.2.4

"""Then, we load our Google Drive as a local folder so we can access the image files.

(Notice we expect you to have already this notebook under your `Colab Notebooks` in a folder called `U-Net-Super-resolution`. Inside that folder you should add the `train` and `test` image folders.)
"""

"""Now we should be able to read the list of **100 training images**."""

from keras.callbacks import EarlyStopping

def create_patches( imgs, num_x_patches, num_y_patches ):
    ''' Create a list of images patches out of a list of images
    Args:
        imgs: list of input images
        num_x_patches: number of patches in the X axis
        num_y_patches: number of patches in the Y axis
        
    Returns:
        list of image patches
    '''
    original_size = imgs[0].shape
    patch_width = original_size[ 0 ] // num_x_patches
    patch_height = original_size[ 1 ] // num_y_patches
    
    patches = []
    for n in range( 0, len( imgs ) ):
        image = imgs[ n ]
        for i in range( 0, num_x_patches ):
            for j in range( 0, num_y_patches ):
                patches.append( image[ i * patch_width : (i+1) * patch_width,
                                      j * patch_height : (j+1) * patch_height ])#.astype(dtype='uint8') ) # All .astype comments can be removed for greater # of files (but ML is slower)
    return patches

import os
from skimage.util import img_as_uint
import tifffile as tif
from tifffile import imread
from matplotlib import pyplot as plt
import tensorflow as tf
import numpy as np

ground_truth_path = './training_data/ground_truth/'
train_scan_path = './training_data/scans/'


load_model_location = "./model/untrained" # POSITION WHERE THE MODEL IS TO SAVED
OUTPUT_DIR = "./model/trained"

# Reading ground truth filenames -----------------------------------------------------------------------------------------------

ground_truth_filenames = [x for x in os.listdir( ground_truth_path ) if x.endswith(".tif")]

print( 'Ground Truth Images loaded: ' + str( len(ground_truth_filenames)) )

# Reading training scan filenames ----------------------------------------------------------------------------------------------

train_scan_filenames = [x for x in os.listdir( train_scan_path ) if x.endswith(".tif")]

print( 'gVXR Scans loaded: ' + str( len(train_scan_filenames)) )

# Loading the model ------------------------------------------------------------------------------------------------------------

model = tf.keras.models.load_model(load_model_location)

#model = tf.keras.models.load_model(OUTPUT_DIR)
model.summary()

# Initial Parameters -----------------------------------------------------------------------------------------------------------

numEpochs = 5 # i.e. number of times the training data set has had the chance to update the model parameters
earlystopper = EarlyStopping(patience=5, verbose=1, restore_best_weights=True) # stops the code if it doesn't get any better
tiffs_per_train = 50 
files = 0 # Do not change value

train_width = 2000 # change depedning on image size
train_height = 2000 # change depedning on image size

test_width = 2000 # change depedning on image size
test_height = 2000 # change depedning on image size

x_patches = 5 # number of patches image is divided into in x direction
y_patches = 5 # number of patches image is divided into in y direction

patch_train_width = int(train_width/x_patches)
patch_train_height = int(train_height/y_patches)

patch_test_width = int(test_width/x_patches)
patch_test_height = int(test_height/y_patches)

# Training the model -----------------------------------------------------------------------------------------------------------


tiff_number = 0
while tiff_number <len(ground_truth_filenames):

    ground_truth_img = [ img_as_uint( imread( ground_truth_path + x ) ) for x in ground_truth_filenames[tiff_number:tiff_number+tiffs_per_train] ]
    train_scan_img = [ img_as_uint( imread( train_scan_path + x ) ) for x in train_scan_filenames[tiff_number:tiff_number+tiffs_per_train] ]

    ground_truth_patches = create_patches(ground_truth_img, x_patches, y_patches)
    train_scan_patches = create_patches(train_scan_img, x_patches, y_patches)

    print("Done 1")
    input_shape = ( patch_train_width, patch_train_height, 1 )
    X_train = [np.reshape(x, input_shape ) for x in ground_truth_patches]
    X_train = np.asarray(X_train)/65535

    print("Done 2")
    # training ground truth
    output_shape = ( patch_train_width, patch_train_height, 1 )
    Y_train = [x/65535 for x in train_scan_patches] # normalize between 0 and 1
    Y_train = [np.reshape( x, output_shape ) for x in Y_train] # Clear images
    Y_train = np.asarray(Y_train)

    print("Done 3")
    history = model.fit(Y_train, X_train, validation_split=0.1, batch_size = 50, #NOTE, THIS MAY NEED CHANGING BASED ON SIZE OF TIFFS
                            epochs=numEpochs, callbacks=[earlystopper])

    tiff_number += tiffs_per_train

print("Done 4")
model.save(OUTPUT_DIR)


"""We can now plot the loss and MAE curves for the training and validation sets."""

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(20, 8))

# summarize history for loss
ax1.plot(np.arange(1, (numEpochs+1), 1), history.history['loss'])
ax1.plot(np.arange(1, (numEpochs+1), 1), history.history['val_loss'])
ax1.set_title('Model Mean-squared Error', fontsize=25)
ax1.set_ylabel('MSE', fontsize=25)
ax1.set_xlabel('Epoch', fontsize=25)
ax1.tick_params(labelsize=18)
ax1.legend(['Training data', 'Validation data'], loc='upper right', fontsize=16)

# summarize history for MAE
ax2.plot([1,2,3,4,5], history.history['mean_absolute_error'])
ax2.plot([1,2,3,4,5], history.history['val_mean_absolute_error'])
ax2.set_title('Model Mean Absolute Error', fontsize=25)
ax2.set_ylabel('MAE', fontsize=25)
ax2.set_xlabel('Epoch', fontsize=25)
ax2.tick_params(labelsize=18)
ax2.legend(['Training data', 'Validation data'], loc='upper right', fontsize=16)

plt.show()
fig.savefig('loss.png', dpi=500)


