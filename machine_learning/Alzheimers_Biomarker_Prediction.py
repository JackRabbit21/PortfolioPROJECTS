#!/usr/bin/env python
# coding: utf-8

# In[18]:


pip install tqdm


# In[35]:


pip install numpy


# In[2]:


pip install tensorflow


# In[1]:


import pandas as pd
import os
from skimage.transform import resize
from skimage.io import imread
import numpy as np
import matplotlib.pyplot as plt


# In[2]:


Categories=['BIOMARKERS','CONTROL','NO BIOMARKERS']
flat_data_arr=[] #input array
target_arr=[] #output array
datadir='C:\\Users\\Kevin\\unecessary stuff\\Desktop\\DATA\\ADNI' 
#path which contains all the categories of images


# In[3]:


for i in Categories:
    
    print(f'loading... category : {i}')
    path=os.path.join(datadir,i)
    for img in os.listdir(path):
        img_array=imread(os.path.join(path,img))
        img_resized=resize(img_array,(150,150,3))
        flat_data_arr.append(img_resized.flatten())
        target_arr.append(Categories.index(i))
    print(f'loaded category:{i} successfully')


# In[4]:


flat_data=np.array(flat_data_arr)
target=np.array(target_arr)
df=pd.DataFrame(flat_data) #dataframe
df['Target']=target
x=df.iloc[:,:-1] #input data 
y=df.iloc[:,-1] #output data


# In[8]:


from sklearn import svm
from sklearn.model_selection import GridSearchCV
param_grid={'C':[0.1,1,10,100],'gamma':[0.0001,0.001,0.1,1],'kernel':['rbf','poly']}
svc=svm.SVC(probability=True)
model=GridSearchCV(svc,param_grid)


# In[7]:


from sklearn.model_selection import train_test_split
x_train,x_test,y_train,y_test=train_test_split(x,y,test_size=0.20,random_state=77,stratify=y)
print('Splitted Successfully')
model.fit(x_train,y_train)
print('The Model is trained well with the given images')
# model.best_params_ contains the best parameters obtained from GridSearchCV


# In[8]:


from sklearn.metrics import accuracy_score

y_pred=model.predict(x_test)
print("The predicted Data is :")
print(y_pred)
print("The actual data is:")
print(np.array(y_test))
print(f"The model is {accuracy_score(y_pred,y_test)*100}% accurate")


# In[10]:


url=input('C:\\Users\\Kevin\\unecessary stuff\\Desktop\\ALZHEIMERS.png')
img=imread(url)
plt.imshow(img)
plt.show()
img_resize=resize(img,(150,150,3))
l=[img_resize.flatten()]
probability=model.predict_proba(l)
for ind,val in enumerate(Categories):
    print(f'{val} = {probability[0][ind]*100}%')
print("The predicted image is having : "+Categories[model.predict(l)[0]])

from matplotlib import pyplot as plt

img =plt.imread("C:\\Users\Kevin\\unecessary stuff\\Desktop\\SCAN Final.png")
plt.imshow(img)
plt.title("Axial Image")
plt.axis("off")
plt.show()

import numpy as np
import matplotlib.pyplot as plt
import cv2
 
%matplotlib inline
 
# Read in the image
image = cv2.imread("C:\\Users\Kevin\\unecessary stuff\\Desktop\\SCAN Final.png")
 
# Change color to RGB (from BGR)
image = cv2.cvtColor(image, cv2.COLOR_BGR2RGB)
 
plt.imshow(image)


pixel_vals = image.reshape((-1,3))
pixel_vals = np.float32(pixel_vals)
criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 100, 0.85)
k = 10
retval, labels, centers = cv2.kmeans(pixel_vals, k, None, criteria, 10, cv2.KMEANS_RANDOM_CENTERS)
centers = np.uint8(centers)
segmented_data = centers[labels.flatten()]
segmented_image = segmented_data.reshape((image.shape))
 
plt.imshow(segmented_image)

import numpy as np
import glob
import matplotlib.pyplot as plt
import skimage.io
import skimage.color
import skimage.filters

image = skimage.io.imread("C:\\Users\\Kevin\\unecessary stuff\\Desktop\\DATA\\ADNI")

fig, ax = plt.subplots()
plt.imshow(image)
plt.show()

# convert the image to grayscale
gray_image = skimage.color.rgb2gray(image)

# blur the image to denoise
blurred_image = skimage.filters.gaussian(gray_image, sigma=1.0)

fig, ax = plt.subplots()
plt.imshow(blurred_image, cmap='gray')
plt.show()

# create a mask based on the threshold
t = 0.8
binary_mask = blurred_image < t

fig, ax = plt.subplots()
plt.imshow(binary_mask, cmap='gray')
plt.show()

# use the binary_mask to select the "interesting" part of the image
selection = np.zeros_like(image)
selection[binary_mask] = image[binary_mask]

fig, ax = plt.subplots()
plt.imshow(selection)
plt.show()

from nipype.interfaces.ants import N4BiasFieldCorrection
n4 = N4BiasFieldCorrection()
n4.inputs.dimension = 3
n4.inputs.input_image = 'structural.nii'
n4.inputs.bspline_fitting_distance = 300
n4.inputs.shrink_factor = 3
n4.inputs.n_iterations = [50,50,30,20]
n4.cmdline

import cv2

# read the image
image = cv2.imread("C:\\Users\\Kevin\\unecessary stuff\\Desktop\\SCAN Final.png")
# convert the image to grayscale format
img_gray = cv2.cvtColor(image, cv2.COLOR_BGR2GRAY)
# apply binary thresholding

ret, thresh = cv2.threshold(img_gray, 150, 255, cv2.THRESH_BINARY)

# visualize the binary image

cv2.imshow('Binary image', thresh)

cv2.waitKey(0)

cv2.imwrite('image_thres1.jpg', thresh)

cv2.destroyAllWindows()

import cv2

# read the image
image = cv2.imread("C:\\Users\Kevin\\unecessary stuff\\Desktop\\SCAN Final.png")

# B, G, R channel splitting
blue, green, red = cv2.split(image)

# detect contours using blue channel and without thresholding
contours1, hierarchy1 = cv2.findContours(image=blue, mode=cv2.RETR_TREE, method=cv2.CHAIN_APPROX_NONE)
 
# draw contours on the original image
image_contour_blue = image.copy()
cv2.drawContours(image=image_contour_blue, contours=contours1, contourIdx=-1, color=(0, 255, 0), thickness=2, lineType=cv2.LINE_AA)

# see the results
cv2.imshow('Contour detection using blue channels only', image_contour_blue)
cv2.waitKey(0)
cv2.imwrite('blue_channel.jpg', image_contour_blue)
cv2.destroyAllWindows()


# detect contours using green channel and without thresholding
contours2, hierarchy2 = cv2.findContours(image=green, mode=cv2.RETR_TREE, method=cv2.CHAIN_APPROX_NONE)

# draw contours on the original image

image_contour_green = image.copy()
cv2.drawContours(image=image_contour_green, contours=contours2, contourIdx=-1, color=(0, 255, 0), thickness=2, lineType=cv2.LINE_AA)

# see the results

cv2.imshow('Contour detection using green channels only', image_contour_green)
cv2.waitKey(0)
cv2.imwrite('green_channel.jpg', image_contour_green)
cv2.destroyAllWindows()

# detect contours using red channel and without thresholding
contours3, hierarchy3 = cv2.findContours(image=red, mode=cv2.RETR_TREE, method=cv2.CHAIN_APPROX_NONE)

# draw contours on the original image
image_contour_red = image.copy()
cv2.drawContours(image=image_contour_red, contours=contours3, contourIdx=-1, color=(0, 255, 0), thickness=2, lineType=cv2.LINE_AA)

# see the results
cv2.imshow('Contour detection using red channels only', image_contour_red)
cv2.waitKey(0)
cv2.imwrite('red_channel.jpg', image_contour_red)
cv2.destroyAllWindows()


# import opencv
import cv2

# Read image
src = cv2.imread("C:\\Users\\Kevin\\unecessary stuff\\Desktop\\SCAN.png", cv2.IMREAD_GRAYSCALE)

# Set threshold and maxValue
thresh = 0
maxValue = 255
 
# Basic threshold example
th, dst = cv2.threshold(src, thresh, maxValue, cv2.THRESH_BINARY);

# Applying Otsu's method setting the flag value into cv.THRESH_OTSU.

# Use a bimodal image as an input.

# Optimal threshold value is determined automatically.
otsu_threshold, image_result = cv2.threshold(

    image, 0, 255, cv2.THRESH_BINARY + cv2.THRESH_OTSU,

)

print("Obtained threshold: ", otsu_threshold)


# Python program to illustrate
# Otsu thresholding type on an image
  
# organizing imports
import cv2         
import numpy as np    
  
# path to input image is specified and
# image is loaded with imread command
image1 = cv2.imread("C:\\Users\\Kevin\\unecessary stuff\\Desktop\\SCAN.png")
  
# cv2.cvtColor is applied over the
# image input with applied parameters
# to convert the image in grayscale
img = cv2.cvtColor(image1, cv2.COLOR_BGR2GRAY)
  
# applying Otsu thresholding
# as an extra flag in binary 
# thresholding     
ret, thresh1 = cv2.threshold(img, 120, 255, cv2.THRESH_BINARY + 
                                            cv2.THRESH_OTSU)     
  
# the window showing output image         
# with the corresponding thresholding         
# techniques applied to the input image    
cv2.imshow('Otsu Threshold', thresh1)         
       
# De-allocate any associated memory usage         
if cv2.waitKey(0) & 0xff == 27:
    cv2.destroyAllWindows()     
    
    
batch_size = 64
epochs = 30
IMG_HEIGHT = 256
IMG_WIDTH = 256
NUM_LABELS = 3

single_classifier.compile(optimizer='adam',
                  loss="categorical_crossentropy",
                  metrics=['accuracy'])

single_classifier.summary()

# Importing the Keras libraries and packages
from keras.models import Sequential
from keras.layers import Conv2D
from keras.layers import MaxPooling2D
from keras.layers import Flatten
from keras.layers import Dense
from keras.callbacks import ModelCheckpoint
import tensorflow as tf
import numpy as np


# Initialising the CNN
classifier = Sequential()
# Step 1 - Convolution
classifier.add(Conv2D(32, (3, 3), input_shape = (64, 64, 3), activation = 'relu'))
# Step 2 - Pooling
classifier.add(MaxPooling2D(pool_size = (2, 2)))
# Adding a second convolutional layer
classifier.add(Conv2D(32, (3, 3), activation = 'relu'))
classifier.add(MaxPooling2D(pool_size = (2, 2)))
# Step 3 - Flattening
classifier.add(Flatten())
# Step 4 - Full connection
classifier.add(Dense(units = 128, activation = 'relu'))
classifier.add(Dense(units = 1, activation = 'sigmoid'))


# Compiling the CNN
classifier.compile(optimizer = 'adam', loss = 'binary_crossentropy', metrics = ['accuracy'])

# Part 2 - Fitting the CNN to the images
from keras.preprocessing.image import ImageDataGenerator
train_datagen = ImageDataGenerator(rescale = 1./255,
shear_range = 0.2,
zoom_range = 0.2,
horizontal_flip = True)
test_datagen = ImageDataGenerator(rescale = 1./255)
training_set = train_datagen.flow_from_directory('C:\\Users\\Kevin\\unecessary stuff\\Desktop\\DATA\\ADNI',
target_size = (64, 64),
batch_size = 32,
class_mode = 'binary')
test_set = test_datagen.flow_from_directory('C:\\Users\\Kevin\\unecessary stuff\\Desktop\\DATA\\ADNI',
target_size = (64, 64),
batch_size = 32,
class_mode = 'binary')
classifier.fit_generator(training_set,
steps_per_epoch = 8000,
epochs = 25,
validation_data = test_set,
validation_steps = 2000)

model = tf.keras.models.Sequential([tf.keras.layers.Dense(10)])
model.compile(tf.keras.optimizers.SGD(), loss='mse')
history = model.fit(np.arange(100).reshape(5, 20), np.zeros(5),
                    epochs=10)
print(history.params)

# check the keys of history object
print(history.history.keys())

# Part 3 - Making new predictions
import numpy as np
from keras.preprocessing import image
test_image = image.load_img('C:\\Users\\Kevin\\unecessary stuff\\Desktop\\SCAN FINAL.png', target_size = (64, 64))
test_image = image.img_to_array(test_image)
test_image = np.expand_dims(test_image, axis = 0)
result = classifier.predict(test_image)
training_set.class_indices
if result[0][0] == 1:
 prediction = 'ALZHEIMERS'
else:
 prediction = 'NO ALZHEIMERS'
