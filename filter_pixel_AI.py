# Description: AI filtering
# Author: genger
# Data: 2023-07-05

import cv2 as cv
import numpy as np
import pandas as pd

import getopt
import sys

## usage
def usage():
    print("")
    print("AI filtering your image.")
    print("uage: python %s -option <argument>" %sys.argv[0])
    print(" -h/--help ")
    print(" --img=<STRING> input image file.")
    print(" --barcodeNum=<INTEGER> num of your barcode.")
    print(" --thre=<INTEGER> threshold for filtering.")

## deal with options
try: 
    opts, args = getopt.getopt(sys.argv[1:], "h", ["help","img=","barcodeNum=","thre="])
except getopt.GetoptError:
    print("ERROR: Get option error.\nYou can contact the author through wechat 13958598285.")
    usage()
    sys.exit(2)

for opt, val in opts:
    if opt in ("-h", "--help"):
        usage()
        sys.exit(1)
    else:
        if opt in ("--img"):
            img_path = val 
        if opt in ("--barcodeNum"):
            num = int(val)
        if opt in ("--thre"):
            threshold = float(val)

try:
    print("Your img is:",img_path)
    print("Your barcodeNum is:",num)
    print("Your thre is:",threshold)
except:
    print("ERROR: Missing option.")
    usage()
    sys.exit(2)
    
#read image
img = cv.imread(img_path, 0)

#filter non-tissue part
blur = cv.GaussianBlur(img, (5,5), 0)
ret3, th3 = cv.threshold(blur, 0, 255, cv.THRESH_BINARY+cv.THRESH_OTSU)
numRows,numCols = th3.shape

pixel_count = 2*num-1
pixel_w = numCols/pixel_count
pixel_h = numRows/pixel_count

treds = np.percentile(blur, threshold)
filter_pixel = ""

#make position
for i in range(1, num+1):
    y = round(2*(i-1)*pixel_h+1)
    for j in range(1, num+1):
        x = round(2*(j-1)*pixel_w+1)
        pixel = th3[y:round(y+pixel_h-1), x:round(x+pixel_w-1)]
        C = np.mean(pixel)
        if C >= treds:
            filter_pixel = filter_pixel+","+str(j)+"x"+str(i)

# save position file
data = pd.DataFrame(filter_pixel.split(",")[1:])
data.to_csv("~/STvis.AI.filter.genger.filterd.pixels.csv", index = False, header = False)

