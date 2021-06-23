#!/bin/bash

echo "this test script needs testing!"
# create input file
printf "DimSize = 30 30 30 \n\
Offset =      0    0    0 \n\
ElementSize = 1   1   1 \n\
Unit =      1 \n\
replaceRange 0 255 1 \n\
Paint   cylinder 0  15 15   20 15 15  10 0; //< point1(on-axis)  point2  radius value=0 water \n\
PaintAdd  sphere 15 15 15   12  2; //< x_centre y_centre z_centre  radius value=2(oil) \n\
replaceRange    3 255  1; //< range toValue, set outside of cylinder as rock=1 \n\
write dump.tif \n\
" > dropInCylinder.mhd

rm -rf dropInCylinder/
AllRunContAngle

echo return $?
