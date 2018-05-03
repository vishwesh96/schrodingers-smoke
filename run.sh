#!/bin/bash
mkdir images 2> /dev/null 
mkdir vol 2> /dev/null
rm -rf ./vol/$1/ 2> /dev/null
rm -rf ./images/$1/ 2> /dev/null
mkdir ./vol/$1
mkdir ./vol/$1/density
mkdir ./vol/$1/color
mkdir ./images/$1 
mkdir ./images/$1/exr
mkdir ./images/$1/png
./Smoke $1 $2
for (( x=0; x<$2; x++ ))
do
	xmlstarlet ed -u '/scene/medium/volume[@name="density"]/string/@value' -v "./vol/$1/density/$x.vol" ./render/hetvol.xml > ./render/temp.xml
	xmlstarlet ed -u '/scene/medium/volume[@name="albedo"]/string/@value' -v "./vol/$1/color/$x.vol" ./render/temp.xml > ./render/temp1.xml
	mitsuba ./render/temp1.xml -o ./images/$1/exr/$x
	convert -colorspace RGB -colorspace sRGB ./images/$1/exr/$x.exr ./images/$1/png/$x.png
done
convert -delay $3 -loop 0 ./images/$1/png/*.png ./images/$1/$1.gif
rm ./render/temp.xml 2> /dev/null
rm ./render/temp1.xml 2> /dev/null