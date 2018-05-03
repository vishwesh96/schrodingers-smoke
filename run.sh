#!/bin/bash
mkdir images 2> /dev/null 
mkdir vol 2> /dev/null
rm -rf ./vol/$1/ 2> /dev/null
rm -rf ./images/$1/ 2> /dev/null
mkdir ./vol/$1
mkdir ./images/$1 
mkdir ./images/$1/exr
mkdir ./images/$1/png
./Smoke $1 $2
for (( x=0; x<$2; x++ ))
do
	xmlstarlet ed -u '/scene/medium/volume/string/@value' -v "./vol/$1/$x.vol" ./render/hetvol.xml > temp.xml
	mitsuba temp.xml -o ./images/$1/exr/$x
	convert -colorspace RGB -colorspace sRGB ./images/$1/exr/$x.exr ./images/$1/png/$x.png
done
convert -delay $3 -loop 0 ./images/$1/png/*.png ./images/$1/$1.gif
rm temp.xml 2> /dev/null