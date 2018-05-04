#!/bin/bash
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/misc/extra/data/guest01/mitsuba-af602c6fd98a/dist
#alias mitsuba='/net/voxel04/misc/extra/data/guest01/mitsuba-af602c6fd98a/dist/mitsuba'
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
	loc1=./vol/$1/density$x.vol
	loc2=./vol/$1/color/$x.vol
	sed "s/density.vol/.\/vol\/$1\/density\/$x.vol/g" ./render/hetvol_$1.xml > ./render/temp_$1.xml
	sed "s/color.vol/.\/vol\/$1\/color\/$x.vol/g" ./render/temp_$1.xml > ./render/temp1_$1.xml
	mitsuba ./render/temp1_$1.xml -o ./images/$1/exr/$x
	convert -colorspace RGB -colorspace sRGB ./images/$1/exr/$x.exr ./images/$1/png/$x.png
done
convert -delay $3 -loop 0 ./images/$1/png/*.png ./images/$1/$1.gif
rm ./render/temp_$1.xml 2> /dev/null
rm ./render/temp1_$1.xml 2> /dev/null
