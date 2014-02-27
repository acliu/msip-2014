# !/bin/bash

boxFname=$1
boxLength=200
sliceNum=1
outputFname=$2

if [ -f args.dat ];
then
rm -f args.dat
fi
touch args.dat
echo $boxFname >> args.dat
echo $boxLength >> args.dat
echo $sliceNum >> args.dat
echo $outputFname >> args.dat

echo $boxFname
./outputslice.x