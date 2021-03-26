#!/bin/bash
endian=$1
if [ $endian == 1 ]
then
    echo "In endian conversion"
    for i in {0..3000..100}
    do
	./wrongendian "output_$i.dat" "outputnew_$i.dat"
    done
fi
if [ $endian == 1 ]
then
    for file in *new*
    do
        ./mmsp2txt $file
    done
else
    for file in *.dat
    do
        ./mmsp2txt $file
    done
fi






		

