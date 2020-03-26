#!/bin/bash
# This method generates two walker files for use with the 
# forcefieldx/matlab/Walkers.m and forcefieldx/matlab/Walkers-Mallory.m
# matlab files. The new-walkers.txt that is generated works with
# Walkers-Mallory, and the old-walkers.txt works with Walkers.m.
# To run the script, use the following command:
# ./prepareWalkers.txt amber-nomeld.5997651.log 8 true
# Where the first argument is the MC-OST log file, the second
# argument is the number of walkers, and the third argument is a 
# boolean indicating if you want to generate the old-walkers.txt file
# or not. 

file=$1
numWalkers=$2
oldWalkers=$3

grep -e "Accept" -e "Reject" ${file} | awk -F " " '{print $1 " " $4}' >> new-walkers.txt
sed -i "s/\[//g" new-walkers.txt
sed -i "s/\]//g" new-walkers.txt
sed -i "s/L=//g" new-walkers.txt
sed -i "s/,//g" new-walkers.txt

if [$oldWalkers]
then
	numWalkers=$((numWalkers-1))
	for i in $(seq 0 $numWalkers)
	do
		grep "${i} " new-walkers.txt | awk -F " " '{print $2}'  >> ${i}.txt
	done

	string=""
	for i in $(seq 0 $numWalkers)
	do
		string="${string} ${i}.txt"
	done

	paste ${string} | column -s $'\t' -t >> old-walkers.txt

	for i in $(seq 0 $numWalkers)
	do
		rm ${i}.txt
	done
fi
