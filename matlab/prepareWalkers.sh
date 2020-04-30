#!/bin/bash
# This method generates two walker files for use with the 
# forcefieldx/matlab/Walkers.m and forcefieldx/matlab/Walkers-Mallory.m
# matlab files. The walkers-mallory.txt that is generated works with
# Walkers-Mallory, and the walkers.txt works with Walkers.m.
# To run the script, use the following command:
# ./prepareWalkers.txt amber-nomeld.5997651.log 8 true
# Where the first argument is the MC-OST log file, the second
# argument is the number of walkers, and the third argument is a 
# boolean indicating if you want to generate the walkers.txt file
# or not. 

file=$1
numWalkers=$2
walkers=$3

grep -e "Accept" -e "Reject" ${file} | awk -F " " '{print $1 " " $4}' >> walkers-mallory.txt
sed -i "s/\[//g" walkers-mallory.txt
sed -i "s/\]//g" walkers-mallory.txt
sed -i "s/L=//g" walkers-mallory.txt
sed -i "s/,//g" walkers-mallory.txt

if $walkers
then
	#Grep the lambda values for each walker into individual files.
	numWalkers=$((numWalkers-1))
	for i in $(seq 0 $numWalkers)
	do
		grep "${i} " walkers-mallory.txt | awk -F " " '{print $2}'  >> ${i}.txt
		lengthArray[$i]=$(wc ${i}.txt | cut -d " " -f 2)
	done

	#Find the walker that sampled the smallest number of lambdas (aka, slowest walker).
	min=0
	for i in ${lengthArray[@]}; do
    		(( $i < min || min == 0)) && min=$i
	done

	#Print out the "min" number of lambdas for each walker (all columns have to be 
	# the same length for the matlab script).
	for i in $(seq 0 $numWalkers)
	do
		tail -n ${min} ${i}.txt >> ${i}-min.txt
	done 

	#Set up a string to concatenate all of the walkers lambda columns into one file.
	string=""
	for i in $(seq 0 $numWalkers)
	do
		string="${string} ${i}-min.txt"
	done

	#Paste the lambda columns into one file.
	paste ${string} | column -s $'\t' -t >> walkers.txt

	#Clean up generated files.
	for i in $(seq 0 $numWalkers)
	do
		rm ${i}.txt
		rm ${i}-min.txt
	done
fi
