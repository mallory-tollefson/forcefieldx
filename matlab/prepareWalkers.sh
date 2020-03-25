#!/bin/bash
file=$1
grep -e "Accept" -e "Reject" ${file} | awk -F " " '{print $1 " " $4}' >> walkers.txt
sed -i "s/\[//g" walkers.txt
sed -i "s/\]//g" walkers.txt
sed -i "s/L=//g" walkers.txt
sed -i "s/,//g" walkers.txt
