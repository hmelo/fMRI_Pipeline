#!/bin/bash

for item in $@;do
	SubNo=${item:5:4}; 
	echo $SubNo
done