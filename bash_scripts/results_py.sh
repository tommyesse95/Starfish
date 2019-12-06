#!/bin/bash

cd 
cd Starfish/discharge 
rm output.txt
rm data_sparc.csv 
cd 
scp -r 192.168.250.68:/data/tommyesse/Starfish/discharge/output.txt Starfish/discharge
scp -r 192.168.250.68:/data/tommyesse/Starfish/discharge/data_sparc.csv Starfish/discharge
