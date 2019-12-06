#!/bin/bash

cd 
cd Starfish/cxt1
rm -r results
cd 
scp -r 192.168.250.68:/data/tommyesse/Starfish/cxt1/results Starfish/cxt1
