#!/bin/bash
NAME=$1
cd /home/dspath/Downloads/flow/flowstar-2.1.0/
./flowstar < /home/dspath/Documents/Dsgit/Iot\ Replan\ Simulations/replan.model  > /dev/null
cp ./outputs/$NAME.plt /home/dspath/Documents/Dsgit/Iot\ Replan\ Simulations/octagon.txt

