#!bin/bash

if [ -z "$(ls -A heat-dat-files)" ]; then
   echo "heat-dat-files is Empty"
else
   rm heat-dat-files/Heat*
fi