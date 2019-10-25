#!/bin/bash

if [ "$(uname)" == "Darwin" ]
then
   echo -framework Accelerate
else
   if [[ ! -z $(pkg-config --cflags --libs lapack 2>&-) ]]
   then
      echo $(pkg-config --cflags --libs lapack)
   else
      #echo -llapack -lblas -lg2c -lm
       echo -llapack
   fi
fi

exit
