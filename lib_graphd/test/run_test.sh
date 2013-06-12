#!/bin/bash

for i in "testGraphReader" "testGraphCreatorFile"  "testGraph" "testWeightedGraph" "testVertexWeightedGraph" "testAdjMatrixGraphWriter" "testAdjMatrixGraphReader" "testDIMACSGraphWriter" "testMetisGraphWriter" "testMetisGraphReader" "testGraphVizGraphWriter" "testGraphProperty" "testGraphUtil" "testGraphEOUtil" "parmetis.sh"
do 
echo "======== Starting $i ==================="
#cd bin
if [ ! -f bin/${i} ] ; then
    echo "ERROR:  bin/${i} does not exist - did it get made?"
else
    ./bin/$i 2>/dev/null
fi
#cd ..
echo "=========End $i=========================="
echo " "
echo " "
sleep 2

done
