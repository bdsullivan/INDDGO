#!/bin/bash

for i in "testDimacsReader" "testGraphReader" "testGraphCreatorFile" "testGraphReaderWriterFactory" "testGraph" "testWeightedGraph" "testMutableGraph" "testWeightedMutableGraph" "testAdjMatrixGraphWriter" "testAdjMatrixGraphReader" "testDIMACSGraphWriter" "testMetisGraphWriter" "testMetisGraphReader" "testGraphVizGraphWriter" "testGraphProperty" "testGraphUtil" "testGraphEOUtil" "parmetis.sh"
do 
echo "======== Starting $i ==================="
#cd bin
./bin/$i 2>/dev/null
#cd ..
echo "=========End $i=========================="
echo " "
echo " "
sleep 2

done
