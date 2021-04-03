#!/bin/bash

voxelWidth=0.046
samplingStep=0.0023
upperBound=1e10
lowerBound=0
K=50

pathToSrc="put your path here"
filename=".h5 file of the reconstructed tree"
measureFilename=".csv file for storing the measure"
dirname_originalVolume="path to the folder containing the original volumes"
dirname_noisyVolume="path to the folder containing the noisy volumes"

for num in {001..015} ; do
for thresholdBelow in 0.005 0.01 0.012 0.014 0.016 0.018 `seq 0.02 0.02 0.20` `seq 0.024 0.004 0.036` ; do
# thresholdBelow is a threshold to select voxels with Frangi vesselness measure larger than this threshold 

dirname_reconstructed="$thresholdBelow/$filename"
python "$pathToSrc/evaluation_metrics.py" "doComputeOverlapMeasure" "$dirname_originalVolume" "$dirname_noisyVolume" "image$num" "$dirname_reconstructed" $voxelWidth $samplingStep $upperBound $lowerBound $K --points="positions" --doOutputHeader --prependHeaderStr="ImageName,ThresholdValue," --prependRowStr="image$num,$thresholdBelow," > "image$num/$thresholdBelow/$measureFilename"

done
done

outputFilename=$measureFilename
# concatenated file is stored in the same folder as the script
cat "./image001/0.005/$measureFilename" | head -n 1 > "./$outputFilename"

for num in {001..015} ; do
for thresholdBelow in 0.005 0.01 0.012 0.014 0.016 0.018 `seq 0.02 0.02 0.20` `seq 0.024 0.004 0.036` ; do

if [ -s "./image$num/$thresholdBelow/$measureFilename" ]; then
cat "./image$num/$thresholdBelow/$measureFilename"  | tail -n+2 >> "./$outputFilename"
else
echo "./image$num/$thresholdBelow/$measureFilename"
fi
done
done

python "$pathToSrc/evaluation_metrics.py" "doAnalyzeOurOverlapMeasureCsv" "." "$outputFilename" $voxelWidth
