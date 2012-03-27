#!/bin/bash

cd tests

# test Coils
echo "test Coils"
pscoils -f ./ferritin.fasta -l F > ./ferritin.Coils
# test PCoils
echo "test PSCoils"
pscoils -p ./ferritin.prof -l F > ./ferritin.PCoils
# test PSCoils
echo "test PSCoils"
pscoils -p ./ferritin.prof -f ./ferritin.fasta -l F > ./ferritin.PSCoils

cd -
