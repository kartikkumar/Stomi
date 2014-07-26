#!/bin/bash

# Set environment variables.
stomi="/home/kartik/app/stomi/bin/stomi"

# Declare array of strings with variants for random walk database generator configuration files.
declare -a randomWalkDatabaseGeneratorVariants=("nominal" "light" "heavy" "sparse" "dense")

# Set template command to populate database with random walk runs.
declare -r randomWalkTemplateCommand="$stomi ../RandomWalkDatabaseGenerator/RUN_VARIANT_randomWalkDatabaseGeneratorSettings.cfg"

# Loop through randoom walk runs and populate database.
for run in {1..15}
do
    for variant in "${randomWalkDatabaseGeneratorVariants[@]}" 
    do
        # Execute command to populate database.
        temporaryCommand="${randomWalkTemplateCommand/RUN/$run}"
        command="${temporaryCommand/VARIANT/$variant}"
        $command

        # Wait for previous process to end before proceeding.
        wait
    done
done
