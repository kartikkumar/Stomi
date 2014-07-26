#!/bin/bash

# Set environment variables.
stomi="/home/kartik/app/stomi/bin/stomi"

# Set template command to populate database with random walk runs.
declare -r randomWalkTemplateCommand="$stomi ../RandomWalkDatabaseGenerator/runID_randomWalkDatabaseGeneratorSettings.cfg"

# Loop through randoom walk runs and populate database.
for runId in {1..15}
do
    # Execute command to populate database.
    command="${randomWalkTemplateCommand/ID/$runId}"
    echo  $command

    # Wait for previous process to end before proceeding.
    wait
done
