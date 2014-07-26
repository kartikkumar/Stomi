#!/bin/bash

# Set environment variables.
stomi="/home/kartik/app/stomi/bin/stomi"

# Set template command to execute random walk runs.
declare -r randomWalkTemplateCommand="$stomi ../RandomWalkSimulator/runID_randomWalkSimulatorSettings.cfg"

# Loop through randoom walk runs and execute.
for runId in {1..15}
do
    # Execute run.
    command="${randomWalkTemplateCommand/ID/$runId}"
    $command

    # Wait for previous process to end before proceeding.
    wait
done
