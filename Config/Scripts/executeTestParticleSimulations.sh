#!/bin/bash

# Set environment variables.
stomi="/home/kartik/app/stomi/bin/stomi"

# Declare array of strings with cases for test particle simulator configuration files.
declare -a testParticleSimulatorCases=("case1_large_sma" "case1" "case2" "case3")

# Set template command to execute test particle cases.
declare -r testParticleTemplateCommand="$stomi ../TestParticleSimulator/CASE_testParticleSimulatorSettings.cfg"

# Loop through test particle simulation cases and execute.
for case in "${testParticleSimulatorCases[@]}" 
do
    # Execute case.
    command="${testParticleTemplateCommand/CASE/$case}"
    $command

    # Wait for previous process to end before proceeding.
    wait
done