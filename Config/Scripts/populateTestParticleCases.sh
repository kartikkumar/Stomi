#!/bin/bash

# Declare array of strings with cases for test particle database generator configuration files.
declare -a testParticleDatabaseGeneratorCases=("case1_large_sma" "case1" "case2" "case3")

# Set template command to populate database with test particle cases.
declare -r testParticleTemplateCommand="stomi ../TestParticleDatabaseGenerator/CASE_testParticleDatabaseGeneratorSettings.cfg"
# declare -r testParticleTemplateCommand="echo $PATH"

# Loop through test particle cases and populate database.
for case in "${testParticleDatabaseGeneratorCases[@]}" 
do
    # Execute command to populate database.
    command="${testParticleTemplateCommand/CASE/$case}"
    $command

    # Wait for previous process to end before proceeding.
    wait
done
