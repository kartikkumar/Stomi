'''   
Copyright (c) 2010-2014, Delft University of Technology
Copyright (c) 2010-2014, K. Kumar (me@kartikkumar.com)
All rights reserved.
See http://bit.ly/12SHPLR for license details.

This script copies the semi-major axis values for test particle Case 1 to Cases 2 and 3.
'''

###################################################################################################
# Set up input deck
###################################################################################################

# Set absolute path to SQLite database with simulation data.
databasePath                = "/Users/kartikkumar/Desktop/mabPaper.sqlite"

# Set test particle case name for Case 1.
testParticleCase1Name       = "circular_equatorial"

# Set test particle case name for Case 2.
testParticleCase2Name       = "eccentric_equatorial"

# Set test particle case name for Case 3.
testParticleCase3Name       = "circular_inclined"

###################################################################################################



'''
    
                            DO NOT EDIT PARAMETERS BEYOND THIS POINT!!!

'''


###################################################################################################
# Import Python packages and modules
###################################################################################################

# Import necessary external packages.
import sqlite3
import time

###################################################################################################


###################################################################################################
# Start timer
###################################################################################################

startTime = time.time()

###################################################################################################


###################################################################################################
# Update database.
###################################################################################################

# Connect to SQLite database.
database = sqlite3.connect(databasePath)

# Enter database and update test particle cases.
with database:
    # Set cursor to scan through database and execute queries.
    cursor = database.cursor()
    
    # Get IDs cases.
    cursor.execute("SELECT testParticleCaseId FROM test_particle_case WHERE "
                   + "testParticleCaseName = \"" + testParticleCase1Name + "\"")
    testParticleCase1Id = cursor.fetchall()[0][0]  

    cursor.execute("SELECT testParticleCaseId FROM test_particle_case WHERE "
                   + "testParticleCaseName = \"" + testParticleCase2Name + "\"")
    testParticleCase2Id = cursor.fetchall()[0][0]  

    cursor.execute("SELECT testParticleCaseId FROM test_particle_case WHERE "
                   + "testParticleCaseName = \"" + testParticleCase3Name + "\"")
    testParticleCase3Id = cursor.fetchall()[0][0]  

    # Get test particle semi-major axis data for each case.
    cursor.execute("SELECT semiMajorAxis FROM test_particle_input WHERE testParticleCaseId = "
                   + str(testParticleCase1Id))
    rawCase1Data = cursor.fetchall()

    cursor.execute("SELECT testParticleSimulationId, semiMajorAxis FROM test_particle_input WHERE "
                   + "testParticleCaseId = " + str(testParticleCase2Id))
    rawCase2Data = cursor.fetchall()

    cursor.execute("SELECT testParticleSimulationId, semiMajorAxis FROM test_particle_input WHERE "
                   + "testParticleCaseId = " + str(testParticleCase3Id))
    rawCase3Data = cursor.fetchall()

    # Copy semi-major axis values for Case 1 to Cases 2 and 3 and update in database.
    for i in xrange(0, len(rawCase1Data)):
        cursor.execute("UPDATE test_particle_input SET semiMajorAxis = " + str(rawCase1Data[i][0])
                       + " WHERE testParticleSimulationId = " + str(rawCase2Data[i][1]))
        cursor.execute("UPDATE test_particle_input SET semiMajorAxis = " + str(rawCase1Data[i][0])
                       + " WHERE testParticleSimulationId = " + str(rawCase3Data[i][1]))        

###################################################################################################
# Finalize timer and print elapsed time.
###################################################################################################

# Finalize timer.
endTime = time.time()

# Print elapsed time for script [s].
print "This script took " + str(endTime - startTime) + "s"

###################################################################################################