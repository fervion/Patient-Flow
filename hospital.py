# -*- coding: utf-8 -*-
"""
Simulation of patients in hospital
Name: Hospital.py
Author: Taige Hou (thou1@partners.org)
"""

"""
Contains classes for the layout of the hospital
"""
from room import *
from patient import *
from copy import deepcopy
from random import random
from pfUtil import *

#Class for defining the parameters of the hospital
class HospitalSettings(object):
    def __init__(self, matchingProfile = [True, True]):
        self.matchingProfile = matchingProfile[:]#matches patients to bedmates based on observed CP status (MRSA, VRE)
        self.currInputs = None
        self.pfSim = None
class Hospital(object):
    elapsedTime = 0
    #imports a layout file containing the layout of the hospital. 
    #A layout is stored as a dictionary matching room ID's to room objects
    def __init__(self):
        self.layout={}
        self.startRoomID = None
        self.endRoomID = None
        self.settings = HospitalSettings()
        self.blockBedMovement = False

        #allPatients is a dictionary mapping ids to patients of all patients that have ever existed (including dead and discharged)
        self.allPatientsInHistory = {}
        #patients that are still in the hospital
        self.patientsInHospital = {}
    #adds list of patients to start room
    def add_patients_to_hospital(self, patientList, t):
        for patient in patientList:
            self.allPatientsInHistory[patient.id] = patient
            self.patientsInHospital[patient.id] = patient
            patient.timeEnteringHospital = t
            patient.timeOfLastAcuityChange = t
        self.layout[self.startRoomID].add_patients(patientList)
        
        self.write_trace("\n".join([str(patient) for patient in patientList]), level = 2)
    #imports Layout from a list of rooms with connections in place,
    def import_layout(self, roomList):
        self.layout = dict([(room.id, room) for room in roomList])
        for room in self.layout.values():
            room.initialize_room()
        
    #validates layout and puts all patients in the start room
    def initialize_hospital(self, initialPatients = []):
        self.startRoomID, self.endRoomID = self.validate_layout(self.layout)
    #initialize certain outputs for this time step (called at start of each time step)
    def initialize_time_step_outputs(self, simOutputs):
        #initialize curr time step outputs
        for disease in DISEASES:
            simOutputs.numTestsGiven[disease].append(0)
            simOutputs.numTestsGivenUnique[disease].append(0)
            simOutputs.numDiffObservedCPDischarge[disease].append(0)
            simOutputs.numFlagsAquired[disease].append(0)
            simOutputs.numFlagsLost[disease].append(0)
        for combo in product((MRSA, VRE),(True, False),(True, False)):
            simOutputs.testResults[combo].append(0)
        for reason in (DEAD, DISCHARGE_HOME, DISCHARGE_FACILITY):
            simOutputs.numPatientsLeaving[reason].append(0)
        for acuity in ACUITY_VALUES:
            simOutputs.numPeopleToAcuityFive[acuity].append(0)
            simOutputs.numPeopleFromAcuityThree[acuity].append(0)
        simOutputs.numBedMovesAcuity.append(0)
        simOutputs.numBedMovesCP.append(0)

    #steps a time interval in the model
    #t is the time step, currtime is a datetime object
    def time_step(self, t, currTime, simOutputs):
        self.write_trace("Time Step: {}".format(t))
        self.initialize_time_step_outputs(simOutputs)
        moves=[]
        simInputs = self.settings.currInputs

        #stop any blocks
        for doBlock, startTime, endTime in simInputs.blockBedMovements:
            if doBlock and endTime.hour == currTime.hour:
                self.blockBedMovement = False
        #start any blocks
        for doBlock, startTime, endTime in simInputs.blockBedMovements:
            if doBlock and startTime.hour == currTime.hour:
                self.blockBedMovement = True
        
        #change acuities and update death/discharge (acuity are only changed if it is time to change acuities)
        self.update_acuity(t, currTime)

        #update patient parameters before time step
        for patient in self.patientsInHospital.values():
            #update observed cp for people who's test results have returned
            indicesToDelete = []

            for i,test in enumerate(patient.waitingTests):
                disease, result, timeToResult = test
                if timeToResult == t:
                    indicesToDelete.append(i)
                    patient.testResults[disease].append(result)
                    self.write_trace("Patient {} received test result for {} as {}".format(patient.id, DISEASE_STRINGS[disease], result, level = 2))
            #delete completed results
            patient.waitingTests = [test for i,test in enumerate(patient.waitingTests) if i not in indicesToDelete]
            #change observed status if applicable
            for disease in (MRSA, VRE):
                if disease == MRSA:
                    testInputs = simInputs.testMRSA
                else:
                    testInputs = simInputs.testVRE
                resultHist = patient.testResults[disease]
                if resultHist:
                    if patient.observedCP[disease]:
                        #check for deflagging
                        numTestsRequired = testInputs.numTestsForDeflag
                        if numTestsRequired > 0 and len(resultHist) >= numTestsRequired and not any(resultHist[-1*numTestsRequired:]):
                            patient.change_observed(disease, False)
                            simOutputs.numFlagsLost[disease][-1]+=1
                    else:
                        if resultHist[-1]:
                            #true result only takes one test to flag
                            patient.change_observed(disease, True)
                            simOutputs.numFlagsAquired[disease][-1]+=1

        #add new people to the hospital
        if isWeekday(currTime):
            distrib = simInputs.numArrivingPatientsMeanWeekday
        else:
            distrib = simInputs.numArrivingPatientsMeanWeekend

        numNewPatients = np.random.poisson(distrib[currTime.hour])

        newPatients = []
        numUnflaggedCarrier = [0,0]
        numFlaggedNonCarrier = [0,0]
        for i in range(numNewPatients):
            gender, service, pstatusObs, acuity = draw_from_distribution(simInputs.arrivingPatientDistribution)
            obsvMRSA, obsvVRE = convertPStatusToCP(pstatusObs)
            trueMRSA =  draw_from_distribution(simInputs.trueCPDistributions[MRSA][obsvMRSA])
            trueVRE = draw_from_distribution(simInputs.trueCPDistributions[VRE][obsvVRE])

            if obsvMRSA != trueMRSA:
                if trueMRSA == True:
                    numUnflaggedCarrier[MRSA]+=1
                else:
                    numFlaggedNonCarrier[MRSA]+=1
            if obsvVRE != trueVRE:
                if trueVRE == True:
                    numUnflaggedCarrier[VRE]+=1
                else:
                    numFlaggedNonCarrier[VRE]+=1
            workUpTime = draw_from_distribution(simInputs.workUpTime[acuity])
            newPatients.append(RealPatient(gender, acuity, service, (trueMRSA,trueVRE),convertPStatusToCP(pstatusObs), workUpTime))
        self.write_trace("Adding Patients to Hospital", level = 1)
        self.add_patients_to_hospital(newPatients, t)
        simOutputs.numNewPatients.append(len(newPatients))
        for disease in DISEASES:
            simOutputs.numUnflaggedCarrierAtAdmit[disease].append(numUnflaggedCarrier[disease])
            simOutputs.numFlaggedNonCarrierAtAdmit[disease].append(numFlaggedNonCarrier[disease])


        #roll for beds being blocked due to staffing
        for room in self.layout.values():
            if isinstance(room, BedRoom):
                for i in range(len(room.beds)):
                    #check to see if beds should be unblocked
                    if room.bedsBlockedStaffing[i] and t>=room.timeToStopBlockingBedsStaffing[i]:
                        room.bedsBlockedStaffing[i] = False
                
                    #bed is not occupied or blocked
                    if not room.beds[i] and not room.bedsBlockedStaffing[i]:
                        #roll for chance to block off bed
                        if random() < simInputs.probClosedBed:
                            room.bedsBlockedStaffing[i] = True
                            lam = 1/float(simInputs.distribClosedBed)
                            room.timeToStopBlockingBedsStaffing[i] = t+int(np.random.exponential(lam)+.5)
                            
                    
        #time step for the room
        for room in self.layout.values():
            moves.extend(room.time_step(self.settings, t))
                        
        #move patients according to move tables
        for move in moves:
            patient, currRoomID, nextRoomID = move
            patientID = patient.id

            #dead or discharged/move to end room
            if nextRoomID is None:
                #remove cleaning patients
                if isinstance(patient, CleaningPatient):
                    self.layout[currRoomID].remove_patient(patient, removeCleaningPatient = True)
                    continue
                if not isinstance(patient, RealPatient):
                    self.layout[currRoomID].remove_patient(patient)
                    continue

                nextRoomID = self.endRoomID             
                #get reason for patient leaving
                simOutputs.numPatientsLeaving[patient.status][-1]+=1
                simOutputs.numPatientsLeavingHour[currTime.isoweekday()][currTime.hour]+=1
                
                #delete patient from inhospital
                del self.patientsInHospital[patientID]
                
            currRoom = self.layout[currRoomID]
            nextRoom = self.layout[nextRoomID]
            self.write_trace("{} {} moved from {} to {}".format(patient.prefix, patient.id, currRoom.name, nextRoom.name), level = 1)
                
            if isinstance(currRoom, BedRoom):
                bedIndex = currRoom.beds.index(patient)
 
            currRoom.move_patient_to(patientID, self.layout[nextRoomID])
            #add cleaning patient in bedrooms
            if isinstance(patient, RealPatient) and isinstance(currRoom, BedRoom):
                currRoom.add_patient(CleaningPatient(t+1), bedIndex)
                
        #resets occupancy numbers
        for room in self.layout.values():
            room.reset_room()

        numBedMovesShuffle = 0
        #do unit shufflings
        if not self.blockBedMovement and simInputs.shufflePeriod > 0 and t%simInputs.shufflePeriod==0:
            numBedMovesShuffle  = self.shuffle_units()

        simOutputs.numBedMovesShuffle.append(numBedMovesShuffle)

        #do transmission
        numNewCasesBackground = np.array([0,0])
        numNewCasesCohorting = np.array([0,0])
        for room in self.layout.values():
            newBackground, newCohorting = room.do_transmission(simInputs, self.write_trace)
            numNewCasesBackground += newBackground
            numNewCasesCohorting += newCohorting
        for disease in DISEASES:
            simOutputs.numNewInfectedBackground[disease].append(numNewCasesBackground[disease])
            simOutputs.numNewInfectedCohorting[disease].append(numNewCasesCohorting[disease])

        #resets occupancy numbers
        for room in self.layout.values():
            room.reset_room()

        for room in self.layout.values():
            room.finalize_room()
    #applies bedmoves to patients in each unit in order to free up closed beds
    def shuffle_units(self):
        def find_matching_patient(unfilledDoubleRooms, singlePatients):
            for unfilledRoom in unfilledDoubleRooms:
                currPatient = unfilledRoom.patients.values()[0]
                for singlePatient in singlePatients:
                    if singlePatient is not currPatient and unfilledRoom.patient_matches(singlePatient, self.settings.matchingProfile) and singlePatient.currRoom.unitName==unfilledRoom.unitName:
                        return singlePatient, unfilledRoom
            return None, None

        print "Shuffling Units"
        bedAllocationRooms = [room for room in self.layout.values() if isinstance(room, BedAllocationRoom)]
        numBedMoves = 0
        for baRoom in bedAllocationRooms:
            #first see if we can accomodate patients with 0 moves
            for waitingPatient in baRoom.patients.values():
                patientMatch = False
                for bedroom in baRoom.nextRooms:
                    if bedroom.patient_matches(waitingPatient, self.settings.matchingProfile):
                        patientMatch = True
                        baRoom.move_patient_to(waitingPatient.id, bedroom)
                        baRoom.reset_room()
                        bedroom.reset_room()
                        break
                        
                if patientMatch:
                    continue
            
            #skip shuffle for this unit if no one's waiting
            moveOccured = False
         
            
            #bed moves that only require a single move
            for waitingPatient in baRoom.patients.values():                
                singlePatients = [] #patients by themselves
                unfilledDoubleRooms = [] #double rooms with single occupancy
                for bedroom in baRoom.nextRooms:
                    if bedroom.occupancy == 1:
                        if not isinstance(bedroom.patients.values()[0], RealPatient):
                            continue
                        singlePatients.append(bedroom.patients.values()[0])
                        if bedroom.capacity == 2:
                            unfilledDoubleRooms.append(bedroom)
                patientToMove, unfilledRoom = find_matching_patient(unfilledDoubleRooms, singlePatients)
   
                if patientToMove:
                    numBedMoves+=1
                    movePatientRoom = patientToMove.currRoom
                    movePatientRoom.move_patient_to(patientToMove.id, unfilledRoom)
                    unfilledRoom.reset_room()
                    movePatientRoom.reset_room()
                    baRoom.move_patient_to(waitingPatient.id, movePatientRoom)
                    movePatientRoom.reset_room()
                    moveOccured=True
                            
                if not moveOccured:
                    break

            #bed moves that require a double move
            for waitingPatient in baRoom.patients.values():
                matchingSingleRoomPatient = None
                for bedroom in baRoom.nextRooms:
                    if bedroom.capacity == 1 and bedroom.occupancy == 1:
                        matchingPatient = bedroom.patients.values()[0]
                        if not isinstance(matchingPatient, RealPatient):
                            continue
                        if waitingPatient.matches(matchingPatient):
                            matchingSingleRoomPatient = matchingPatient
                            break
                if not matchingSingleRoomPatient:
                    continue
                for bedroom in baRoom.nextRooms:
                    numBlockedBeds = len([bed for bed in bedroom.bedsBlockedStaffing if bed])
                    if bedroom.capacity == 2 and bedroom.occupancy == 1 and numBlockedBeds == 0:
                        originalOccupent = bedroom.patients.values()[0]
                        if not isinstance(originalOccupent, RealPatient):
                            continue

                        numBedMoves += 2
                        singleRoom = matchingSingleRoomPatient.currRoom
                        singleRoom.move_patient_to(matchingSingleRoomPatient.id, bedroom)
                        singleRoom.reset_room()

                        bedroom.move_patient_to(originalOccupent.id, singleRoom)
                        bedroom.reset_room()
                        singleRoom.reset_room()

                        baRoom.move_patient_to(waitingPatient.id, bedroom)
                        bedroom.reset_room()
                        break


        return numBedMoves
    #draws a new acuity value for that person
    def draw_new_acuity(self, patient, t):
        simInputs = self.settings.currInputs
        fromAcuity = patient.acuity
        timeOnAcuity = t - patient.timeOfLastAcuityChange #time spent at this acuity
        #use threshold parameters
        if timeOnAcuity >= simInputs.acuityChangeThresholdTimes[fromAcuity]:
            distribution = simInputs.acuityChangeThresholdDistribution[fromAcuity]
        else: #use hourly table
            distribution = simInputs.acuityChangeDistribution[fromAcuity][timeOnAcuity]
            
        toAcuity = draw_from_distribution(distribution)
        return toAcuity
        
    #updates acuity for each person.  Also takes care of death and discharge
    def update_acuity(self, t, currTime):
        for patient in self.patientsInHospital.values():
            #freeze acuity for patients with a bed assigned
            if patient.realPatient.delayPatient:
                continue
            
            simInputs = self.settings.currInputs

            rollForLeave = True

            if simInputs.blockDischarge and simInputs.blockDischargeTimeSteps>=(t-patient.timeEnteringHospital):
                rollForLeave = False
                
            if isWeekday(currTime):
                probLeave = simInputs.probLeaveWeekday[patient.acuity][currTime.hour]
            else:
                probLeave = simInputs.probLeaveWeekend[patient.acuity][currTime.hour]


            #check if patient leaves hospital
            if rollForLeave and random() < probLeave:
                #get method of leaving
                if isWeekday(currTime):
                    distribution = simInputs.leaveDistributionsWeekday[patient.acuity][currTime.hour]
                else:
                    distribution = simInputs.leaveDistributionsWeekend[patient.acuity][currTime.hour]
                toAcuity = draw_from_distribution(distribution)
            else:
                toAcuity = self.draw_new_acuity(patient, t)

            timeSpentInHospital = t - patient.timeEnteringHospital
            isTimeToChangeAcuity = (timeSpentInHospital%simInputs.lengthOfStayToChangeAcuity)==0
            leavingHospital = toAcuity == ACUITY_DISCHARGE_HOME or toAcuity == ACUITY_DISCHARGE_FACILITY or toAcuity == ACUITY_DEATH 
            if leavingHospital:
                for disease in DISEASES:
                    if patient.observedCP[disease]!=patient.originalObservedCP[disease]:
                        self.settings.pfSim.outputs.numDiffObservedCPDischarge[disease][-1]+=1
            
            #change acuity (could be death/discharge)
            if isTimeToChangeAcuity or leavingHospital:
                if patient.acuity != toAcuity:
                    self.write_trace("Patient {} changed acuity to {}".format(patient.id, toAcuity),level = 1)
                patient.change_acuity(toAcuity, t)
            
    #returns the start room id from a layout
    def get_start_room_id_from_layout(self,layout):
        for roomID, room in layout.iteritems():
            if type(room).__name__ == "StartRoom":
                return roomID
        return None

    #returns the end room id from a layout
    def get_end_room_id_from_layout(self,layout):
        for roomID, room in layout.iteritems():
            if type(room).__name__ == "EndRoom":
                return roomID
        return None
                
    """
    Checks the layout is valid.
    Raises error if not valid.
    All layouts must contain one start room
    All Bed Allocation Rooms must be connected to Bed rooms
    All Bed Rooms must be owned by exactly one Bed Allocation Room
    """
    def validate_layout(self, layout):
        #All layouts must have exactly one start and end room
        startRoomID = self.get_start_room_id_from_layout(layout)
        if startRoomID is None:
            raise InvalidLayoutError

        endRoomID = self.get_end_room_id_from_layout(layout)
        if endRoomID is None:
            raise InvalidLayoutError
        
        allConnections = {}
        reverseConnections = dict([(rID, []) for rID in layout.keys()])
        for roomID, room in layout.iteritems():
            allConnections[roomID] = [r.id for r in room.nextRooms]
            for nextRoom in room.nextRooms:
                reverseConnections[nextRoom.id].append(roomID)

        #All Bed rooms must be owned by exacty one bed allocation room
        bedRoomIDs = [room.id for room in layout.values() if type(room).__name__=="BedRoom"]
        for bRoomID in bedRoomIDs:
            if len(reverseConnections[bRoomID])!= 1:
                raise InvalidLayoutError
            if type(layout[reverseConnections[bRoomID][0]]).__name__ != "BedAllocationRoom":
                raise InvalidLayoutError
        
        #Bed allocation Rooms can only connect to Bed rooms
        bedAllocationRooms = [room for room in layout.values() if type(room).__name__=="BedAllocationRoom"]
        for baRoom in bedAllocationRooms:
            for room in baRoom.nextRooms:
                if type(room).__name__ != "BedRoom":
                    raise InvalidLayoutError

        return startRoomID, endRoomID


    #writes events to trace file
    def write_trace(self, event, level = 0):
        if self.settings.currInputs.saveTraceFile:
            with open(self.settings.currInputs.traceFilePath, 'a') as fTrace:
                fTrace.write("\n".join(["\t"*level+line for line in event.split("\n")])+"\n")
        
    #prints lists of rooms with patients in each room
    def print_hospital(self):
        for roomID, room in self.layout.iteritems():
            print "="*15
            print room
            for pid, patient in room.patients.iteritems():
                print patient
            print "="*15+"\n"

if __name__=="__main__":
    hospital=Hospital()
    sRoom = StartRoom()
    eRoom = EndRoom()
    sRoom.add_connections(eRoom.id)

    hospital.import_layout((sRoom, eRoom))
