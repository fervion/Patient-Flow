
# -*- coding: utf-8 -*-
"""
Simulation of patients in hospital
Name: Room.py
Author: Taige Hou (thou1@partners.org)
"""

"""
Contains classes for the room types that patients transition between
"""
from patient import *
from pfUtil import *
import numpy as np
from itertools import chain, product
import random
from random import sample
from time import time


"""
Base class for any type of room
"""
class Room(object):
    #unique id for each Room object
    idCounter=0
    @staticmethod
    def reset_counter():
        Room.idCounter = 0
    def __init__(self, capacity,name="",guiCoords=(0,0), writeTrace = None, pfSim = None):
        self.id = Room.idCounter
        Room.idCounter += 1

        #A Dictionary matching patient ID's to patients in the room
        self.patients = {}

        #maxmimum number of people allowed in room (-1 means unlimited)
        self.capacity = capacity

        #only used if gui is enabled
        self.guiCoords = guiCoords[:]
        self.nextRooms = []

        #displayed name of the room
        self.name=name
        
        #number going to be in this room after moves occur
        self.occupancy = 0
        
        #function to write trace file
        self.writeTrace = writeTrace

        #outputs
        self.pfSim = pfSim
    def __str__(self):
        return "{roomtype} {roomid}:{roomname}\n\tConnections:{conn}".format(roomtype=type(self).__name__,
                                                                    roomname = self.name,
                                                                  roomid=self.id,
                                                                  conn=[nr.id for nr in self.nextRooms])
    def add_connections(self, nextRooms):
        #A list containing the next rooms this node is connected to
        try:
            self.nextRooms.extend(nextRooms)
        except TypeError:
            self.nextRooms.append(nextRooms)
    def get_num_patients(self):
        return len([p for p in self.patients.values() if isinstance(p,RealPatient)])
    def add_patient(self, p):
        self.patients[p.id] = p
        p.currRoom = self
        p.timeStepLastRoomChange = self.pfSim.currTimeStep
    def add_patients(self, pList):
        [self.add_patient(p) for p in pList]
    #removes patient from room and returns patient or None if nonexistant p can be patient ID or patient
    def remove_patient(self, p):
        patientToRemove = None
        if isinstance(p,BasePerson):
            patientToRemove = p
        else:
            try:
                patientToRemove = self.patients[p]
            except KeyError:
                patientToRemove = None
            
        try:
            patientToRemove.currRoom = None
            pID = patientToRemove.id
        except AttributeError:
            return None
        else:
            return self.patients.pop(pID, None)
    #moves patient from this room to given room
    def move_patient_to(self, pID, room):
        patientToAdd=self.remove_patient(pID)

        if patientToAdd:
            room.add_patient(patientToAdd)
    #This is called after rooms are created and connections made
    def initialize_room(self):
        pass
    #this is called to calculate transmission (either background or cohorting)
    def do_transmission(self, simInputs, writeTrace):
        transmProbBackground = simInputs.backgroundTransmissionProb
        transmProbCohorting = simInputs.cohortingTransmissionProb
        numNewCasesBackground = [0,0]
        numNewCasesCohorting = [0,0]
        for patient in self.patients.values():
            if not isinstance(patient, RealPatient):
                continue
            for disease in (MRSA, VRE):
                oldStatus = patient.trueCP[disease]
                if random.random() < transmProbBackground[disease]:
                    patient.become_infected(disease)
                if oldStatus != patient.trueCP[disease]:
                    writeTrace("Patient {} became infected with {}".format(patient.id, DISEASE_STRINGS[disease]), level = 1)
                    numNewCasesBackground[disease]+=1
        return numNewCasesBackground, numNewCasesCohorting

    #This should be only once called after each timestep moves have occured and is the cleanup for each room
    def finalize_room(self):
        self.reset_room()
    #resets room occupancy and status and can be called multiple times in the time step
    def reset_room(self):
        self.occupancy = len(self.patients.keys())
        
"""
All patients start in the model from the start room
"""
class StartRoom(Room):
    def __init__(self, name="", guiCoords=(0,0), writeTrace = None, pfSim = None):
        super(StartRoom, self).__init__(-1, name=name, guiCoords=guiCoords, writeTrace = writeTrace, pfSim = pfSim)
    #time step steps one time interval and returns the moves that need to happen if any
    def time_step(self, hospitalSettings, currTimeStep):
        moves=[]
        #each move is represented by a tuple of (patientID, currRoomID, nextRoomID)
        for pID,patient in self.patients.iteritems():
            if patient.status == DEAD or patient.status == DISCHARGE_HOME or patient.status == DISCHARGE_FACILITY:
                nextRoomID = None
            else:
                nextRoomID = self.nextRooms[0].id
            moves.append((patient, self.id, nextRoomID))
            if nextRoomID:
                #give patients a test if they were to recieve test on admission
                testMRSA = hospitalSettings.currInputs.testMRSA
                testVRE = hospitalSettings.currInputs.testVRE
                #make sure we are testing real patient
                patientToTest = patient
                if isinstance(patient, PhantomPatient):
                    patientToTest = patient.realPatient
		    assert isinstance(patientToTest, RealPatient)
				
                testPatient = False
                if testMRSA.testOnAdmit:
                    testPatient = True
                if testMRSA.testOnHist and patientToTest.has_observed_MRSA():
                    testPatient = True
                if testPatient:
                    patientToTest.do_test(MRSA, testMRSA, currTimeStep, hospitalSettings.pfSim.outputs)
                    self.writeTrace("Patient {} Tested for MRSA at Entry".format(patientToTest.id, level = 2))
                testPatient = False
                if testVRE.testOnAdmit:
                    testPatient = True
                if testVRE.testOnHist and patientToTest.has_observed_VRE():
                    testPatient = True
                if testPatient:
                    patientToTest.do_test(VRE, testVRE, currTimeStep, hospitalSettings.pfSim.outputs)
                    self.writeTrace("Patient {} Tested for VRE at Entry".format(patientToTest.id, level = 2))
                
        return moves

    #draws room on canvas
    def draw(self, canvas):
        canvas.AddRectangle( self.guiCoords, (5,10), LineColor="Black", FillColor=None)

"""
Exit point for patients in the model
"""
class EndRoom(Room):
    def __init__(self, name="", guiCoords=(0,0), writeTrace = None, pfSim = None):
        super(EndRoom, self).__init__(-1, name=name,guiCoords=guiCoords, writeTrace = writeTrace, pfSim = pfSim)
        self.lengthsOfStay = [] #array of patient los times
    def time_step(self, hospitalSettings, currTimeStep):
        for patient in self.patients.values():
            self.lengthsOfStay.append(currTimeStep - patient.timeEnteringHospital)
            self.pfSim.outputs.timeSpentExposed.append(patient.timeStepsExposedToInfection)
        self.patients.clear()
        return []

"""
Decides which BedAllocationRooms to put patient in based on acuity
WardAllocationRooms can only be connected to BedAllocationRooms
"""
class WardAllocationRoom(Room):
    def __init__(self, name="", guiCoords=(0,0), writeTrace = None, pfSim = None):
        super(WardAllocationRoom,self).__init__(-1, name=name, guiCoords=guiCoords, writeTrace = writeTrace, pfSim = pfSim)
    def time_step(self, hospitalSettings, currTimeStep):
        moves = []
        #looks for a matching BedAllocationRoom for the patient based on acuity
        for pID, patient in self.patients.iteritems():
            if patient.status == DEAD or patient.status == DISCHARGE_HOME  or patient.status == DISCHARGE_FACILITY:
                moves.append((patient, self.id, None))
                continue
            for baRoom in self.nextRooms:
                if patient.acuity in baRoom.allowedAcuity and patient.service in baRoom.allowedService:
                    moves.append((patient, self.id, baRoom.id))
                    break
        return moves
    def initialize_room(self):
        for baRoom in self.nextRooms:
            baRoom.wardAllocationRoom = self
            for bRoom in baRoom.nextRooms:
                bRoom.wardAllocationRoom = self
            
        

"""
Decides which rooms to put patients in, patients loop in this room until a suitable room is found
BedAllocationRooms can only be connected to BedRooms
BedRooms take the acuity of whichever allocation room it is connected to
BedAllocationRooms have a specific service and acuity
"""
class BedAllocationRoom(Room):
    def __init__(self, allowedAcuity, allowedService, name="", guiCoords=(0,0), writeTrace = None, pfSim = None):
        super(BedAllocationRoom,self).__init__(-1, name=name, guiCoords=guiCoords, writeTrace = writeTrace, pfSim = pfSim)
        self.allowedAcuity = allowedAcuity[:]
        self.allowedService = allowedService[:]
        self.wardAllocationRoom = None #reference to the ward allocation room
        self.units = {} #dictionary of unit names to beds in that unit
    def time_step(self, hospitalSettings, currTimeStep):
        moves = []

        #sort patients by real/phantom status and then time of arrival
        if self.pfSim.inputs.prioritizeRealPatientQueueing:
            patients = sorted(self.patients.values(), key = lambda p:(0 if isinstance(p,RealPatient) else 1,p.timeStepLastRoomChange))
        else:
            patients = sorted(self.patients.values(), key = lambda p:p.timeStepLastRoomChange)
        #looks in each room connected to it and looks for match
        for patient in patients:
            if patient.status == DEAD or patient.status == DISCHARGE_HOME  or patient.status == DISCHARGE_FACILITY:
                moves.append((patient, self.id, None))
                continue
            

                
            #if acuity or service mismatch move patient back to ward allocation room
            if patient.acuity not in self.allowedAcuity or patient.service not in self.allowedService:
                moves.append((patient, self.id, self.wardAllocationRoom.id))
                continue

            #don't move patients with a bed already assigned
            if patient.realPatient.delayPatient:
                continue

            #if we are restricting bed movement dont allow phantoms (transfer patients) to move
            if self.pfSim.hospital.blockBedMovement and isinstance(patient, PhantomPatient):
                continue

            #don't allow patients into rooms if they havent finished work up time
            if isinstance(patient, RealPatient) and patient.workUpTime > 0:
                continue                
            for room in sample(self.nextRooms, len(self.nextRooms)):
                #dont allow patients to be put into the same room as themselves
                if patient.realPatient in room.patients.values():
                    continue
                
                if room.patient_matches(patient, hospitalSettings.matchingProfile):
                    #found match
                    #roll for delay to bed
                    if patient.realPatient.hasBeenAllocated:
                        dist = hospitalSettings.currInputs.transferBedDelayDist
                    else:
                        dist = hospitalSettings.currInputs.initialBedDelayDist

                    timeStepToEndDelay = currTimeStep + draw_from_distribution(dist)
                    
                    #Create Delay Patient
                    delayPatient = DelayPatient(patient.realPatient, timeStepToEndDelay)
                    #Put Delay patient in Room
                    room.add_patient(delayPatient)
                    
                    if isinstance(patient, QueueingPatient):
                        #remove queueing patient if room found
                        if patient.reason == REASON_ACUITY: #patient moved due to acuity change
                            self.pfSim.outputs.numBedMovesAcuity[-1]+=1
                            day = self.pfSim.getCurrDateTime().isoweekday()
                            self.pfSim.outputs.numBedMovesAcuityByDay[day]+=1
                        elif patient.reason == REASON_CP:
                            self.pfSim.outputs.numBedMovesCP[-1]+=1
                        realPatient = patient.realPatient
                        moves.append((patient, self.id, None))
                        realPatient.queueingPatient = None
                    
                    """
                    moves.append((patient, self.id, room.id))
                    """
                    room.occupancy += 1
                    room.currRoomGender = patient.gender
                    room.currCPStatus[MRSA] = room.currCPStatus[MRSA] or patient.observedCP[MRSA]
		    room.currCPStatus[VRE] = room.currCPStatus[VRE] or patient.observedCP[VRE]

                    if patient.acuity == 5:
                        if isinstance(patient, PhantomPatient):
                            self.pfSim.outputs.numPeopleToAcuityFive[patient.realPatient.currRoom.allowedAcuity[0]][-1]+=1
                        else:
                            self.pfSim.outputs.numPeopleToAcuityFive[5][-1]+=1
                    break
        return moves
    def __str__(self):
        return super(BedAllocationRoom, self).__str__()+"\n\tAllowed Acuity:{acuity}\n\tAllowed Service:{service}\n\tTime to Initial Bed Assignment:{tiba}"\
               .format(acuity=self.allowedAcuity, service = self.allowedService, tiba = self.get_time_to_initial_bed_assignment())
    #Returns mean and std dev of time to initial bed assignment for this allocation room
    def get_time_to_initial_bed_assignment(self):
 
        times = []
        for bedRoom in self.get_bed_rooms():
            times.extend(bedRoom.timesToInitialBedAssignment)
        times = np.array(times)
        
        if len(times):
            return times.mean(), times.std()
        else:
            return None, None

    #returns num closed beds controlled by this allocation as a tuple of (gender, mrsa, vre, both)  each bed can only be in one category of tuples
    def get_num_closed_beds_by_category(self):
        bedRooms = list(chain(*self.units.values()))
        unfilledBedRooms = [room for room in bedRooms if room.occupancy < room.capacity and room.occupancy > 0]
        waitingPatients = self.patients.values()

        bins = product((1, 0), (1, 0), (1,0)) #tuples of (gender?, mrsa?, vre?) (1,1,0) would be closed bed due to gender and mrsa
        closedBedDict = dict([(combo, 0) for combo in bins])
        numClosed = 0
        numClosedStaffing = 0
        for bedRoom in bedRooms:
            currRoomBedsClosedStaffing = len([bed for bed in bedRoom.bedsBlockedStaffing if bed])
            numClosedStaffing += currRoomBedsClosedStaffing

            currRoomNumClosed = currRoomBedsClosedStaffing
            #unfilled double room
            if bedRoom.occupancy < bedRoom.capacity and bedRoom.occupancy > 0:
                genderMismatch = False
                MRSAMismatch = False
                VREMismatch = False
                bedClosed = True

                for patient in waitingPatients:
                    patientGenderMismatch = False
                    patientMRSAMismatch = False
                    patientVREMismatch = False
                    if patient.gender != bedRoom.currRoomGender:
                        genderMismatch = True
                        patientGenderMismatch = True
                    if patient.observedCP[MRSA] != bedRoom.currCPStatus[MRSA]:
                        MRSAMismatch = True
                        patientMRSAMismatch = True
                    if patient.observedCP[VRE] != bedRoom.currCPStatus[VRE]:
                        VREMismatch = True
                        patientVREMismatch = True
                    if not any((patientGenderMismatch, patientMRSAMismatch, patientVREMismatch)):
                        bedClosed = False
                #if there are no patients waiting bed is not closed
                if not waitingPatients:
                    bedClosed = False
                if bedClosed:
                    currRoomNumClosed = max(1, currRoomBedsClosedStaffing)
                    closedBedDict[(int(genderMismatch), int(MRSAMismatch), int(VREMismatch))]+=1
            numClosed += currRoomNumClosed

        return closedBedDict, numClosed, numClosedStaffing
                
        
    def get_bed_rooms(self):
        bedRooms = []
        for unitName, unit in self.units.iteritems():
            bedRooms.extend(unit)
        return bedRooms
    def initialize_room(self):
        for bRoom in self.nextRooms:
            bRoom.allowedAcuity = self.allowedAcuity
            bRoom.allowedService = self.allowedService
            bRoom.bedAllocationRoom = self
            try:
                self.units[bRoom.unitName].append(bRoom)
            except KeyError:
                self.units[bRoom.unitName] = [bRoom]
    def finalize_room(self):
        super(BedAllocationRoom, self).finalize_room()
        realPatients = [p for p in self.patients.values() if isinstance(p, RealPatient)]

        for p in realPatients:
            if not p.hasBeenAllocated:
                p.timeToInitialBedAssignment+=1
                #subtract one hour from work up time
                if p.workUpTime > 0:
                    p.workUpTime -= 1
"""
Bed Room for patients to stay in
"""
class BedRoom(Room):
    def __init__(self, capacity, name="", unitName="", guiCoords=(0,0), writeTrace = None, pfSim = None):
        super(BedRoom,self).__init__(capacity, name=name, guiCoords=guiCoords, writeTrace = writeTrace, pfSim = pfSim)
        self.currRoomGender = None #room gender shifts based on whos occupying room
        self.currCPStatus = [None, None] #currCPStatus is a list of booleans (hasMRSA, hasVRE) of the CP status of the occupants of the room
        self.allowedAcuity = None #This is set by the allocation room that owns it
        self.allowedService = None #This is set by the allocation room that owns it
        self.wardAllocationRoom = None #reference to the ward allocation room
        self.bedAllocationRoom = None #reference to the bed allocation room
        self.beds = [None for i in range(capacity)] #list of beds with patients in them
        self.bedsBlockedStaffing = [False for i in range(capacity)] #Blocking status of beds due to staffing
        self.timeToStopBlockingBedsStaffing = [None for i in range(capacity)] #Time to stop blocking each bed due to staffing
        self.unitName = unitName #name of the unit this bed is in
        self.timesToInitialBedAssignment = []
    def time_step(self, hospitalSettings, currTimeStep):
        moves = []

        #create Queueing patient if acuity or cp mismatch
        for pID, patient in self.patients.iteritems():
            if isinstance(patient, CleaningPatient):
                if currTimeStep == patient.timeStepToEndCleaning:
                    moves.append((patient, self.id, None))
                continue
            if patient.status == DEAD or patient.status == DISCHARGE_HOME  or patient.status == DISCHARGE_FACILITY:
                moves.append((patient, self.id, None))
                continue
            #check if is time to end delay
            if isinstance(patient, DelayPatient):
                if currTimeStep >= patient.timeStepToEndDelay:
                    realPatient = patient.realPatient
                    
                    moves.append((realPatient, realPatient.currRoom.id, self.id))

		    #do mrsa/vre test for entry into icu if required
		    if 5 in self.allowedAcuity:
                        testMRSA = hospitalSettings.currInputs.testMRSA
                        testVRE = hospitalSettings.currInputs.testVRE
                        if testMRSA.testOnICU:
                            realPatient.do_test(MRSA, testMRSA, currTimeStep, hospitalSettings.pfSim.outputs)
                            self.writeTrace("Patient {} Tested for MRSA at ICU".format(realPatient.id, level = 2))
                        if testVRE.testOnICU:
                            realPatient.do_test(VRE, testVRE, currTimeStep, hospitalSettings.pfSim.outputs)
                            self.writeTrace("Patient {} Tested for VRE at ICU".format(realPatient.id, level = 2))
                    
            if isinstance(patient, RealPatient):
                #if patient has delayPatient do not create phantom
                if not patient.delayPatient:
                    if not patient.queueingPatient:
                        if patient.acuity not in self.allowedAcuity:
                            queueingPatient = QueueingPatient(patient, reason = REASON_ACUITY)

                            self.wardAllocationRoom.add_patient(queueingPatient)

                            if 3 in self.allowedAcuity:
                                    self.pfSim.outputs.numPeopleFromAcuityThree[patient.acuity][-1]+=1
                            continue
                        if patient.observedCP != self.currCPStatus:
                            queueingPatient = QueueingPatient(patient, reason = REASON_CP)
                            self.bedAllocationRoom.add_patient(queueingPatient)
                            continue
                    #delete phantom if acuity and cp matches again
                    else:
                        if patient.acuity in self.allowedAcuity and patient.observedCP == self.currCPStatus:
                            patient.queueingPatient.currRoom.remove_patient(patient.queueingPatient)
                            patient.queueingPatient = None
                
        return moves
    def add_patient(self,p, bedIndex = None):

        #move patient to bed and destroy phantom if bed found
        if (isinstance(p, QueueingPatient)):
            if p.reason == REASON_ACUITY: #patient moved due to acuity change
                self.pfSim.outputs.numBedMovesAcuity[-1]+=1
                day = self.pfSim.getCurrDateTime().isoweekday()
                self.pfSim.outputs.numBedMovesAcuityByDay[day]+=1
            elif p.reason == REASON_CP:
                self.pfSim.outputs.numBedMovesCP[-1]+=1
            realPatient = p.realPatient
            self.remove_patient(p)
            realPatient.currRoom.move_patient_to(realPatient.id, self)
            realPatient.queueingPatient = None
            return
        
        if isinstance(p, RealPatient):
            #first time patient has ever been in a bed
            if not p.hasBeenAllocated:
                p.hasBeenAllocated = True
                self.timesToInitialBedAssignment.append(p.timeToInitialBedAssignment)
                currTime = self.pfSim.getCurrDateTime()
                self.pfSim.outputs.numPatientsBedAdmissionHour[currTime.isoweekday()][currTime.hour]+=1
            #delete delay patient
            if p.delayPatient:
                p.delayPatient.currRoom.remove_patient(p.delayPatient)
                p.delayPatient = None


            
        super(BedRoom,self).add_patient(p)

        
        if bedIndex is not None:
            self.beds[bedIndex] = p
        else:
            for i in range(len(self.beds)):
                if not self.beds[i] and not self.bedsBlockedStaffing[i]:
                    self.beds[i] = p
                    break
                


        self.reset_room()
    def remove_patient(self,p, removeCleaningPatient = False):
        #Cannot remove cleaning Patient unless specified
        if isinstance(p, CleaningPatient) and not removeCleaningPatient:
            return None

        removedPatient = super(BedRoom,self).remove_patient(p)

        for i in range(len(self.beds)):
            if self.beds[i] is removedPatient:
                self.beds[i] = None
                break
        self.reset_room()
        return removedPatient
    #gets the matching of pstatus for patients.  Can be either 0,1,2,3 (0 means there is no mismatch)
    def get_pstatus_mismatch(self):
        realPatients = [p for p in self.patients.values() if isinstance(p, RealPatient)]
        if len(realPatients) <= 1:
            return PSTATUS_NONE
        trueMRSA = [p.trueCP[0] for p in realPatients]
        trueVRE = [p.trueCP[1] for p in realPatients]
        matchMRSA = all([status == trueMRSA[0] for status in trueMRSA])
        matchVRE = all([status == trueVRE[0] for status in trueVRE])
        if matchMRSA and matchVRE:
            return PSTATUS_NONE
        elif matchVRE and not matchMRSA:
            return PSTATUS_MRSA
        elif matchMRSA and not matchVRE:
            return PSTATUS_VRE
        else:
            return PSTATUS_BOTH
    #returns number of acuity mismatches by high/low
    def num_acuity_mismatch(self):
        realPatients = [p for p in self.patients.values() if isinstance(p, RealPatient)]
        lowMismatches = [p for p in realPatients if p.acuity < max(self.allowedAcuity)]
        highMismatches = [p for p in realPatients if p.acuity > min(self.allowedAcuity)]
        return len(lowMismatches), len(highMismatches)
    #returns true if the patient can be put in the room
    def patient_matches(self, patient, matchingProfile):
        numBlockedBeds = len([bed for bed in self.bedsBlockedStaffing if bed])
        if self.occupancy >= self.capacity - numBlockedBeds:
            return False
        if self.currRoomGender is not None and self.currRoomGender != patient.gender:
            return False
        if matchingProfile[MRSA] and self.currCPStatus[MRSA] is not None and self.currCPStatus[MRSA] != patient.observedCP[MRSA]:
            return False
        if matchingProfile[VRE] and self.currCPStatus[VRE] is not None and self.currCPStatus[VRE] != patient.observedCP[VRE]:
            return False
        return True
    #this is called to calculate transmission (either background or cohorting)
    def do_transmission(self, simInputs, writeTrace):
        transmProbBackground = simInputs.backgroundTransmissionProb
        transmProbCohorting = simInputs.cohortingTransmissionProb
        numNewCasesBackground = [0,0]
        numNewCasesCohorting = [0,0]  
        roomTrueCPStatus = [False, False]
        roomUninfectedPeople = [[],[]]
        realPatients = [p for p in self.patients.values() if isinstance(p, RealPatient)]
        for patient in realPatients:
            for disease in (MRSA, VRE):
                if patient.trueCP[disease]:
                    roomTrueCPStatus[disease] = True
                else:
                    roomUninfectedPeople[disease].append(patient)
        for disease in (MRSA, VRE):          
            for uninfected in roomUninfectedPeople[disease]:
                becomeInfected = False
                becomeInfectedCause = [False, False] #background, cohorting
                oldStatus = uninfected.trueCP[disease]
                if random.random() < transmProbBackground[disease]:
                    becomeInfected = True
                    becomeInfectedCause[0] = True
                if roomTrueCPStatus[disease]:
                    if random.random() < transmProbCohorting[disease]:
                        becomeInfected = True
                        becomeInfectedCause[1] = True
                if becomeInfected:
                    uninfected.become_infected(disease)
                if oldStatus != uninfected.trueCP[disease]:                        
                    writeTrace("Patient {} became infected with {}".format(uninfected.id, DISEASE_STRINGS[disease]), level = 1)
                    if random.choice([index for index, result in enumerate(becomeInfectedCause) if result]) == 0: 
                        numNewCasesBackground[disease]+=1
                    else:
                        numNewCasesCohorting[disease]+=1
        return numNewCasesBackground, numNewCasesCohorting
    def reset_room(self):
        super(BedRoom,self).reset_room()
        allPatients = [p for p in self.patients.values() if isinstance(p, BasePatient)]
        if not allPatients:
            self.currRoomGender = None
            self.currCPStatus = [None, None]
	else:
            self.currRoomGender = allPatients[0].gender
	    self.currCPStatus[MRSA] = any([patient.observedCP[MRSA] for patient in allPatients])
	    self.currCPStatus[VRE] = any([patient.observedCP[VRE] for patient in allPatients])
	    
    def finalize_room(self):
        super(BedRoom, self).finalize_room()
        allPatients = [p for p in self.patients.values() if isinstance(p, RealPatient)]
        MRSAStatus = any([p.has_MRSA() for p in allPatients])
        VREStatus = any([p.has_VRE() for p in allPatients])

        for p in allPatients:
            MRSAExposure = False
            VREExposure = False
            if MRSAStatus and not p.has_MRSA():
                MRSAExposure = True
            if VREStatus and not p.has_VRE():
                VREExposure = True
            if MRSAExposure and VREExposure:
                p.timeStepsExposedToInfection[EXPOSURE_BOTH]+=1
            if MRSAExposure:
                p.timeStepsExposedToInfection[EXPOSURE_MRSA]+=1
            if VREExposure:
                p.timeStepsExposedToInfection[EXPOSURE_VRE]+=1
 

    def __str__(self):
        return super(BedRoom, self).__str__()+"\n\tCurr Room Gender:{roomgender}\n\tMRSA, VRE:{cp}\n\tAllowed Acuity:{acuity}"\
               .format(roomgender=self.currRoomGender,
                       cp=self.currCPStatus,
                       acuity=self.allowedAcuity)
    
if __name__=='__main__':
    for i in range(10):
        r=StartRoom()
    print r

