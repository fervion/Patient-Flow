# -*- coding: utf-8 -*-
"""
Graphical Interface for CEPAC Cluster tool Library
Name: CEPACClusterLibGui.py
Author: Taige Hou (thou1@partners.org)
"""


"""
Contains classes for running the patient flow simulation
"""
from hospital import *
from patient import *
from room import *
from itertools import product
from xlrd import open_workbook, xldate_as_tuple
import numpy as np
import datetime as dt
from datetime import timedelta, datetime
from os import path
from math import log, exp

class InvalidHospitalError(Exception):
    pass

class InputSheetError(Exception):
    #define error codes and messages
    def __init__(self, code):
        self.code = code
    def error_message(self):
        codes = {1:"Error while reading acuity table"}
        return codes.get(self.code,"Input Error")

#Function that calculates geometric mean and geometric std deviation
def gmean(s):
    mean = exp(np.mean([log(a) for a in s]))
    std = exp(np.mean([log(a/mean)**2 for a in s])**(0.5))
    return mean, std


#stores all the inputs for the simulation
class Inputs(object):
    PATIENT_TYPES = product(GENDERS,SERVICE_VALUES, PSTATUS_VALUES,ACUITY_VALUES)
    #class for test characteristics and strategy
    class Test(object):
        def __init__(self, sensitivity, specificity, tat, testOnICU, testOnHist,  testOnAdmit, numTests):
            self.sensitivity = sensitivity
            self.specificity = specificity
            self.tat = tat
            self.testOnICU = testOnICU
            self.testOnHist = testOnHist
            self.numTestsForDeflag = numTests
            self.testOnAdmit = testOnAdmit
        def __str__(self):
            return str(vars(self))
    def __init__(self):
        self.inputFilePath = None

        #run parameters
        #datetime object for when model starts
        self.startTime = None
        self.saveEveryN = None
        self.shufflePeriod = None
        self.saveTraceFile = None
        self.traceFilePath = None
        self.useFixedSeed = None
        self.runInPeriod = None #number of steps in run-in period.  The outputs are reset at the end of the run in period.
        self.prioritizeRealPatientQueueing = None #prioritize real patient before phantom patient in queue
        self.blockBedMovements = []
        self.blockBedMovement = None #enable block bed movement period
        self.blockBedMovementStartTime = None
        self.blockBedMovementEndTime = None
        self.blockDischarge = None #enable blocking of dishcarges for first few time steps since arrival
        self.blockDischargeTimeSteps = None #number of time steps since arrival to block for
        
        #new patient inputs
        self.numArrivingPatientsMeanWeekday = None #lambda of poisson distribution for number of patients arriving per time step stratified by time of day
        self.numArrivingPatientsMeanWeekend = None
        self.arrivingPatientDistribution = dict([(pType, 0) for pType in Inputs.PATIENT_TYPES]) #dictionary which matches (gender, service, pstatus, acuity) to frequency (note this is not normalized)
        self.trueCPDistributions = {MRSA:{True:{}, False:{}}, VRE:{True:{}, False:{}}} #distribution of True CP based on observed for incoming patients
        self.probLeaveWeekday = dict([(fromAcuity, {}) for fromAcuity in ACUITY_VALUES])
        self.probLeaveWeekend = dict([(fromAcuity, {}) for fromAcuity in ACUITY_VALUES])
        
        #movement inputs
        self.probClosedBed = None #Prob of closing bed for staffing reasons
        self.distribClosedBed = None #lambda of exponential distribution for duration of bed closure
        self.initialBedDelayDist = {} #Dict of distribution to draw delay to arrival in bed
        self.transferBedDelayDist = {}
        self.workUpTime = dict([(acuity, {}) for acuity in ACUITY_VALUES]) #Dict of diststribution of work up times based on acuity
        
        #acuity change table
        #times to change acuity.  Patient will draw acuity if the time they have spent in the hospital is a multiple of this value
        self.lengthOfStayToChangeAcuity = None
        #dict of threshold times (times after which to use threshold distribution)
        self.acuityChangeThresholdTimes = {}
        self.acuityChangeThresholdDistribution = dict([(fromAcuity, {}) for fromAcuity in ACUITY_VALUES])
        #this is a dict[fromAcuity][hoursAtAcuity][toAcuity] to count
        self.acuityChangeDistribution = dict([(fromAcuity, {}) for fromAcuity in ACUITY_VALUES])
        #distribution of people leaving the hospital strat by acuity, time, and leave type a[fromAcuity][hour][endAcuity]
        self.leaveDistributionsWeekday = dict([(fromAcuity, {}) for fromAcuity in ACUITY_VALUES])
        self.leaveDistributionsWeekend = dict([(fromAcuity, {}) for fromAcuity in ACUITY_VALUES])
        
        #transmission characteristics
        #probabilities of acquiring MRSA or vre
        self.backgroundTransmissionProb = {}
        self.cohortingTransmissionProb = {}

        #testing inputs
        self.testMRSA = None
        self.testVRE = None
        
class OutputStatistics(object):
    def __init__(self):
        self.totalNumBeds = None #number of beds in hospital
        self.totalNumSemiPrivateBeds = None #number of semiprivate beds in hospital
        self.numPatientsInHospital = [] # number of people in the hospital per time step
        self.numPatientsInHospitalAcuity =  dict([(acuity, []) for acuity in ACUITY_VALUES])# number of people in hospital per time step by acuity
        self.numNewPatients = [] #number of new patients per time step
        self.numPatientsLeaving = {DEAD:[], DISCHARGE_HOME:[], DISCHARGE_FACILITY:[]} #number of people leaving the hospital at each time step stratified by reason
        self.numPatientsLeavingHour = dict([(day,dict([(hour, 0) for hour in range(24)]))for day in range(1,8)]) #total number of people leaving the hospital over the period stratified by day of week and hour
        self.numPatientsBedAdmissionHour = dict([(day,dict([(hour, 0) for hour in range(24)]))for day in range(1,8)]) #number of patients physically admitted into a bed that hour(does not count transfers)

        self.numClosedBeds = [] #number of closed beds per time step
        self.numClosedBedsByCategory = [] #num closed beds by category per time step
        self.numClosedBedsCP = [] #num closed beds due to CP per time step
        self.numClosedBedsStaffing = [] #num closed beds for staffing reasons per time step
        self.numBedsCP = dict([(pstatus, [])for pstatus in PSTATUS_VALUES]) #number of beds by cp per time step
        self.operationalCapacity = [] #operational capacity by time step
        self.timesToIBA = [] #Array of people (not stratified by time) Times to initial Bed assignment
        self.lengthsOfStay = [] #Array of people (not stratified by time)
        self.numMismatchCPRooms = [] #time stratifed array of number of mismatched rooms (true cp wise) array is stratified by pstatus
        self.numMismatchCPBeds = [] #time stratifed array of number of mismatched beds (true cp wise) array is stratified by pstatus
        self.numMismatchAcuityRooms = []
        self.numMismatchAcuityPatients = []
        self.numMismatchAcuityBedsHigh = [] #patient is higher than room acuity
        self.numMismatchAcuityBedsLow = [] #patient is lower than room acuity
        self.numBedMovesShuffle = [] #number of bed moves by shuffling
        self.numBedMovesAcuity = [] #number of bed moves per time step from acuity mismatch
        self.numBedMovesAcuityByDay = dict([(day, 0) for day in range(1,8)]) #number of bed moves due to acuity mismatch by day of week
        self.numBedMovesCP = [] #number of bed moves per time step from cp mismatch
        self.numClosedBedsGender = None
        self.numClosedBedsMRSA = None
        self.numClosedBedsVRE = None
        self.numClosedBedsBoth = None

        self.numQueuing = [] #number of patients that are queueing for a room (either in ward or bed allocation room) can be real or phantom
        self.numQueuingReal = []
        self.numQueuingPhantom = []

        self.numBeddedPatients = []#number of patients in a bed per time step
        self.bedCapacity = []# number of bedded patients/num avail beds per time step 

        self.numUnflaggedCarrierAtAdmit = [[],[]] #number of people with discordant between true and observed (unflagged carrier) at admit per time step stratified by mrsa and vre. 
        self.numFlaggedNonCarrierAtAdmit = [[],[]] #number of people with discordant between true and observed (flagged non carrier) at admit per time step stratified by mrsa and vre. 
        self.numFlagsAquired = {MRSA:[], VRE:[]}#number of cp flags aquired per time step
        self.numFlagsLost = {MRSA:[], VRE:[]}#number of cp flags deflagged per time step
        
        self.numDiffObservedCPDischarge = [[],[]] #number of people who have different observed cp status at discharge than entry (strat by disease) per time step
        self.numNewInfectedBackground = [[],[]] #number of new infections per time step stratified by disease (incidence)
        self.numNewInfectedCohorting = [[],[]] #number of new infections per time step stratified by disease (incidence)
        self.numInfected = [[],[]] #number of people infected per time step (prev)
        self.timeSpentExposed = []#time spent exposed to infection per person

        self.numTestsGiven = {MRSA:[], VRE:[]}#number of tests given per time step strat by disease
        self.numTestsGivenUnique = {MRSA:[], VRE:[]} #number of tests given to people who have never been tested previously per time step
        self.testResults = dict([(combo, []) for combo in product((MRSA, VRE),(True, False),(True, False))])#number of Test results per time step as dict of a[disease][trueStatus][testResult] 

        #temp outputs
        self.numPeopleToAcuityFive = dict([(acuity, []) for acuity in ACUITY_VALUES])
        self.numPeopleFromAcuityThree = dict([(acuity, []) for acuity in ACUITY_VALUES])

class Sim(object):
    def __init__(self):
        #time delta object defaults to 1 hour
        self.timeStep = timedelta(hours = 1)
        
        self.hospital = None
        self.timeStepsRemaining = -1
        self.currTimeStep = 0
        self.inputs = Inputs()
        self.outputs = OutputStatistics()

        #stored layout stores the inputs for the last hospital layout imported
        self.storedLayout = None
        self.initialPatients = []
        
        #set flags
        self.PAUSE_SIM = False #used to interupt run of model
        self.SIM_INITIALIZED = False #set to true when the hospitial is initialized
        self.INITIAL_PATIENTS_LOADED = False
        self.INPUTS_LOADED = False
    #Creates hospital object
    def create_hospital(self):
        self.hospital = Hospital()
    #resets hospital to blank state and reloads map
    def reset_hospital(self):
        BasePerson.reset_counter()
        Room.reset_counter()
        self.hospital = Hospital()
        self.timeStepsRemaining = -1
        self.currTimeStep = 0
        self.inputs = Inputs()
        self.outputs = OutputStatistics()
        self.initialPatients = []
        self.PAUSE_SIM = False #used to interupt run of model
        self.SIM_INITIALIZED = False #set to true when the hospitial is initialized
        self.INITIAL_PATIENTS_LOADED = False
        self.INPUTS_LOADED = False

        self.import_layout(*self.storedLayout)
        

    """
    imports and create initial patients
    must have already created or imported hospital
    importType specifies what type of import we are doing and args specifies additional information
    import type can be
    'l' list of tuples of patient characteristics
    """
    def import_initial_patients(self, importType, arg):
        BasePerson.reset_counter()
        if not self.hospital:
            raise InvalidHospitalError

        #list import expect a list of tuple of patient characteristics (gender,)
        if importType == 'l':
           self.initialPatients = [RealPatient(*p) for p in arg]
        else:
            raise

    """
    imports and create hospital layout
    must have already created or imported hospital
    importType specifies what type of import we are doing and args specifies additional information
    import type can be
    'l' list of tuples of rooms with (unique identifier, room type, roomArgs, connections)
        unique id is any identifier to specify the room
        room type is the type of room e.g. Start, End
        room Args is a tuple of values used to initialize the room
        connections is a tuple of unique ID's
    """
    def import_layout(self, importType, arg):
        Room.reset_counter()
        if not self.hospital:
            raise InvalidHospitalError

        self.storedLayout = (importType, arg)
        
        self.outputs.totalNumBeds = 0
        self.outputs.totalNumSemiPrivateBeds = 0
        
        #list import expect a list of tuple of patient characteristics (gender,)
        if importType == 'l':
            roomDict = {}
            for rName, rType, rArgs, rConnections in arg:
                if rType+"Room" in [cls.__name__ for cls in Room.__subclasses__()]:
                    
                    rkwArgs = {'writeTrace': self.hospital.write_trace, 'pfSim': self}
                    roomDict[rName] = globals()[rType+"Room"](*rArgs, **rkwArgs)
                    if isinstance(roomDict[rName], BedRoom):
                        if roomDict[rName].capacity >1:
                            self.outputs.totalNumSemiPrivateBeds+=roomDict[rName].capacity
                        self.outputs.totalNumBeds+=roomDict[rName].capacity

                else:
                    raise NameError

            for rName, rType, rArgs, rConnections in arg:
                for connection in rConnections:
                    roomDict[rName].add_connections(roomDict[connection])
            
            self.hospital.import_layout(roomDict.values())
            
        else:
            raise

        self.initialize_run()
    #reads in the input file
    def load_input_file(self, inputFilePath):
        self.inputs.inputFilePath = inputFilePath
        wb = open_workbook(inputFilePath)
        
        #get simulation parameters
        simSheet = wb.sheet_by_name('Simulation')
        startDate = xldate_as_tuple(simSheet.cell(3,2).value, wb.datemode)
        startTime = xldate_as_tuple(simSheet.cell(3,3).value, wb.datemode)
        
        self.inputs.startTime = datetime(startDate[0],startDate[1], startDate[2],startTime[3], startTime[4])
        self.inputs.saveEveryN = int(simSheet.cell(3,7).value)
        self.inputs.shufflePeriod = int(simSheet.cell(5,2).value)
        self.inputs.saveTraceFile = simSheet.cell(5,7).value=="yes"
        if self.inputs.saveTraceFile:
            self.inputs.traceFilePath = path.splitext(path.basename(inputFilePath))[0]+".pfTrace"
            #create trace file
            with open(self.inputs.traceFilePath, 'w') as fTrace:
                pass
        self.inputs.useFixedSeed = simSheet.cell(7,2).value=="yes"
        if self.inputs.useFixedSeed:#set the seed
            random.seed(10)
            np.random.seed(10)
        else:
            random.seed(None)
            np.random.seed(None)
        self.inputs.numStepsToRun = int(simSheet.cell(7,7).value)
        self.inputs.runInPeriod = int(simSheet.cell(9,7).value)
        self.inputs.prioritizeRealPatientQueueing = simSheet.cell(12,2).value=="yes"

        for i in range(5):
            doBlock = simSheet.cell(16+i,2).value=="yes"
            startTime = dt.time(xldate_as_tuple(simSheet.cell(16+i,3).value, wb.datemode)[3])
            endTime = dt.time(xldate_as_tuple(simSheet.cell(16+i,4).value, wb.datemode)[3])
            self.inputs.blockBedMovements.append((doBlock, startTime, endTime))

        self.inputs.blockDischarge = simSheet.cell(23, 2).value=="yes"
        self.inputs.blockDischargeTimeSteps = int(simSheet.cell(24,2).value)
        
        #get new patient distributions
        patientSheet = wb.sheet_by_name('Patients')
        self.inputs.numArrivingPatientsMeanWeekday = patientSheet.col_values(14,7,31)
        self.inputs.numArrivingPatientsMeanWeekend = patientSheet.col_values(17,7,31)

        colStartNum = 14
        rowStartNum = 37
        for i, fromAcuity in enumerate(ACUITY_VALUES):
            for j in range(24):
                self.inputs.probLeaveWeekday[fromAcuity][j] = patientSheet.cell(rowStartNum+j, colStartNum+i*13).value
                self.inputs.probLeaveWeekend[fromAcuity][j] = patientSheet.cell(rowStartNum+j, colStartNum+i*13+2).value
                self.inputs.leaveDistributionsWeekday[fromAcuity][j] = {}
                self.inputs.leaveDistributionsWeekend[fromAcuity][j] = {}
                for k,toAcuity in enumerate(ACUITY_VALUES_LEAVE):
                    self.inputs.leaveDistributionsWeekday[fromAcuity][j][toAcuity] = patientSheet.cell(rowStartNum+j, colStartNum+i*13+k+4).value
                    self.inputs.leaveDistributionsWeekend[fromAcuity][j][toAcuity] = patientSheet.cell(rowStartNum+j, colStartNum+i*13+k+8).value

        distribRows = range(7,23)+range(25,41)
        distribCols = ((0,5),(6,11))
        for row in distribRows:
            for colSet in distribCols:
                gender, service, pstatus, acuity, count = patientSheet.row_values(row,colSet[0],colSet[1])
                if gender == "M":
                    gender = MALE
                else:
                    gender = FEMALE
                if service == "Medical":
                    service = MEDICINE
                else:
                    service = SURGERY
                pstatus = int(pstatus)
                acuity = int(acuity)
                self.inputs.arrivingPatientDistribution[(gender, service, pstatus, acuity)] = count


        cpDistributionCells =((5,21,(True, True), MRSA), (6,21,(True, False),MRSA), (8,21,(False, True),MRSA),(9,21,(False, False), MRSA),(5,25,(True, True), VRE),(6,25,(True, False), VRE),(8,25,(False, True), VRE),(9,25,(False, False), VRE))
        for row, col, status, condition in cpDistributionCells:
            self.inputs.trueCPDistributions[condition][status[0]][status[1]] = float(patientSheet.cell(row,col).value)


        #get movement sheet
        movementSheet = wb.sheet_by_name('Movement')
                
        self.inputs.probClosedBed = float(movementSheet.cell(3,1).value)
        self.inputs.distribClosedBed = float(movementSheet.cell(5, 1).value)


        i = 0
        while True:
            try:
                self.inputs.initialBedDelayDist[int(movementSheet.cell(4+i,7).value)] = float(movementSheet.cell(4+i,8).value)
            except (ValueError, IndexError):
                break
            i+=1

        i = 0
        while True:
            try:
                self.inputs.transferBedDelayDist[int(movementSheet.cell(4+i,10).value)] = float(movementSheet.cell(4+i,11).value)
            except (ValueError, IndexError):
                break
            i+=1

        for j,acuity in enumerate(ACUITY_VALUES):
            i = 0
            while True:
                try:
                    self.inputs.workUpTime[acuity][int(movementSheet.cell(3+i,14).value)] = float(movementSheet.cell(3+i,15+j).value)
                except (ValueError, IndexError):
                    break
                i+=1
        
        #get acuity change sheet
        acuitySheet = wb.sheet_by_name('Acuity Change')

        #when to change acuity
        self.inputs.lengthOfStayToChangeAcuity = int(acuitySheet.cell(1,12).value)
        #acuity change threshold distribution
        self.inputs.acuityChangeThresholdTimes = dict([(ACUITY_VALUES[i],int(t)) for i,t in enumerate(acuitySheet.col_values(2,6,10))])
        colStartNum = 3
        for i,fromAcuity in enumerate(ACUITY_VALUES):
            rowNum = 6 + i
            frequencyTable = acuitySheet.row_values(rowNum,colStartNum, colStartNum + len(ACUITY_VALUES))
            for j,toAcuity in enumerate(ACUITY_VALUES):
                self.inputs.acuityChangeThresholdDistribution[fromAcuity][toAcuity] = frequencyTable[j]
        
        #acuity change table
        col = 0
        countStartRow = 14
        while True:
            try:
                fromAcuity, toAcuity = [int(n) for n in acuitySheet.col_values(col+1, 11,13)]
            except IndexError:
                break
            try:
                hours = [int(hour) for hour in acuitySheet.col_values(col, countStartRow, countStartRow+ self.inputs.acuityChangeThresholdTimes[fromAcuity])]
            except ValueError:
                raise InputSheetError(1)
                
            counts = acuitySheet.col_values(col+1, countStartRow, countStartRow+ self.inputs.acuityChangeThresholdTimes[fromAcuity])

            for i, hour in enumerate(hours):
                try:
                    self.inputs.acuityChangeDistribution[fromAcuity][hour][toAcuity] = counts[i]
                except KeyError:
                    self.inputs.acuityChangeDistribution[fromAcuity][hour] = {toAcuity:counts[i]}
    
            col += 3

        #transmission characteristics
        transmSheet = wb.sheet_by_name('Transmission')
        self.inputs.backgroundTransmissionProb[MRSA] = float(transmSheet.cell(3,1).value)
        self.inputs.backgroundTransmissionProb[VRE] = float(transmSheet.cell(4,1).value) 
        self.inputs.cohortingTransmissionProb[MRSA] = float(transmSheet.cell(7,1).value)
        self.inputs.cohortingTransmissionProb[VRE] = float(transmSheet.cell(8,1).value)

        #testing inputs
        testingSheet = wb.sheet_by_name('Testing')
        sen, spe, tat = testingSheet.row_values(4, 1, 4)
        testICU, testHist, testAdmit, numTest  = testingSheet.col_values(3,7,11)
        self.inputs.testMRSA = Inputs.Test(float(sen), float(spe), int(tat),testICU =="yes",testHist == "yes", testAdmit == "yes", int(numTest))   

        sen, spe, tat = testingSheet.row_values(4, 7, 10)
        testICU, testHist, testAdmit, numTest = testingSheet.col_values(9,7,11)
        self.inputs.testVRE = Inputs.Test(float(sen), float(spe), int(tat),testICU =="yes",testHist == "yes", testAdmit == "yes", int(numTest))   

        self.hospital.settings.currInputs = self.inputs
        self.hospital.settings.pfSim = self
        self.INPUTS_LOADED = True
    #Saves an output file
    def save_output_file(self, outputFilePath):
        #aggregate overall stats
        self.outputs.timesToIBA = []
        for room in self.hospital.layout.values():
            if isinstance(room, BedRoom):
                self.outputs.timesToIBA.extend(room.timesToInitialBedAssignment)
        def writeAvg(outFile, array, label):
            outFile.write("\n{}\t{}\t{}\t{}".format(label, np.mean(array),np.median(array), np.std(array)))
        with open(outputFilePath, 'w') as fOut:
            fOut.write("Output File for {}".format(self.inputs.inputFilePath))
            fOut.write("\nModel run for {} steps".format(self.currTimeStep))
            fOut.write("\nOutputs tracked for {} steps".format(self.currTimeStep-self.inputs.runInPeriod))
            fOut.write("\n")
            #overall statistics
            fOut.write("\n\tTotal\tSemi-Private")
            fOut.write("\nNumber of Beds\t{}\t{}".format(self.outputs.totalNumBeds, self.outputs.totalNumSemiPrivateBeds))
            fOut.write("\n\tNone\tMRSA Only\tVRE Only\tBoth")
            fOut.write("\nNumber of Filled Beds by CP per time step\t{}\t{}\t{}\t{}".format(np.mean(self.outputs.numBedsCP[PSTATUS_NONE]),
                                                                                     np.mean(self.outputs.numBedsCP[PSTATUS_MRSA]),
                                                                                     np.mean(self.outputs.numBedsCP[PSTATUS_VRE]),
                                                                                     np.mean(self.outputs.numBedsCP[PSTATUS_BOTH]),
                                                                                     ))
            fOut.write("\n\tMean\tMedian\tStd Dev")
            fOut.write("\nTotal Number of Patients (Incident)\t{}".format(sum(self.outputs.numNewPatients)))
            writeAvg(fOut, self.outputs.numNewPatients, "Number of Arriving Patients per time step")
            writeAvg(fOut,self.outputs.numPatientsInHospital, "Number of patients in hospital per time step")
            fOut.write("\n\t"+"\t".join([str(acuity) for acuity in ACUITY_VALUES]))
            fOut.write("\nNumber of patients in hospital by Acuity per time step\t"+"\t".join([str(np.mean((self.outputs.numPatientsInHospitalAcuity[acuity])))for acuity in ACUITY_VALUES]))
            fOut.write("\n\tHome\tFacility\tDeath\tTotal")
            fOut.write("\nNumber of People Leaving per time step\t{}\t{}\t{}\t{}".format(np.mean(self.outputs.numPatientsLeaving[DISCHARGE_HOME]),
                                                                                         np.mean(self.outputs.numPatientsLeaving[DISCHARGE_FACILITY]),
                                                                                         np.mean(self.outputs.numPatientsLeaving[DEAD]),
                                                                                         np.mean(sum([np.array(self.outputs.numPatientsLeaving[reason]) for reason in (DEAD, DISCHARGE_HOME, DISCHARGE_FACILITY)]))))
            fOut.write("\nTotal Number of Bed Admissions")
            fOut.write("\n\tMon\tTue\tWed\tThu\tFri\tSat\tSun")
            for hour in range(24):
                fOut.write("\n{}".format(hour))
                for day in range(1,8):
                    fOut.write("\t{}".format(self.outputs.numPatientsBedAdmissionHour[day][hour]))
                    
            fOut.write("\nTotal Number of People Leaving")
            fOut.write("\n\tMon\tTue\tWed\tThu\tFri\tSat\tSun")
            for hour in range(24):
                fOut.write("\n{}".format(hour))
                for day in range(1,8):
                    fOut.write("\t{}".format(self.outputs.numPatientsLeavingHour[day][hour]))

            writeAvg(fOut,self.outputs.lengthsOfStay, "Length of stay per patient per time step")

            fOut.write("\nLength of stay\tgeometric mean\tstd dev)")
            fOut.write("\n\t{}\t{}".format(*gmean(self.outputs.lengthsOfStay)))
            fOut.write("\nLength of stay distribution (hours)")
            bins = list(range(400))+[1000]
            hist, edges = np.histogram(self.outputs.lengthsOfStay, bins)
            fOut.write("\n"+"\t".join([str(val) for val in edges]))
            fOut.write("\n"+"\t".join([str(val) for val in hist]))
            fOut.write("\n\tTotal\tWaiting for first bed\tWaiting for transfer")
            fOut.write("\nNumber of People Queuing per time step\t{}\t{}\t{}".format(np.mean(self.outputs.numQueuing),
                                                                                 np.mean(self.outputs.numQueuingReal),
                                                                                 np.mean(self.outputs.numQueuingPhantom)))
            fOut.write("\n\tMean\tMedian\tStd Dev")
            writeAvg(fOut, self.outputs.bedCapacity, "Licensed Bed Capacity per time step")
            writeAvg(fOut, self.outputs.operationalCapacity, "Operational Capacity per time step")
            writeAvg(fOut, self.outputs.numClosedBeds, "Number of Closed Beds per time step")
            writeAvg(fOut, self.outputs.numClosedBedsCP, "Number of Closed Beds per time step due to CP")
            writeAvg(fOut, self.outputs.numClosedBedsStaffing, "Number of Closed Beds per time step due to Staffing")
            writeAvg(fOut, self.outputs.timesToIBA, "Time to Initial Bed Assignment per patient")
            fOut.write("\nTime to Initial Bed Assignment\tgeometric mean\tstd dev)")
            fOut.write("\n\t{}\t{}".format(*gmean(self.outputs.timesToIBA)))
                       
            fOut.write("\nDistribution of time to Initial Bed Assignment (hrs)")
            bins = range(25)+[1000]
            fOut.write("\n"+"\t".join([str(i) for i in range(25)])+"\t24+")
            dist, b = np.histogram(self.outputs.timesToIBA, bins=bins)
            fOut.write("\n"+"\t".join([str(i) for i in dist]))
            writeAvg(fOut, self.outputs.numBedMovesAcuity, "Number of Bed Moves due to Acuity per time step")

            fOut.write("\nNumber of Bed Moves due to Acuity by Day of Week")
            fOut.write("\n\tMon\tTue\tWed\tThu\tFri\tSat\tSun\n")
            for day in range(1,8):
                fOut.write("\t{}".format(self.outputs.numBedMovesAcuityByDay[day]))
            
            writeAvg(fOut, self.outputs.numBedMovesCP, "Number of Bed Moves due to CP per time step")
            writeAvg(fOut, self.outputs.numBedMovesShuffle, "Number of Bed Moves due to Shuffling per time step")


            writeAvg(fOut, self.outputs.numMismatchAcuityRooms, "Number of Rooms with Acuity Mismatch per time step")
            writeAvg(fOut, self.outputs.numMismatchAcuityBedsLow, "Number of Beds with Acuity Mismatch per time step (Low looking for High)")
            writeAvg(fOut, self.outputs.numMismatchAcuityBedsHigh, "Number of Beds with Acuity Mismatch per time step (High looking for Low)")
            fOut.write("\n\tMRSA Only\tVRE Only\tBoth\tTotal")
            mismatches = {}
            for category in (PSTATUS_MRSA, PSTATUS_VRE, PSTATUS_BOTH):
                mismatches[category] = [mismatch[category] for mismatch in self.outputs.numMismatchCPRooms]
            fOut.write("\nPStatus Mismatched Rooms (rooms per time step)\t{}\t{}\t{}\t{}".format(np.mean(mismatches[PSTATUS_MRSA]),
                                                                                                np.mean(mismatches[PSTATUS_VRE]),
                                                                                                np.mean(mismatches[PSTATUS_BOTH]),
                                                                                                np.mean(mismatches[PSTATUS_MRSA])+np.mean(mismatches[PSTATUS_VRE])+np.mean(mismatches[PSTATUS_BOTH])))

            mismatches = {}
            for category in (PSTATUS_MRSA, PSTATUS_VRE, PSTATUS_BOTH):
                mismatches[category] = [mismatch[category] for mismatch in self.outputs.numMismatchCPBeds]
            fOut.write("\nPStatus Mismatched Beds (beds per time step)\t{}\t{}\t{}\t{}".format(np.mean(mismatches[PSTATUS_MRSA]),
                                                                                                np.mean(mismatches[PSTATUS_VRE]),
                                                                                                np.mean(mismatches[PSTATUS_BOTH]),
                                                                                                np.mean(mismatches[PSTATUS_MRSA])+np.mean(mismatches[PSTATUS_VRE])+np.mean(mismatches[PSTATUS_BOTH])))


            fOut.write("\n\tMRSA\tVRE")
            fOut.write("\nNumber Truly Colonized (prevalence) per time step\t{}\t{}".format(np.mean(self.outputs.numInfected[MRSA]), np.mean(self.outputs.numInfected[VRE])))
            fOut.write("\nNumber Newly Truly Colonized(incidence) due to Background per time step\t{}\t{}".format(np.mean(self.outputs.numNewInfectedBackground[MRSA]), np.mean(self.outputs.numNewInfectedBackground[VRE])))
            fOut.write("\nNumber Newly Truly Colonized(incidence) due to false Cohorting per time step\t{}\t{}".format(np.mean(self.outputs.numNewInfectedCohorting[MRSA]), np.mean(self.outputs.numNewInfectedCohorting[VRE])))

            MRSAExposures = [exposure[EXPOSURE_MRSA] for exposure in self.outputs.timeSpentExposed if exposure[EXPOSURE_MRSA]]
            VREExposures = [exposure[EXPOSURE_VRE] for exposure in self.outputs.timeSpentExposed if exposure[EXPOSURE_VRE]]
            bothExposures = [exposure[EXPOSURE_BOTH] for exposure in self.outputs.timeSpentExposed if exposure[EXPOSURE_BOTH]]
            numLeaving = len(self.outputs.timeSpentExposed)
            
            fOut.write("\n\tMRSA\tVRE\tBoth")
            fOut.write("\nProportion of Patients Ever Exposed to Colonization\t{}\t{}\t{}".format(len(MRSAExposures)/float(numLeaving),len(VREExposures)/float(numLeaving),len(bothExposures)/float(numLeaving)))
            fOut.write("\nTotal Patient Hours Exposed to Colonization\t{}\t{}\t{}".format(sum(MRSAExposures),sum(VREExposures),sum(bothExposures)))
            fOut.write("\nAvg Length of Exposure to Colonization\t{}\t{}\t{}".format(np.mean(MRSAExposures),np.mean(VREExposures),np.mean(bothExposures)))
            
            fOut.write("\n\tMRSA\tVRE")
            fOut.write("\nNumber of Flags Aquired\t{}\t{}".format(sum(self.outputs.numFlagsAquired[MRSA]), sum(self.outputs.numFlagsAquired[VRE])))
            fOut.write("\nNumber of Flags Lost\t{}\t{}".format(sum(self.outputs.numFlagsLost[MRSA]), sum(self.outputs.numFlagsLost[VRE])))
            fOut.write("\nTotal Number Observed Non-Colonized/Truly Colonized at Admit\t{}\t{}".format(sum(self.outputs.numUnflaggedCarrierAtAdmit[MRSA]), sum(self.outputs.numUnflaggedCarrierAtAdmit[VRE])))
            fOut.write("\nProportion Observed Non-Colonized/Truly Colonized at Admit\t{}\t{}".format(sum(self.outputs.numUnflaggedCarrierAtAdmit[MRSA])/float(sum(self.outputs.numNewPatients)), sum(self.outputs.numUnflaggedCarrierAtAdmit[VRE])/float(sum(self.outputs.numNewPatients))))
            fOut.write("\nTotal Number Observed Colonized/Truly Non-Colonized at Admit\t{}\t{}".format(sum(self.outputs.numFlaggedNonCarrierAtAdmit[MRSA]), sum(self.outputs.numFlaggedNonCarrierAtAdmit[VRE])))
            fOut.write("\nProportion Observed Colonized/Truly Non-Colonized at Admit\t{}\t{}".format(sum(self.outputs.numFlaggedNonCarrierAtAdmit[MRSA])/float(sum(self.outputs.numNewPatients)), sum(self.outputs.numFlaggedNonCarrierAtAdmit[VRE])/float(sum(self.outputs.numNewPatients))))

            fOut.write("\nNumber Observed Non-Colonized at Admit/Observed Colonized at Discharge\t{}\t{}".format(sum(self.outputs.numDiffObservedCPDischarge[MRSA]),sum(self.outputs.numDiffObservedCPDischarge[VRE])))

            writeAvg(fOut, self.outputs.numTestsGiven[MRSA], "Number of MRSA Tests Given per time step")
            writeAvg(fOut, self.outputs.numTestsGiven[VRE], "Number of VRE Tests Given per time step")
            writeAvg(fOut, self.outputs.numTestsGivenUnique[MRSA], "Number of MRSA Tests Given (never tested) per time step")
            writeAvg(fOut, self.outputs.numTestsGivenUnique[VRE], "Number of VRE Tests Given (never tested) per time step")
            
            fOut.write("\nTest Results\tMRSA Test Result\t\tVRE Test Result")
            fOut.write("\nTrue Status\tTrue\tFalse\tTrue\tFalse")
            fOut.write("\nTrue\t{}\t{}\t{}\t{}".format(sum(self.outputs.testResults[(MRSA,True,True)]),sum(self.outputs.testResults[(MRSA,True,False)]),sum(self.outputs.testResults[(VRE,True,True)]),sum(self.outputs.testResults[(VRE,True,False)])))
            fOut.write("\nFalse\t{}\t{}\t{}\t{}".format(sum(self.outputs.testResults[(MRSA,False,True)]),sum(self.outputs.testResults[(MRSA,False,False)]),sum(self.outputs.testResults[(VRE,False,True)]),sum(self.outputs.testResults[(VRE,False,False)])))
          
            #time based statistics
            numTimeSteps = len(self.outputs.numPatientsInHospital)
            fOut.write("\n")
            fOut.write("\nHourly Statistics")
            for t in range(numTimeSteps):
                fOut.write("\nTime Step {}".format(t))
                fOut.write("\n\tDate\t{}".format(self.inputs.startTime + t * self.timeStep))
                fOut.write("\n\tNumber of Patients\t{}".format(self.outputs.numPatientsInHospital[t]))
                fOut.write("\n\tNumber of New Patients\t{}".format(self.outputs.numNewPatients[t]))
                fOut.write("\n\t\tHome\tFacility\tDeath\tTotal")
                fOut.write("\n\tNumber of Patients Leaving\t{}\t{}\t{}\t{}".format(self.outputs.numPatientsLeaving[DISCHARGE_HOME][t], self.outputs.numPatientsLeaving[DISCHARGE_FACILITY][t],self.outputs.numPatientsLeaving[DEAD][t],self.outputs.numPatientsLeaving[DISCHARGE_HOME][t]+ self.outputs.numPatientsLeaving[DISCHARGE_FACILITY][t]+self.outputs.numPatientsLeaving[DEAD][t]))
                fOut.write("\n\t\tTotal\tWaiting for first bed\tWaiting for transfer")
                fOut.write("\n\tNumber of People Queuing\t{}\t{}\t{}".format(self.outputs.numQueuing[t],
                                                                                 self.outputs.numQueuingReal[t],
                                                                                 self.outputs.numQueuingPhantom[t]))
                
                fOut.write("\n\tLiscenced Bed Capacity\t{}".format(self.outputs.bedCapacity[t]))
                fOut.write("\n\tOperational Capacity\t{}".format(self.outputs.operationalCapacity[t]))

                fOut.write("\n\t\tNone\tMRSA Only\tVRE Only\tBoth")
                fOut.write("\n\tNumber of Filled Beds by CP\t{}\t{}\t{}\t{}".format(self.outputs.numBedsCP[PSTATUS_NONE][t],
                                                                                     self.outputs.numBedsCP[PSTATUS_MRSA][t],
                                                                                     self.outputs.numBedsCP[PSTATUS_VRE][t],
                                                                                     self.outputs.numBedsCP[PSTATUS_BOTH][t],
                                                                                     ))

                
                fOut.write("\n\tNumber Closed Beds\t{}".format(self.outputs.numClosedBeds[t]))
                fOut.write("\n\tNumber Closed Beds CP\t{}".format(self.outputs.numClosedBedsCP[t]))
                fOut.write("\n\tNumber Closed Beds Staffing\t{}".format(self.outputs.numClosedBedsStaffing[t]))
                fOut.write("\n\tNumber of Bed Moves due to Acuity\t{}".format(self.outputs.numBedMovesAcuity[t]))
                fOut.write("\n\tNumber of Bed Moves due to CP\t{}".format(self.outputs.numBedMovesCP[t]))
                fOut.write("\n\tNumber of Bed Moves due to Shuffling\t{}".format(self.outputs.numBedMovesShuffle[t]))
                fOut.write("\n\tNumber of Mismatched Rooms Due to Acuity\t{}".format(self.outputs.numMismatchAcuityRooms[t]))
                fOut.write("\n\tNumber of Mismatched Beds Due to Acuity\t{}".format(self.outputs.numMismatchAcuityPatients[t]))
                fOut.write("\n\tNumber of Mismatched Beds Due to Acuity (Low looking for High)\t{}".format(self.outputs.numMismatchAcuityBedsLow[t]))
                fOut.write("\n\tNumber of Mismatched Beds Due to Acuity (High looking for Low)\t{}".format(self.outputs.numMismatchAcuityBedsHigh[t]))
                
                fOut.write("\n\t\tMRSA Only\tVRE Only\tBoth\tTotal")
                mismatches = self.outputs.numMismatchCPRooms[t]
                fOut.write("\n\tNumber of Mismatched Rooms By True CP\t{}\t{}\t{}\t{}".format(mismatches[PSTATUS_MRSA], mismatches[PSTATUS_VRE], mismatches[PSTATUS_BOTH], sum(mismatches[1:])))
                fOut.write("\n\t\tMRSA\tVRE")
                fOut.write("\n\tNumber of Flags Aquired\t{}\t{}".format(self.outputs.numFlagsAquired[MRSA][t],self.outputs.numFlagsAquired[VRE][t]))
                fOut.write("\n\tNumber of Flags Lost\t{}\t{}".format(self.outputs.numFlagsLost[MRSA][t],self.outputs.numFlagsLost[VRE][t]))
                fOut.write("\n\tTotal Number Observed Non-Colonized/Truly Colonized at Admit\t{}\t{}".format(self.outputs.numUnflaggedCarrierAtAdmit[MRSA][t],self.outputs.numUnflaggedCarrierAtAdmit[VRE][t]))
                fOut.write("\n\tTotal Number Observed Colonized/Truly Non-Colonized at Admit\t{}\t{}".format(self.outputs.numFlaggedNonCarrierAtAdmit[MRSA][t],self.outputs.numFlaggedNonCarrierAtAdmit[VRE][t]))
                fOut.write("\n\tNumber Observed Non-Colonized at Admit/Observed Colonized at Discharge\t{}\t{}".format(self.outputs.numDiffObservedCPDischarge[MRSA][t],self.outputs.numDiffObservedCPDischarge[VRE][t]))


                fOut.write("\n\tNumber of New Infections Background\t{}\t{}".format(self.outputs.numNewInfectedBackground[MRSA][t], self.outputs.numNewInfectedBackground[VRE][t]))
                fOut.write("\n\tNumber of New Infections Cohorting\t{}\t{}".format(self.outputs.numNewInfectedCohorting[MRSA][t], self.outputs.numNewInfectedCohorting[VRE][t]))
                fOut.write("\n\tNumber Infected\t{}\t{}".format(self.outputs.numInfected[MRSA][t], self.outputs.numInfected[VRE][t]))

                fOut.write("\n\t\tMRSA\tVRE")
                fOut.write("\n\tNumber of Tests Given\t{}\t{}".format(self.outputs.numTestsGiven[MRSA][t], self.outputs.numTestsGiven[VRE][t]))

                fOut.write("\n\tTest Results\tMRSA Test Result\t\tVRE Test Result")
                fOut.write("\n\tTrue Status\tTrue\tFalse\tTrue\tFalse")
                fOut.write("\n\tTrue\t{}\t{}\t{}\t{}".format(self.outputs.testResults[(MRSA,True,True)][t],self.outputs.testResults[(MRSA,True,False)][t],self.outputs.testResults[(VRE,True,True)][t],self.outputs.testResults[(VRE,True,False)][t]))
                fOut.write("\n\tFalse\t{}\t{}\t{}\t{}".format(self.outputs.testResults[(MRSA,False,True)][t],self.outputs.testResults[(MRSA,False,False)][t],self.outputs.testResults[(VRE,False,True)][t],self.outputs.testResults[(VRE,False,False)][t]))
    #Saves a summary file 
    def save_summary_file(self, summaryFilePath, outputsDict):        
        with open(summaryFilePath, 'w') as fOut:
            fOut.write("\tMean Operational Capacity per time step(mean)\tMean Operational Capacity per time step(geo mean)\tMean Operational Capacity per time step(std dev)")
            fOut.write("\tLicensed Bed Capacity per time step(mean)\tLicensed Bed Capacity per time step(geo mean\tLicensed Bed Capacity per time step(std dev)")
            fOut.write("\tLOS per patient(mean)\tLOS per patient(geo mean)\tLOS per patient(std dev)\tLOS per patient(geom std dev)")
            fOut.write("\tTime to IBA per patient(mean)\tTime to IBA per patient(geo mean)\tTime to IBA per patient(std dev)\tTime to IBA per patient(geom std dev)")
            fOut.write("\tNum People Leaving per time step\tNum People Leaving to Home\tNum People Leaving to Facility\tNum People Leaving to Death")
            fOut.write("\tNumber of Mismatched Beds Due to Acuity (total)")
            fOut.write("\tNumber of Mismatched Beds Due to Acuity (Low looking for High)")
            fOut.write("\tNumber of Mismatched Beds Due to Acuity (High looking for Low)")
            fOut.write("\tPStatus Mismatched Rooms (rooms per time step)\tPStatus Mismatch MRSA Only\tPStatus Mismatch VRE Only\tPstatus Mismatch Both")
            fOut.write("\tNum New MRSA cases due to Background per time step\tNum New VRE cases due to Background per time step\tNum New MRSA cases due to False Cohorting per time step\tNum New VRE Cases due to False Cohorting per time step")
            fOut.write("\tClosed Beds(mean)\tClosed Beds(std dev)")
            fOut.write("\tClosed Beds due to non CP/acuity (mean)\tClosed Beds due to non CP/acuity (std dev)")
            
            fOut.write("\tBed Moves Due to Acuity (mean)\tBed Moves Due to Acuity (std dev)")
            fOut.write("\tBed Moves Due to Acuity by Day")
            fOut.write("\tMon\tTue\tWed\tThu\tFri\tSat\tSun")
            fOut.write("\tPatient Hours By Acuity")
            fOut.write("\t"+"\t".join([str(acuity) for acuity in ACUITY_VALUES]))

            for fileName, output in outputsDict.iteritems():
                fOut.write("\n{}".format(fileName))
                fOut.write("\t{}\t{}\t{}".format(np.mean(output.operationalCapacity),gmean(output.operationalCapacity)[0], np.std(output.operationalCapacity)))
                fOut.write("\t{}\t{}\t{}".format(np.mean(output.bedCapacity),gmean(output.bedCapacity)[0], np.std(output.bedCapacity)))
                fOut.write("\t{}\t{}\t{}\t{}".format(np.mean(output.lengthsOfStay),gmean(output.lengthsOfStay)[0], np.std(output.lengthsOfStay), gmean(output.lengthsOfStay)[1]))
                fOut.write("\t{}\t{}\t{}\t{}".format(np.mean(output.timesToIBA),gmean(output.timesToIBA)[0], np.std(output.timesToIBA), gmean(output.timesToIBA)[1]))
                
                fOut.write("\t{}\t{}\t{}\t{}".format(np.mean(sum([np.array(output.numPatientsLeaving[reason]) for reason in (DEAD, DISCHARGE_HOME, DISCHARGE_FACILITY)])),
                                                     np.mean(output.numPatientsLeaving[DISCHARGE_HOME]),
                                                     np.mean(output.numPatientsLeaving[DISCHARGE_FACILITY]),
                                                     np.mean(output.numPatientsLeaving[DEAD])))

                fOut.write("\t{}".format(np.mean(output.numMismatchAcuityPatients)))
                fOut.write("\t{}".format(np.mean(output.numMismatchAcuityBedsLow)))
                fOut.write("\t{}".format(np.mean(output.numMismatchAcuityBedsHigh)))
                mismatches = {}
                for category in (PSTATUS_MRSA, PSTATUS_VRE, PSTATUS_BOTH):
                    mismatches[category] = [mismatch[category] for mismatch in output.numMismatchCPRooms]
                fOut.write("\t{}\t{}\t{}\t{}".format(np.mean(mismatches[PSTATUS_MRSA])+np.mean(mismatches[PSTATUS_VRE])+np.mean(mismatches[PSTATUS_BOTH]),
                                                     np.mean(mismatches[PSTATUS_MRSA]),
                                                     np.mean(mismatches[PSTATUS_VRE]),
                                                     np.mean(mismatches[PSTATUS_BOTH]),))
                fOut.write("\t{}\t{}\t{}\t{}".format(np.mean(output.numNewInfectedBackground[MRSA]),
                                                     np.mean(output.numNewInfectedBackground[VRE]),
                                                     np.mean(output.numNewInfectedCohorting[MRSA]),
                                                     np.mean(output.numNewInfectedCohorting[VRE])))


                fOut.write("\t{}\t{}".format(np.mean(output.numClosedBeds), np.std(output.numClosedBeds)))
                fOut.write("\t{}\t{}".format(np.mean(output.numClosedBedsStaffing), np.std(output.numClosedBedsStaffing)))
                fOut.write("\t{}\t{}".format(np.mean(output.numBedMovesAcuity), np.std(output.numBedMovesAcuity)))
                fOut.write("\t")
                for day in range(1,8):
                    fOut.write("\t{}".format(output.numBedMovesAcuityByDay[day]))
                fOut.write("\t")
                for acuity in ACUITY_VALUES:
                    fOut.write("\t{}".format(np.mean(output.numPatientsInHospitalAcuity[acuity])))

    def load_initial_patients(self):
        if not self.SIM_INITIALIZED:
            raise
        self.INITIAL_PATIENTS_LOADED = True
        self.hospital.add_patients_to_hospital(self.initialPatients, self.currTimeStep)
    def initialize_run(self):
        self.currTimeStep = 0
        self.hospital.initialize_hospital()
        self.SIM_INITIALIZED = True

    #must be called after start_run
    def resume_run(self, timeSteps = -1):
        if not self.SIM_INITIALIZED:
            return False
        self.timeStepsRemaining = timeSteps
        #if timeSteps is set to -1 run for unlimited time
        while self.timeStepsRemaining:
            if self.PAUSE_SIM:
                break
            self.run_single_time_step()

        return True
    def run_single_time_step(self):
        if not self.SIM_INITIALIZED:
            return False
        if self.inputs.runInPeriod == self.currTimeStep:
            self.outputs = OutputStatistics()
        self.hospital.time_step(self.currTimeStep, self.getCurrDateTime(), self.outputs)
        self.currTimeStep += 1
        #self.print_time_step()
        self.accumulate_statistics()
        print self.get_statistics_text()

        if self.currTimeStep and not self.currTimeStep%self.inputs.saveEveryN:
            self.save_output_file("defaultSave.pfout")
        if self.timeStepsRemaining > 0:
            self.timeStepsRemaining -= 1

        return True
    def changeHospitalCPMatching(self, matchOnCP):
        matching = [matchOnCP, matchOnCP]
        self.hospital.settings.matchingProfile = matching
    #gets the current time as datetime object
    def getCurrDateTime(self):
        return self.inputs.startTime + self.currTimeStep * self.timeStep
    #calculates the statistics for that time step
    def accumulate_statistics(self):
        self.outputs.numPatientsInHospital.append(len(self.hospital.patientsInHospital))
        for acuity in ACUITY_VALUES:
            self.outputs.numPatientsInHospitalAcuity[acuity].append(0)
        for p in self.hospital.patientsInHospital.values():
            self.outputs.numPatientsInHospitalAcuity[p.acuity][-1]+=1
        numClosedBeds = 0
        numBeds = 0
        numInBeds = 0
        numCPMismatchRooms = [0 for status in PSTATUS_VALUES]
        numCPMismatchBeds = [0 for status in PSTATUS_VALUES]
        numAcuityMismatchPatients = 0
        numAcuityMismatchRooms = 0
        numAcuityMismatchBedsHigh = 0
        numAcuityMismatchBedsLow = 0
        numClosedByCategory = dict([(combo, 0) for combo in product((1,0),(1,0),(1,0))])
        numClosedBedsStaffing = 0
        numClosedCP = 0
        for pstatus in PSTATUS_VALUES:
            self.outputs.numBedsCP[pstatus].append(0)
        self.outputs.numClosedBedsGender = 0
        self.outputs.numClosedBedsMRSA = 0
        self.outputs.numClosedBedsVRE = 0
        self.outputs.numClosedBedsBoth = 0
        self.outputs.numQueuing.append(0)
        self.outputs.numQueuingReal.append(0)
        self.outputs.numQueuingPhantom.append(0)
        numInfected = [0,0]
        #loop through current patients
        for patient in self.hospital.patientsInHospital.values():
            for disease in DISEASES:
                numInfected[disease]+= int(patient.trueCP[disease])
        for disease in DISEASES:
            self.outputs.numInfected[disease].append(numInfected[disease])
        for room in self.hospital.layout.values():
            if isinstance(room, BedAllocationRoom) or isinstance(room, WardAllocationRoom):
                self.outputs.numQueuing[-1]+= len(room.patients)
                for patient in room.patients.values():
                    if isinstance(patient, PhantomPatient):
                        self.outputs.numQueuingPhantom[-1]+=1
                    else:
                        self.outputs.numQueuingReal[-1]+=1
            if isinstance(room, BedAllocationRoom):
                closedBedDict, numClosed, numClosedStaffing = room.get_num_closed_beds_by_category()
                numClosedBeds+=numClosed
                numClosedBedsStaffing += numClosedStaffing
                for combo, num in closedBedDict.iteritems():
                    numClosedByCategory[combo]+=num
                    if combo[0] == 1:
                        self.outputs.numClosedBedsGender += num
                    if combo[1] ==1:
                        self.outputs.numClosedBedsMRSA += num
                    if combo[2] == 1:
                        self.outputs.numClosedBedsVRE += num
                    if combo[1] == 1 and combo[2] ==1:
                        self.outputs.numClosedBedsBoth += num
                    if combo[1]==1 or combo[2]==1:
                        numClosedCP+=num
            if isinstance(room, BedRoom):
                numBeds+=room.capacity
                if room.capacity > 1:
                    mismatchStatus = room.get_pstatus_mismatch()
                    numCPMismatchRooms[mismatchStatus]+=1
                    numCPMismatchBeds[mismatchStatus]+=room.capacity

                    numLowMismatch, numHighMismatch = room.num_acuity_mismatch()
                    if numLowMismatch+numHighMismatch:
                        numAcuityMismatchRooms+=1
                    numAcuityMismatchPatients+=numLowMismatch+numHighMismatch
                        
                    numAcuityMismatchBedsLow+=numLowMismatch
                    numAcuityMismatchBedsHigh+=numHighMismatch

                        
                numInBeds+= room.get_num_patients()
                for patient in room.patients.values():
                    if not isinstance(patient, RealPatient):
                        continue
                    self.outputs.numBedsCP[patient.get_observed_pstatus()][-1]+=1
                    
            if isinstance(room, EndRoom):
                self.outputs.lengthsOfStay = room.lengthsOfStay
                
        self.outputs.numMismatchCPRooms.append(numCPMismatchRooms)
        self.outputs.numMismatchCPBeds.append(numCPMismatchBeds)
        self.outputs.numMismatchAcuityRooms.append(numAcuityMismatchRooms)
        self.outputs.numMismatchAcuityPatients.append(numAcuityMismatchPatients)
        self.outputs.numMismatchAcuityBedsLow.append(numAcuityMismatchBedsLow)
        self.outputs.numMismatchAcuityBedsHigh.append(numAcuityMismatchBedsHigh)
        #operational capacity only counts num closed due to cp or other, and not gender/acuity
        self.outputs.operationalCapacity.append(numInBeds/float(numBeds - numClosedCP - numClosedStaffing))
        self.outputs.numBeddedPatients.append(numInBeds)
        self.outputs.bedCapacity.append(numInBeds/float(numBeds))
        self.outputs.numClosedBeds.append(numClosedBeds)
        self.outputs.numClosedBedsByCategory.append(numClosedByCategory)
        self.outputs.numClosedBedsCP.append(numClosedCP)
        self.outputs.numClosedBedsStaffing.append(numClosedBedsStaffing)
    def get_statistics_text(self):
        statistics = ""
        statistics += "\n"+ "#"*20
        statistics += "\n"+  "Time Step {}".format(self.currTimeStep)
        statistics += "\n"+  "Date {}".format(self.getCurrDateTime())
        timesToIBA = []
        for room in self.hospital.layout.values():
            if isinstance(room, BedRoom):
                timesToIBA.extend(room.timesToInitialBedAssignment)
        
        timesToIBA = np.array(timesToIBA)
        numClosedBeds = np.array(self.outputs.numClosedBeds)
        if len(timesToIBA):
            timeToIBA = timesToIBA.mean(), timesToIBA.std()
        else:
            timeToIBA = (None, None)
        statistics += "\n"+  "Number of People: {}".format(len(self.hospital.patientsInHospital))
        statistics += "\n"+  "Time to Initial Bed Assignment: {}({})({})".format(timeToIBA[0],timeToIBA[1], np.median(timesToIBA))
        statistics += "\n"+  "Number Closed Beds Total: {}".format(numClosedBeds[-1])
        statistics += "\n"+  "\tNumber Closed Beds Gender: {}".format(self.outputs.numClosedBedsGender)
        statistics += "\n"+  "\tNumber Closed Beds MRSA: {}".format(self.outputs.numClosedBedsMRSA)
        statistics += "\n"+  "\tNumber Closed Beds VRE: {}".format(self.outputs.numClosedBedsVRE)
        statistics += "\n"+  "\tNumber Closed Beds MRSA and VRE: {}".format(self.outputs.numClosedBedsBoth)
        statistics += "\n"+  "\n"+  "Operational Capacity: {}".format(self.outputs.operationalCapacity[-1])
        statistics += "\n"+  "Num Rooms with mismatched Pstatus: {}".format(sum(self.outputs.numMismatchCPRooms[-1][1:]))
        statistics += "\n"+  "Num Closed Beds per time step: {} ({})".format(numClosedBeds.mean(), numClosedBeds.std())
        statistics += "\n"+  "\n"*2
        return statistics
    def print_time_step(self):
        print "#"*20
        print "Time {}".format(self.currTimeStep)
        self.hospital.print_hospital()
        print "#"*20+"\n"*2
        
if __name__ == "__main__":
    sim = Sim()

    sim.create_hospital()
    
    initialPatients = ((MALE,1,(True, False)),
                       (FEMALE,3,(True,True)),
                       (MALE,2,(False,False)))
    sim.import_initial_patients('l',initialPatients)

    layout = (("s","Start",(),('wa',)),
              ("wa","WardAllocation",(),('ba',)),
              ("ba","BedAllocation",((2,),),("b1","b2","b3")),
              ("b1","Bed",(1,),("e")),
              ("b2","Bed",(2,),("e")),
              ("b3","Bed",(2,),("e")),
              ("e","End",(),()))
    sim.import_layout('l', layout)

    sim.start_run(6)


    
