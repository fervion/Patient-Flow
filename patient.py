from copy import copy
from abc import ABCMeta, abstractmethod
import random
"""
Contains the class for managing people that flow through the model
"""

MALE="Male"
FEMALE="Female"
GENDERS = (MALE, FEMALE)
MRSA=0
VRE=1
DISEASES = MRSA, VRE
DISEASE_STRINGS = ["MRSA", "VRE"]

PSTATUS_NONE, PSTATUS_MRSA, PSTATUS_VRE, PSTATUS_BOTH = (0,1,2,3)
EXPOSURE_MRSA, EXPOSURE_VRE, EXPOSURE_BOTH = (0,1,2)
PSTATUS_VALUES = (PSTATUS_NONE, PSTATUS_MRSA, PSTATUS_VRE, PSTATUS_BOTH)
SURGERY="s"
MEDICINE="m"
REASON_ACUITY, REASON_CP = (0,1)
ACTIVE, DEAD, DISCHARGE_HOME, DISCHARGE_FACILITY = range(4)
ACUITY_VALUES = (2,3,4,5)
ACUITY_VALUES_LEAVE = (6,7,8)
ACUITY_VALUES_GENERAL = (2,3,4,5,6,7,8) #acuity values including 6-discharge 7-death
ACUITY_DISCHARGE_HOME = 6
ACUITY_DISCHARGE_FACILITY = 7
ACUITY_DEATH = 8

ACUITY_INDICES = range(len(ACUITY_VALUES))
SERVICE_VALUES = (SURGERY, MEDICINE) #allowed service values

#convert a pstatus to a tuple of [MRSA, VRE] booleans
def convertPStatusToCP(pStatus):
    hasMRSA = pStatus == PSTATUS_MRSA or pStatus == PSTATUS_BOTH
    hasVRE = pStatus == PSTATUS_VRE or pStatus == PSTATUS_BOTH
    return (hasMRSA, hasVRE)

#The base class for all patients/persons/place holders in the hospital.  Should never call this class directly.
class BasePerson(object):
    __metaclass__ = ABCMeta
    
    #unique id for each Patient object
    idCounter=0

    @staticmethod
    def reset_counter():
        BasePerson.idCounter = 0

    def __init__(self):
        self.id = BasePerson.idCounter
        BasePerson.idCounter += 1
        self.currRoom = None
    def __str__(self,printID = None):
        if printID is None:
            printID = self.id
        return self.prefix + " {pid}\n\tCurr Room:{roomName}\n\tCurr Room ID:{roomID}"\
        .format(roomName = self.currRoom.name,
                roomID = self.currRoom.id,
                pid=printID,)
   
#The abstract class for all patients (real or phantom) in the hospital.  Should never instance this class directly.
class BasePatient(BasePerson):
    __metaclass__ = ABCMeta

    def __init__(self, gender, acuity, service, trueCPStatus, observedCPStatus):
        super(BasePatient,self).__init__()
        self.gender = gender
        self.acuity = acuity
        self.service = service

        #cp status is a tuple of booleans with (hasMRSA, hasVRE)
        self.trueCP = list(trueCPStatus)
        self.observedCP = list(observedCPStatus)
        self.originalObservedCP = list(observedCPStatus)
        self.CPHistory = list(observedCPStatus)
        

        self.timeStepLastRoomChange = None #time that they last moved rooms of any kind

        self.timeEnteringHospital = -1

        self.status = ACTIVE #patient can either be active, dead, or discharge
        
    def has_MRSA(self):
        return self.trueCP[MRSA]
    def has_VRE(self):
        return self.trueCP[VRE]
    def has_observed_MRSA(self):
        return self.observedCP[MRSA]
    def has_observed_VRE(self):
        return self.observedCP[VRE]
    def __str__(self,printID = None):
        if printID is None:
            printID = self.id
        return "{gender} {pid}\n\tCurr Room:{roomName}\n\tCurr Room ID:{roomID}\n\tRoomAcuity:{roomAcuity}\n\tAcuity:{acuity}\n\tService:{service}\n\tTrue MRSA:{trueMRSA}\n\tTrue VRE:{trueVRE}\n\tObsv MRSA:{obsvMRSA}\n\tObsv Vre:{obsvVRE}"\
        .format(gender=self.gender,
                roomName = self.currRoom.name,
                roomID = self.currRoom.id,
                roomAcuity = getattr(self.currRoom, "allowedAcuity", "None"),
                pid=printID,
                acuity=self.acuity,
                service=self.service,
                trueMRSA=self.has_MRSA(),
                trueVRE=self.has_VRE(),
                obsvMRSA=self.has_observed_MRSA(),
                obsvVRE=self.has_observed_VRE())
    
    #checks if this patient matches another patient based on gender and observed cp status
    def matches(self, otherPatient):            
        if isinstance(self, CleaningPatient) or isinstance(otherPatient, CleaningPatient):
            return True
        if self.gender != otherPatient.gender:
            return False
        if self.observedCP != otherPatient.observedCP:
            return False
        return True
#The class for all patients physically in the hospital
class RealPatient(BasePatient): 
    def __init__(self, gender, acuity, service, trueCPStatus, observedCPStatus, workUpTime):
        super(RealPatient,self).__init__(gender, acuity, service, trueCPStatus, observedCPStatus)

        #test results in history stratified by disease
        self.testResults = [[],[]]
        #list of tests that are waiting.  stored as tuple of (disease, result, timeToGetResult)
        self.waitingTests = []
        #if patient has ever been tested strat by disease
        self.hasBeenTested = [False, False]
        self.timeStepsExposedToInfection = [0,0,0] #time spent as negative falsly cohorted with positive roommate (MRSA, VRE, BOTH) Not mutally exclusive
        
        self.timeOfLastAcuityChange = -1
        self.timeOfDeathDischarge = -1

        self.realPatient = self
        self.queueingPatient = None
        self.delayPatient = None
        
        self.isIncident = True #patient can either be incident or prevalent

        self.hasBeenAllocated = False #true if patient has ever been allocated to a bed
        self.timeToInitialBedAssignment = 0 #time spent in acuity/ser allocation rooms before being allocated first bed
        self.workUpTime = workUpTime #time spent queueing before patient can be assigned a bed
        self.prefix = ""
        
    #change acuity.  If new acuity is death or discharge sets status to dead or discharged
    def change_acuity(self, newAcuity, time):
        if newAcuity == ACUITY_DISCHARGE_HOME:
            self.status = DISCHARGE_HOME
            self.process_death_discharge(time)
            return
        if newAcuity == ACUITY_DISCHARGE_FACILITY:
            self.status = DISCHARGE_FACILITY
            self.process_death_discharge(time)
            return
        if newAcuity == ACUITY_DEATH:
            self.status = DEAD
            self.process_death_discharge(time)
            return
        if self.acuity!=newAcuity:
            timeOfLastAcuityChange = time
            self.acuity = newAcuity
            if self.queueingPatient:
                self.queueingPatient.acuity = newAcuity
            if self.delayPatient:
		self.delayPatient.acuity = newAcuity
    #change observed status
    def change_observed(self, disease, newStatus):
        self.observedCP[disease] = newStatus
        self.CPHistory[disease] = self.CPHistory[disease] or newStatus
        if self.queueingPatient:
            self.queueingPatient.change_observed(disease, newStatus)
	if self.delayPatient:
            self.delayPatient.change_observed(disease, newStatus)
    def process_death_discharge(self, time):
        self.timeOfDeathDischarge = time
        #remove phantoms if exists
        if self.queueingPatient:
            self.queueingPatient.currRoom.remove_patient(self.queueingPatient)
            self.queueingPatient = None
        if self.delayPatient:
	    self.delayPatient.currRoom.remove_patient(self.delayPatient)
            self.delayPatient = None
    #person becomes infected with MRSA or VRE through background transmission or false -cohorting
    def become_infected(self, disease):
        self.trueCP[disease] = True
        if self.queueingPatient:
            self.queueingPatient.trueCP[disease] = True
        if self.delayPatient:
            self.delayPatient.trueCP[disease] = True
    def get_observed_pstatus(self):
        obsMrsa = self.has_observed_MRSA()
        obsVre = self.has_observed_VRE()
        if obsMrsa and obsVre:
            return PSTATUS_BOTH
        if obsMrsa:
            return PSTATUS_MRSA
        if obsVre:
            return PSTATUS_VRE
        return PSTATUS_NONE
    #perform testing on patient for MRSA/VRE based on testing inputs
    def do_test(self, disease, testInputs, currTime, simOutputs):
        if self.trueCP[disease]:
            prob = testInputs.sensitivity
        else:
            prob = testInputs.specificity
        result = self.trueCP[disease]
        if random.random() > prob:
            result = not result
        self.waitingTests.append((disease, result, currTime+ testInputs.tat))
        simOutputs.numTestsGiven[disease][-1]+=1
        if not self.hasBeenTested[disease]:
            simOutputs.numTestsGivenUnique[disease][-1]+=1
        simOutputs.testResults[(disease,self.trueCP[disease],result)][-1]+=1
        self.hasBeenTested[disease] = True

        return result


#Phantom Patient is used to wait in queue for patients that are in a bed but need to look for a bed move
class PhantomPatient(BasePatient):
    def __init__(self, patient):
        assert isinstance(patient, RealPatient)    
        super(PhantomPatient,self).__init__(patient.gender,patient.acuity, patient.service, patient.trueCP, patient.observedCP)
        #The RealPatient that the PhantomPatientMirrors
        self.realPatient = patient
        self.prefix = "Phantom"
    def __str__(self):
        return self.prefix + " " +super(PhantomPatient,self).__str__(self.realPatient.id)
    def change_observed(self, disease, newStatus):
        self.observedCP[disease] = newStatus
        self.CPHistory[disease] = self.CPHistory[disease] or newStatus

#Queueing Patient is used to wait in queue for patients that are in a bed but need to look for a bed move
class QueueingPatient(PhantomPatient):
    def __init__(self, patient, reason):
        super(QueueingPatient,self).__init__(patient)
        self.prefix = "Queueing"
        self.reason = reason
        patient.queueingPatient = self

#Delay Patient is used to put a placeholder in a bed to cause a delay for the patient reaching the bed
class DelayPatient(PhantomPatient):
    def __init__(self, patient, timeStepToEndDelay):
        super(DelayPatient,self).__init__(patient)
        self.timeStepToEndDelay = timeStepToEndDelay
        self.prefix = "Delay"
        patient.delayPatient = self
        
#Cleaning patient is used to block off beds that are being cleaned after a patient has left
class CleaningPatient(BasePerson):
    def __init__(self,timeStepToEndCleaning):
        super(CleaningPatient,self).__init__()
        self.timeStepToEndCleaning = timeStepToEndCleaning
        self.prefix = "Cleaning"
    def __str__(self):
        return "Cleaning "+super(CleaningPatient,self).__str__(self.id)
if __name__ == '__main__':
    patients=[]
    for i in range(1000):
        patients.append(RealPatient(MALE))
    print patients[99]
    Patient.reset_counter()
    p = Patient(MALE)
    p = Patient(FEMALE)
    print p
    for i in range(1000):
        patients.append(Patient(MALE))
    print patients[-1]
