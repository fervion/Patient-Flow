# -*- coding: utf-8 -*-
"""
Simulation of patients in hospital
Name: hospitalBedSimulation.py
Author: Taige Hou (thou1@partners.org)
"""

from gui import *
from itertools import product
from random import choice, seed
def parse_patient_beds(bedFile):
    with open(bedFile) as fBed:
        bedData = [line.strip() for line in fBed.readlines()[1:]]
    beds = []
    unitAcuityServiceTable = {} #matches unit id to acuity service
    services = ("Surgical", "Medical", "Observation")
    
    for bedRow in bedData:
        cells = bedRow.split("\t")
        unitID = cells[4]
        acuity = cells[6]
        service = cells[5]
        roomID = cells[2]
        bedID =  cells[3]
        virtual = cells[11]
        if service not in services:
            continue
        if virtual != "N":
            continue
        if service and acuity:
            acuity = int(acuity)
            unitAcuityServiceTable[unitID]=(acuity, service)
            beds.append((unitID, roomID, bedID))
    return unitAcuityServiceTable, beds
def get_bedAllocation_id(acuity, service):
    if acuity == 4:
        return "ba4"
    if acuity == 2:
        return "ba2"
    serviceKey = {"Surgical":"s","Medical":"m", "Observation":"o"}
    return "ba{0}{1}".format(acuity, serviceKey[service])
if __name__== '__main__':
    unitTable, bedTable = parse_patient_beds("PatientBeds.txt")
    sim = Sim()
    sim.create_hospital()
    
    rooms = {} #database of rooms of lists of (unitid, capacity)
    for bed in bedTable:
        unitID, roomID, bedID = bed
        try:
            if rooms[roomID][1]>=2:
                #make second room
                roomID = roomID + "_2"
                #print rooms[roomID], roomID
            rooms[roomID][1]+=1

        except KeyError:
            rooms[roomID] = [unitID, 1]
    CPStatus = [(True, True), (True, False), (False, True), (False, False)]
    allocation = (SURGERY, MEDICINE)
    gender = (MALE, FEMALE)

    num_patients = 60
    numWithCP = 5
    
    initialPatients = [(choice(gender),4,MEDICINE, CPStatus[3], CPStatus[3]) for i in range(num_patients-numWithCP*3)]
    initialPatients.extend([(choice(gender),4,MEDICINE, CPStatus[0], CPStatus[0]) for i in range(numWithCP)])
    initialPatients.extend([(choice(gender),4,MEDICINE, CPStatus[1], CPStatus[1]) for i in range(numWithCP)])
    initialPatients.extend([(choice(gender),4,MEDICINE, CPStatus[2], CPStatus[2]) for i in range(numWithCP)])

    sim.import_initial_patients('l',initialPatients)
    

    
    layout = [("s","Start",("Start",),('wa',)),
              ("wa","WardAllocation",("Acu/Ser Alloc",),('ba5s','ba5m','ba4','ba3s','ba3m','ba2')),
              ("e","End",("End",),())]

    bedAllocations = (("ba5s", "BedAllocation", ((5,),(SURGERY,),"Sur ICU"), []),
                      ("ba5m", "BedAllocation", ((5,),(MEDICINE,),"Med ICU"), []),
                      ("ba4", "BedAllocation", ((4,),(MEDICINE,),"Step Down"), []),
                      ("ba3s", "BedAllocation", ((3,),(SURGERY,),"Sur General"), []),
                      ("ba3m", "BedAllocation", ((3,),(MEDICINE,),"Med General"), []),
                      ("ba2", "BedAllocation", ((2,),(SURGERY,MEDICINE),"Observation"), []),
                      )

    bedRooms = []
    count = 0
    for room in rooms:
        roomID = room
        unitID, capacity = rooms[room]
        bedRooms.append((roomID,"Bed", (capacity,roomID,unitID),()))
        acuity, service = unitTable[unitID]
        baID = get_bedAllocation_id(acuity, service)
        for baRoom in bedAllocations:
            if baRoom[0] == baID:
                baRoom[3].append(roomID)
                break
        
        count+=1

    layout.extend(bedAllocations)
    layout.extend(bedRooms)
    
    sim.import_layout('l', layout)
    app=wx.App()

    frame = MainWindow(None, "Patient Flow", sim)
    frame.Show()
    app.MainLoop()

