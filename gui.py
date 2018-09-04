# -*- coding: utf-8 -*-
"""
Simulation of patients in hospital
Name: gui.py
Author: Taige Hou (thou1@partners.org)
"""

import wx
from wx.lib.floatcanvas import NavCanvas, FloatCanvas
from patientFlowSim import *
from collections import deque
from time import sleep
from numpy import mean
import os, threading
from os import path
from glob import glob
from random import choice, seed
from copy import deepcopy
#Dictionary mapping room type to how to draw the room
#Tuple of (width, height, linewidth, linecolor, fillcolor)
WIDTH, HEIGHT, LINEWIDTH, LINECOLOR, FILLCOLOR = range(5)
GENERAL_SETTINGS = {'label_size':8,
                    }
ROOM_SETTINGS= {'room_mapping':{StartRoom: (40,50,2,"Gold",None),
                    EndRoom: (40,50,2,"Firebrick",None),
                    WardAllocationRoom: (40,50,2,"Black",None),
                    BedAllocationRoom: (40,50,2,"Black",None),
                    BedRoom: (27,20,2,"Black",None),},
                'name_text_size':3.8,
                'room_data_text_size':3,
                'double_room_line_width':2,
                'double_room_line_color':"Black",
                'unit_room_line_width':2,
                'unit_room_line_color':"Light Steel Blue",
                'unit_room_border':(10, 8),
                'max_display_patients':9,
                }
PATIENT_SETTINGS = {'radius': 4,
                    'line_width': 2,
                    'line_color_female': "Orange Red",
                    'line_color_male': "Sky Blue",
                    'line_color_cleaning': "Black",
                    'fill_color_phantom': "Grey",
                    'fill_color_cleaning': "Green",
                    'fill_color': None,
                    'num_per_row': 3,
                    'id_text_size': 4,
                    'acuity_text_size': 2.5,
                    }
class RunThread(threading.Thread):
    def __init__(self, mainWindow, n):
        threading.Thread.__init__(self)
        self.mainWindow = mainWindow
        self.n = n
    def run(self):
        self.mainWindow.OnStepMany(None, self.n)
        
class MainWindow(wx.Frame):
    elapsedTime = 0
    timeStepTime = 0

            
    def __init__(self, parent, title, sim):
        wx.Frame.__init__(self, parent, title=title, size=(950,800))

        self.sim = sim
        self.runThread = None
        layout =self.sim.hospital.layout
        startRoomID, endRoomID = sim.hospital.validate_layout(layout)
        self.auto_layout(layout, layout[startRoomID], layout[endRoomID])

        #Menu bar
        MenuBar = wx.MenuBar()
        
        file_menu = wx.Menu()
        item = file_menu.Append(-1, "&Load Input File", "Load Inputs")
        self.Bind(wx.EVT_MENU, self.OnLoadInputFile, item)
        item = file_menu.Append(-1, "Load Initial &Patients","Load")
        self.Bind(wx.EVT_MENU, self.OnLoadPatients, item)
        item = file_menu.Append(-1, "Run &Batch","Batch")
        self.Bind(wx.EVT_MENU, self.OnRunBatch, item)
        item = file_menu.Append(-1, "Write Output File","Write Output")
        self.Bind(wx.EVT_MENU, self.OnWriteOutput, item)

        item = file_menu.Append(-1, "&SavePNG","Save the current image as a PNG")
        self.Bind(wx.EVT_MENU, self.OnSavePNG, item)
        item = file_menu.Append(-1, "&Close","Close this frame")
        self.Bind(wx.EVT_MENU, self.OnQuit, item)
        MenuBar.Append(file_menu, "&File")

        settings_menu = wx.Menu()
        self.matchCPMI = settings_menu.Append(-1, "&Match on CP status","Match patients based on observed CP status", kind = wx.ITEM_CHECK)
        self.matchCPMI.Check(True)
        self.Bind(wx.EVT_MENU, self.OnMatchCP, self.matchCPMI)
        MenuBar.Append(settings_menu, "&Settings")

        self.SetMenuBar(MenuBar)
        self.CreateStatusBar()
        
        # Add the Canvas
        NC = NavCanvas.NavCanvas(self,
                                Debug = 0,
                                BackgroundColor = "White")
        self.canvas = NC.Canvas

        #tool bar
        tb = NC.ToolBar
        tb.AddSeparator()

        Button = wx.Button(tb, wx.ID_ANY, "Step")
        tb.AddControl(Button)
        Button.Bind(wx.EVT_BUTTON,self.OnStep)

        Button = wx.Button(tb, wx.ID_ANY, "Step 10")
        tb.AddControl(Button)
        Button.Bind(wx.EVT_BUTTON,self.OnStepMany)

        Button = wx.Button(tb, wx.ID_ANY, "Step n")
        tb.AddControl(Button)
        Button.Bind(wx.EVT_BUTTON,self.OnStepN)

        self.stepNTc = wx.TextCtrl(tb, -1, size=(60,-1))
        tb.AddControl(self.stepNTc)
        
        tb.Realize()
        
        self.patientGroup = FloatCanvas.Group()
        self.canvas.AddObject(self.patientGroup)
        self.labelGroup = FloatCanvas.Group()
        self.canvas.AddObject(self.labelGroup)
        self.draw_layout()

        #call zoom to fit seems to only work as an event
        wx.PostEvent(NC.ZoomButton, wx.PyCommandEvent(wx.EVT_BUTTON.typeId))
    def draw_layout(self):
        wx.GetApp().Yield(True)
        Canvas = self.canvas
        Canvas.InitAll()
        for room in self.sim.hospital.layout.values():
            settings = ROOM_SETTINGS['room_mapping'][type(room)]
            roomWidth = settings[WIDTH]
            
            if isinstance(room, BedRoom):
                roomWidth = settings[WIDTH]*(2-1/float(room.capacity))
                #add line for double rooms
                if room.capacity == 2:
                    self.canvas.AddLine(((room.guiCoords[0],room.guiCoords[1]-settings[HEIGHT]/2.0),(room.guiCoords[0],room.guiCoords[1]+settings[HEIGHT]/2.0)),LineWidth = ROOM_SETTINGS['double_room_line_width'],LineColor=ROOM_SETTINGS['double_room_line_color'])    

            #draw unit if bed allocation room
            if isinstance(room, BedAllocationRoom):
                for unitName, unit in room.units.items():
                    border = ROOM_SETTINGS['unit_room_border']
                    left = min([bRoom.guiCoords[0]-self.get_room_dim(bRoom)[WIDTH]/2.0 for bRoom in unit])
                    right = max([bRoom.guiCoords[0]+self.get_room_dim(bRoom)[WIDTH]/2.0 for bRoom in unit])
                    bottom = min([bRoom.guiCoords[1]-self.get_room_dim(bRoom)[HEIGHT]/2.0 for bRoom in unit])
                    top = max([bRoom.guiCoords[1]+self.get_room_dim(bRoom)[HEIGHT]/2.0 for bRoom in unit])
                    dim = (right - left+border[0]*2, top - bottom+2*border[1])
                    self.canvas.AddRectangle((left-border[0], bottom-border[1]), dim, LineWidth=ROOM_SETTINGS['unit_room_line_width'], LineColor=ROOM_SETTINGS['unit_room_line_color'])
                    #unit name
                    Canvas.AddScaledText(unitName,((right+left)/2.0, bottom-border[1]),Position="tc", Size=ROOM_SETTINGS['name_text_size'], Weight=wx.BOLD)
                    #draw arrow
                    Canvas.AddArrowLine([(room.guiCoords[0]+roomWidth/2.0, room.guiCoords[1]),(left-border[0], (top+bottom)/2.0)])
            #draw room
            self.canvas.AddRectangle((room.guiCoords[0]-roomWidth/2.0,room.guiCoords[1]-settings[HEIGHT]/2.0), (roomWidth,settings[HEIGHT]), LineWidth = settings[LINEWIDTH], LineColor = settings[LINECOLOR], FillColor = settings[FILLCOLOR])

             #room name
            Canvas.AddScaledText(str(room.name),(room.guiCoords[0], room.guiCoords[1]-settings[HEIGHT]/2.0), Position="tc", Size=ROOM_SETTINGS['name_text_size'], Weight=wx.BOLD)
            #room allowed acuity, service
            roomString = ""
            try:
                roomString += "Acu={};".format(str(",".join([str(acu) for acu in room.allowedAcuity])))
            except AttributeError:
                pass
            try:
                roomString += "Ser={}".format(str(",".join([str(ser) for ser in room.allowedService])))
            except AttributeError:
                pass

            Canvas.AddScaledText(roomString, (room.guiCoords[0], room.guiCoords[1]+settings[HEIGHT]/2.0), Position="bc", Size=ROOM_SETTINGS['room_data_text_size'], Weight=wx.BOLD)
           
            #draw arrows to each of connecting rooms
            arrowStart = (room.guiCoords[0]+roomWidth/2.0, room.guiCoords[1])
            for childRoom in room.nextRooms:
                if isinstance(childRoom, BedRoom):
                    continue
                childRoomWidth = ROOM_SETTINGS['room_mapping'][type(childRoom)][WIDTH]
                if isinstance(childRoom, BedRoom):
                    childRoomWidth = ROOM_SETTINGS['room_mapping'][type(childRoom)][WIDTH]*(2-1/float(childRoom.capacity))
                arrowEnd = (childRoom.guiCoords[0]-childRoomWidth/2.0,childRoom.guiCoords[1])
                Canvas.AddArrowLine([arrowStart,arrowEnd])
        #C=Canvas.AddCircle((0,0),7,InForeground =True)
        #C.Bind(FloatCanvas.EVT_FC_LEFT_DOWN,self.ClickLeft)

        Canvas.ZoomToBB()

    #updates labels at each time step
    def update_labels(self):
        Canvas = self.canvas
        try:
            Canvas.RemoveObject(self.labelGroup)
        except ValueError:
            pass
        self.labelGroup = FloatCanvas.Group()
        for room in self.sim.hospital.layout.values():
            #draw output for allocation room
            if isinstance(room, BedAllocationRoom):
                settings = ROOM_SETTINGS['room_mapping'][type(room)]
                timeToIBA = room.get_time_to_initial_bed_assignment()
                numClosedBeds = room.get_num_closed_beds_by_category()[1]
                self.labelGroup.AddObject(FloatCanvas.ScaledText("t to IBA = {}({})".format(timeToIBA[0],timeToIBA[1]), (room.guiCoords[0], room.guiCoords[1]-settings[HEIGHT]/2.0 - ROOM_SETTINGS['name_text_size']), Position = "tc", Size=ROOM_SETTINGS['name_text_size']))
                self.labelGroup.AddObject(FloatCanvas.ScaledText("Num Closed Beds = {}".format(numClosedBeds), (room.guiCoords[0], room.guiCoords[1]-settings[HEIGHT]/2.0 - 2*ROOM_SETTINGS['name_text_size']), Position = "tc", Size=ROOM_SETTINGS['name_text_size']))
                
        Canvas.AddObject(self.labelGroup)
        Canvas.Draw()
    #draws all patients in the hospital
    def draw_patients(self):
        Canvas = self.canvas
        try:
            Canvas.RemoveObject(self.patientGroup)
        except ValueError:
            pass

        self.patientGroup = FloatCanvas.Group()

        #loop through all rooms and draw patients in each room
        for room in self.sim.hospital.layout.values():
            settings = PATIENT_SETTINGS
            roomSettings = ROOM_SETTINGS['room_mapping'][type(room)]
            roomCoords = room.guiCoords
            patientIDs = room.patients.keys()
            numPatients = len(patientIDs)
            roomWidth, roomHeight = self.get_room_dim(room)
            subRoomWidth = roomSettings[WIDTH]
                
            x_num = min(settings['num_per_row'], numPatients)
            y_num = (numPatients-1)/settings['num_per_row']+1
            if isinstance(room, BedRoom):
                subRoomWidth = roomWidth/float(room.capacity)
                x_num, y_num = 1,1

                #draw blocked beds
                for bed in range(len(room.beds)):
                    if room.bedsBlockedStaffing[bed]:
                        self.draw_blocked_bed(room, bed)

            avgBorderSpacing = 2
            x_spacing = (subRoomWidth-2*x_num*(settings['radius']+avgBorderSpacing))/float(numPatients+1)+2*(settings['radius']+avgBorderSpacing)
            y_spacing = (roomHeight-2*y_num*(settings['radius']+avgBorderSpacing))/float(numPatients+1)+2*(settings['radius']+avgBorderSpacing)

            if numPatients > ROOM_SETTINGS['max_display_patients']:
                from itertools import product
                patientTypes = tuple(product((MALE, FEMALE),(True, False), (True, False)))
                patientTypesDict = dict([(pType, 0) for pType in patientTypes])
                for patient in room.patients.values():
                    
                    pType = (patient.gender, patient.has_observed_MRSA(), patient.has_observed_VRE())
                    patientTypesDict[pType]+=1
                #draw generic male and female icon
                for n,pType in enumerate(patientTypes):
                    num = patientTypesDict[pType]
                    i,j = n%2, n/2
                    x_num = 2
                    y_num = 4
                    x_spacing = (roomWidth-2*x_num*(settings['radius']+avgBorderSpacing))/float(8+1)+2*(settings['radius']+avgBorderSpacing)+avgBorderSpacing*3
                    y_spacing = (roomHeight-2*y_num*(settings['radius']+avgBorderSpacing))/float(8+1)+2*(settings['radius']+avgBorderSpacing)
                    self.patientGroup.AddObject(self.draw_generic_patient(pType[0],num,(pType[1],pType[2]),(roomCoords[0]+(i-(x_num-1)/2.0)*x_spacing,roomCoords[1]+(j-(y_num-1)/2.0)*y_spacing)))

            else:
                for n in range(numPatients):
                    patient = room.patients[patientIDs[n]]
                    i, j = n%settings['num_per_row'],n/settings['num_per_row']
                    #draw patient in correct bed
                    if isinstance(room, BedRoom):
                        bedNum = room.beds.index(patient)
                        roomCoords = (room.guiCoords[0] - roomWidth/2.0+subRoomWidth*(bedNum)+subRoomWidth/2.0, room.guiCoords[1])
                        i, j = 0, 0
                    self.patientGroup.AddObject(self.draw_patient(patient,(roomCoords[0]+(i-(x_num-1)/2.0)*x_spacing,roomCoords[1]+(j-(y_num-1)/2.0)*y_spacing)))
        Canvas.AddObject(self.patientGroup)                    
        Canvas.Draw()
    #draws blocked bed
    def draw_blocked_bed(self, room, bed):
        settings = ROOM_SETTINGS['room_mapping'][BedRoom]
        roomWidth = settings[WIDTH]*(2-1/float(room.capacity))
        bedWidth = roomWidth/float(room.capacity)
        bedCoord = (room.guiCoords[0] - roomWidth/2.0+bedWidth*bed, room.guiCoords[1]-settings[HEIGHT]/2.0)
        self.patientGroup.AddObject(FloatCanvas.Line(((bedCoord[0],bedCoord[1]),
                                     (bedCoord[0]+bedWidth,bedCoord[1]+settings[HEIGHT])),
                                    LineColor = 'Grey',
                                    LineWidth = ROOM_SETTINGS['double_room_line_width']))
        self.patientGroup.AddObject(FloatCanvas.Line(((bedCoord[0],bedCoord[1]+settings[HEIGHT]),
                                     (bedCoord[0]+bedWidth,bedCoord[1])),
                                    LineColor = 'Grey',
                                    LineWidth = ROOM_SETTINGS['double_room_line_width']))
    #draws generic empty patient
    def draw_generic_patient(self, gender, number, observedCP, patientCoords):
        settings = PATIENT_SETTINGS
        group = FloatCanvas.Group()
        
        if gender == MALE:
            lineColor = settings['line_color_male']
        else:
            lineColor = settings['line_color_female']

        MRSAStr = observedCP[0] and "M+" or "M-"
        VREStr = observedCP[1] and "V+" or "V-"
        group.AddObject(FloatCanvas.Circle(patientCoords, settings['radius']*2, LineWidth=settings['line_width'], LineColor=lineColor, FillStyle="Solid"))
        group.AddObject(FloatCanvas.ScaledText("x {0}".format(number), patientCoords, Position="cc", Size = settings['id_text_size'], Weight = wx.BOLD))
        group.AddObject(FloatCanvas.ScaledText(str(MRSAStr), (patientCoords[0]-settings['radius'],patientCoords[1]+settings['radius']), Position="cc", Size = settings['acuity_text_size'], Weight = wx.BOLD))
        group.AddObject(FloatCanvas.ScaledText(str(VREStr), (patientCoords[0]+settings['radius'],patientCoords[1]+settings['radius']), Position="cc", Size = settings['acuity_text_size'], Weight = wx.BOLD))
            
        return group
    #draws patient at specified coord
    def draw_patient(self, patient, patientCoords):
        settings = PATIENT_SETTINGS
        group = FloatCanvas.Group()
        if isinstance(patient, CleaningPatient):
            fillColor = settings['fill_color_cleaning']
            lineWidth = settings['line_width']
            lineColor = settings['line_color_cleaning']
        else:
            if patient.gender == MALE:
                lineColor = settings['line_color_male']
            else:
                lineColor = settings['line_color_female']

            patID = patient.id
            fillColor = None
            lineWidth = patient.acuity
            if isinstance(patient, PhantomPatient):
                fillColor = settings['fill_color_phantom']
                patID = patient.realPatient.id
            
        group.AddObject(FloatCanvas.Circle(patientCoords, settings['radius']*2, LineWidth=lineWidth, LineColor=lineColor, FillColor=fillColor, FillStyle="Solid"))

        if isinstance(patient, BasePatient):
            #id num
            group.AddObject(FloatCanvas.ScaledText(str(patID), patientCoords, Position="cc", Size = settings['id_text_size'], Weight = wx.BOLD))
            #acuity num
            group.AddObject(FloatCanvas.ScaledText(str(patient.acuity), (patientCoords[0]-settings['radius'],patientCoords[1]-settings['radius']), Position="cc", Size = settings['acuity_text_size'], Weight = wx.BOLD))
            #service
            group.AddObject(FloatCanvas.ScaledText(str(patient.service), (patientCoords[0]+settings['radius'],patientCoords[1]-settings['radius']), Position="cc", Size = settings['acuity_text_size'], Weight = wx.BOLD))
            #MRSA VRE
            MRSAStr = patient.has_observed_MRSA() and "M+" or "M-"
            VREStr = patient.has_observed_VRE() and "V+" or "V-"
            group.AddObject(FloatCanvas.ScaledText(str(MRSAStr), (patientCoords[0]-settings['radius'],patientCoords[1]+settings['radius']), Position="cc", Size = settings['acuity_text_size'], Weight = wx.BOLD))
            group.AddObject(FloatCanvas.ScaledText(str(VREStr), (patientCoords[0]+settings['radius'],patientCoords[1]+settings['radius']), Position="cc", Size = settings['acuity_text_size'], Weight = wx.BOLD))
        else:
            group.AddObject(FloatCanvas.ScaledText("C", patientCoords, Position="cc", Size = settings['id_text_size'], Weight = wx.BOLD))

        return group
    #gets the true (width, height) of the room.
    def get_room_dim(self, room):
        roomSettings = ROOM_SETTINGS['room_mapping'][type(room)]
        width = roomSettings[WIDTH]
        height = roomSettings[HEIGHT]
        if isinstance(room, BedRoom):
            width = roomSettings[WIDTH]*(2-1/float(room.capacity))

        return (width, height)
    #automatically determines coordinates for each room in the hospital
    def auto_layout(self,layout, startRoom, endRoom):
        #horizontal and vertical spacing
        hSpacing = 55
        vSpacing = 70

        #set start room to (0,0)
        startRoom.guiCoords = (0,0)
        
        currNodes = deque(layout.itervalues())
        numVertChildDict = {}
        while currNodes:
            currRoom = currNodes.pop()
            if not currRoom.nextRooms:
                numVertChildDict[currRoom] = 1
                continue
            if isinstance(currRoom, BedAllocationRoom):
                numVertChildDict[currRoom] = len(currRoom.units)
                continue
            
            numVertNeeded = 0
            try:
                for childRoom in currRoom.nextRooms:
                    numVertNeeded += numVertChildDict[childRoom]
            except KeyError:
                currNodes.appendleft(currRoom)
            else:
                numVertChildDict[currRoom] = numVertNeeded

    
        #bredth first layout level by level
        currNodes = deque([startRoom])
        ParentNodes = {startRoom: None}
        maxHDistance = 0
        while currNodes:
            currRoom = currNodes.popleft()
            currIndex = 0
            #layout unit if this room is bed allocation room 
            if isinstance(currRoom, BedAllocationRoom):
                numUnits = len(currRoom.units.keys())

                for i in range(numUnits):
                    unitY = currRoom.guiCoords[1]+((numUnits-1)/2.0 - i)* vSpacing
                    currUnitX = currRoom.guiCoords[0]+hSpacing
                    prevUnitWidth = 0
                    roomSpacing =10
                    for bedRoom in currRoom.units.values()[i]:
                        currUnitX += prevUnitWidth/2.0+self.get_room_dim(bedRoom)[WIDTH]/2.0+roomSpacing
                        prevUnitWidth = self.get_room_dim(bedRoom)[WIDTH]
                        
                        bedRoom.guiCoords = (currUnitX, unitY)
                        maxHDistance = max(maxHDistance, bedRoom.guiCoords[0])
                        ParentNodes[bedRoom] = currRoom
                        currNodes.append(bedRoom)
                    
            else:
                n = len(currRoom.nextRooms)
                
                for i in range(n):
                    childRoom = currRoom.nextRooms[i]
                    if childRoom in ParentNodes:
                        continue
                    indexRange = range(currIndex, currIndex+numVertChildDict[childRoom])
                    indexMean = mean(indexRange)
                    currIndex += numVertChildDict[childRoom]
                    childRoom.guiCoords = (currRoom.guiCoords[0]+hSpacing, currRoom.guiCoords[1]+((numVertChildDict[currRoom]-1)/2.0-indexMean)*vSpacing)
                    maxHDistance = max(maxHDistance, childRoom.guiCoords[0])
                    ParentNodes[childRoom] = currRoom
                    currNodes.append(childRoom)

        endRoom.guiCoords = (maxHDistance+hSpacing, 0)
        
    def OnQuit(self,event):
        self.Close(True)
    def OnMatchCP(self,event):
        self.sim.changeHospitalCPMatching(self.matchCPMI.IsChecked())
    def OnLoadInputFile(self, event):
        dlg = wx.FileDialog(self, "Choose an input file", os.getcwd(),"","*.xlsx",wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.sim.load_input_file(path.join(dlg.GetDirectory(), dlg.GetFilename()))
    def OnRunBatch(self, event):
        batchOutputs = {}
        dlg = wx.DirDialog(self, "Choose a folder with runs", os.getcwd())
        if dlg.ShowModal() == wx.ID_OK:
            runPath = dlg.GetPath()
        for inputfile in glob(path.join(runPath, "*.xlsx")):
            #ignore temp files
            if path.split(inputfile)[1].startswith("~"):
                continue

            self.sim.reset_hospital()
            layout =self.sim.hospital.layout
            startRoomID, endRoomID = self.sim.hospital.validate_layout(layout)
            self.auto_layout(layout, layout[startRoomID], layout[endRoomID])
            self.sim.load_input_file(inputfile)
            self.OnStepMany(None, numSteps = self.sim.inputs.numStepsToRun)
            self.sim.save_output_file(path.splitext(inputfile)[0]+".pfout")
            batchOutputs[inputfile] = deepcopy(self.sim.outputs)
            self.sim.save_summary_file(path.join(path.split(inputfile)[0],"summaryFile.pfout"),batchOutputs)
    def OnWriteOutput(self, event):
        dlg = wx.FileDialog(self, "Save an output file", os.getcwd(),"","*.pfout",wx.OPEN)
        if dlg.ShowModal() == wx.ID_OK:
            self.sim.save_output_file(path.join(dlg.GetDirectory(), dlg.GetFilename()))

    def OnLoadPatients(self, event):
        if self.sim.INITIAL_PATIENTS_LOADED:
            dlg = wx.MessageDialog(self, "Patients Already Loaded", "Warning", wx.OK|wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
        else:
            self.sim.load_initial_patients()
            self.draw_patients()
    def OnStep(self,event, updateDrawing = True):
        if not self.sim.INPUTS_LOADED:
            dlg = wx.MessageDialog(self, "No Input Sheet Loaded", "Error", wx.OK|wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
            return

        isInit = self.sim.run_single_time_step()

        if isInit:
            if updateDrawing:
                self.draw_patients()
                self.update_labels()
        else:
            dlg = wx.MessageDialog(self, "Hospital not initialized", "Error", wx.OK|wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
            return
    def OnStepMany(self, event, numSteps = 10):
        if not self.sim.INPUTS_LOADED:
            dlg = wx.MessageDialog(self, "No Input Sheet Loaded", "Error", wx.OK|wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
            return
        for i in range(numSteps):
            self.OnStep(None, False)
        self.draw_patients()
        self.update_labels()
    def OnStepN(self, event):
        try:
            n = int(self.stepNTc.GetValue())
        except ValueError:
            dlg = wx.MessageDialog(self, "Invalid Value of N", "Error", wx.OK|wx.ICON_ERROR)
            dlg.ShowModal()
            dlg.Destroy()
        else:
            self.OnStepMany(None,n)
            """
            if not self.runThread:
                self.runThread = RunThread(self, n)
            self.runThread.start()
            """
    def OnSavePNG(self, event=None):
        import os
        dlg = wx.FileDialog(
            self, message="Save file as ...", defaultDir=os.getcwd(), 
            defaultFile="", wildcard="*.png", style=wx.SAVE
            )
        if dlg.ShowModal() == wx.ID_OK:
            path = dlg.GetPath()
            if not(path[-4:].lower() == ".png"):
                path = path+".png"
            self.canvas.SaveAsImage(path)            
            

if __name__ == '__main__':
    sim = Sim()

    sim.create_hospital()

    layout = [("s","Start",("Start",),('wa',)),
              ("wa","WardAllocation",("Acu/Ser Alloc",),('sicu','micu','mstepdown','sgen','mgen','smobs')),
              ("e","End",("End",),())]
    bedAllocations = (("sicu", "BedAllocation", ((5,),(SURGERY,),"Sur ICU"), ('sicu1', 'sicu2', 'sicu3', 'sicu4')),
                      ("micu", "BedAllocation", ((5,),(MEDICINE,),"Med ICU"), ('micu1', 'micu2')),
                      ("mstepdown", "BedAllocation", ((4,),(MEDICINE,),"Step Down"), ('mstepdown1', 'mstepdown2')),
                      ("sgen", "BedAllocation", ((3,),(SURGERY,),"Sur General"), ('sgen1', 'sgen2')),
                      ("mgen", "BedAllocation", ((3,),(MEDICINE,),"Med General"), ('mgen1', 'mgen2')),
                      ("smobs", "BedAllocation", ((2,),(SURGERY,MEDICINE),"Observation"), ('smobs1', 'smobs2')),
                      )

    beds = [(bedA[0]+"1", "Bed", (1,bedA[2][2]+"1","unit 1"),()) for bedA in bedAllocations]
    beds.extend([(bedA[0]+"2", "Bed", (2,bedA[2][2]+"2","unit 1"),()) for bedA in bedAllocations])
    beds.extend([("sicu3","Bed",(2,"Surgery ICU3", "unit 2"),())])
    beds.extend([("sicu4","Bed",(2,"Surgery ICU4", "unit 2"),())])
    layout.extend(bedAllocations)
    layout.extend(beds)
    sim.import_layout('l', layout)

    app=wx.App()

    frame = MainWindow(None, "Patient Flow", sim)
    frame.Show()

    app.MainLoop()
