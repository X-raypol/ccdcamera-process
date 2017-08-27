#add Information to the header of a fits file

import numpy as np
from astropy.io import fits
import sys
import os.path


def addHeaderToFile():
    #get the path to the current directory
    path = os.getcwd()            
    #inicialize most of what will be added to a fits file
    addHdr = [["Date","",""],
           ["Time","",""],
           ["Frames","",""],
           ["Exposure","","seconds"],
           ["Temp","","degrees C (of detector)"],
           ["ROI","","pixels (Region of Interest)"],
           ["ADC","","kilohertz (Analog-Digital Conversion)"],
           ["ExpTime","","seconds (Experiment Time, as in run time)"],
           ["Voltage","","volts"],
           ["Current","","milliamps (Current Average)"],
           ["CurRange","","milliamps (Current Range)"],
           ["PolAngle","","degrees (Polarization Angle)"],
           ["Anode","",""],
           ["MirPos","","millimeters (Mirror Position)"],
           ["MirRot","","degrees (Mirror Rotation)"],
           ["ApatTran","","centimeters (Apature Translation)"],
           ["GratTran","","centimeters (Grating Translation)"],
           ["GratRot","","degrees (Grating Rotation)"],
           ["CamTran","","centimeters (Camera Translation)"],
           ["CamVert","","centimeters (Camera Vertical)"],
           ["DetMirRo","","degrees (Detector Mirror Rotation)"],
           ["GratVert","","centimeters (Grating Verticle)"],
           ["Diode1Av","","(Diode 1 Average)"],
           ["Diode2Av","","(Diode 2 Average)"]]

    #get what the user wants
    option = input("What are you adding? \n1:adding both header and stats\n2:adding stats only\n3:adding header only\n4:adding indevidual tag/value pairs\nOption--> ")

    #adding both header data and stats
    if int(option) == 1:
        txt = input("What is the name of the text file with additional header information made by the labview code --> ")
        stats = input("What is the name of the stats file from that day (make sure it is in the same directory as this file)--> ")
        if os.path.isfile(path+"\\"+txt) and os.path.isfile(path+"\\"+stats):
            #open the file with header information in memory
            infile = open(txt,"r")
            #read in the contents of the file
            hdrData = infile.read()
            #close file in memory
            infile.close()
            #split the string on the comma into a list
            hdrDataArray = hdrData.split(", ")

            #loop through this new list
            for element in hdrDataArray:
                #for each element, headerEntry is a list of the form ["tag",data]
                headerEntry = element.split(":")
                #add the data to the correct place in addhdr
                if "Date" in headerEntry:
                    addHdr[0][1] = headerEntry[1][1:]
                elif "Time" in headerEntry:
                    hours = 0
                    if "PM" in headerEntry[3] and int(headerEntry[1]) != 12:
                        hours = 12
                    addHdr[1][1] = str(int(headerEntry[1])+hours)+":"+headerEntry[2]+":"+headerEntry[3][:2]
                    #startTime is necesary for reading in the stats
                    if int(headerEntry[3][:2]) == 57 or int(headerEntry[3][:2]) == 58 or int(headerEntry[3][:2]) == 59:
                        startTime1 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":"+str(int(headerEntry[3][:2])-2)
                        startTime2 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":"+str(int(headerEntry[3][:2])-1)
                        startTime3 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":"+str(int(headerEntry[3][:2]))
                    elif int(headerEntry[3][:2]) < 10 and int(headerEntry[3][:2]) > 7:
                        startTime1 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":0"+str(int(headerEntry[3][:2])-2)
                        startTime2 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":0"+str(int(headerEntry[3][:2])-1)
                        startTime3 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":0"+str(int(headerEntry[3][:2]))
                    elif int(headerEntry[3][:2]) <= 7:
                        startTime1 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":0"+str(int(headerEntry[3][:2])+2)
                        startTime2 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":0"+str(int(headerEntry[3][:2])+1)
                        startTime3 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":0"+str(int(headerEntry[3][:2]))
                    else:
                        startTime1 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":"+str(int(headerEntry[3][:2]))
                        startTime2 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":"+str(int(headerEntry[3][:2])+1)
                        startTime3 = str(int(headerEntry[1])+hours)+":"+str(headerEntry[2])+":"+str(int(headerEntry[3][:2])+2)
                elif "Frames" in headerEntry:
                    addHdr[2][1] = int(headerEntry[1])
                elif "Exposure" in headerEntry:
                    addHdr[3][1] = float(headerEntry[1])
                elif "Temp" in headerEntry:
                    addHdr[4][1] = int(headerEntry[1])
                elif "ROI" in headerEntry:
                    a = headerEntry[1].split()
                    addHdr[5][1] = "x1:"+a[0]+", x2:"+a[1]+", y1:"+a[2]+", y2:"+a[3]
                elif "ADC" in headerEntry:
                    addHdr[6][1] = int(headerEntry[1]) 
                elif "ExpTime" in headerEntry:
                    addHdr[7][1] = float(headerEntry[1])
                    #stats file is written to every 3 seconds, so the number of lines that have to do with the
                    #experiment is the length of the experiment divided by 3
                    numStatsLines = round(float(headerEntry[1]))//3


            #open stats file in memory
            infile = open(stats,"r")
            #read it into a list of the different lines of the file
            statsdata = infile.readlines()
            #close the file in memory
            infile.close()

            startIndex = 0
            #loop through the lines of the stats data looking for the start time
            for i in range(len(statsdata)):
                if (startTime1 in statsdata[i] or startTime2 in statsdata[i] or startTime3 in statsdata[i]) and startIndex == 0:
                    startIndex = i
            #endIndex is the start time + the length of the experiment
            endIndex = startIndex + numStatsLines

            #if the start index doesn't change, something is probably wrong
            if startIndex == 0:
                print("Something went wrong when reading the stats file, one of the times "+startTime1+" / "+startTime2+" / "+startTime3+" / "
                      "was not found in the stats file,\nplease check that you have labled the stats file "
                      "correctly, or simply open today's stats file to check to see if the start time is present")

            #if the end index is larger than the last line in the file, something is probably wrong
            elif endIndex > len(statsdata):
                print("Something went wrong when reading the stats file, apparently the important stats information "
                      "starts on line",startIndex,"and ends on line",endIndex,"yet there are only",len(statsdata),
                      "lines in the file.\n please check that you have labled the stats file "
                      "correctly, or simply open today's stats file to see if there is an obvious problem")
                
            else:
                experimentData = []
                #cut out the block of data where the experiment took place and split up each row
                for j in statsdata[startIndex:endIndex]:
                    experimentData.append(j.split())
                #make final numpy array
                statsArray = np.array(experimentData)
                #get non changing values from the center of the data block
                statsArrayMiddle = len(statsArray)//2

                #add the different values from the stats file to the header array
                addHdr[8][1] = float(statsArray[statsArrayMiddle][2]) #voltage shouldn't change so get a value from the middle
                addHdr[9][1] = statsArray[:,3].astype(float).mean() #Current average
                addHdr[10][1] = statsArray[:,3].astype(float).max()-statsArray[:,3].astype(float).min() #current range
                #same as voltage these values should not have changed over the course of the experiment
                for i in range(11,22):
                    addHdr[i][1] = float(statsArray[statsArrayMiddle][i-7])
                addHdr[22][1] = statsArray[:,15].astype(float).mean() #Diode 1 Average
                addHdr[23][1] = statsArray[:,16].astype(float).mean() #Diode 2 Average

            return addHdr 

        elif not os.path.isfile(path+"\\"+txt):
            print("The file "+txt+" was not found in the directory "+path)
            sys.exit()

        elif not os.path.isfile(path+"\\"+stats):
            print("The file "+stats+" was not found in the directory "+path)
            sys.exit()

        else:
            print("something went wrong reading the files")
            sys.exit()

    


    elif int(option) == 2:
        stats = input("What is the name of the stats file from that day (make sure it is in the same directory as this file)--> ")
        if os.path.isfile(path+"\\"+stats):
            starttime = input("What time would you like to start reading the stats file from? (hh:mm:ss) -> ")
            exptime = input("how long was the experiment (in seconds please) -> ")
            numStatsLines = int(exptime)//3
            timelist = starttime.split(":")
            if int(timelist[2]) == 57 or int(timelist[2]) == 58 or int(timelist[2]) == 59:
                startTime1 = str(int(timelist[0]))+":"+str(timelist[1])+":"+str(int(timelist[2])-2)
                startTime2 = str(int(timelist[0]))+":"+str(timelist[1])+":"+str(int(timelist[2])-1)
                startTime3 = str(int(timelist[0]))+":"+str(timelist[1])+":"+str(int(timelist[2]))
            elif int(timelist[2]) < 10 and int(timelist[2]) > 7:
                startTime1 = str(int(timelist[0]))+":"+str(timelist[1])+":0"+str(int(timelist[2])-2)
                startTime2 = str(int(timelist[0]))+":"+str(timelist[1])+":0"+str(int(timelist[2])-1)
                startTime3 = str(int(timelist[0]))+":"+str(timelist[1])+":0"+str(int(timelist[2]))
            elif int(timelist[2]) <= 7:
                startTime1 = str(int(timelist[0]))+":"+str(timelist[1])+":0"+str(int(timelist[2])+2)
                startTime2 = str(int(timelist[0]))+":"+str(timelist[1])+":0"+str(int(timelist[2])+1)
                startTime3 = str(int(timelist[0]))+":"+str(timelist[1])+":0"+str(int(timelist[2]))
            else:
                startTime1 = str(int(timelist[0]))+":"+str(timelist[1])+":"+str(int(timelist[2])+2)
                startTime2 = str(int(timelist[0]))+":"+str(timelist[1])+":"+str(int(timelist[2])+1)
                startTime3 = str(int(timelist[0]))+":"+str(timelist[1])+":"+str(int(timelist[2]))


            #open stats file in memory
            infile = open(stats,"r")
            #read it into a list of the different lines of the file
            statsdata = infile.readlines()
            #close the file in memory
            infile.close()

            startIndex = 0
            #loop through the lines of the stats data looking for the start time
            for i in range(len(statsdata)):
                if (startTime1 in statsdata[i] or startTime2 in statsdata[i] or startTime3 in statsdata[i]) and startIndex == 0:
                    startIndex = i
            #endIndex is the start time + the length of the experiment
            endIndex = startIndex + numStatsLines

            #if the start index doesn't change, something is probably wrong
            if startIndex == 0:
                print("Something went wrong when reading the stats file, one of the times "+startTime1+" / "+startTime2+" / "+startTime3+" / "
                      "was not found in the stats file,\nplease check that you have labled the stats file "
                      "correctly, or simply open today's stats file to check to see if the start time is present")

            #if the end index is larger than the last line in the file, something is probably wrong
            elif endIndex > len(statsdata):
                print("Something went wrong when reading the stats file, apparently the important stats information "
                      "starts on line",startIndex,"and ends on line",endIndex,"yet there are only",len(statsdata),
                      "lines in the file.\n please check that you have labled the stats file "
                      "correctly, or simply open today's stats file to see if there is an obvious problem")
                
            else:
                experimentData = []
                #cut out the block of data where the experiment took place and split up each row
                for j in statsdata[startIndex:endIndex]:
                    experimentData.append(j.split())
                #make final numpy array
                statsArray = np.array(experimentData)
                #get non changing values from the center of the data block
                statsArrayMiddle = len(statsArray)//2

                #add the different values from the stats file to the header array
                addHdr[8][1] = float(statsArray[statsArrayMiddle][2]) #voltage shouldn't change so get a value from the middle
                addHdr[9][1] = statsArray[:,3].astype(float).mean() #Current average
                addHdr[10][1] = statsArray[:,3].astype(float).max()-statsArray[:,3].astype(float).min() #current range
                #same as voltage these values should not have changed over the course of the experiment
                for i in range(11,22):
                    addHdr[i][1] = float(statsArray[statsArrayMiddle][i-7])
                addHdr[22][1] = statsArray[:,15].astype(float).mean() #Diode 1 Average
                addHdr[23][1] = statsArray[:,16].astype(float).mean() #Diode 2 Average

        else:
            print("The file "+stats+" was not found in the directory "+path)
            sys.exit()

        return addHdr[7:]

        
    elif int(option) == 3:
        txt = input("What is the name of the text file with additional header information made by the labview code --> ")
        if os.path.isfile(path+"\\"+txt):
            #open the file with header information in memory
            infile = open(txt,"r")
            #read in the contents of the file
            hdrData = infile.read()
            #close file in memory
            infile.close()
            #split the string on the comma into a list
            hdrDataArray = hdrData.split(", ")

            #loop through this new list
            for element in hdrDataArray:
                #for each element, headerEntry is a list of the form ["tag",data]
                headerEntry = element.split(":")
                #add the data to the correct place in addhdr
                if "Date" in headerEntry:
                    addHdr[0][1] = headerEntry[1][1:]
                elif "Time" in headerEntry:
                    hours = 0
                    if "PM" in headerEntry[3] and int(headerEntry[1]) != 12:
                        hours = 12
                    addHdr[1][1] = str(int(headerEntry[1])+hours)+":"+headerEntry[2]+":"+headerEntry[3][:2]
                elif "Frames" in headerEntry:
                    addHdr[2][1] = int(headerEntry[1])
                elif "Exposure" in headerEntry:
                    addHdr[3][1] = float(headerEntry[1])
                elif "Temp" in headerEntry:
                    addHdr[4][1] = int(headerEntry[1])
                elif "ROI" in headerEntry:
                    a = headerEntry[1].split()
                    addHdr[5][1] = "x1:"+a[0]+", x2:"+a[1]+", y1:"+a[2]+", y2:"+a[3]
                elif "ADC" in headerEntry:
                    addHdr[6][1] = int(headerEntry[1]) 
                elif "ExpTime" in headerEntry:
                    addHdr[7][1] = float(headerEntry[1])
        else:
            print("The file "+txt+" was not found in the directory "+path)
            sys.exit()

        return addHdr[:7]

    elif int(option) == 4:
        return 1
        


def main():
    file = input("Name of the file you are adding header information to?\n(Note: please include extentions for all files) -> ")
    #open the fits file in memory
    hdulist = fits.open(file, mode='update')
    #get the header
    header = hdulist[0].header
    #get what needs to be updated
    newInfo = addHeaderToFile()

    #if they chose to indevidually hard code header entries
    if newInfo == 1:
        print("For each entry type the header tag (must be 8 or fewer characters), then the value, and a comment (optiona), seperated by spaces")
        option = "y"
        while option == "y":
            entry = input("-> ")
            entryList = entry.split()
            if len(entryList) == 2:
                header[entryList[0]] = entryList[1]
            else:
                header[entryList[0]] = (entryList[1],entry[entry.find(entryList[2]):])
            option = input("Enter another? -> ")
            
    #otherwise add whatever was in the documents
    else:
        for key_value in newInfo:
                header[key_value[0]] = (key_value[1],key_value[2])
    #update the fits header
    hdulist.flush()
    hdulist.close()
    print("Header updated")

if __name__ == "__main__":
    main()





            
