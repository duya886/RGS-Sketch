#USing this code, we turn the downloaded csv files into binary dat files 
import csv
import os
import time

def divideDataSet(lines,outputDir,pieceNum):
    totalLineNum=len(lines)
    pieceLineNum=int(totalLineNum/pieceNum)

    for i in range(pieceNum):
        outputFilePath=os.path.join(outputDir,"%02d.txt"%i)

        if i<pieceNum-1:
            outputFile=open(outputFilePath,"w")
            outputFile.writelines(lines[pieceLineNum*i:pieceLineNum*(i+1)])
            outputFile.close()
        else:
            outputFile=open(outputFilePath,"w")
            outputFile.writelines(lines[pieceLineNum*i:])
            outputFile.close()

def processCosOrMulCsvFile(filePaths,outputDir,prodcutid_as_flowID=True):
    lines=[]
    for filePath in filePaths:
        f=open(filePath)
        reader=csv.reader(f)

        lineNum = -1
        for row in reader:
            lineNum+=1
            if lineNum==0:
                continue
            else:
                userid=row[7]
                prodcutid=row[2]

                flowId =prodcutid
                elementId=userid

                if not prodcutid_as_flowID:
                    flowId=userid
                    elementId=prodcutid

                lines.append(str(flowId)+" "+str(elementId)+"\n");

            if lineNum % 1000000 == 0:
                print("have read %d lines in %s" % (lineNum, filePath))

        f.close()
        print("finish reading %s" % (filePath))

    divideDataSet(lines,outputDir,10)


def getDataNum(data):
    newData = data.split(" ")  # newData should have two elements, the first is flow id, and the second is the element, both are strings

    SrcDstNum=(int(newData[1])<<32)+int(newData[0])
    return SrcDstNum

def changeDataSetToBinary(dirName,fileName,outPutDirName):
    startTime=time.time()

    numOfRecords=0
    txtFileName = os.path.join(dirName, fileName)
    txtFile = open(txtFileName,'r')

    outPutDataFileName = os.path.join(outPutDirName,fileName.replace(".txt",".dat"))
    outPutDataFile = open(outPutDataFileName,"wb")

    data = txtFile.readline()

    while data != "":
        data = data.strip()
        if data == "":
            continue
        numOfRecords+=1

        dataNum = getDataNum(data)
        key_bin=dataNum.to_bytes(8,byteorder='little',signed=False)
        outPutDataFile.write(key_bin)

        data = txtFile.readline()
        if numOfRecords % 1000000 == 0:
            currentTime=time.time()
            print("processing file %s, have processed a total of %d lines, cost %.2f seconds"%(txtFileName,numOfRecords,currentTime-startTime))

    txtFile.close()
    outPutDataFile.close()

if __name__=="__main__":


    # rawDirName=r".\Cos2\raw"
    # fileNames=["2019-Oct.csv","2019-Nov.csv","2019-Dec.csv","2020-Jan.csv","2020-Feb.csv"]
    #
    # filePaths=[]
    #
    # for fileName in fileNames:
    #     filePath=os.path.join(rawDirName,fileName)
    #     filePaths.append(filePath)
    #
    # outputDir = r".\Cos2"
    # processCosOrMulCsvFile(filePaths,outputDir)
    # for i in range(10):
    #     fileName="%02d.txt"%i
    #     changeDataSetToBinary(outputDir,fileName,outputDir)


    rawDirName=r".\Mul\raw"
    fileNames=["2020-Apr.csv"]
    filePaths=[]
    for fileName in fileNames:
        filePath=os.path.join(rawDirName,fileName)
        filePaths.append(filePath)

    outputDir = r".\Mul"
    processCosOrMulCsvFile(filePaths,outputDir)
    for i in range(10):
        fileName="%02d.txt"%i
        changeDataSetToBinary(outputDir,fileName,outputDir)