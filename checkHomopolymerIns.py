import os
import pandas as pd

workingDirectory = os.path.dirname(__file__)
allBedFiles = ["rad30_sorted_insertions_54nuc.bed"]             #WT_sorted_insertions_12nuc.bed
output = "projectFiles/54Nuc_and_HP_Ins/rad30_insertions_HP_54Nuc.bed"

dataBedFile = {
        "Chromosome" : [],
        "Start(0)" : [],
        "Start(1)" : [],
        "Nucleotides" : [],
        "Insertion" : [], 
        "Strand" : []

}
newDf = {
       "HP Count" : [],
        "Type" : []
}

myDfBed = pd.DataFrame(dataBedFile) #make dictionary into dataframe


def parseLine(bedLine, dataframe):
       
        cleanSplitLine = bedLine.split()

        newRow = {"Chromosome" : cleanSplitLine[0], "Start(0)" : cleanSplitLine[1], "Start(1)" : 
        cleanSplitLine[2], "Nucleotides" : cleanSplitLine[3], 
        "Insertion" : cleanSplitLine[4], "Strand" : cleanSplitLine[5]}
     
        newRowDf = pd.DataFrame([newRow], columns = dataframe.columns) #make the new row in dataframe format
        dataframe = pd.concat([dataframe, newRowDf], ignore_index = True) #combine the two dataframes
     
        return dataframe

def fileIntoDf(bedFileNum, myDfBed):
        with open(os.path.join(workingDirectory, allBedFiles[bedFileNum]), "r") as file:
                myDfBed = pd.DataFrame(dataBedFile) #make dictionary into dataframe

                for line in file: #iterate through each line

                        if line.strip(): #if line not empty:
                                
                                myDfBed = parseLine(line, myDfBed)

        return myDfBed

def reverse(dna):
        return dna[::-1]

def findHomopolyRight(itemNum, start):
        maximum = 0
        homoType = "NA"
        
        sequence = parsedDF.loc[itemNum, "Nucleotides"]  #Grab sequence
        # if parsedDF.loc[itemNum, "Strand"] == "-":
        #         sequence = reverse(sequence)     
        base = start

        while base < start + 1:

                HPcount = 1 #start with one because minimum is one

                while base + 1 < len(sequence) and sequence[base] == sequence[base + 1]: #iterate while not at the end / is the same nucleotide
                        HPcount += 1
                        base += 1 

                if HPcount > maximum:   #check maxiumm
                        maximum = HPcount
                        homoType = sequence[start]

                base += 1


        return (maximum, homoType)

def findHomopolyLeft(itemNum, start):
        maximum = 0
        homoType = "NA"
        
        sequence = parsedDF.loc[itemNum, "Nucleotides"]  #Grab sequence
        # if parsedDF.loc[itemNum, "Strand"] == "-":
        #         sequence = reverse(sequence)     
        base = start - 1 

        while base > start - 2:

                HPcount = 1 #start with one because minimum is one

                while base - 1 > 0 and sequence[base] == sequence[base - 1]: #iterate while not at the end / is the same nucleotide
                        HPcount += 1
                        base -= 1 

                if HPcount > maximum:   #check maxiumm
                        maximum = HPcount
                        homoType = sequence[start - 1]

                base -= 1


        return (maximum, homoType)




parsedDF = fileIntoDf(0, myDfBed).copy() 

for items in range(len(parsedDF)):
      start =  int(54/2)
      dataRight = findHomopolyRight(items, start)
      dataLeft = findHomopolyLeft(items, start)


      if dataRight[1] == dataLeft[1]:
        hpcountTotal = dataRight[0] + dataLeft[0]
        if hpcountTotal > 3:
               newDf["Type"].append(dataRight[1])
        else:
               newDf["Type"].append("NA")

      elif dataLeft[0] > dataRight[0]:
        
        hpcountTotal = dataLeft[0]
        if hpcountTotal > 3:
                newDf["Type"].append(dataLeft[1])        
        else:
               newDf["Type"].append("NA")

      else:
        hpcountTotal = dataRight[0]    
        if hpcountTotal > 3:
               newDf["Type"].append(dataRight[1])
        else:
               newDf["Type"].append("NA")



      newDf["HP Count"].append(hpcountTotal)    #append to add onto dctionary   


homoPolyDF = pd.DataFrame(newDf)

finalDF = pd.DataFrame.join(parsedDF, homoPolyDF)

finalDF.to_csv(output, sep = '\t', header = False, index = False)