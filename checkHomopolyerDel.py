import os
import pandas as pd

workingDirectory = os.path.dirname(__file__)
allBedFiles = ["WT_sorted_deletions_25nuc.bed"]             #WT_sorted_deletions_12nuc.bed
output = "projectFiles/25Nuc_and_HP_Del/WT_Deletions_HP_25Nuc.bed"

dataBedFile = {
        "Chromosome" : [],
        "Start(0)" : [],
        "Start(1)" : [],
        "Nucleotides" : [],
        "Deletion" : [], 
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
        "Deletion" : cleanSplitLine[4], "Strand" : cleanSplitLine[5]}
     
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

def findHomopoly(itemNum):
        maximum = 0
        homoType = "NA"
        
        sequence = parsedDF.loc[itemNum, "Nucleotides"]  #Grab sequence
        if parsedDF.loc[itemNum, "Strand"] == "-":
                sequence = reverse(sequence)     
        base = 0

        while base < 2:

                HPcount = 1 #start with one because minimum is one

                while base + 1 < len(sequence) and sequence[base] == sequence[base + 1]: #iterate while not at the end / is the same nucleotide
                        HPcount += 1
                        base += 1 

                if HPcount > maximum:   #check maxiumm
                        maximum = HPcount
                        homoType = sequence[1]

                base += 1

        if maximum > 3:
                return (maximum, homoType)
        else:
                return (maximum, "NA")





parsedDF = fileIntoDf(0, myDfBed).copy() 

for items in range(len(parsedDF)):
      data = findHomopoly(items)
      newDf["HP Count"].append(data[0])    #append to add onto dctionary
      newDf["Type"].append(data[1])

homoPolyDF = pd.DataFrame(newDf)

finalDF = pd.DataFrame.join(parsedDF, homoPolyDF)

finalDF.to_csv(output, sep = '\t', header = False, index = False)


