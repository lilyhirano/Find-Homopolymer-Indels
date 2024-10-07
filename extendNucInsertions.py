import os
import pandas as pd

workingDirectory = os.path.dirname(__file__)
bedFile = "rad30_sorted_insertions_Trinuc.bed"                  # rad16_sorted_insertions_Trinuc       WT_sorted_insertions_Trinuc
allChrTxt = ["chr1.txt", "chr2.txt", "chr3.txt", "chr4.txt", "chr5.txt", "chr6.txt", "chr7.txt", "chr8.txt",
             "chr9.txt", "chr10.txt", "chr11.txt", "chr12.txt", "chr13.txt", "chr14.txt", "chr15.txt",
             "chr16.txt"]
roman = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]
output = "projectFiles/rad30_sorted_insertions_54nuc.bed"

desiredNucLen = 54  #must be EVEN

dataBedFile = {
        "Chromosome" : [],
        "Start(0)" : [],
        "Start(1)" : [],
        "Nucleotides" : [],
        "Insertion" : [],
        "Strand" : []

}

revComplementDict  = {
        'A' : 'T',
        'G' : 'C',
        'T' : 'A',
        'C' : 'G',
        '*' : '*'
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
       
def fileIntoDf(myDfBed):
        with open(os.path.join(workingDirectory, bedFile), "r") as file:
                myDfBed = pd.DataFrame(dataBedFile) #make dictionary into dataframe

                for line in file: #iterate through each line

                        if line.strip(): #if line not empty:
                                
                                myDfBed = parseLine(line, myDfBed)

        return myDfBed


def checkIfValid5(start, desiredNucLen):
        if start < desiredNucLen / 2:
                return False
        else:
                return True
def checkIfValid3(start, end, desiredNucLen):
        if desiredNucLen / 2 > end - start:
                return False
        else:
                return True

def grabSomeNuc5(start, desiredNucLen):
        grabbedSequence = ""
        x = start
        while x > 0:
                grabbedSequence += dna[start - x]
                x-= 1
        numStars = int(desiredNucLen/2) - start
        
        stars = ""
        for x in range(numStars):
                stars += "*"

        grabbedSequence = stars + grabbedSequence
        
        return grabbedSequence

def grabSomeNuc3(endOfDna, start, desiredNucLen):
        grabbedSequence = ""
        x = 0
        while x < endOfDna - start  + 1:
                grabbedSequence += dna[start + x]
                x+= 1

        numStars =  -1 *(endOfDna- start - int(desiredNucLen/2) + 1)
        
        stars = ""
        for x in range(numStars):
                stars += "*"
        grabbedSequence = grabbedSequence + stars

        return grabbedSequence

def revComplement(dna):
        dna2 = ""
        for base in dna: 
                dna2 += revComplementDict[base]
        return dna2[::-1]

def grabMoreNuc3(start, dna, desiredNucLen):
        grabbedSequence = ""
        x = 0
        while x < int(desiredNucLen/2) : #(num nucs wanted/ 2)
                grabbedSequence += dna[start + x]
                x+= 1
        return grabbedSequence

def grabMoreNuc5(start, dna, desiredNucLen):
        grabbedSequence = ""
        x = int(desiredNucLen/2) #num nucs wanted /2
        while x > 0:
                grabbedSequence += dna[start - x]
                x-= 1
        return grabbedSequence


insertionDF = fileIntoDf(myDfBed).copy() #copy dataframe into variable

x = 0
while x < 16:
        with open(os.path.join(workingDirectory, allChrTxt[x]), "r") as chrFile:
                y = 0

                dna = chrFile.read().replace('\n', '')   #turn txt to string
                end = len(dna) - 1

                filteredDf = insertionDF[insertionDF["Chromosome"] == ("chr%s"  %roman[x])]

                if filteredDf.empty:
                        x+= 1
                        continue
                else:
                        rangeIndex = len(filteredDf.index)
                        x += 1

                        for y in filteredDf.index:
                                
                                start = int(filteredDf.at[y, "Start(0)"])  #row by column
                                
                                if(insertionDF.at[y,"Strand"] == ("-")):
                                      
                                        if checkIfValid3(start, end, desiredNucLen) and checkIfValid5(start, desiredNucLen):
                                                sequence5 = grabMoreNuc5(start, dna, desiredNucLen)
                                                sequence3 = grabMoreNuc3(start, dna, desiredNucLen)
                                        elif checkIfValid3(start,end,desiredNucLen) and not checkIfValid5(start, desiredNucLen):
                                                sequence5 = grabSomeNuc5(start, desiredNucLen)
                                                sequence3 = grabMoreNuc3(start, dna, desiredNucLen)
                                        elif checkIfValid5(start, desiredNucLen) and not checkIfValid3(start,end,desiredNucLen):
                                                sequence5 = grabMoreNuc5(start, dna, desiredNucLen)
                                                sequence3 = grabSomeNuc3(end, start, desiredNucLen)

                                        sequence = sequence5 +sequence3
                                        sequence = revComplement(sequence)
                                        insertionDF.at[y, "Nucleotides"] = sequence    
                                else:
                                        
                                        if checkIfValid3(start, end, desiredNucLen) and checkIfValid5(start, desiredNucLen):
                                                sequence5 = grabMoreNuc5(start, dna, desiredNucLen)
                                                sequence3 = grabMoreNuc3(start, dna, desiredNucLen)
                                        elif checkIfValid3(start,end,desiredNucLen) and not checkIfValid5(start, desiredNucLen):
                                                sequence5 = grabSomeNuc5(start, desiredNucLen)
                                                sequence3 = grabMoreNuc3(start, dna, desiredNucLen)
                                        elif checkIfValid5(start, desiredNucLen) and not checkIfValid3(start,end,desiredNucLen):
                                                sequence5 = grabMoreNuc5(start, dna, desiredNucLen)
                                                sequence3 = grabSomeNuc3(end, start, desiredNucLen)

                                        sequence = sequence5 +sequence3
                                        insertionDF.at[y, "Nucleotides"] = sequence
                                                
insertionDF.to_csv(output, sep = '\t', header = False, index = False)