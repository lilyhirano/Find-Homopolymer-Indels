import os
import pandas as pd

workingDirectory = os.path.dirname(__file__)
bedFile = "WT_sorted_deletions_Trinuc.bed"                  #rad30_sorted_deletions_Trinuc.bed        WT_sorted_deletions_Trinuc.bed
allChrTxt = ["chr1.txt", "chr2.txt", "chr3.txt", "chr4.txt", "chr5.txt", "chr6.txt", "chr7.txt", "chr8.txt",
             "chr9.txt", "chr10.txt", "chr11.txt", "chr12.txt", "chr13.txt", "chr14.txt", "chr15.txt",
             "chr16.txt"]
roman = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]
output = "projectFiles/WT_sorted_deletions_25nuc.bed"
x = 0


dataBedFile = {
        "Chromosome" : [],
        "Start(0)" : [],
        "Start(1)" : [],
        "Nucleotides" : [],
        "Deletion" : [],
        "Strand" : []

}

revComplementDict  = {
        'A' : 'T',
        'G' : 'C',
        'T' : 'A',
        'C' : 'G'
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
       
def fileIntoDf(x, myDfBed):
        with open(os.path.join(workingDirectory, bedFile), "r") as file:
                myDfBed = pd.DataFrame(dataBedFile) #make dictionary into dataframe

                for line in file: #iterate through each line

                        if line.strip(): #if line not empty:
                                
                                myDfBed = parseLine(line, myDfBed)

        return myDfBed


def checkIfValid(start, end):
        if start > end - 23:  #nucs wanted -2
                return False
        else:
                return True

def grabSomeNuc(endOfDna, start):
        grabbedSequence = ""
        numNucLeft = endOfDna - start
        x = -1
        while x < numNucLeft: #num nucs wanted - 2
                grabbedSequence += dna[start + x]
                x+= 1
        return grabbedSequence

def revComplement(dna):
        dna2 = ""
        for base in dna: 
                dna2 += revComplementDict[base]
        return dna2[::-1]

def grabMoreNuc(start, dna):
        grabbedSequence = ""
        x = -1
        while x < 24: #num nucs wanted - 1
                grabbedSequence += dna[start + x]
                x+= 1
        return grabbedSequence



deletionDF = fileIntoDf(0, myDfBed).copy() #copy dataframe into variable


while x < 16:
        with open(os.path.join(workingDirectory, allChrTxt[x]), "r") as chrFile:
                y = 0

                dna = chrFile.read().replace('\n', '')   #turn txt to string
                end = len(dna) - 1

                filteredDf = deletionDF[deletionDF["Chromosome"] == ("chr%s"  %roman[x])]

                if filteredDf.empty:
                        x+= 1
                        continue
                else:
                        rangeIndex = len(filteredDf.index)
                        x += 1

                        for y in filteredDf.index:
                                
                                start = int(filteredDf.at[y, "Start(0)"])  #row by column
                                
                                if(deletionDF.at[y,"Strand"] == ("-")):
                                      
                                        if checkIfValid(start, end):
                                                sequence = grabMoreNuc(start, dna)
                                        else:
                                                sequence =grabSomeNuc(end, start)

                                        sequence = revComplement(sequence)
                                        deletionDF.at[y, "Nucleotides"] = sequence    
                                else:
                                        
                                        if checkIfValid(start, end):
                                                sequence = grabMoreNuc(start, dna)
                                        else:
                                                sequence =grabSomeNuc(end, start)


                                        deletionDF.at[y, "Nucleotides"] = sequence
                                                
deletionDF.to_csv(output, sep = '\t', header = False, index = False)