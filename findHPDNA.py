import os
import numpy as np
import pandas as pd

workingDirectory = os.path.dirname(__file__)
allChrTxt = ["chr1.txt", "chr2.txt", "chr3.txt", "chr4.txt", "chr5.txt", "chr6.txt", "chr7.txt", "chr8.txt",
             "chr9.txt", "chr10.txt", "chr11.txt", "chr12.txt", "chr13.txt", "chr14.txt", "chr15.txt",
             "chr16.txt"]
roman = ["I", "II", "III", "IV", "V", "VI", "VII", "VIII", "IX", "X", "XI", "XII", "XIII", "XIV", "XV", "XVI"]
output = "genomeHPcount.txt"


HPLengthsAll = np.zeros(45, dtype =int)
data = (HPLengthsAll, 0)


def countHomopoly(dna, lenDNA, data):

        TempHPLengthsAll = np.zeros(45, dtype = int)
        base = 0
        HPcount = 0

        while base < lenDNA:

                HPlength = 1            #start with one because minimum is one

                while base + 1 < lenDNA and dna[base] == dna[base + 1]: #iterate while not at the end / is the same nucleotide
                        HPlength += 1
                        base += 1 

                if HPlength > 3:
                      HPcount+= 1
                
                TempHPLengthsAll[HPlength - 1] += 1
                
                base += 1


        finalHPLengthsAll = np.add(data[0], TempHPLengthsAll) #add to previous values
        finalHPCount = data[1] + HPcount

        return (finalHPLengthsAll, finalHPCount)


for x in range(len(allChrTxt)):
    with open(os.path.join(workingDirectory, allChrTxt[x]), "r") as file:
        dna = file.read().replace('\n', '')
        dnaLen = len(dna)

        data = countHomopoly(dna, dnaLen, data)


df = pd.DataFrame(data[0]) #make the array from the tuple into a dataframe

with open(os.path.join(workingDirectory, output), "w") as w:
    df_str = df.to_string(header = False, index = False)
    w.write("HP Count at Each Length Starting at 1: \n")
    w.write(df_str)
    w.write("\n\n\nHP Count: ")
    w.write(str(data[1]))

