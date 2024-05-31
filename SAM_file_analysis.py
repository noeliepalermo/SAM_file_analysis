#!/usr/bin/env python3
# -*- coding: utf-8 -*-
 #Project done for the Master 1 in Bioinformatics of Montpellier - UE HAI724I Système
# Our bioinformatic project consisted in the processing and analysis of the data resulting from the alignment of the sequencing data stored in a SAM format file.
# __authors__ = ("Fatima-Zahra ABANI", "Noëlie PALERMO")
#__contact__ = ("abanifatimazahra@gmail.com>,"palermo.n@live.fr")
#__version__ = "0.0.1"
#__date__ = "16/10/2022"
#__licence__ ="This program is free software: you can redistribute it and/or modify
    ##it under the terms of the GNU General Public License as published by
    ##the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
    ##This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of
    ##MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
    ##You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>."
# OPTION LIST:
    ##-h or --help : help information
    ##-i or --input: input file (.sam)
    ##-o or --output: output name files (.txt)

#Synopsis:
    ##SamReader.py -h or --help # launch the help.
    ##SamReader.py -i or --input <file> # Launch SamReader to analyze a samtools file (.sam) and print the result in the terminal
    ##SamReader.py -i or --input <file> -o or --output <name> # Launch SamReader to analyze a samtools file (.sam) and print the result in the file called <name>

############### IMPORT MODULES ###############
import os, sys, re, pathlib

############### FUNCTIONS TO :

###Each function can work with a SAM file input in the sam function

## 1/ Check, the input file, it's must be in the directory of your terminal

def sam():
    #Ask the user to input their SAM file
    fichier = input("Please input the SAM file: ")
    #Lecture of the SAM file : the input must be a file name
    if os.path.isdir(fichier):
        print("ERROR: the input is not for a file")
        exit()
    else:
        #Verification of the format of the SAM file
        extension = fichier.split(".") #Stock the extension of the file
        if len(extension)>1:
            if extension[-1] == "sam" or extension[-1] == "SAM":
                return(fichier)
            else:
                print("ERROR: the file is not in the SAM format")
                exit()
            #Verification if the file is empty
            if os.path.getsize(fichier) == 0:
                print("ERROR: the SAM file is empty!") #They must be no blank space in the begining of the file !
                exit()

## 2/ Read and stock, the input SAM file in a dictionnary

def read(fichier):
    fichiers = open(fichier)
    #Creation of the dictionnary for the SAM datas
    data = {}
    header = []
    for lines in fichiers:
        if lines[0] == "@": #Stock headers in a list
            header.append(lines[0:])
        else: #Remove headers
            #Split the elements of the SAM file into columns
            lignes = lines.split("\t")
            readName = lignes[0]
            #Check if the name of the reads have special characters
            if readName[-2:] == ".1" or readName[-2:] == ".2" or readName[-2:] == "_1" or readName[-2:] == "_2" or readName[-2:] == "/1" or readName[-2:] == "/2":
                readName = readName[:len(readName)-2]
            #Check if double reads in the file
            if readName in data.keys(): #The read's name are keys of the dictionnary
                data[readName].append(lignes[1:]) #Create a second list in the dictionnary coresponding to the mate of the read
            #Create a dictionnary
            else:
            #For this script the informations are used: l[0] is the FLAG, l[1] is the reference, l[2] is the reference's position l[3] is the quality of mapping, l[4] is the CIGAR informations
            #others will not be analysed
                data[readName] = [lignes[1:]]
    return(header, data)

## 3/Write return results of analyses functions in SAM file
#Parameters are the dictionnary, headers of the SAM file, the name of the input file, and the name of the function

def writeSAM(data, header, fichier, nom):
    #Creation of a variable with: the name of the input SAM file_the name of the analysis step.sam
    name = fichier.rstrip(".sam")+"_"+nom+".sam"
    with open(name, "w", encoding="utf-8") as line:
        #Add the header in the begining of the file
        for h in header:
            line.write(h)
        for values in data.items():
            for v in range(len(values[1])):
                #Separate each elements of the dictionnary to create a SAM file
                line.write(values[0]+"\t"+"\t".join(values[1][v]))

## 4/ Analyse

    ### 4/a Counting number of reads in a dictionnary

def count(data):
    counter = 0
    #Count the number of reads for each keys of a dictionnary and return the sum
    for read in data.values():
        counter = counter+len(read)
    return(counter)

    ### 4/b Verification of the "Reference" annotation (reads are mapped or not)
    #If the reference = "*" or the position reference = 0, the reads are not mapped

def reference(data):
    #Creation of two dictionnary to stock the referenced and non referenced reads
    readReferenced = {}
    noReferencedRead ={}
    for read in data.items():
        #Check the references annotations for the reads and it's mate and return the two dictionnaries
        if read[1][0][1] == "*" or read[1][0][2] == "0" or read[1][1][1] == "*" or read[1][1][2] == "*": 
            noReferencedRead[read[0]] = read[1]
        else:
            readReferenced[read[0]] = read[1]
    return(readReferenced, noReferencedRead)

 ### 4/c Quality control
    #If a reads is strictly inferior to the limit value of the quality control it is dismiss. 
    #Users input they value and the default value is 20.

def quality(data):
    #Creation of two dictionnary to stock the reads with a superior or equal value to the quality and to contain the reads with a quality value strictly inferior to the quality
    readQuality = {}
    noQualityRead = {}
    #Default quality value
    qualite = 20
    #The user can input their own quality value.
    ql = input("The default quality value of mapping is 20. Do you want to change the value ? \n y/n ")
    if ql == 'y' or ql == 'Y':
        qualite = int(input("Input your quality value: "))
    #The default quality value of mapping still be 20 
    elif ql == 'y' or ql == 'N': 
        qualite = 20
    #Comparison between the given/default quality value and the quality of mapping (colum 3) of each reads.
    for q in data.items():
        #If a reads is superiror or equal to the quality value it is stock in the quality dictionnary. 
        if int(q[1][0][3]) <= qualite or int(q[1][1][3]) <= qualite:
            noQualityRead[q[0]] = q[1]
        else:
            readQuality[q[0]] = q[1] 
    return(qualite, readQuality, noQualityRead)

### 4/d Flag control

def samFlag(data):
    #Creation of a list with SAM Flags (https://www.samformat.info/sam-format-flag)
    #For more explanaition of the SAM Flags signification, please take a look at the Read_me.txt or the upper internet link
    decimal = [2048, 1024, 512, 256, 128, 64, 32, 16, 8, 4, 2 ,1]
    #Creation of the dictionnary of will contains paired reads or not
    readPaired = {}
    noPairedRead = {}
    for read in data.items():
        decimalRead = 2048
        decimalMate = 2048
        #Recuperation of the Flag references in the SAM file for the read and it's mate
        flagRead = read[1][0][0]
        flagMate = read[1][1][0]
        for n in range(len(decimal)):
            #If the decimal in the list is inferior to the Flag reference of the read, it will be substract
            #The difference will be stock in a variable before being again compare to the next decimal in the list
            if decimal[n] <= int(flagRead):
                flagRead = int(flagRead)-decimal[n]
                decimalRead = decimal[n]
                #If the decimal is equal to 2, then we can compare the Flag reference of the mate to the decimals of the list
                if decimalRead == 2:
                    for m in range(len(decimal)):
                        if decimal[m] <= int(flagMate):
                            flagMate = int(flagMate )-decimal[m]
                            decimalMate = decimal[m]
                            #If the decimal is again a 2, then reads are properly paired and will be stock in a new dictionnary
                        if decimalMate == 2:
                            readPaired[read[0]] = read[1]
                            break
                        #If the decimal of the mate is equal to only a 1 (read paired) or to 8 (the mate is unmapped) or to 4 (the read is unmapped)
                        #Or to 256 (the read is not in the primary alignement) or 1024 (the read is PCR or optical duplicate) or to 2048 (the read is a suplementary alignement)
                        #Then reads are not paired and will be stock in a new dictionanry
                        elif decimalMate == 1 or decimalMate == 8 or decimalMate == 4 or decimalMate == 256 or decimalMate == 1024 or decimalMate == 2048 :
                            noPairedRead[read[0]] = read[1]
                    break
                #If the decimal of the read is equal to only a 1 or to 8 or to 4 or to 256 or to 1024 or to 2048
                #Then reads are not paired and will be stock in a new dictionanry
                elif decimalRead == 1 or decimalRead == 8 or decimalRead == 4 or decimalMate == 256 or decimalMate == 1024 or decimalMate == 2048 :
                    noPairedRead[read[0]] = read[1]
    return(readPaired, noPairedRead)

### 4/d CIGAR control

def cigar(data):
    #Control of the CIGAR quality if it's equal to a number superior to 100M, the reads are aligned
    readAligned = {}
    readNoAligned = {}
    #Pattern of the CIGAR must be a number strictly superior to 0
    pattern = "[1-9]{1}\d{1,4}M"
    for read in data.items():
        #Match between the pattern and the read CIGAR and the mate CIGAR
        resultRead = re.match(pattern, read[1][0][4])
        resultMate = re.match(pattern, read[1][1][4])
        if resultRead:
            if resultMate:
                readAligned[read[0]] = read[1]
            else:
                readNoAligned[read[0]] = read[1]
        else:
            readNoAligned[read[0]] = read[1]
    return(readAligned, readNoAligned)

### 5/ Summary of the analysis

def resume(counter, readReferenced, qualite, readQuality, readPaired, readAligned):
    #Summary of the number of reads which passed the refererence, quality, paired and aligned controls
    print("Resume of the SAM file analysis: ")
    print("Total number of reads: ", counter)
    print("Number of referenced reads: ", count(readReferenced))
    print("Number of reads with a quality superior to ", int(qualite), ": ", count(readQuality))
    print("Number of paired reads: ", count(readPaired))
    print("Number of aligned reads: ", count(readAligned))

#### Main function ####     
def main(argv):
    #Function sam for read the file
    fichier = sam()
    #Function to create the dictionnary of the sam file
    header, dictio = read(fichier) 
    #Function to count reads inside a dictionnary
    counter = count(dictio)
    #Function to check if the reads are referenced
    readReferenced, noReferencedRead = reference(dictio)
    #Function to control the quality of the reads
    qualite, readQuality, noQualityRead = quality(dictio)
    #Function to control the FLAG references of the reads
    readPaired, noPairedRead = samFlag(dictio)
    #Function to control the CIGAR of the reads
    readAligned, readNoAligned =cigar(dictio)
    #Function to print in the terminal the resume of analysis
    result = resume(counter, readReferenced, qualite, readQuality, readPaired, readAligned)
    #Function to write the results of functions into SAM file
    writeSAM(readReferenced, header, fichier, "referenced")
    writeSAM(noReferencedRead, header, fichier, "no_referenced")
    writeSAM(readQuality, header, fichier, "quality"+str(qualite))
    writeSAM(noQualityRead, header, fichier, "no_quality"+str(qualite))
    writeSAM(readPaired, header, fichier, "paired")
    writeSAM(noPairedRead, header, fichier, "no_paired")
    writeSAM(readAligned, header, fichier, "aligned")
    writeSAM(readNoAligned, header, fichier, "no_aligned")
    

############### LAUNCH THE SCRIPT ###############
if __name__ == "__main__":
    main(sys.argv[1:])
