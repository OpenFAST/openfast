import sys
from shutil import move

def tempFileName(inputFile):
    return inputFile+".tmp"

def setupWorkingFiles(inputFile):
    fin = open(inputFile, "r")
    fout = open(tempFileName(inputFile), "w")
    return fin, fout

def moveTempToPermanent(inputFile):
    move(tempFileName(inputFile), inputFile)

def commentLineContaining(targetWord, inputFile):
    fin, fout = setupWorkingFiles(inputFile)
    for line in fin:
        if targetWord.lower() in line.lower():
            line = "!" + line
        fout.write(line)
    moveTempToPermanent(inputFile)

def swapWordInFile(oldWord, newWord, inputFile):
    # capitalization matters
    fin, fout = setupWorkingFiles(inputFile)
    for i,line in enumerate(fin):
        if oldWord in line:
            line = line.replace(oldWord, newWord)
        fout.write(line)
    moveTempToPermanent(inputFile)

if __name__=="__main__":
    if sys.argv[1] == "privateToPublic":
        commentLineContaining("private", sys.argv[2])
    elif sys.argv[1] == "swapWordInFile":
        swapWordInFile(sys.argv[2], sys.argv[3], sys.argv[4])
