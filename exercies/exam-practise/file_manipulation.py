import os

file2_path = "exam-practise/file2.txt"
file1_path = "exam-practise/file.txt"

filemode = 'a' if os.path.exists(file2_path) else "w+"

f = open(file1_path,'r+')
f2 = open(file2_path,filemode)
tempLines = []

lines = f.readlines()
f.close()

with open(file1_path,"w+") as tempFile:
    for line in lines:
        # if line == "5000"
        replaced = line.replace('5000','1')
        # replaced = line.replace('1','5000')
        tempLines.append(replaced)
        tempFile.write(replaced)


# f2.write(tempLines)
print(tempLines)
for line in tempLines:
    f2.write(line)

f2.close()