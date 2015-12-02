
import random as rdm 
fp = open("dnaData","w")
dnaUnits = ['a','c','g','u']
dnaDataSize = input("Enter the Size of DNA sequence:")
dnaTestDataSize = input("Enter the Size of DNA Test Sequence:")

fp.write(str(dnaDataSize)+'\n')
fp.write(str(dnaTestDataSize) +'\n')

for i in range(0,dnaDataSize):
	fp.write(str(dnaUnits[rdm.randint(0,3)])+' ')

fp.write('\n')


for i in range(0,dnaTestDataSize):
	fp.write(str(dnaUnits[rdm.randint(0,3)])+' ')

fp.write('\n')