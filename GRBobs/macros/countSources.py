#!/usr/bin/env python
#import os
import sys

class source:
    id        = 0
    name      = 'default'
    generated = 0
    generated = 0
    detected  = 0
    fileName  = 'default.txt'
    def dump(self):
        print '%4i, %30s, %6i, %6i in %10s '% (self.id,self.name,self.generated,self.detected,self.fileName)
    def save(self,outfile):
        outfile.writelines('%4i, %30s, %6i, %6i in %10s \n'% (self.id,self.name,self.generated,self.detected,self.fileName))

    def append(self,sourceList, initial):
        if(self.name.find(initial)>=0):
            added=0
#            N=1
            for s in sourceList:
                if (s.name.find(initial)>=0):
                    s.name=initial
                    s.generated+=self.generated
                    s.detected+=self.detected
                    s.fileName='more than one'
                    s.id+=1
#                    N=s.id
                    added=1
            if(added==0):
                self.id=1
                sourceList.append(self)
#            print 'Sources with ',initial,' : ', N

def SortbyID(a,b):
    return cmp(a.id,b.id)
def SortbyName(a,b):
    return cmp(a.name,b.name)
def SortbyFileName(a,b):
    return cmp(a.fileName,b.fileName)
def SortbyDetected(a,b):
    return cmp(a.detected,b.detected)
def SortbyGenerated(a,b):
    return cmp(a.generated,b.generated)



def countSources(fileName='list2.txt',append=0,initial='GRB'):
    f=file(fileName,'r')
    source_list=[]
    L=f.readlines()
    for l in L:
        if ( l.find('*')>=0 ):
            break
        line = l.split(',')
        s=source()        
        s.id        = 1
        s.name      = line[1]
        s.generated = int(line[2])
        lastTwo     = line[3].split('in')
        s.detected  = int(lastTwo[0])
        s.fileName  = lastTwo[1]
        source_list.append(s)
    f.close()
    source_list.sort(SortbyName)   
    New_source_list=[]
    for s in source_list:
        s.append(New_source_list,initial)
    New_source_list.sort(SortbyName)
    if(append==1):
        outputFile = file('list3.txt','a')
    else: 
        outputFile = file('list3.txt','w')
    s=New_source_list[0]
    s.save(outputFile)
    s.dump()
    return s.detected

    
if __name__ == '__main__':

    if sys.argv[2:]:
        print ' File name: ',sys.argv[1],' initials: ',sys.argv[2]
        foo = countSources(sys.argv[1],0,sys.argv[2])

    else:
        N= countSources(sys.argv[1],0," GRB_")
        N+= countSources(sys.argv[1],1," GRBtemplate")
        N+= countSources(sys.argv[1],1," pbh")
        N+= countSources(sys.argv[1],1," DC2_")
        N+= countSources(sys.argv[1],1," PSR_")
        N+= countSources(sys.argv[1],1," moon_")
        N+= countSources(sys.argv[1],1," sun_")
        N+= countSources(sys.argv[1],1," SolarFlare")
        N+= countSources(sys.argv[1],1," lcc2")
        N+= countSources(sys.argv[1],1," Plerion")
        N+= countSources(sys.argv[1],1," SC_blazar")
        N+= countSources(sys.argv[1],1," SC1_extragal")
        N+= countSources(sys.argv[1],1," SC1_GALPROP")
        N+= countSources(sys.argv[1],1," Galactic")
        N+= countSources(sys.argv[1],1," SNR")
        N+= countSources(sys.argv[1],1," Galaxies_")
        N+= countSources(sys.argv[1],1," Giommi_BLLac")
        N+= countSources(sys.argv[1],1," Giommi_FSRQ")
        N+= countSources(sys.argv[1],1," OB")
        N+= countSources(sys.argv[1],1," J")
        print N,19279683-N
        
