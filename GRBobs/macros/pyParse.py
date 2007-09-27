#!/usr/bin/env python
import os
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

    def append(self,sourceList):
        for s in sourceList:
            if s.name==self.name:
                s.generated+=self.generated
                s.detected+=self.detected
                s.fileName='more than one'
                return
        sourceList.append(self)
                
        


        
def SortbyID(a,b):
    return cmp(a.id,b.id)
def SortbyFileName(a,b):
    return cmp(a.fileName,b.fileName)
def SortbyName(a,b):
    return cmp(a.Name,b.Name)
def SortbyDetected(a,b):
    return cmp(a.detected,b.detected)
def SortbyGenerated(a,b):
    return cmp(a.generated,b.generated)


            
def pyParse(dirname,baseName, sort=0):
    
    if not os.path.exists(dirname):
        print 'pyParse: ',dirname,' not exists'
        return 0 
    file_list=os.listdir(dirname)
    new_file_list=[]
    for file_name in file_list:
        if '_srcIds.txt' in file_name and baseName in file_name:
            new_file_list.append(file_name)
            
    print new_file_list
    outputFile = file('MergedIdSource.txt','w')
    outputFile.writelines('mc_src_id source generated detected\n')

    sourcelist=[]
    onlygrbs=[]
    for file_name in new_file_list:
        f=file(dirname+'/'+file_name,'r')
        L=f.readlines()
        outputFile.writelines(f.readlines())
        for l in L:
            line = l.split()
            s=source()
            s.id        = int(line[0])
            s.name      = line[1]
            s.generated = int(line[2])
            s.detected  = int(line[3])
            s.fileName  = file_name
            sourcelist.append(s)
	    if (s.name.find('GRB')>=0):
		s.dump()
		onlygrbs.append(s)

	f.close()
    outputFile.close()

    if(sort==0):
        print '*************************    SORT BY ID:  **********************************'
        sourcelist.sort(SortbyID)
    elif(sort==1):
        print '***********************    SORT BY FILE NAME *******************************'
        sourcelist.sort(SortbyFileName)
    elif(sort==2):
        print '***********************    SORT BY GENERATED *******************************'
        sourcelist.sort(SortbyGenerated)
    elif(sort==3):
        print '***********************    SORT BY DETECTED  *******************************'
        sourcelist.sort(SortbyDetected)
    elif(sort==4):
        print '***********************    SORT BY NAME  *******************************'
        sourcelist.sort(SortbyName)
	

    onlygrbs.sort(SortbyID)
    outputFile = file('GRB_list.txt','w')
    for s in onlygrbs:
        s.save(outputFile)
        s.dump()
    outputFile.close()

    
    print '**************************************************'


    outputFile = file('SourceList1.txt','w')
    for s in sourcelist:
	s.save(outputFile)
    #s.dump()    
    outputFile.close()
    print '**************************************************'
    Newsourcelist=[]
    sourcelist.sort(SortbyID)
    for s in sourcelist:
	s.append(Newsourcelist)
    print '*******    SUMMED SOURCES (Sorted by Generated) ********'
    #    Newsourcelist.sort(SortbyID)
    
    Newsourcelist.sort(SortbyGenerated)
    outputFile = file('SortedSourceList.txt','w')
    for s in Newsourcelist:
	s.save(outputFile)
    #s.dump()
    outputFile.close()
    print '**************************************************'

    totalGenerated=0
    totalDetected=0

    moreThan1p =0
    moreThan5p =0
    moreThan10p =0
    moreThan100p =0
    moreThan1000p =0
    moreThan10000p =0
    
    for s in Newsourcelist:
        totalGenerated += s.generated
        totalDetected  += s.detected
        if(s.detected>0): moreThan1p+=1
        if(s.detected>=5): moreThan5p+=1
        if(s.detected>=10): moreThan10p+=1
        if(s.detected>=100): moreThan100p+=1
        if(s.detected>=1000): moreThan1000p+=1
        if(s.detected>=10000): moreThan10000p+=1
        
    print 'Number of Generated photons: ',    totalGenerated,' Detected : ', totalDetected
    print 'Total Number of sources: ', len(Newsourcelist)
    print 'Number of sources (>1     ph): ',moreThan1p
    print 'Number of sources (>5     ph): ',moreThan5p
    print 'Number of sources (>10    ph): ',moreThan10p
    print 'Number of sources (>100   ph): ',moreThan100p
    print 'Number of sources (>1000  ph): ',moreThan1000p
    print 'Number of sources (>10000 ph): ',moreThan10000p

if __name__ == '__main__':
    if sys.argv[2:]:
        foo = pyParse(sys.argv[1],sys.argv[2])
    elif sys.argv[1:]:
        foo = pyParse(sys.argv[1])

