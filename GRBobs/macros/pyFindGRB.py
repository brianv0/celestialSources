#!/usr/bin/env python2.3
import os
import re
import sys
#import time    
import datetime
import pyParse
import math

from datetime import timedelta
#DataDir='/users/omodei/Documents/DATA/SC1/output_GRBobs_Diffuse_61221_HANDOFF_FT2'
#output_only_GRB_061128_HANDOFF'
#XMLFile='/users/omodei/Documents/DATA/SC1/skymodel/grb/GRBobs_user_library_SC1_v1.xml'
#DataDir="../data/output_OnlyGRBs_061108"
#XMLDir="../data/xml/"
#baseName='all'

class grb:
    def __init__(self,mc_src_id ):
        self.mc_src_id  = mc_src_id
	self.mc_src_num = int(mc_src_id.split('GRB_')[1])
        self.grb_name = 'GRB010101000'
        self.grbdate=datetime.datetime(2001,1,1,0,0,0) 
        self.MET  =  0.0    
        self.l    =  0.0
        self.b    =  0.0
        self.theta=  0.0
        self.phi  =  0.0
        #parameters from the XML library:        
        self.time = 0.0
        self.duration  = 0.0
        self.peakFlux  = 0.0
        self.redshift  = 0.0
        self.alpha  = 0.0
        self.beta  = 0.0
        self.Ep  = 0.0
        self.MinExtractedPhotonEnergy  = 0.0
        self.Essc_Esyn =0.0
        self.Fssc_Fsyn =0.0
        self.GBMFlag =0.0
        self.Nph =0.0
        self.Dt_extended=0.0
        self.Duration_extended=0.0
        self.CutOff=0.0
        # get from the list srcId file:
        self.mc_id=0
        self.generated=0
        self.detected=0
        self.file=''
    def computeGRBName(self):
        refdate = datetime.datetime(2008, 1, 1,0,0,0)
        metdate = datetime.datetime(2001, 1, 1,0,0,0)
        dt     = timedelta(seconds=self.time)
        self.grbdate=refdate+dt
        tt = (self.grbdate).timetuple()
        print tt
        yr = tt[0]
        mo = tt[1]
        da = tt[2]
        hr = tt[3]
        mi = tt[4]
        se = tt[5]
#        self.grbdate=datetime.datetime(yr, mo,da,hr,mi,se)
        self.MET=(self.grbdate-metdate).days*86400+(self.grbdate-metdate).seconds
       
        yrs=str(yr)[2:]
        if(mo<10):
            mos='0'+str(mo)
        else:
            mos=str(mo)
        if(da<10):
            das='0'+str(da)
        else:
            das=str(da)
    
        rem = (se+60*mi+3600*hr)*1000/86400
        if(rem<10):
            rems='00'+str(rem)
        elif(rem<100):
            rems='0'+str(rem)
        else:
            rems=str(rem)
            
        self.grb_name = yrs+mos+das+rems
        
    def CheckDEFFile(self):
        fileName = DataDir+'/GRBOBS_'+self.grb_name+'.DEF'
        if os.path.exists(fileName): 
            return 1
        else:
            print fileName+' not found'
            return 0
	
    def ParseParFile(self):
        fileName = DataDir+'/GRBOBS_'+self.grb_name+'_PAR.txt'

        if not os.path.exists(fileName):
            print fileName+' not found'
            return 0
        
        f=file(fileName,'r')
        f.readline()
        l=f.readline()
        params=l.split()
	
        if(params[0]!=self.grb_name):  # is this enougth?
            print '***************** WARNING!!! *************' 
        self.l=float(params[3])
        self.b=float(params[4])
        self.theta=float(params[5])
        self.phi=float(params[6])
        return 1
    def ParseListFile(self):
        fileName = 'SourceList1.txt'
        f=file(fileName,'r')
        L=f.readlines()
        for l in L:
            if l.find(self.mc_src_id.strip())>0: 
                vars = l.split()    
                self.mc_id= vars[0]
		print vars[2],' - ',vars[3]
                self.generated = int(vars[2].split(',')[0])
                self.detected = int(vars[3].split(',')[0])
		print self.generated,' - ',self.detected
                self.file = vars[5].strip()
                return
        
#        print params

	
    def Print(self):
        print 100*'_'
        print self.mc_src_id,' alias GRB',self.grb_name
        print 'Date: ', self.grbdate,' MET: ',self.MET,' time since sim. start: ',self.time
        print 'Duration: ', self.duration,' PeakFlux: ',self.peakFlux,' z: ',self.redshift, 'alpha: ',self.alpha,' beta: ', self.beta ,'ep: ', self.Ep , self.MinExtractedPhotonEnergy
        print 'l: ',self.l,' b: ', self.b ,'theta: ', self.theta , 'phi ',self.phi
        print 'ssc: ',self.Essc_Esyn ,  self.Fssc_Fsyn , self.GBMFlag
        print 'mc Id: ',self.mc_id , ' generated photons: ', self.generated,'detected photons: ',self.detected,' in file: ',self.file

    def Dump(self, fileOut):
	fileOut.write('<tr>')
	list = [str(self.mc_src_id),str(self.grb_name),str(self.MET),str(self.time),str(self.duration),str( \
            self.peakFlux),str(self.redshift),str(self.alpha),str(self.beta),str(self.Ep),str( \
            self.l),str(self.b),str(self.theta),str(self.phi),str(self.Fssc_Fsyn),str(self.generated),str(self.detected)]
	for i in list:
	    line='<td style="vertical-align: top;">'+i+'<br></td>'
	    fileOut.write(line)
	    print(line)
	fileOut.write('</tr>')

def HTML_Header(fileOut):
    fileOut.write('<html>')
    fileOut.write('<head>')
    fileOut.write('<meta content="text/html; charset=ISO-8859-1"')
    fileOut.write('http-equiv="content-type">')
    fileOut.write('<title></title>')
    fileOut.write('</head>')
    fileOut.write('<body>')
    fileOut.write('<table style="text-align: left; width: 100%;" border="1">')
    fileOut.write('<tbody>')
    fileOut.write('<tr>')
    list = ['mc_src_id','grb_name','MET','time','duration','peakFlux','redshift',
	    'alpha','beta','Ep','l','b','theta','phi','externalComponent','generated','detected']
    for i in list:
	line='<td style="vertical-align: top;">'+i+'<br></td>'
	fileOut.write(line)
	print(line)
    fileOut.write('</tr>')

    
    
    
def HTML_Tail(fileOut):
    fileOut.write('</tbody>')
    fileOut.write('</table>')
    fileOut.write('<br>')
    fileOut.write('<br>')
    fileOut.write('</body>')
    fileOut.write('</html>')
    
#    def SetRootTree(self, tree):
#        tree.SetBranchAddress('src_id',self.mc_src_id,'src_id\I')
        

def Parse_xml_file(xmlLib, src_name):
    try:
        ''' values = re.search('(?<=abc)def', 'abcdef') '''
        values = re.search('(?<=<source name=\" %s).*(?=\">)' %(src_name) ,xmlLib)
        return values.groups()
    except AttributeError:
        return None;
        
def ReadLib(xml_file):
     f=file(xml_file,'r')
     L=f.readlines()
     xmlLib=''
     for l in L:
         xmlLib+=l 
     return xmlLib

def ReadLib1(xml_file):
    
     f=file(xml_file,'r')
     L=f.readlines()
     grbList=[]
     params=''
     for l in L:
         b=1

         try:
             name   = re.search('(?<=<source name=\").*(?=\">)' ,l)
             if(len(name.group())>0):
                 mc_src_name = name.group()
                 print 'Name:', mc_src_name
             b=1
         except AttributeError:
             b=0
         try:
             params = re.search('(?<= <SpectrumClass name="GRBobsmanager" params=").*(?=\"/>)' ,l).group()
             if(len(params)>0):
                 params=params.split(' , ')
                 _grb=grb(mc_src_name)     
                 _grb.time = int(params[0])
                 _grb.duration =float(params[1])
                 _grb.peakFlux =float(params[2])
                 _grb.redshift =float(params[3])
                 _grb.alpha =float(params[4])
                 _grb.beta =float(params[5])
                 _grb.Ep =float(params[6])
                 _grb.MinExtractedPhotonEnergy =float(params[7])
                 _grb.Essc_Esyn=float(params[8])
                 _grb.Fssc_Fsyn=float(params[9])
                 _grb.GBMFlag=float(params[10])
                 _grb.Nph=float(params[11])
                 _grb.Dt_extended=float(params[12])
                 _grb.Duration_extended=float(params[13])
                 _grb.CutOff=float(params[14])
                 _grb.computeGRBName()
                 print "Alias:",_grb.grb_name 
                 b=1
                                  
         except AttributeError:
             b=0
             
         if (b==1):
             grbList.append(_grb)
     return grbList

#Call with the xml file as input
    
if __name__ == '__main__':
    sort=0
    if sys.argv[3:]:
	DataDir=sys.argv[1]
	xml_file=sys.argv[2]
	baseName=sys.argv[3]
    else:
	print ' --------------------------------------------------'
	print ' To run this macro, you need a datadir, with all the GRB*par and def files'
	print ' and a xml file, to cross check the sources'
	print ' and a basename defined'
	print ' pyParse datadir xml_file basename'
	sys.exit()
	    
    pyParse.pyParse(DataDir,baseName, sort)
    grb_list = ReadLib1(xml_file)
    print len(grb_list)
    burst_ok=0;
    burst_gbm=0;
    burst_batse=0;
    burst_tgbm=0;
    burst_tbatse=0;

    fullSky=4*math.pi

    thresholds=[0.2,0.3,0.35,0.5, 0.75, 1] # 1.024 s accumulation time
    fovs=[0.5,8/fullSky,9/fullSky]
    detected=[]
    t = len(thresholds)
    f = len(fovs)
    i=0
    while i < t*f:
        detected.append(0)        
        i+=1    
	
    fout = file('All_src_List.html','w')
    fout_detected = file('Detected_List.html','w')
    import ROOT
    from array import array
    tree=ROOT.TTree('GrbCatalog','GRB Catalog')
    print 'imported root'

##################################################
    #mc_src_id = array('s',[0])
    mc_src_num = array('i',[0])
    MET  =  array('f',[0]) 
    l    =  array('f',[0])
    b    =  array('f',[0])
    theta=  array('f',[0])
    phi  =  array('f',[0])
    duration  = array('f',[0])
    peakFlux  = array('f',[0])
    redshift  = array('f',[0])
    alpha  = array('f',[0])
    beta  = array('f',[0])
    Ep  = array('f',[0])
    Essc_Esyn =array('f',[0])
    Fssc_Fsyn =array('f',[0])
    GBMFlag =array('f',[0])
    Nph =array('f',[0])
    Dt_extended=array('f',[0])
    Duration_extended=array('f',[0])
    CutOff=array('f',[0])

    generated_ph = array('i',[0])
    detected_ph  = array('i',[0])

#    tree.Branch('mc_src_id',mc_src_id,'mc_src_id/S')
    tree.Branch('MET',MET,'MET/F')
    tree.Branch('l',l,'l/F')
    tree.Branch('b',b,'b/F')
    tree.Branch('theta',theta,'theta/F')
    tree.Branch('phi',phi,'phi/F')
    tree.Branch('duration',duration,'duration/F')
    tree.Branch('peakFlux',peakFlux,'peakFlux/F')
    tree.Branch('redshift',redshift,'redshift/F')
    tree.Branch('alpha  ',alpha,'alpha/F')
    tree.Branch('beta  ',beta,'beta/F')
    tree.Branch('Ep  ',Ep,'Ep/F')
    tree.Branch('Essc_Esyn',Essc_Esyn,'Essc_Esyn/F')
    tree.Branch('Fssc_Fsyn',Fssc_Fsyn,'Fssc_Fsyn/F')
    tree.Branch('Nph',Nph,'Nph/F')
    tree.Branch('Dt_extended',Dt_extended,'Dt_extended/F')
    tree.Branch('Duration_extended',Duration_extended,'Duration_extended/F')
    tree.Branch('CutOff',CutOff,'CutOff/F')
    tree.Branch('generated',generated_ph,'generated/I')
    tree.Branch('detected',detected_ph,'detected/I')
    print '**************************************** fava **************************************************'
    
    HTML_Header(fout_detected)
    HTML_Header(fout)

    for grb in grb_list:
        check = grb.ParseParFile()
        check *= grb.CheckDEFFile()
        if check>0:
            theta_rad = (90.0 - float(grb.theta))/180.0*math.pi
            grb.ParseListFile()
            grb.Print()
            grb.Dump(fout)
	    if(grb.detected>0):
		grb.Dump(fout_detected)
            burst_ok+=1 
	    ############

	    mc_src_num[0]=grb.mc_src_num
	    l[0]=grb.l
	    b[0]=grb.b
	    theta[0] = grb.theta
	    phi[0]   = grb.phi
	    alpha[0] = grb.alpha
	    beta[0]  = grb.beta
	    Ep[0]    = grb.Ep
	    peakFlux[0] = grb.peakFlux
	    generated_ph[0] = grb.generated
	    detected_ph[0]  = grb.detected
	    Essc_Esyn[0]    = grb.Essc_Esyn
	    Fssc_Fsyn[0]    = grb.Fssc_Fsyn

	    tree.Fill()	
	    ############
            ti=0
            fi=0
            i=0

	    

            for fov in fovs:
                for th in thresholds:
                    angle=(1.0 - 2 * fov)
                    if(grb.peakFlux>th and math.cos(theta_rad)>angle): 
                        detected[i]+=1
                        
                    i+=1

    HTML_Tail(fout_detected)
    HTML_Tail(fout)

    fout.flush()
    fout.close()
    fout_detected.flush()
    fout_detected.close()
    rootf=ROOT.TFile('GRBCatalog.root','RECREATE')
    tree.Write()
    
    print'Successfully Done',burst_ok,' out of',len(grb_list)
    i=0
    for fov in fovs:    
        for th in thresholds:
            print' Threshold: ',th,'(ph/cm^2/s), fov: ',fov,' : ',detected[i],' GRB/yr, assuming SAA (14%:)',detected[i]*(1-0.14)
            i+=1
