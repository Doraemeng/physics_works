#!/usr/bin/python

#
# A code for reading pgs output files (.lhco) or LHE (.lhe) files
# into event lists that can be manipulated in further analysis 
# (the main program comes at the end)
#


#***************************************** 
#*******   included libraries   ********** 
#***************************************** 

import sys
import os
from numpy import *        #... need to install numpy and check that you can do this command without error ...


#***************************************** 
#********   class definitions  *********** 
#***************************************** 

class run:
	
	#
	#... a "run" is a collection of events
	#

	def __init__(self,filename,fmt):
		#... PGS output file ...
		if fmt=='.lhco':
			fn=filename+fmt
			self.readinr(fn);

		#... MadEvent output file ...
		elif fmt=='.lhe':
			fn=filename+fmt
			self.readinl(fn);

		else:
			print ''
			print 'unsupported format...'
			print ''


	#
	#... readinr: creates a list of PGS events ...
	#

	def readinr(self, filename):
		self.eventlist=[];
		tt=open(filename,'r');
		tt.readline();
		att=tt.readlines();
		for i in range(len(att)):
			x=(att[i]).split()
			if (len(x) == 3 and x[0] == '0'):
				en=x[1];
				#print en
				j=1;
				rlist=[];				
				y=(att[i+1]).split();
				while len(y)>3:
					rlist.append(y);
					j+=1;
					k=i+j;
					if k<len(att):
						y=(att[k]).split();
					else:
						y=[];
				self.eventlist.append(revent(en,rlist));

		tt.close();


	#
	#... readinl: creates a list of LHE events ...
        #  


	def readinl(self,filename):
		self.eventlist=[];
		tt=open(filename,'r');

		s=tt.readline();
		while s!='<event>\n' :
			s=tt.readline()

	       	att=tt.readlines()
	       	for i in range(len(att)):
	       		x=(att[i]).split()
	       		if len(x)==6:
	       			en=x[1]
				wt=x[3]
				sc=x[4]
	       			j=1
	       			rlist=[]		
				y=(att[i+1]).split()
			       	while len(y)>3:
			       		rlist.append(y)
			       		j+=1
			       		k=i+j
			       		if k<len(att):
			       			y=(att[k]).split()
			       		else:
			       			y=[]
				self.eventlist.append(lhevent(en,wt,sc,rlist))
		tt.close();

class revent(run):

#
# ... an "revent" is a PGS event (list of PGS objects)
#		

	def __init__(self,evtnumber,llist):
		self.enum=evtnumber;
		self.objlist=[];
		for line in llist:
			#print line
			self.objlist.append(robj(line));
		self.objnum=len(self.objlist)	



	#... truthInvMass is a method of lhevent that looks for two particles that come from the same parent
	#... (specified by the caller) and calculates their invariant mass ...

	def truthInvMass(self,pid1,pid2):
		ob1list=[];
		ob2list=[];
		for obj in self.objlist:
			if (obj.typ==pid1 and obj.ntrk == '1.0'):
				ob1list.append(obj);
			if (obj.typ==pid2 and obj.ntrk == '-1.0'):
				ob2list.append(obj);
		if (len(ob1list)==0 or len(ob2list)==0) : 
			#print 'No such pair in this event...';
			#return -555;
			pass

		#else:
		if (len(ob1list)!=0 and len(ob2list)!=0) :
			for ob1 in ob1list:
				evtpt1 = array(ob1.pt, dtype = float64)
				evteta1 = array(ob1.eta, dtype = float64)
				evtphi1 = array(ob1.phi, dtype = float64)
				p41 = [evtpt1 * cosh(evteta1), evtpt1 * sin(evtphi1), evtpt1 * cos(evtphi1), evtpt1 * sinh(evteta1)]
			for ob2 in ob2list:
			    	evtpt2 = array(ob2.pt, dtype = float64)
				evteta2 = array(ob2.eta, dtype = float64)
				evtphi2 = array(ob2.phi, dtype = float64)
				p42 = [evtpt2 * cosh(evteta2), evtpt2 * sin(evtphi2), evtpt2 * cos(evtphi2), evtpt2 * sinh(evteta2)]
			minv=sqrt((p41[0]+p42[0])**2 - (p41[1] + p42[1])**2- (p41[2] + p42[2])**2 - (p41[3] + p42[3])**2)
			return minv;
			

		#else:
		if (len(ob1list)!=0 and len(ob2list)!=0) :
			for ob1 in ob1list:
				evtpt1 = array(ob1.pt, dtype = float64)
				evteta1 = array(ob1.eta, dtype = float64)
				evtphi1 = array(ob1.phi, dtype = float64)
				p41 = [evtpt1 * cosh(evteta1), evtpt1 * sin(evtphi1), evtpt1 * cos(evtphi1), evtpt1 * sinh(evteta1)]
			for ob2 in ob2list:
			    	evtpt2 = array(ob2.pt, dtype = float64)
				evteta2 = array(ob2.eta, dtype = float64)
				evtphi2 = array(ob2.phi, dtype = float64)
				p42 = [evtpt2 * cosh(evteta2), evtpt2 * sin(evtphi2), evtpt2 * cos(evtphi2), evtpt2 * sinh(evteta2)]
			minv=sqrt((p41[0]+p42[0])**2 - (p41[1] + p42[1])**2- (p41[2] + p42[2])**2 - (p41[3] + p42[3])**2)
			return minv;

	def eta_l(self,pid1,pid2):
		etaPlus = 0
		etaMinus = 0
		for obj in self.objlist:
			if (obj.typ==pid1 and obj.ntrk == '1.0'):
				etaPlus = obj.eta
		for obj in self.objlist:
			if (obj.typ==pid1 and obj.ntrk == '-1.0'):
				etaMinus = obj.eta
		if (etaPlus != 0 and etaMinus != 0):
			return etaPlus
		else:
			etaPlus = 0
			return etaPlus
		
	def phi_l(self,pid1,pid2):
		phiPlus = 0
		phiMinus = 0
		for obj in self.objlist:
			if (obj.typ==pid1 and obj.ntrk == '1.0'):
				phiPlus = obj.phi
		for obj in self.objlist:
			if (obj.typ==pid1 and obj.ntrk == '-1.0'):
				phiMinus = obj.eta
		if (phiPlus != 0 and phiMinus != 0):		
			return phiPlus
		else:
			phiPlus = 0
			return phiPlus
		

			
	def eta_j(self):
		j_pt_cut = 30;
		jet_mass = 0;
		eta = 0;
		for obj in self.objlist:
			if (obj.typ == '4' and obj.jmas > jet_mass and float(obj.pt) > j_pt_cut and float(obj.btag) != 0.0):
				jet_mass = obj.jmas	
				eta = obj.eta
				phi = obj.phi
		#if eta != 0:
		return eta
		
	def phi_j(self):
		j_pt_cut = 30;
		jet_mass = 0;
		phi = 0;
		for obj in self.objlist:
			if (obj.typ == '4' and obj.jmas > jet_mass and float(obj.pt) > j_pt_cut and float(obj.btag) != 0.0):
				jet_mass = obj.jmas	
				phi = obj.phi
		#if eta != 0:
		return phi
	
	def pt_j(self):
		j_pt_cut = 30;
		jet_mass = 0;
		pt = 0;
		for obj in self.objlist:
			if (obj.typ == '4' and obj.jmas > jet_mass and float(obj.pt) > j_pt_cut and float(obj.btag) != 0.0):
				jet_mass = obj.jmas	
				pt = obj.pt
		#if eta != 0:
		return pt





class lhevent(run):

#
# ... an "lhevent" is an LHE event (list of LHE objects, which are particles)
#

        def __init__(self,evtnumber,weight,scale,llist):
                self.enum=evtnumber;
		self.scale=scale;
		self.weight=weight;
                self.objlist=[];
                for line in llist:
                        self.objlist.append(lhobj(line));
                self.objnum=len(self.objlist)



	#... truthInvMass is a method of lhevent that looks for two particles that come from the same parent
	#... (specified by the caller) and calculates their invariant mass ...

	def truthInvMass(self,pid1,pid2):
		ob1list=[];
		ob2list=[];
		for obj in self.objlist:
			if obj.pdg==pid1:
				ob1list.append(obj);
			if obj.pdg==pid2:
				ob2list.append(obj);
		if (len(ob1list)==0 or len(ob2list)==0) :
			print 'No such pair in this event...';
			return -666;
		else:
			for ob1 in ob1list:
				p1=ob1.m1;
				p2=ob1.m2;
				for ob2 in ob2list:
					if (ob2.m1==p1 and ob2.m2==p2):
						fv1=array(ob1.fv,dtype=float64)
						fv2=array(ob2.fv,dtype=float64)
						f=fv1+fv2
						minv=sqrt(-f[0]*f[0]-f[1]*f[1]-f[2]*f[2]+f[3]*f[3])
						return minv;
			

class robj(revent):

#
# an "robj" is a PGS object (i.e., jet, lepton, photon,...)
#

	def __init__(self,str):
		self.num=str[0]   #The definitions of these fields can be found at J.Thaler's LHC Olympics wiki under "Data File Format"
		self.typ=str[1]
		self.eta=str[2]
                self.phi=str[3]
                self.pt=str[4]
                self.jmas=str[5]
                self.ntrk=str[6]
                self.btag=str[7]
                self.hadem=str[8]
		
	def printobj(self):
		print self.num + '\t' + self.typ + '\t' + self.eta + '\t' + self.phi + '\t' + self.pt + '\t' + self.jmas + '\t' + self.ntrk + '\t' + self.btag + '\t' + self.hadem

class lhobj(lhevent):

#
# an "lhobj" is an LHE object (a particle from the hard process)
#
 
 	def __init__(self,str):
		self.pdg=str[0];   #The definitions of these fields can be found on pg. 5 of hep-ph/0609017v1
                self.stat=str[1];
                self.m1=str[2];
                self.m2=str[3];
                self.c1=str[4];
                self.c2=str[5];
                self.fv=str[6:10];  #four-vector defined as a list
                self.mass=str[10];
                self.lt=str[11];
                self.spin=str[12];
	def printobj(self):
                print self.pdg," ", self.stat," ", self.m1," ", self.m2," ", self.c1," ", self.c2," ",self.fv," ", self.mass," ", self.lt," ", self.spin



#*****************************************
#**********   main program   *************
#***************************************** 

if len(sys.argv)!=2 :
	print 'Usage: ./analysis.py filename.lhco or filename.lhe'
else:	
	fname=sys.argv[1]
	fileName, fileExtension = os.path.splitext(fname)


	#... loads all events into a run object ...
	r2=run(fileName,fileExtension)

	#.... go through each event, calculate and print bbbar invt. mass ...
	count = 0
	for event in r2.eventlist:
		s = event.eta_j()
		t=event.truthInvMass('2','2')
		u = event.phi_j()
		v = event.pt_j()
		x = event.phi_l('2', '2')
		y = event.eta_l('2', '2')
		#w = float(event.phi_j()) - float(x)
		r = float(event.eta_j()) - float(y)
		
		if (float(s) != 0 and float(x) != 0):
			print u
			
	#print count
		
		

			
		

#****************************************
