#include <iostream>
#include <stdio.h>     
#include <stdlib.h>

#include "TTree.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "Math/Vector3D.h"

#include "WCSimRootOptions.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootEvent.hh"

using std::endl;
using std::cerr;
using std::cout;
using ROOT::Math::XYZVector;

// Simple example of reading a generated Root file
int pmt_plotting(const char *filename="../wcsim.root") 
{


	// Open the file
	TFile * file = new TFile(filename,"read");
	if (!file->IsOpen()){
		cout << "Error, could not open input file: " << filename << endl;
		return -1;
	}

	// Get the a pointer to the tree from the file
	TTree *tree = (TTree*)file->Get("wcsimT");
	tree->SetBranchStatus("wcsimrootevent2", 0);
	tree->SetBranchStatus("wcsimrootevent_OD", 0);
	tree->SetBranchStatus("wcsimrootevent", 1);

	// Get the number of events
	const long nevent = tree->GetEntries();

	// Create a WCSimRootEvent to put stuff from the tree in
	WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();

	// Set the branch address for reading from the tree
	TBranch *branch = tree->GetBranch("wcsimrootevent");
	branch->SetAddress(&wcsimrootsuperevent);

	// Force deletion to prevent memory leak 
	tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

	// Geometry tree - only need 1 "event"
	TTree *geotree = (TTree*)file->Get("wcsimGeoT");
	WCSimRootGeom *geo = 0; 
	geotree->SetBranchAddress("wcsimrootgeom", &geo);
	if (geotree->GetEntries() == 0) {
		exit(9);
	}
	geotree->GetEntry(0);
	int npmt = geo->GetWCNumPMT();
	for(int i = 0; i < npmt; i++){
		WCSimRootPMT pmt = geo->GetPMT(i);
		XYZVector pmt_pos (pmt.GetPosition(0), pmt.GetPosition(1), pmt.GetPosition(2)); 
		if (pmt.GetCylLoc() == 0){ //Get all PMTs on the TOP CAP
			cout << "TOP CAP PMTS HAVE COORDINATES (" << pmt_pos.X() << ", " << pmt_pos.Y() << ", " << pmt_pos.Z() << ")" << endl;
			
		}
		if (pmt.GetCylLoc() == 1){ //Get all PMTs on the WALLS
			cout << "WALL PMTS HAVE COORDINATES (" << pmt_pos.X() << ", " << pmt_pos.Y() << ", " << pmt_pos.Z() << ")" << endl;
			
		}
		if (pmt.GetCylLoc() == 2){ //Get all PMTs on the BOTTOM CAP
			cout << "BOTTOM CAP PMTS HAVE COORDINATES (" << pmt_pos.X() << ", " << pmt_pos.Y() << ", " << pmt_pos.Z() << ")" << endl;
			
		}
	}//end loop over pmts

	WCSimRootTrigger* trig; 

	return 0;
}
