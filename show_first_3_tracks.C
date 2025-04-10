

#include <Math/Vector3Dfwd.h>
#include <TH1.h>
#include <cmath>
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
#include "TPolyLine3D.h"
#include "TH3D.h"
#include "Math/Vector3D.h"

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using ROOT::Math::XYZVector;

void show_first_3_tracks(const char* filename = "../wcsim.root"){
	TFile* file = new TFile(filename, "read");

	if (!file || file->IsOpen() == false) {
		cerr << "Error: could not open file " << filename << endl;
	}
	file->ls();
	
	TTree* tree = (TTree*)file->Get("wcsimT");
	tree->SetBranchStatus("wcsimrootevent2", 0);
	tree->SetBranchStatus("wcsimrootevent_OD", 0);
	tree->SetBranchStatus("wcsimrootevent", 1);

	cout << "There are " << tree->GetEntries() << " events" << endl; 
	
	TTree* geotree = (TTree*) file->Get("wcsimGeoT");
	WCSimRootGeom* geo = 0;
	geotree->SetBranchAddress("wcsimrootgeom", &geo);
	geotree->GetEntry(0);
	
	// Go over the PMTs, get maximum and minimum coordinates to get barrel size
	int nPMTs = geo->GetWCNumPMT();

	cout << "Geometry has " << nPMTs << " PMTs" << endl;
	double maximum_x = 0, maximum_y = 0, maximum_z = 0;
	double minimum_x = 0, minimum_y = 0, minimum_z = 0;

	for (int iPMT = 0; iPMT < nPMTs; iPMT++){
		WCSimRootPMT pmt = geo->GetPMT(iPMT);
		double pmt_x = pmt.GetPosition(0);
		double pmt_z = pmt.GetPosition(1);
		double pmt_y = pmt.GetPosition(2);

		if(pmt_x > maximum_x)maximum_x = pmt_x;
		if(pmt_y > maximum_y)maximum_y = pmt_y;
		if(pmt_z > maximum_z)maximum_z = pmt_z;
		if(pmt_x < minimum_x)minimum_x = pmt_x;
		if(pmt_y < minimum_y)minimum_y = pmt_y;
		if(pmt_z < minimum_z)minimum_z = pmt_z;
	}//end loop over PMTs
	maximum_x += 20;
	maximum_y += 20;
	maximum_z += 20;		//set approximate WCTE barrel dimensions
	
	minimum_x -= 20;
	minimum_y -= 20;
	minimum_z -= 20;
	
	double WCTE_radius;
	if(maximum_x > maximum_y) WCTE_radius = maximum_x;
	else double WCTE_radius = maximum_y;
	double WCTE_height = maximum_z;

	cout << "WCTE radius is " << WCTE_radius << " cm and WCTE height is " << WCTE_height << " cm" << endl;

	WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent(); //Initialising the wcsimsuperevent - the event started in WCSim
	TBranch* branch = tree->GetBranch("wcsimrootevent");
	branch->SetAddress(&wcsimrootsuperevent);

	tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

	WCSimRootTrigger* trig;





	TCanvas* vis_canvas = new TCanvas("", "", 1800, 900);
	for(long i = 0; i < tree->GetEntries(); i++){
		cout << "Event " << i << " starting " << endl;
		tree->GetEntry(i);
		if (!wcsimrootsuperevent) cout << "WCSIM super event is null" << endl;

		trig = wcsimrootsuperevent->GetTrigger(0);

		int ntracks = trig->GetNtrack_slots();
		cout << "There are " << ntracks << " tracks in event " << i << endl;
		
		for(int itrack = 0; itrack < 3; itrack++){

			WCSimRootTrack* track = dynamic_cast<WCSimRootTrack*> (trig->GetTracks()->At(itrack));
			double particle_E = track->GetE() / 1000;	
			int p_type = track->GetIpnu();
			//int parent = track->GetParenttype();
			XYZVector track_start;
			XYZVector track_dir;
			XYZVector track_end;
			XYZVector vtx_pos;
			vector<double> b_times;
			vector<float> b_energy;
			track_start.SetCoordinates(track->GetStart(0), track->GetStart(2), track->GetStart(1));
			track_dir.SetCoordinates(track->GetDir(0), track->GetDir(2), track->GetDir(1));
			track_end.SetCoordinates(track->GetStop(0), track->GetStop(2), track->GetStop(1));
			vtx_pos.SetCoordinates(trig->GetVtx(0), trig->GetVtx(2), trig->GetVtx(1));
			cout << endl;
			cout << "The particle type is " << p_type << endl; 
			cout << "The particle track starts at (" << track_start.X() << ", " << track_start.Y() << ", " << track_start.Z() << ")" << endl;
			cout << "The vtx of event starts at (" << vtx_pos.X() << ", " << vtx_pos.Y() << ", " << vtx_pos.Z() << ")" << endl;
			cout << "Direction of the particle is (" << track_dir.X() << ", " << track_dir.Y() << ", " << track_dir.Z() << ")" << endl;
			cout << "Theta of the particle is " << track_dir.Theta() << endl;
			cout << "The particle track stops at (" << track_end.X() << ", " << track_end.Y() << ", " << track_end.Z() << ")" << endl;
			cout << endl;



		}//end loop over tracks




	}// end loop over events
	
	file->Close();
} // end file
