

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

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"


using std::endl;
using std::cerr;
using std::cout;


void analyse_muon_tracks(const char* filename = "../wcsim.root"){
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
		double pmt_y = pmt.GetPosition(1);
		double pmt_z = pmt.GetPosition(2);

		if(pmt_x > maximum_x){
			maximum_x = pmt_x;
			if(pmt_y > maximum_y){
				maximum_y = pmt_y;
				if(pmt_z > maximum_z){
					maximum_z = pmt_z;
				}
			}
		}

		if(pmt_x < minimum_x){
			minimum_x = pmt_x;
			if(pmt_y < minimum_y){
				minimum_y = pmt_y;
				if(pmt_z < minimum_z){
					minimum_z = pmt_z;
				}
			}
		}
	}//end loop over PMTs
	maximum_x += 20;
	maximum_y += 20;
	maximum_z += 20;
	
	minimum_x -= 20;
	minimum_y -= 20;
	minimum_z -= 20;

	double WCTE_radius = sqrt(maximum_x * maximum_x + maximum_y* maximum_y);
	double WCTE_height = 2 * maximum_z;

	cout << "WCTE radius is " << WCTE_radius << " cm and WCTE height is " << WCTE_height << " cm" << endl;

	WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent(); //Initialising the wcsimsuperevent - the event started in WCSim
	TBranch* branch = tree->GetBranch("wcsimrootevent");
	branch->SetAddress(&wcsimrootsuperevent);

	tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

	WCSimRootTrigger* trig;



	TH3D* visualize_tracks = new TH3D("visualize_tracks", "Event visualization;x;y;z", 50, -250, 250, 50, -250, 250, 50, -250, 250);  


	TCanvas* vis_canvas = new TCanvas("", "", 1800, 900);
	for(long i = 0; i < 1; i++){
		tree->GetEntry(i);
		if (!wcsimrootsuperevent) cout << "WCSIM super event is null" << endl;

		trig = wcsimrootsuperevent->GetTrigger(0);

		int ntracks = trig->GetNtrack_slots();
		cout << "There are " << ntracks << " tracks in event " << i << endl;

		TPolyLine3D* vis_track = new TPolyLine3D(ntracks + 1);
		vis_track->SetLineColor(kBlue);
		for(int itrack = 0; itrack < ntracks; itrack++){
			WCSimRootTrack* track = dynamic_cast<WCSimRootTrack*> (wcsimrootsuperevent->GetTrigger(0)->GetTracks()->At(itrack));
			double particle_E = track->GetE() / 1000;	
			int p_type = track->GetIpnu();
			int parent = track->GetParenttype();
			bool first_event = true;
			if(p_type == 13) {
				if(first_event){
					double_t start_x = track->GetStart(0);
					double_t start_y = track->GetStart(1);
					double_t start_z = track->GetStart(2);
					vis_track->SetPoint(itrack, start_x, start_y, start_z);
					first_event = false;
				} //end first event

				double_t end_x = track->GetStop(0);
				double_t end_y = track->GetStop(1);
				double_t end_z = track->GetStop(2);
				vis_track->SetPoint(itrack + 1, end_x, end_y, end_z);

				cout << "Muon energy (track " << itrack << ") is " << particle_E << " GeV" << endl;

			}// end if condition



		}//end loop over tracks

		visualize_tracks->Fill(1, 1, 1);
		visualize_tracks->Draw("BOX");
		vis_track->Draw();



		file->Close();
	}// end loop over events
} // end file
