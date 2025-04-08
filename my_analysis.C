

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


void my_analysis(const char* filename = "../wcsim.root"){
	TFile* file = new TFile(filename, "read");

	if (!file || file->IsOpen() == false) {
		cerr << "Error: could not open file " << filename << endl;
	}
	file->ls();
	
	TTree* tree = (TTree*)file->Get("wcsimT");
	tree->SetBranchStatus("wcsimrootevent2", 0);
	tree->SetBranchStatus("wcsimrootevent_OD", 0);
	tree->SetBranchStatus("wcsimrootevent", 1);
	tree->Print();

	cout << "There are " << tree->GetEntries() << " events" << endl; 
	
	TTree* geotree = (TTree*) file->Get("wcsimGeoT");
	WCSimRootGeom* geo = 0;
	geotree->SetBranchAddress("wcsimrootgeom", &geo);
	geotree->GetEntry(0);
	
	




	WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent(); //Initialising the wcsimsuperevent - the event started in WCSim
	TBranch* branch = tree->GetBranch("wcsimrootevent");
	branch->SetAddress(&wcsimrootsuperevent);

	tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);
	
	WCSimRootTrigger* trig;

	geotree->Print();

	TTree *opttree = (TTree*)file->Get("wcsimRootOptionsT");
	WCSimRootOptions *opt = 0; 
	opttree->SetBranchAddress("wcsimrootoptions", &opt);
	cout << "Options tree has: " << opttree->GetEntries() << " entries (1 expected)" << endl;
	if (opttree->GetEntries() == 0) {
		exit(9);
	}
	opttree->GetEntry(0);
	opt->Print();
	
	cout << "WCTE barrel has radius " << geo->GetWCCylRadius() << " and height " << 2 * geo->GetWCCylLength() << endl;

	auto nPMT = geo->GetWCNumPMT();
	cout << "Detector has " << nPMT << " pmts" << endl;
	TH1D* pmt_charge = new TH1D("pmt_charge", "total charge collected by PMTs by ID number; pmt_id; charge", nPMT, 0, nPMT); 
	TH1D* pmt_type_charge = new TH1D("pmt_type_charge", "0 = top cap, 1 = wall, 2 = bottom cap; position; charge", 3, 0, 3);
	TH1D* mPMT_charge = new TH1D("mPMT_charge", "Total charge collected by an mpmt; mpmt_id; charge", nPMT/19, 0, nPMT/19 );
	//TH3D* visualize_tracks = new TH3D("visualize_tracks", "Event visualization;x;y;z", 50, -250, 250, 50, -250, 250, 50, -250, 250);  
	TH1D* particle_types = new TH1D("particle types", "Types of particles in simulation; MC PID# ; counts", 1000, 0, 1000);
	TH1D* parent_types = new TH1D("parent types", "Types of particles parents in simulation; MC PID# ; counts", 1000, 0, 1000);
	TH1D* electron_energies = new TH1D("electron_energies", "Energy of electrons; Energy[GeV]; count", 500, 0, 2); 
	TH1D* muon_energies = new TH1D("moun_energies", "Energy of muons; Energy[GeV]; count", 1000, 0, 200); 
	cout << "File has " << tree->GetEntries() << " entries" << endl;


	TCanvas* vis_canvas = new TCanvas("", "", 1800, 900);
	for(long i = 0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		if (!wcsimrootsuperevent) cout << "WCSIM super event is null" << endl;
		
		trig = wcsimrootsuperevent->GetTrigger(0);

		int ntracks = trig->GetNtrack_slots();
		cout << "There are " << ntracks << " tracks in event " << i << endl;

		//TPolyLine3D* vis_track = new TPolyLine3D(ntracks + 1);
		//vis_track->SetLineColor(kBlue);
		for(int itrack = 0; itrack < ntracks; itrack++){
			WCSimRootTrack* track = dynamic_cast<WCSimRootTrack*> (wcsimrootsuperevent->GetTrigger(0)->GetTracks()->At(itrack));
			double particle_E = track->GetE() / 1000;	
			int p_type = track->GetIpnu();
			int parent = track->GetParenttype();
			cout << "Track " << itrack << " is a particle type " << p_type << " with parent particle " << track->GetParenttype() << endl;
			//if (p_type != 13) cout << "WARNING, Different particle than muon in track!" << endl;
			particle_types->Fill(p_type);
			parent_types->Fill(parent);
			if(p_type == 11) electron_energies->Fill(particle_E);
			if(p_type == 13 && itrack == 0) muon_energies->Fill(particle_E);
			//double_t start_x = track->GetStart(0);
			//double_t start_y = track->GetStart(1);
			//double_t start_z = track->GetStart(2);
			//vis_track->SetPoint(itrack, start_x, start_y, start_z);
			

			//double_t end_x = track->GetStop(0);
			//double_t end_y = track->GetStop(1);
			//double_t end_z = track->GetStop(2);
			
			//cout << "Track point start is: (" << start_x << ", " << start_y << ", " << start_z << ") and end point is: (" << end_x << ", " << end_y << ", " << end_z << ")" << endl;
			//vis_track->SetPoint(itrack + 1, end_x, end_y, end_z);

		}//end loop over tracks
		//cout << "fill test TH3D" << endl;
		//visualize_tracks->Fill(1, 1, 1);
		//visualize_tracks->Draw("BOX");
		//vis_track->Draw();
		//cout << "Box has been drawn" << endl;

		for(int j = 0; j < trig-> GetNcherenkovdigihits(); j++){

			WCSimRootCherenkovDigiHit* hit = (WCSimRootCherenkovDigiHit*) trig->GetCherenkovDigiHits()->At(j);
			double charge = hit->GetQ();
			int tubeID = hit->GetTubeId();
			pmt_charge->Fill(tubeID, charge);
			
			WCSimRootPMT pmt =  geo->GetPMT(tubeID - 1);
			int location = pmt.GetCylLoc();
			//cout << "PMT " << tubeID - 1 << " is at X = " << pmt.GetPosition(0) << " Y = " << pmt.GetPosition(1) << " Z = " << pmt.GetPosition(2) << endl;
		
			pmt_type_charge->Fill(location, charge);
			
			mPMT_charge->Fill(pmt.GetmPMTNo(), charge);
			
			hit->GetmPMTId();


		}//loop over subevents end
		//cout << "Superevent has " << wcsimrootsuperevent->GetNumberOfEvents() << " events" << endl;


	}//end of loop over events
	TCanvas* can = new TCanvas("", "", 1800, 900);
	pmt_charge->Draw("hist");
	can->Print("plots/pmt_charge.pdf");
	can->Print("plots/pmt_charge.png");

	pmt_type_charge->Draw("hist");
	can->Print("plots/pmt_type_charge.pdf");

	mPMT_charge->Draw("hist");
	can->Print("plots/mpmt_charge.pdf");

	particle_types->GetXaxis()->SetRangeUser(0, 30);
	particle_types->Draw();
	can->Print("plots/particle_types.pdf");
	
	parent_types->Draw();
	can->Print("plots/parent_types.pdf");

	muon_energies->Draw();
	can->Print("plots/muon_energies.pdf");

	electron_energies->GetXaxis()->SetRangeUser(0, 0.1);
	electron_energies->Draw();
	can->Print("plots/electron_energies.pdf");
	delete can;


	file->Close();
} 
