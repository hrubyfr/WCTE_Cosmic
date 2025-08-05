

#include <TH1.h>
#include <cmath>
#include <iostream>
#include <stdio.h>     
#include <stdlib.h>
#include <string>

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
	TString help = filename;
	TString plots_folder = "plots/" + help.Remove(0, 3);
	if (gSystem->AccessPathName(plots_folder)) gSystem->Exec("mkdir -p " + plots_folder);
	

	TFile* file = new TFile(filename, "read");

	if (!file || file->IsOpen() == false) {
		cerr << "Error: could not open file " << filename << endl;
	}
	file->ls();

	gSystem->cd(plots_folder);
	
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
	TH1D* muon_energies = new TH1D("muon_energies", "Energy of muons; Energy[GeV]; count", 1000, 0, 200); 
	TH1D* creator_processes = new TH1D("creator_processes", "Which creator process created an electron?;ID;count", 1000, 0, 1000); 
	TH1D *h_hit_pmts = new TH1D("hit_pmts", "number of hit PMTs in low energy events;hit pmts;count", 50, 0, 1000);
	TH1D *h_hit_pmts_2 = new TH1D("hit_pmts_2", "number of hit PMTs in events with higher energy;hit pmts;count", 50, 0, 1000);
	cout << "File has " << tree->GetEntries() << " entries" << endl;


	TCanvas* vis_canvas = new TCanvas("", "", 1800, 900);
	for(long i = 0; i < tree->GetEntries(); i++){
		tree->GetEntry(i);
		if (!wcsimrootsuperevent) cout << "WCSIM super event is null" << endl;
		
		trig = wcsimrootsuperevent->GetTrigger(0);

		int ntracks = trig->GetNtrack_slots();
		// cout << "There are " << ntracks << " tracks in event " << i << endl;
		
		double event_energy;
		//TPolyLine3D* vis_track = new TPolyLine3D(ntracks + 1);
		//vis_track->SetLineColor(kBlue);
		for(int itrack = 0; itrack < ntracks; itrack++){
			WCSimRootTrack* track = dynamic_cast<WCSimRootTrack*> (wcsimrootsuperevent->GetTrigger(0)->GetTracks()->At(itrack));
			double particle_E = track->GetE() / 1000;	
			int p_type = track->GetIpnu();
			int parent = track->GetParenttype();
			// cout << "Track " << itrack << " is a particle type " << p_type << " with parent particle " << track->GetParenttype() << endl;
			//if (p_type != 13) cout << "WARNING, Different particle than muon in track!" << endl;
			particle_types->Fill(p_type);
			parent_types->Fill(parent);
			if(p_type == 11) electron_energies->Fill(particle_E);
			if(p_type == 13 && itrack == 0) {muon_energies->Fill(particle_E); event_energy = particle_E;}
			/* if(p_type == 13) cout << "Start and stop coordinates of muon in simulation: Start - (" << track->GetStart(0) << ", " << track->GetStart(1) << ", " << track->GetStart(2) << ") Stop - (" <<track->GetStop(0) << ", " << track->GetStop(1) << ", " << track->GetStop(2) << "), With energy " << track->GetE()/1000 << " GeV" << endl;  
			if(p_type == 11 && "muMinusCaptureAtRest" == track->GetCreatorProcessName()) cout << "Start and stop coordinates of electron from muon capture at rest: Start - (" << track->GetStart(0) << ", " << track->GetStart(1) << ", " << track->GetStart(2) << ") Stop - (" <<track->GetStop(0) << ", " << track->GetStop(1) << ", " << track->GetStop(2) << "), With energy " << track->GetE()/1000 << " GeV" << endl;   */
			double track_x = track->GetStart(0) - track->GetStop(0), track_y = track->GetStart(1) - track->GetStop(1), track_z = track->GetStart(2) - track->GetStop(2);
			double track_distance = std::sqrt(track_x * track_x + track_y * track_y + track_z * track_z);
			double track_dir_x = track->GetDir(0);
			double track_dir_y = track->GetDir(1);
			double track_dir_z = track->GetDir(2);
			/* 
			if(p_type == 13 && itrack > 1) cout << "Direction of the muon track: (" << track_dir_x << ", " << track_dir_y << ", " << track_dir_z << ") \t Muon track distance: " << track_distance << endl;  
			if(p_type == 11 && "muMinusCaptureAtRest" == track->GetCreatorProcessName() && track->GetE() > 1) cout << "Direction of the electron: (" << track_dir_x << ", " << track_dir_y << ", " << track_dir_z << ")" << endl; */
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
		TString hits_folder = "events_pmt_hits/";
		gSystem->Exec("mkdir -p " + hits_folder);
		
		TString h_name = "WCSim cosmic muon event " + std::to_string(i) + " Energy " + std::to_string(event_energy) + " GeV" + ";pmt_ID;charge";
		TH1D *h_event_charge = new TH1D("h_event_charge", h_name, nPMT, 0, nPMT);
		int n_hit_pmts = 0; 
		for(int j = 0; j < trig-> GetNcherenkovdigihits(); j++){

			WCSimRootCherenkovDigiHit* hit = (WCSimRootCherenkovDigiHit*) trig->GetCherenkovDigiHits()->At(j);
			double charge = hit->GetQ();
			int tubeID = hit->GetTubeId();
			pmt_charge->Fill(tubeID, charge);

			h_event_charge->Fill(tubeID, charge);

			WCSimRootPMT pmt =  geo->GetPMT(tubeID - 1);
			int location = pmt.GetCylLoc();
			//cout << "PMT " << tubeID - 1 << " is at X = " << pmt.GetPosition(0) << " Y = " << pmt.GetPosition(1) << " Z = " << pmt.GetPosition(2) << endl;

			pmt_type_charge->Fill(location, charge);

			mPMT_charge->Fill(pmt.GetmPMTNo(), charge);

			hit->GetmPMTId();
			if (charge > 5) n_hit_pmts++;


		}//loop over subevents end
		 //cout << "Superevent has " << wcsimrootsuperevent->GetNumberOfEvents() << " events" << endl;
		if (event_energy < 1.5) h_hit_pmts->Fill(n_hit_pmts);
		else h_hit_pmts_2->Fill(n_hit_pmts);

		if (gSystem->AccessPathName(hits_folder)){
		TCanvas* c_event_charge = new TCanvas("", "", 1800, 900);
		h_event_charge->Draw("hist");
		c_event_charge->Print(hits_folder + Form("hit_%i.png", i));
		delete c_event_charge;
		delete h_event_charge;
		}

	}//end of loop over events
	TCanvas* can = new TCanvas("", "", 1800, 900);
	pmt_charge->Draw("hist");
	can->Print("pmt_charge.pdf");
	can->Print("pmt_charge.png");

	pmt_type_charge->Draw("hist");
	can->Print("pmt_type_charge.pdf");

	mPMT_charge->Draw("hist");
	can->Print("mpmt_charge.pdf");

	particle_types->GetXaxis()->SetRangeUser(0, 30);
	particle_types->Draw();
	can->Print("particle_types.pdf");

	parent_types->Draw();
	can->Print("parent_types.pdf");

	muon_energies->Draw();
	can->Print("muon_energies.pdf");
	h_hit_pmts_2->Scale(1/h_hit_pmts_2->Integral());

	h_hit_pmts->Scale(1/h_hit_pmts->Integral());
	h_hit_pmts->Draw("hist");
	h_hit_pmts_2->SetLineColor(kRed);
	h_hit_pmts_2->Draw("hist same");
	can->Print("h_hit_PMTS.png");
	can->Print("h_hit_PMTS.pdf");

	electron_energies->GetXaxis()->SetRangeUser(0, 2);
	electron_energies->GetYaxis()->SetRangeUser(0, 50);
	electron_energies->Draw();
	can->Print("electron_energies.pdf");
	delete can;


	file->Close();
} 
