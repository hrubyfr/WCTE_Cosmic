

#include <Math/Vector3Dfwd.h>
#include <TH1.h>
#include <TH3.h>
#include <TMath.h>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <stdio.h>     
#include <stdlib.h>
#include <system_error>

#include "TLegend.h"
#include "TTree.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "TH2D.h"
#include "Math/Vector3D.h"

#include "WCSimRootEvent.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootOptions.hh"

using std::vector;
using std::endl;
using std::cerr;
using std::cout;
using ROOT::Math::XYZVector;

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
	
	long nentries = tree->GetEntries();
	cout << "There are " << nentries << " events" << endl; 
	
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
	
	double WCTE_radius = 0;
	if(maximum_x >= maximum_y) WCTE_radius = maximum_x;
	else double WCTE_radius = maximum_y;
	double WCTE_height = maximum_z;

	cout << "WCTE radius is " << WCTE_radius << " cm and WCTE height is " << WCTE_height << " cm" << endl;

	WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent(); //Initialising the wcsimsuperevent - the event started in WCSim
	TBranch* branch = tree->GetBranch("wcsimrootevent");
	branch->SetAddress(&wcsimrootsuperevent);

	tree->GetBranch("wcsimrootevent")->SetAutoDelete(kTRUE);

	WCSimRootTrigger* trig;

	// Initialise histograms
	TH1D* h_theta = new TH1D("h_theta", "Entry angle; #theta [rad]; counts", 90, 0, M_PI_2);
	TH1D* h_cos_theta = new TH1D("h_cos_theta", "Entry angle; cos(#theta); counts", 90, 0, 1);
	TH1D* h_cos_theta2 = new TH1D("h_cos_theta2", "cos_theta; cos(#theta); A * cos^{2}(#theta)", 91, 0, 1);

	TH2D* h2_angles = new TH2D("h2-angles", "Incident particle angle; #phi [rad]; #theta [rad]", 180, -M_PI, M_PI, 90, 0, M_PI/2);
	TH1D* entry_energies = new TH1D("entry_energies", "Entry energies; E[GeV]; counts", 1000, 0, 300);
	TH1D* exit_energies = new TH1D("exit_energies", "Exit energies; E[GeV]; counts", 1000, 0, 300);
	TH1D* h_energy = new TH1D("h_energy", "Energy lost in detector; E_{in} - E_{out}[GeV]; counts", 100, 0, 5);

	TH2D* top_cap_entries = new TH2D("top_cap_entries", "Muons entries through the top cap of the detector; x[cm]; y[cm]", 100, minimum_x, maximum_x, 100, minimum_y, maximum_y);
	TH2D* bottom_cap_exits = new TH2D("bottom_cap_exits", "Muons exits through the bottom cap of the detector; x[cm]; y[cm]", 100, minimum_x, maximum_x, 100, minimum_y, maximum_y);
	TH2D* wall_entries = new TH2D("wall_entries", "Muons entries through the wall of the detector; #phi [rad]; z[cm]", 100, -M_PI, M_PI, 
			100, -WCTE_height, WCTE_height);
	TH2D* wall_exits = new TH2D("wall_exits", "Muons exits through the wall of the detector; #phi [rad]; z[cm]", 100, -M_PI, M_PI, 100, -WCTE_height, WCTE_height);
	TH1D* nmuon_tracks = new TH1D("nmuon_tracks", "number of muon tracks per event", 10, 0, 10);
	//TCanvas* vis_canvas = new TCanvas("", "", 1800, 900);
	for(long i = 0; i < nentries; i++){
		tree->GetEntry(i);
		if (!wcsimrootsuperevent) cout << "WCSIM super event is null" << endl;

		cout << endl;
		cout << "################### EVENT " << i << " ############################" << endl;
		cout << endl;
		trig = wcsimrootsuperevent->GetTrigger(0);

		int ntracks = trig->GetNtrack_slots();
		cout << "There are " << ntracks << " tracks in event " << i << endl;

		int n_muons = 0;
		for(int itrack = 2; itrack < ntracks; itrack++){ //Start looping over tracks from 2 - first two tracks are the parent particle producing the mion, and the second is propagating the muon to the WCSim world
			WCSimRootTrack* track = dynamic_cast<WCSimRootTrack*> (wcsimrootsuperevent->GetTrigger(0)->GetTracks()->At(itrack));
			double particle_E = track->GetE() / 1000;	
			int p_type = track->GetIpnu();
			//int parent = track->GetParenttype();
			XYZVector track_start;
			XYZVector track_dir;
			XYZVector track_end;
			XYZVector vtx_pos;
			vector<double> b_times;
			vector<float> b_energy;
			vector<int> b_type;
			if(p_type == 13) {
				bool boundary_check = false;
				track_start.SetCoordinates(track->GetStart(0), track->GetStart(2), track->GetStart(1));
				track_dir.SetCoordinates(track->GetDir(0), track->GetDir(2), track->GetDir(1));
				track_end.SetCoordinates(track->GetStop(0), track->GetStop(2), track->GetStop(1));
				vtx_pos.SetCoordinates(trig->GetVtx(0), trig->GetVtx(2), trig->GetVtx(1));
				b_type = track->GetBoundaryTypes();
				b_times = track->GetBoundaryTimes();
				b_energy = track->GetBoundaryKEs();

				h2_angles->Fill(track_dir.Phi(), -(track_dir.Theta() - M_PI));
				h_theta->Fill(-(track_dir.Theta() - M_PI));
				h_cos_theta->Fill(cos(-(track_dir.Theta() - M_PI)));
				cout << "Angles: " << track_dir.Phi() << " phi, " << - (track_dir.Theta() - M_PI) << " theta" << endl;
				for(unsigned int iKE = 0; iKE < b_times.size(); iKE++){
					cout << "Boundary " << iKE << " of type " << b_type[iKE] << " was crossed at " << b_times[iKE] <<" ns, has kinetic energy " << b_energy[iKE] << endl;					
					cout << "Difference in times with the following boundary crossing is " << b_times[iKE] - b_times[iKE+1] << endl;
					cout << "Following boundary crossing has an energy " << b_energy[iKE+1] << endl;
					if ((b_times[iKE+1] - b_times[iKE]) > 1.0){
						boundary_check = true;
						entry_energies->Fill(b_energy[iKE]/1000);
						exit_energies->Fill(b_energy[iKE+1]/1000);
						h_energy->Fill((b_energy[iKE] - b_energy[iKE+1])/1000.0);
						break;
					}//end if condition checking boundary entry / exit
				}//end loop over boundaries
				if(!boundary_check){
					cout << "##########################################################" << endl;
					cout << "EVENT " << i << " DIDNT GET A PROPER BOUNDARY ENERGY CHECK"  << endl; 
					cout << "##########################################################" << endl;
				}	
				 //cout << "There are " << b_times.size() << " boundary times" << endl;
				 //cout << "The muon starts at (" << track_start.X() << ", " << track_start.Y() << ", " << track_start.Z() << ")" << endl;
				 //cout << "The vtx starts at (" << vtx_pos.X() << ", " << vtx_pos.Y() << ", " << vtx_pos.Z() << ")" << endl;

				 //cout << "Direction of the muon is (" << track_dir.X() << ", " << track_dir.Y() << ", " << track_dir.Z() << ")" << endl;
				 //cout << "Theta of the muon is " << track_dir.Theta() << endl;
				 //cout << "Energy of the muon is " << particle_E << " GeV" << endl;
				 cout << "The muon stops at (" << track_end.X() << ", " << track_end.Y() << ", " << track_end.Z() << ")" << endl;
				if (sqrt(track_end.X() * track_end.X() + track_end.Y() * track_end.Y()) < WCTE_radius && track_end.Z() >= -WCTE_height){
					cout << "TRACK END INSIDE BARREL" << endl;
					return;
				}
				XYZVector muon_pos = track_start;
				while (!(muon_pos.Z() <= maximum_z)){
					muon_pos.SetCoordinates(muon_pos.X() + track_dir.X(),
							muon_pos.Y() + track_dir.Y(),
							muon_pos.Z() + track_dir.Z());

				}//Propagate muon track position through the while loop, until it intersects the WCTE barrel top
				if(sqrt(muon_pos.X() * muon_pos.X() + muon_pos.Y() * muon_pos.Y()) < WCTE_radius){
					top_cap_entries->Fill(muon_pos.X(), muon_pos.Y());
				}
				else {
					while(!(sqrt(muon_pos.X() * muon_pos.X() + muon_pos.Y() * muon_pos.Y()) < WCTE_radius)){
						muon_pos.SetCoordinates(muon_pos.X() + track_dir.X(),
								muon_pos.Y() + track_dir.Y(),
								muon_pos.Z() + track_dir.Z());
					}
					wall_entries->Fill(muon_pos.Phi(), muon_pos.Z());
				}//end else


				bool exit_through_bottom = false;
				while (sqrt(muon_pos.X() * muon_pos.X() + muon_pos.Y() * muon_pos.Y()) < WCTE_radius){	
					muon_pos.SetCoordinates(muon_pos.X() + track_dir.X(),
							muon_pos.Y() + track_dir.Y(),
							muon_pos.Z() + track_dir.Z());
					if(muon_pos.Z() < -WCTE_height){
						bottom_cap_exits->Fill(muon_pos.X(), muon_pos.Y());
						exit_through_bottom = true;
						break;
					}
				}//Propagate muon track position through the while loop, until it intersects the WCTE barrel sides or bottom
				if (!exit_through_bottom){
					wall_exits->Fill(muon_pos.Phi(), muon_pos.Z());
				}
				
				n_muons++; //increment number of muon tracks in the event
			}// end if condition



		}//end loop over tracks
		nmuon_tracks->Fill(n_muons);



	}// end loop over events
	TCanvas* can = new TCanvas("", "", 1800, 900);
	nmuon_tracks->Draw("hist");
	delete can;

	TCanvas* ee_can = new TCanvas("", "", 1800, 900);
	ee_can->Divide(2, 2);
	ee_can->cd(1);
	top_cap_entries->Draw("colz");
	ee_can->cd(2);
	wall_entries->Draw("colz");
	ee_can->cd(3);
	bottom_cap_exits->Draw("colz");
	ee_can->cd(4);
	wall_exits->Draw("colz");
	ee_can->Print("plots/ee_can.png");
	ee_can->Print("plots/ee_can.pdf");
	delete ee_can;

	TCanvas* c_en_ex = new TCanvas("", "", 1800, 900);
	c_en_ex->Divide(2);
	c_en_ex->cd(1);
	entry_energies->Draw("hist");
	c_en_ex->cd(2);
	exit_energies->Draw("hist");
	c_en_ex->Print("plots/c_en_ex.png");
	c_en_ex->Print("plots/c_en_ex.pdf");
	delete c_en_ex;
	
	TCanvas* c_angles = new TCanvas("", "", 1800, 900);
	h2_angles->Draw("colz");
	c_angles->Print("plots/c_angles.png");
	c_angles->Print("plots/c_angles.pdf");
	delete c_angles;

	TCanvas* c_energy = new TCanvas("", "", 1800, 900);
	h_energy->Draw("hist");
	c_energy->Print("plots/c_energy.png");
	c_energy->Print("plots/c_energy.pdf");
	delete c_energy;

	TCanvas* c_theta = new TCanvas("", "", 1800, 900);
	h_theta->Draw("hist");
	c_theta->Print("plots/c_theta.png");
	c_theta->Print("plots/c_theta.pdf");
	delete c_theta;
	
	TCanvas* c_cos_theta = new TCanvas("", "", 1800, 900);
	gStyle->SetOptStat(0);
	
	h_cos_theta->Draw("hist");

	for(double c_t = 0; c_t <= 1.0; c_t += (1.0/90)){
		h_cos_theta2->Fill(c_t, c_t * c_t);	
	}
	h_cos_theta2->SetLineColor(kRed);
	h_cos_theta2->Scale(h_cos_theta->GetMaximum());
	h_cos_theta2->Draw("hist same");
	TLegend* leg = new TLegend(0.1, 0.8, 0.3, 0.9);
	leg->AddEntry(h_cos_theta, "Data entries", "l");
	leg->AddEntry(h_cos_theta2, "Acos^{2}(#theta)", "l");
	leg->Draw("same");
	c_cos_theta->Print("plots/c_cos_theta.png");
	c_cos_theta->Print("plots/c_cos_theta.pdf");

	delete c_cos_theta;
	//file->Close();
} // end file
