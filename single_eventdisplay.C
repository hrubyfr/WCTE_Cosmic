#include <Rtypes.h>
#include <iostream>
#include <ostream>
#include <stdio.h>     
#include <stdlib.h>

#include <TChain.h>
#include "TTree.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include <TH2D.h>
#include <TMarker.h>

#include "WCSimRootOptions.hh"
#include "WCSimRootGeom.hh"
#include "WCSimRootEvent.hh"


void SetDarkStyle(){
	gStyle->SetPalette(1);

	gStyle->SetCanvasColor(kBlack);
	gStyle->SetPadColor(kBlack);

	gStyle->SetLegendFillColor(kBlack);
	gStyle->SetStatColor(kBlack);
	gStyle->SetStatTextColor(kWhite);

	gStyle->SetTitleTextColor(kWhite);
	gStyle->SetTitleColor(kWhite);

	gStyle->SetLabelColor(kWhite);

	gStyle->SetAxisColor(kWhite);
	gStyle->SetFrameLineColor(kWhite);
	gStyle->SetGridColor(kWhite);

	gStyle->SetPadTickX(1);
	gStyle->SetPadTickY(1);
}

void MakeWhiteAxis(TH1D *h1){
	h1->GetXaxis()->SetAxisColor(kWhite);
	h1->GetXaxis()->SetLabelColor(kWhite);
	h1->GetXaxis()->SetTitleColor(kWhite);

	h1->GetYaxis()->SetAxisColor(kWhite);
	h1->GetYaxis()->SetAxisColor(kWhite);
	h1->GetYaxis()->SetLabelColor(kWhite);

}
void MakeWhiteAxes(TH2D *h2){
	h2->GetXaxis()->SetAxisColor(kWhite);
	h2->GetXaxis()->SetLabelColor(kWhite);
	h2->GetXaxis()->SetTitleColor(kWhite);

	h2->GetYaxis()->SetAxisColor(kWhite);
	h2->GetYaxis()->SetAxisColor(kWhite);
	h2->GetYaxis()->SetLabelColor(kWhite);
	
	h2->GetZaxis()->SetTitleColor(kWhite);
	h2->GetZaxis()->SetLabelColor(kWhite);
	h2->GetZaxis()->SetTitleColor(kWhite);
}
int single_eventdisplay(const char * fname = "../wcsim.root", const int evt_nr = 0)
{
	const int susp_event = evt_nr;
	gStyle->SetOptStat(0);

	TString help = fname;
	TString plots_folder = "plots/" + help.Remove(0, 3);
	if (gSystem->AccessPathName(plots_folder)) gSystem->Exec("mkdir -p " + plots_folder);
	TChain *t = new TChain("wcsimT");
	t->Add(fname);
	std::string single_file_name = t->GetFile()->GetName();
	// Get the first file for geometry
	TFile *f = TFile::Open(single_file_name.c_str());
	if (!f->IsOpen()){
		std::cout << "Error, could not open input file: " << single_file_name << std::endl;
		return -1;
	}

	gSystem->cd(plots_folder);

	std::string prefix = fname;
	if (prefix.find_last_of("/")!=std::string::npos) 
	{   
		prefix = prefix.substr(prefix.find_last_of("/")+1);
	}
	if (prefix.find_last_of("[")!=std::string::npos) 
	{   
		prefix = prefix.substr(0,prefix.find_last_of("["));
	}
	if (prefix.find_last_of("*")!=std::string::npos) 
	{   
		prefix = prefix.substr(0,prefix.find_last_of("*"));
	}
	if (prefix.find_last_of(".root")!=std::string::npos)
	{
		prefix = prefix.substr(0,prefix.find_last_of(".root")-4);
	}
	while (prefix.find_last_of("\\")!=std::string::npos)
	{
		prefix.erase(prefix.find_last_of("\\"),1);
	}
	if (char(prefix.back())!='_') prefix += "_";
	std::cout<<"prefix = "<<prefix<<std::endl;

	WCSimRootEvent* wcsimrootsuperevent = new WCSimRootEvent();
	t->SetBranchAddress("wcsimrootevent",&wcsimrootsuperevent);

	WCSimRootTrigger* wcsimrootevent;
	// Get vertex and beam direction from first event
	t->GetEntry(0);
	wcsimrootevent=wcsimrootsuperevent->GetTrigger(0);
	TVector3 vtx(wcsimrootevent->GetVtx(0),wcsimrootevent->GetVtx(1),wcsimrootevent->GetVtx(2));
	TVector3 BeamDir(((WCSimRootTrack*)wcsimrootevent->GetTracks()->At(0))->GetDir(0),((WCSimRootTrack*)wcsimrootevent->GetTracks()->At(0))->GetDir(1),((WCSimRootTrack*)wcsimrootevent->GetTracks()->At(0))->GetDir(2));
	std::cout<<"BeamDir = "<<BeamDir.x()<<" "<<BeamDir.y()<<" "<<BeamDir.z()<<std::endl;

	// Geometry tree - only need 1 "event"
	WCSimRootGeom *geo = 0;
	TTree *geotree = (TTree*)f->Get("wcsimGeoT");
	geotree->SetBranchAddress("wcsimrootgeom", &geo);
	std::cout << "Geotree has " << geotree->GetEntries() << " entries" << std::endl;
	if (geotree->GetEntries() == 0) {
		exit(9);
	}
	geotree->GetEntry(0);
	int nPMTs_type0=geo->GetWCNumPMT();
	std::cout << "geo has " << nPMTs_type0 << " PMTs" << std::endl;
	std::vector<std::vector<double>> pmt_pos(nPMTs_type0);
	std::vector<TVector3> pmt_posT(nPMTs_type0);
	std::vector<std::vector<double>> pmt_dir(nPMTs_type0);
	std::vector<double> pmt_ang(nPMTs_type0);
	std::vector<double> pmt_tof(nPMTs_type0);
	double vg = 2.20027795333758801e8*100/1.e9; // rough speed of light in water in cm/ns
	double max_r = 0, max_z = 0;
	for (int i=0;i<nPMTs_type0;i++) 
	{
		WCSimRootPMT pmt;
		pmt = geo->GetPMT(i);
		std::vector<double> pos(3);
		std::vector<double> dir(3);
		for(int j=0;j<3;j++){
			pos[j] = pmt.GetPosition(j);
			dir[j] = pmt.GetOrientation(j);
		}
		pmt_pos[i] = pos;
		pmt_dir[i] = dir;

		TVector3 pmtpos(pos[0],pos[1],pos[2]);
		pmt_posT[i] = pmtpos;
		pmt_ang[i] = BeamDir.Angle(pmtpos-vtx)*180/TMath::Pi();

		pmt_tof[i] = (pmtpos-vtx).Mag()/vg;

		// y-axis is vertical
		if (max_z<fabs(pos[1])) max_z=fabs(pos[1]);
		if (max_r<sqrt(pos[0]*pos[0]+pos[2]*pos[2]))
			if (fabs(pmt.GetOrientation(1))>0.5) max_r = sqrt(pos[0]*pos[0]+pos[2]*pos[2]);
	}

	double barrelCut = max_z-10;
	TH2D* hist_event_display = new TH2D("Charges","WCSim event display",250,-TMath::Pi()*max_r,TMath::Pi()*max_r,250,-max_z-2*max_r,max_z+2*max_r);
	std::vector<std::vector<double>> eventDiplayXY;
	for (int i=0;i<nPMTs_type0;i++)
	{
		// rotation for event display
		double y = -pmt_pos.at(i).at(0);
		double x = pmt_pos.at(i).at(2);
		double z = pmt_pos.at(i).at(1);
		std::vector<double> pmtXY;
		if (fabs(z)<barrelCut) // barrel
		{
			double th = atan2(y,x);
			pmtXY.push_back(-max_r*th);
			pmtXY.push_back(z);
		}
		else if (z>barrelCut) //top
		{
			pmtXY.push_back(-y);
			pmtXY.push_back(max_z+max_r-x);
		}
		else //bot
		{
			pmtXY.push_back(-y);
			pmtXY.push_back(-max_z-max_r+x);
		}
		eventDiplayXY.push_back(pmtXY);
	}

	std::vector<double> pmt_hit(nPMTs_type0,0.);
	int count1pc = t->GetEntries()/100;
	if (count1pc==0) count1pc=1;
	std::cout<<"Running "<<susp_event<<"-th event of "<<t->GetEntries()<<" events"<<std::endl;

	delete wcsimrootsuperevent;
	wcsimrootsuperevent = 0;  // EXTREMELY IMPORTANT

	t->GetEntry(susp_event); // Load a suspicious event - looks like decay

	wcsimrootevent = wcsimrootsuperevent->GetTrigger(0);

	std::vector<double> triggerInfo = wcsimrootevent->GetTriggerInfo();
	double triggerShift=0, triggerTime=0;
	if(wcsimrootevent->GetTriggerType()!=kTriggerNoTrig && triggerInfo.size()>=3)
	{
		triggerShift = triggerInfo[1];
		triggerTime = triggerInfo[2];
	}

	int nhits = wcsimrootevent->GetNcherenkovdigihits(); 

	// Fill digi hit histogram
	for (int i=0; i< nhits ; i++)
	{
		WCSimRootCherenkovDigiHit* wcsimrootcherenkovdigihit = (WCSimRootCherenkovDigiHit*) (wcsimrootevent->GetCherenkovDigiHits())->At(i);
		int tubeNumber     = wcsimrootcherenkovdigihit->GetTubeId()-1;
		double peForTube      = wcsimrootcherenkovdigihit->GetQ();
		double time = wcsimrootcherenkovdigihit->GetT()+triggerTime-triggerShift;

		pmt_hit[tubeNumber] += peForTube;

	}

	// Fill true hit histgram
	int ncherenkovhits     = wcsimrootevent->GetNcherenkovhits();
	TClonesArray *timeArray = wcsimrootevent->GetCherenkovHitTimes();
	for (int i=0; i< ncherenkovhits ; i++)
	{
		WCSimRootCherenkovHit *wcsimrootcherenkovhit = (WCSimRootCherenkovHit*) (wcsimrootevent->GetCherenkovHits())->At(i);
		int tubeNumber     = wcsimrootcherenkovhit->GetTubeID()-1;
		int timeArrayIndex = wcsimrootcherenkovhit->GetTotalPe(0);
		int peForTube      = wcsimrootcherenkovhit->GetTotalPe(1);
		for (int idx = timeArrayIndex; idx<timeArrayIndex+peForTube; idx++)
		{
			WCSimRootCherenkovHitTime * cht = (WCSimRootCherenkovHitTime*) timeArray->At(idx);
			TVector3 endPos(cht->GetPhotonEndPos(0)/10.,cht->GetPhotonEndPos(1)/10.,cht->GetPhotonEndPos(2)/10.); // mm to cm
			TVector3 endDir(cht->GetPhotonEndDir(0),cht->GetPhotonEndDir(1),cht->GetPhotonEndDir(2));
		}
	}

	for (int i=0;i<nPMTs_type0;i++)
	{
		hist_event_display->Fill(eventDiplayXY.at(i).at(0),eventDiplayXY.at(i).at(1),pmt_hit[i]);

	}

	SetDarkStyle();
	MakeWhiteAxes(hist_event_display);
	TCanvas* c1 = new TCanvas();
	//hist_event_display->GetXaxis()->SetRangeUser(0, 150);
	//hist_event_display->GetYaxis()->SetRangeUser(-350, -200);
	c1->SetRightMargin(0.15);
	hist_event_display->GetZaxis()->SetTitle("charge");
	hist_event_display->Draw("colz");
	double vtx_y = -vtx.x(), vtx_x = vtx.z(), vtx_z = vtx.y();
	// Extrapolate beam target point on the other side of the tank
	double target_y = vtx_y - BeamDir.x(), target_x = vtx_x + BeamDir.z(), target_z = vtx_z + BeamDir.y();
	while (sqrt(target_x*target_x+target_y*target_y)<max_r && fabs(target_z)<max_z)
	{
		target_y += -BeamDir.x(); target_x += BeamDir.z(); target_z += BeamDir.y();
	}
	double evtx, evty;
	if (fabs(vtx_z)<barrelCut)
	{
		double th = atan2(vtx_y,vtx_x);

		evtx = -max_r*th;
		evty = vtx_z;
	}
	else if (vtx_z>barrelCut)
	{
		evtx = -vtx_y;
		evty = max_z+max_r-vtx_x;
	}
	else
	{
		evtx = -vtx_y;
		evty = -max_z-max_r+vtx_x;
	}
	std::cout << "vtx x: " << vtx_x<< " vtx_y: " << vtx_y << " vtx_z: " << vtx_z << std::endl; 
	std::cout << "evtx: " << evtx << " evty: " << evty << std::endl; 
	TMarker m1(evtx,evty,29);
	m1.SetMarkerColor(kRed);
	m1.Draw();
	if (fabs(target_z)<barrelCut)
	{
		double th = atan2(target_y,target_x);

		evtx = -max_r*th;
		evty = target_z;
	}
	else if (target_z>barrelCut)
	{
		evtx = -target_y;
		evty = max_z+max_r-target_x;
	}
	else
	{
		evtx = -target_y;
		evty = -max_z-max_r+target_x;
	}
	std::cout << "beam direction is: " << BeamDir.x() << " " << BeamDir.y() << " " << BeamDir.z() << std::endl;
	std::cout << "target is: " << target_x << " " << target_y << " " << target_z << std::endl;
	std::cout << "evtx: " << evtx << " evty: " << evty << std::endl; 
	TMarker m2(evtx,evty,29);
	m2.SetMarkerColor(kWhite);
	m2.Draw();
	c1->SaveAs(Form("Event_display_evt%i.pdf",evt_nr));


	delete c1;
	f->Close();
	t->Reset();

	return 0;
}//end function
