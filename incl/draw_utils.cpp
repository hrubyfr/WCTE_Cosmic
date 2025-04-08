#include "draw_utils.h"
#include <Rtypes.h>
#include <TStyle.h>
#include <iostream>

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
