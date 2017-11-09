#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
#include <vector>
#include <unistd.h>
#include <getopt.h>
#include <cstring>

#include "TFile.h"
#include "TTree.h"
#include "TH1F.h"
#include "TString.h"
#include "TNtuple.h"
#include "TChain.h"
#include "TMath.h"
#include "TSystem.h"
#include "TMatrixT.h"
#include "TRandom.h"
#include "TError.h"
#include "TCanvas.h"
#include "TH2F.h"
#include "TGraph.h"

#include "prob.h"
#include "params.h"
#include "SBNconfig.h"
#include "SBNchi.h"
#include "SBNdet.h"
#include "SBNspec.h"
#include "SBNosc.h"
#include "SBgeN.h"
#include "SBNfit.h"
#include "SBNfit3pN.h"

#define no_argument 0
#define required_argument 1
#define optional_argument 2

using namespace sbn;


/*************************************************************
 *************************************************************
 *		BEGIN example.cxx
 ************************************************************
 ************************************************************/
int main(int argc, char* argv[])
{




	std::string xml = "default.xml";
	int iarg = 0;
	opterr=1;
	int index; 
	int test_mode=0;
	/*************************************************************
	 *************************************************************
	 *		Command Line Argument Reading
	 ************************************************************
	 ************************************************************/

	const struct option longopts[] = 
	{
		{"xml", 		required_argument, 	0, 'x'},
		{"test",		required_argument,	0, 't'},
		{0,			no_argument, 		0,  0},
	};


	while(iarg != -1)
	{
		iarg = getopt_long(argc,argv, "x:t:", longopts, &index);

		switch(iarg)
		{
			case 'x':
				xml = optarg;
				break;
			case 't':
				test_mode = strtof(optarg,NULL);
				break;
			case '?':
			case 'h':
				std::cout<<"Allowed arguments:"<<std::endl;
				std::cout<<"\t-x\t--xml\t\tInput .xml file for SBNconfig"<<std::endl;
				return 0;
		}

	}


	/*************************************************************
	 *************************************************************
	 *		Example 3: Sterile neutrino, 3+1 sensitivity
	 ************************************************************
	 ************************************************************/

	//Time for a more model dependant example
	//using just SBNspec does not give you a huge amount of power. The idea is to make your own model dependant classes from this
	//such as here, SBNosc, a SBNspec precoded to do oscillation phyiscs (at SBL)

	//Load up your background, uboone was wierdly scaled in the precomp so fix here
	double uboonepot=0.5;
	SBNspec bkg_spec("../../data/precomp/SBN_BKG", xml);
	bkg_spec.Scale("uBooNE", uboonepot);
	bkg_spec.compressVector();


	//again create a SBNchi from this spectra
	TMatrixD Mempty(bkg_spec.fullVec.size(), bkg_spec.fullVec.size());

	SBNchi test_chi(bkg_spec, Mempty);
	test_chi.setStatOnly(true);

	//Now create a oscillation spectra, constructed the same.
	SBNosc oscSig("../../data/precomp/SBN_BKG",xml);

	//Say we want to just look at apperance mode (as its easiest to plot 2d!)
	oscSig.Scale("uBooNE",uboonepot);		
	oscSig.setBothMode();

	//so varying over all Dela M and sin^2 2 theta
	double um4min = -4.0;
	double um4max = 0;
	double um4step =0.05;

	//Want to contour plot sensitivity eventually so som standard root
	TCanvas *c1 = new TCanvas("c1","c1",600,400);
	TH2F *hcontz = new TH2F("hcontz","MicroBooNE 3+1 90\% C.L ",(um4max-um4min)/um4step,um4min,um4max, 100,-2,2);
	hcontz->GetXaxis()->SetTitle("U_{#mu 4}^2");
	hcontz->GetYaxis()->SetTitle("#Delta m^{2}_{41} (eV^{2})");

	for(double m = -2.00; m <2.00; m=m+0.04){
		for(double um4 = um4max ; um4 >= um4min; um4 = um4 - um4step){

			//always work in proper UPMNS elements!!
			double uei = 0.0;
			double umi = pow(10,um4);

			//This is where you can set up 3+N
			double imn[3] = {sqrt(pow(10,m)),0,0};
			double iue[3] = {uei,0,0};
			double ium[3] = {umi,0,0};
			double iph[3] = {0,0,0};

			//construct a signalModel
			neutrinoModel signalModel(imn,iue,ium,iph);
			signalModel.numsterile = 1; //this isnt really necessary as it can tell from imn, but nice for reading

			//And load thus model into our spectra. At this point its comuted all the necessary mass-splittins and which frequencies they are
			oscSig.load_model(signalModel);

			//And apply this oscillaion! Adding to it the bkgSpec that it was initilised with.
			std::vector<double> ans = oscSig.Oscillate();

			
			//std::cout<<"HERE: "<<ans.size()<<" oscSig.num_bins_total_compressed: "<<oscSig.num_bins_total_compressed<<std::endl;

			//Then calculate a chu^2
			double tchi=test_chi.CalcChi(ans); 

			std::cout<<"Dm^2: "<<m<<" um4^2: "<<2*um4<<" chi^2: "<<tchi<<std::endl;
			//and save wherever you like , this si just a quick hodge podge example
			int bb= 1+floor((-um4min-um4)/um4step+0.00001); 
			bb = hcontz->GetXaxis()->FindBin(um4);
			//std::cout<<"Bin: "<<bb<<"/"<<hcontz->GetNbinsX()<<std::endl; 
			hcontz->SetBinContent(bb, 1+floor(-(-2.00-m)/0.04+0.00001), tchi);


		}
	}
	Double_t contours[1];
	contours[0] = 1.64;
	hcontz->SetContour(1, contours);

	c1->cd();

	hcontz->Draw("CONT3");
	TFile * ff = new TFile("example_3.root","RECREATE");
	ff->cd();
	
	c1->Write();
	ff->Close();

	return 0;
}
