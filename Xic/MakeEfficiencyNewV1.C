//efficiency study with TTree MCXic

#include <iostream>
#include <string>
#include <fstream>

#include <TH2.h>
#include <TH3.h>
#include <TH1.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>

using namespace std;

void MakeEfficiencyV5(){

	const int nset = 1;
	//const string setname[nset] = {"AOD235_LHC20MCj7all"};
	const string setname[nset] = {"LHC18MCpythia6"};
	const int nColor[3] = {1, 2, 4};

	//const int npt = 7;
	//const double ptbin[npt+1] = {1, 2, 3, 4, 5, 6, 8, 12};

	const int npt = 10;
	const double ptbin[npt+1] = {0., 1., 2., 3., 4., 5., 6., 8., 12., 16., 20};

	const int nrap = 6;
	const double rapbin[nrap+1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6};

	const int nrap2 = 10;
	const double rapbin2[nrap2+1] = {0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0};

	bool bFIRST = true;

	double centbin[101];
	for (int ii=0; ii<101; ii++) centbin[ii] = ii;

	double spdbin[201];
	for (int ii=0; ii<201; ii++) spdbin[ii] = ii;

	double v0mbin[351];
	for (int ii=0; ii<351; ii++) v0mbin[ii] = ii;

	TFile *infile = new TFile("mc_weight_zvtx.root","read");

	TF1 *fW_zvtx[5];
	for (int ii=0; ii<5; ii++){
		fW_zvtx[ii] = (TF1*)infile->Get(Form("fvtxz_ratio_to_mc_cent%d",ii));
	}

	/*
	TFile *infileWT = new TFile("outfile_MCWeight.root","read");
	TH1D *hWT[3];
	for (int ii=1; ii<3; ii++){
		hWT[ii] = (TH1D*)infileWT->Get(Form("hWeight_cent%d",ii));
	}
	*/

	TFile *infileWT = new TFile("/home/shlim/Work/Model/Pythia8/Xic/outfile_WT.root","read");
	TF1 *fmb = (TF1*)infileWT->Get("fmb");
	TH1D *fmc;

	if ( setname[0]=="AOD235_LHC20MCj7all" ){
		fmc = (TH1D*)infileWT->Get("hXicPt_ratio_1");
	}else if ( setname[0]=="LHC18MCpythia6" ){
		fmc = (TH1D*)infileWT->Get("hXicPt_ratio_0");
	}else{
		fmc = (TH1D*)infileWT->Get("hXicPt_ratio_1");
	}

	fmc->Print();

	TH2D *h2dGen[nset];
	TH3D *h3dGenSPDV[nset];

	TH1D *h1dGenCharm[nset];
	TH1D *h1dGenBottom[nset];

	TH1D *h1dGenWTCharm[nset];
	TH1D *h1dGenWTBottom[nset];

	TH2D *h2dGenSPDVCharm[nset];
	TH2D *h2dGenSPDVBottom[nset];

	TH2D *h2dGenV0MVCharm[nset];
	TH2D *h2dGenV0MVBottom[nset];
	TH2D *h2dGenV0MPCharm[nset];
	TH2D *h2dGenV0MPBottom[nset];

	TH2D *h2dGenVtxZCharm[nset];
	TH2D *h2dGenVtxZBottom[nset];
	TH2D *h2dGenVtxZWTCharm[nset][5];
	TH2D *h2dGenVtxZWTBottom[nset][5];

	TH2D *h2dReco[nset];

	TH2D *h2dRecoSPDVCharm[nset];
	TH2D *h2dRecoSPDVBottom[nset];

	TH2D *h2dRecoMassCharm[nset];
	TH2D *h2dRecoMassBottom[nset];

	TH2D *h2dRecoPtMassCharm[nset];
	TH2D *h2dRecoPtMassBottom[nset];

	TH2D *h2dRecoPtMassWTCharm[nset];
	TH2D *h2dRecoPtMassWTBottom[nset];

	TH2D *h2dRecoV0MVCharm[nset];
	TH2D *h2dRecoV0MVBottom[nset];
	TH2D *h2dRecoV0MPCharm[nset];
	TH2D *h2dRecoV0MPBottom[nset];

	TH2D *h2dRecoVtxZCharm[nset];
	TH2D *h2dRecoVtxZBottom[nset];
	TH2D *h2dRecoVtxZWTCharm[nset][5];
	TH2D *h2dRecoVtxZWTBottom[nset][5];

	TH2D *h2dRecoSPDVCut1[nset];
	TH2D *h2dRecoSPDVCut2[nset];
	TH2D *h2dRecoSPDVCut3[nset];
	TH2D *h2dRecoSPDVCut4[nset];
	TH2D *h2dRecoSPDVCut5[nset];

	TH1D *h1dGen_pt[nset];
	TH1D *h1dGen_ce[nset];
	TH1D *h1dGen_spdce[nset];

	TH1D *h1dReco_pt[nset];
	TH1D *h1dReco_ce[nset];
	TH1D *h1dReco_spdce[nset];

	TH1D *h1dEff_pt[nset];
	TH1D *h1dEff_ce[nset];
	TH1D *h1dEff_spdce[nset];

	TH2D *h2dReco_mass_pt[nset][3];
	TH2D *h2dReco_truth_mass_pt[nset];

	TH2D *hRPM_un[nset][16];

	float pTe; float echarge;
	float pTv; float vcharge;
	float Massv; float cosoa; float In_Mass; float Pt;
	float nSigmaTOF; float nSigmaTPC; float TPCCluster; float ITSCluster; float TPCPIDCluster;
	float CascDecayLength; float V0DecayLength; float DCABachToPrimVertex;
	float V0CosineOfPoiningAngleXi; float DCAV0ToPrimVertex; float DCAPosToPrimVertex; float DCANegToPrimVertex;
	float e_minmass; float e_minmass_ss; float phi; float erap; float Xirap;
	float pionTPCCluster; float protonTPCCluster; float b_pionTPCCluster;
	float e_crossedratio; float e_findable; float pion_crossedratio; float pion_findable;
	float proton_crossedratio; float proton_findable; float bpion_crossedratio; float bpion_findable;
	float XiCosineOfPoiningAngle; float pTpion; float pTproton; float pTbach;
	float MassLambda; float MassAntiLambda;

	float mcpTe; float mcecharge;
	float mcpTv; float mcvcharge;
	float mcpteXi; float mcptXic0; float mcc_flag; float mcb_flag; float mcXib; float XibeXi; float mcXibMass;
	float mcYXic0;

	float t2_V0MCentrality, t2_V0MValue;
	float t2_SPDCentrality, t2_SPDValue;
	float t2_VtxZ;

	float t3_Xic0_pT, t3_Xic0_rap, t3_c_flag, t3_b_flag; 
	float t3_VtxZ, t3_SPDValue, t3_V0MValue;
	float t3_SPDCentrality, t3_V0MCentrality;

	ifstream flist;
	char fname[300];

	TH2D *hXic0pT_V0DL = new TH2D("hXic0pT_V0DL","",10,0,10,100,0,100);
	TH2D *hXic0pT_XiDL = new TH2D("hXic0pT_XiDL","",10,0,10,100,0,100);

	for (int iset=0; iset<nset; iset++){

		float NumGen = 0.0;

		//h2dGen[iset] = new TH2D(Form("h2dGen_set%d",iset),"",npt,ptbin,100,centbin);
		//h3dGenSPDV[iset] = new TH3D(Form("h3dGenSPDV_set%d",iset),"",npt,ptbin,nrap2,rapbin2,200,spdbin);

		h1dGenCharm[iset] = new TH1D(Form("h1dGenCharm_set%d",iset),"",100,0,20);
		h1dGenBottom[iset] = new TH1D(Form("h1dGenBottom_set%d",iset),"",100,0,20);

		h1dGenWTCharm[iset] = new TH1D(Form("h1dGenWTCharm_set%d",iset),"",100,0,20);
		h1dGenWTBottom[iset] = new TH1D(Form("h1dGenWTBottom_set%d",iset),"",100,0,20);

		h2dGenSPDVCharm[iset] = new TH2D(Form("h2dGenSPDVCharm_set%d",iset),"",npt,ptbin,200,0,200);
		h2dGenSPDVBottom[iset] = new TH2D(Form("h2dGenSPDVBottom_set%d",iset),"",npt,ptbin,200,0,200);

		h2dGenV0MVCharm[iset] = new TH2D(Form("h2dGenV0MVCharm_set%d",iset),"",npt,ptbin,360,0,360);
		h2dGenV0MVBottom[iset] = new TH2D(Form("h2dGenV0MVBottom_set%d",iset),"",npt,ptbin,360,0,360);
		h2dGenV0MPCharm[iset] = new TH2D(Form("h2dGenV0MPCharm_set%d",iset),"",npt,ptbin,100,0,100);
		h2dGenV0MPBottom[iset] = new TH2D(Form("h2dGenV0MPBottom_set%d",iset),"",npt,ptbin,100,0,100);

		h2dGenVtxZCharm[iset] = new TH2D(Form("h2dGenVtxZCharm_set%d",iset),"",npt,ptbin,20,-10,10);
		h2dGenVtxZBottom[iset] = new TH2D(Form("h2dGenVtxZBottom_set%d",iset),"",npt,ptbin,20,-10,10);
		for (int ii=0; ii<5; ii++){
			h2dGenVtxZWTCharm[iset][ii] = new TH2D(Form("h2dGenVtxZWTCharm_set%d_%d",iset,ii),"",npt,ptbin,20,-10,10);
			h2dGenVtxZWTBottom[iset][ii] = new TH2D(Form("h2dGenVtxZWTBottom_set%d_%d",iset,ii),"",npt,ptbin,20,-10,10);
		}

		h2dReco[iset] = new TH2D(Form("h2dReco_set%d",iset),"",npt,ptbin,100,centbin);

		h2dRecoMassCharm[iset] = new TH2D(Form("h2dRecoMassCharm_set%d",iset),"",npt,ptbin,40,1.0,3.0);
		h2dRecoMassBottom[iset] = new TH2D(Form("h2dRecoMassBottom_set%d",iset),"",npt,ptbin,40,1.0,3.0);

		h2dRecoPtMassCharm[iset] = new TH2D(Form("h2dRecoPtMassCharm_set%d",iset),"",npt,ptbin,40,1.0,3.0);
		h2dRecoPtMassBottom[iset] = new TH2D(Form("h2dRecoPtMassBottom_set%d",iset),"",npt,ptbin,40,1.0,3.0);

		h2dRecoPtMassWTCharm[iset] = new TH2D(Form("h2dRecoPtMassWTCharm_set%d",iset),"",npt,ptbin,40,1.0,3.0);
		h2dRecoPtMassWTBottom[iset] = new TH2D(Form("h2dRecoPtMassWTBottom_set%d",iset),"",npt,ptbin,40,1.0,3.0);

		h2dRecoSPDVCharm[iset] = new TH2D(Form("h2dRecoSPDVCharm_set%d",iset),"",npt,ptbin,200,0,200);
		h2dRecoSPDVBottom[iset] = new TH2D(Form("h2dRecoSPDVBottom_set%d",iset),"",npt,ptbin,200,0,200);

		h2dRecoV0MVCharm[iset] = new TH2D(Form("h2dRecoV0MVCharm_set%d",iset),"",npt,ptbin,360,0,360);
		h2dRecoV0MVBottom[iset] = new TH2D(Form("h2dRecoV0MVBottom_set%d",iset),"",npt,ptbin,360,0,360);
		h2dRecoV0MPCharm[iset] = new TH2D(Form("h2dRecoV0MPCharm_set%d",iset),"",npt,ptbin,100,0,100);
		h2dRecoV0MPBottom[iset] = new TH2D(Form("h2dRecoV0MPBottom_set%d",iset),"",npt,ptbin,100,0,100);

		h2dRecoVtxZCharm[iset] = new TH2D(Form("h2dRecoVtxZCharm_set%d",iset),"",npt,ptbin,20,-10,10);
		h2dRecoVtxZBottom[iset] = new TH2D(Form("h2dRecoVtxZBottom_set%d",iset),"",npt,ptbin,20,-10,10);
		for (int ii=0; ii<5; ii++){
			h2dRecoVtxZWTCharm[iset][ii] = new TH2D(Form("h2dRecoVtxZWTCharm_set%d_%d",iset,ii),"",npt,ptbin,20,-10,10);
			h2dRecoVtxZWTBottom[iset][ii] = new TH2D(Form("h2dRecoVtxZWTBottom_set%d_%d",iset,ii),"",npt,ptbin,20,-10,10);
		}

		h2dRecoSPDVCut1[iset] = new TH2D(Form("h2dRecoSPDVCut1_set%d",iset),"",npt,ptbin,200,0,200);
		h2dRecoSPDVCut2[iset] = new TH2D(Form("h2dRecoSPDVCut2_set%d",iset),"",npt,ptbin,200,0,200);
		h2dRecoSPDVCut3[iset] = new TH2D(Form("h2dRecoSPDVCut3_set%d",iset),"",npt,ptbin,200,0,200);
		h2dRecoSPDVCut4[iset] = new TH2D(Form("h2dRecoSPDVCut4_set%d",iset),"",npt,ptbin,200,0,200);
		h2dRecoSPDVCut5[iset] = new TH2D(Form("h2dRecoSPDVCut5_set%d",iset),"",npt,ptbin,200,0,200);

		h2dReco_truth_mass_pt[iset] = new TH2D(Form("h2dReco_truth_mass_pt_set%d",iset),"",npt,ptbin,30,1.0,4.0);
		for (int ichg=0; ichg<3; ichg++){
			h2dReco_mass_pt[iset][ichg] = new TH2D(Form("h2dReco_mass_pt_set%d_chg%d",iset,ichg),"",npt,ptbin,30,1.0,4.0);
		}

		for (int ii=0; ii<16; ii++){
			hRPM_un[iset][ii] = new TH2D(Form("hRPM_un_set%d_mass%02d",iset,ii),"",npt,ptbin,npt,ptbin);
		}

		sprintf(fname,"list_%s.lst",setname[iset].c_str());
		flist.open(fname);

		while ( flist >> fname ){

			TFile *infile = new TFile(fname,"read");
			cout << "OPEN: " << fname << endl;

			TDirectoryFile *tdf = (TDirectoryFile*)infile->Get("PWG3_D2H_Xic02eXipp13TeV_HM");

			if ( !tdf ){
				cout << "CAN NOT Get Directory!!" << endl;
				infile->Close();
				delete infile;
				continue;
			}

			TTree *T0 = (TTree*)tdf->Get("eXiTree");
			TTree *T1 = (TTree*)tdf->Get("MCTree");
			TTree *T2 = (TTree*)tdf->Get("EventTree");
			
			TTree *T3 = (TTree*)tdf->Get("MCXicTree");

			if ( !T0 || !T1 || !T2 ){
				cout << "CAN NOT Get TTree!!" << endl;
				infile->Close();
				delete infile;
				continue;
			}

			/*
			TH2D *_h2dGen = (TH2D*)infile->Get("histogram")->FindObject("hTrueXic0");
			TH3D *_h3dGenSPDV = (TH3D*)infile->Get("histogram")->FindObject("hTrueXic0SPD");

			TH2D *_hXic0PtFromCharm1SPD = (TH2D*)infile->Get("histogram")->FindObject("hXic0PtFromCharm1SPD");
			TH2D *_hXic0PtFromCharm2SPD = (TH2D*)infile->Get("histogram")->FindObject("hXic0PtFromCharm2SPD");

			TH2D *_hXic0PtFromBottom1SPD = (TH2D*)infile->Get("histogram")->FindObject("hXic0PtFromBottom1SPD");
			TH2D *_hXic0PtFromBottom2SPD = (TH2D*)infile->Get("histogram")->FindObject("hXic0PtFromBottom2SPD");

			if ( !_h2dGen || !_h3dGenSPDV ){
				infile->Close();
				delete infile;
				continue;
			}

			NumGen += _h2dGen->Integral();

			h2dGen[iset]->Add(_h2dGen);
			h3dGenSPDV[iset]->Add(_h3dGenSPDV);

			h2dGenSPDVCharm[iset]->Add(_hXic0PtFromCharm1SPD);
			h2dGenSPDVCharm[iset]->Add(_hXic0PtFromCharm2SPD);

			h2dGenSPDVBottom[iset]->Add(_hXic0PtFromBottom1SPD);
			h2dGenSPDVBottom[iset]->Add(_hXic0PtFromBottom2SPD);
			*/

			T0->SetBranchAddress("pTe",&pTe);
			T0->SetBranchAddress("echarge",&echarge);
			T0->SetBranchAddress("TOFnSigma",&nSigmaTOF);
			T0->SetBranchAddress("TPCnSigma",&nSigmaTPC);
			T0->SetBranchAddress("TPCPID",&TPCPIDCluster);
			T0->SetBranchAddress("ITS",&ITSCluster);
			T0->SetBranchAddress("e_crossedrows",&e_crossedratio);
			T0->SetBranchAddress("e_findable",&e_findable);
			T0->SetBranchAddress("phi",&phi);
			T0->SetBranchAddress("erap",&erap);
			T0->SetBranchAddress("e_minmass",&e_minmass);
			T0->SetBranchAddress("e_minmass_ss",&e_minmass_ss);
			T0->SetBranchAddress("pTv",&pTv);
			T0->SetBranchAddress("vcharge",&vcharge);
			T0->SetBranchAddress("Massv",&Massv);
			T0->SetBranchAddress("MassLambda",&MassLambda);
			T0->SetBranchAddress("MassAntiLambda",&MassAntiLambda);  //mod
			T0->SetBranchAddress("V0DecayLength",&V0DecayLength);
			T0->SetBranchAddress("CascDecayLength",&CascDecayLength);
			T0->SetBranchAddress("DCABachToPrimVertex",&DCABachToPrimVertex);
			T0->SetBranchAddress("DCAV0NegToPrimVertex",&DCANegToPrimVertex);
			T0->SetBranchAddress("DCAV0PosToPrimVertex",&DCAPosToPrimVertex);
			T0->SetBranchAddress("V0CosineOfPoiningAngleXi",&V0CosineOfPoiningAngleXi);
			T0->SetBranchAddress("XiCosineOfPoiningAngle",&XiCosineOfPoiningAngle);  //modify
			T0->SetBranchAddress("DCAV0ToPrimVertex",&DCAV0ToPrimVertex);   ///new
			T0->SetBranchAddress("Xirap",&Xirap);  //mod
			T0->SetBranchAddress("pion_crossedrows",&pion_crossedratio);
			T0->SetBranchAddress("pion_findable",&pion_findable);
			T0->SetBranchAddress("proton_crossedrows",&proton_crossedratio);
			T0->SetBranchAddress("proton_findable",&proton_findable);
			T0->SetBranchAddress("bpion_crossedratio",&bpion_crossedratio);
			T0->SetBranchAddress("bpion_findable",&bpion_findable);
			T0->SetBranchAddress("pTpion",&pTpion);
			T0->SetBranchAddress("pTproton",&pTproton);
			T0->SetBranchAddress("pTbach",&pTbach);
			T0->SetBranchAddress("cosoa",&cosoa);  //mod
			T0->SetBranchAddress("In_Mass",&In_Mass);
			T0->SetBranchAddress("eXiPt",&Pt);

			T1->SetBranchAddress("mcpTe",&mcpTe);
			T1->SetBranchAddress("mcecharge",&mcecharge);
			T1->SetBranchAddress("mcpTv",&mcpTv);
			T1->SetBranchAddress("mcvcharge",&mcvcharge);
			T1->SetBranchAddress("mcpTXic0",&mcptXic0);
			T1->SetBranchAddress("mcpTeXi",&mcpteXi);
			//T1->SetBranchAddress("mcYXic0",&mcYXic0);
			T1->SetBranchAddress("c_flag",&mcc_flag);
			T1->SetBranchAddress("b_flag",&mcb_flag);
			T1->SetBranchAddress("mcpTXib",&mcXib);
			T1->SetBranchAddress("mceXipTb",&XibeXi);
			//T1->SetBranchAddress("mcXibMass",&mcXibMass);

			T2->SetBranchAddress("fCentrality",&t2_V0MCentrality);
			T2->SetBranchAddress("fCentralSPD",&t2_SPDCentrality);
			T2->SetBranchAddress("fNSPDTracklets",&t2_SPDValue);
			T2->SetBranchAddress("fNV0M",&t2_V0MValue);
			T2->SetBranchAddress("fVtxZ",&t2_VtxZ);

			T3->SetBranchAddress("Xic0_pT",&t3_Xic0_pT);
			T3->SetBranchAddress("Xic0_rap",&t3_Xic0_rap);
			T3->SetBranchAddress("c_flag",&t3_c_flag);
			T3->SetBranchAddress("b_flag",&t3_b_flag);
			T3->SetBranchAddress("VtxZ",&t3_VtxZ);
			T3->SetBranchAddress("NSPDTracklets",&t3_SPDValue);
			T3->SetBranchAddress("NV0M",&t3_V0MValue);
			T3->SetBranchAddress("Centrality",&t3_V0MCentrality);

			//continue;

			int nentries0 = T0->GetEntries();
			int nentries1 = T1->GetEntries();
			int nentries2 = T2->GetEntries();
			int nentries3 = T3->GetEntries();

			for (int ien=0; ien<nentries3; ien++){

				T3->GetEntry(ien);

				if ( fabs(t3_VtxZ)>10.0 ) continue;
				if ( fabs(t3_Xic0_rap)>0.5 ) continue;

				float wt = fmb->Eval(t3_Xic0_pT)/fmc->GetBinContent(fmc->FindBin(t3_Xic0_pT));

				if ( t3_c_flag>0.5 ){

					h1dGenCharm[iset]->Fill(t3_Xic0_pT);
					h1dGenWTCharm[iset]->Fill(t3_Xic0_pT, wt);

					h2dGenSPDVCharm[iset]->Fill(t3_Xic0_pT, t3_SPDValue);
					h2dGenV0MVCharm[iset]->Fill(t3_Xic0_pT, t3_V0MValue);
					h2dGenV0MPCharm[iset]->Fill(t3_Xic0_pT, t3_V0MCentrality);
					h2dGenVtxZCharm[iset]->Fill(t3_Xic0_pT, t3_VtxZ);
					for (int ii=0; ii<5; ii++){
						float wt = fW_zvtx[ii]->Eval(t3_VtxZ);
						h2dGenVtxZWTCharm[iset][ii]->Fill(t3_Xic0_pT, t3_VtxZ, wt);
					}
				}else if ( t3_b_flag>0.5 ){

					h1dGenBottom[iset]->Fill(t3_Xic0_pT);
					h1dGenWTBottom[iset]->Fill(t3_Xic0_pT, wt);

					h2dGenSPDVBottom[iset]->Fill(t3_Xic0_pT, t3_SPDValue);
					h2dGenV0MVBottom[iset]->Fill(t3_Xic0_pT, t3_V0MValue);
					h2dGenV0MPBottom[iset]->Fill(t3_Xic0_pT, t3_V0MCentrality);
					h2dGenVtxZBottom[iset]->Fill(t3_Xic0_pT, t3_VtxZ);
					for (int ii=0; ii<5; ii++){
						float wt = fW_zvtx[ii]->Eval(t3_VtxZ);
						h2dGenVtxZWTBottom[iset][ii]->Fill(t3_Xic0_pT, t3_VtxZ, wt);
					}
				}

			}//

			if ( nentries0!=nentries1 || nentries0!=nentries2 ){
				cout << "inconsistent entries " << nentries0 << " " << nentries1 << " " << nentries2 << endl;
			}

			for (int ien=0; ien<nentries0; ien++){
				T0->GetEntry(ien);
				T1->GetEntry(ien);
				T2->GetEntry(ien);

				if ( fabs(t2_VtxZ)>10.0 ) continue;

				if ( fabs(Massv-1.32171)>0.008 ) continue;  //Xi mass tolerance
				if ( In_Mass<1.3 ) continue;  //pair low limit
				if ( fabs(pTe)>900 ) continue;  //tmp tree reject
				//if ( mcptXic0<0 ) continue;

				float S_e_nsigma_cut	= -3.9 + (1.17*pTe) - (0.094*pTe*pTe);
				float VT_e_nsigma_cut	= -3.5 + (1.15*pTe) - (0.090*pTe*pTe);
				float T_e_nsigma_cut	= -3.7 + (1.17*pTe) - (0.094*pTe*pTe); ///need to modify
				float L_e_nsigma_cut	= -4.1 + (1.17*pTe) - (0.094*pTe*pTe); ///need to modify
				float VL_e_nsigma_cut	= -4.3 + (1.17*pTe) - (0.094*pTe*pTe);

				if ( pTe>=5.0 ){
					S_e_nsigma_cut	= -3.9+(1.17*5)-(0.094*25);
					VT_e_nsigma_cut	= -3.5+(1.15*5)-(0.090*25);
					T_e_nsigma_cut	= -3.7+(1.17*5)-(0.094*25);
					L_e_nsigma_cut	= -4.1+(1.17*5)-(0.094*25);
					VL_e_nsigma_cut	= -4.3+(1.17*5)-(0.094*25);
				}    

				Bool_t Xi_Topology_Stand_flag = kFALSE; if(V0DecayLength>2.67 && CascDecayLength>0.50 && DCABachToPrimVertex>0.07
						&& DCANegToPrimVertex>0.095 && DCAPosToPrimVertex>0.095 && V0CosineOfPoiningAngleXi>0.983 && XiCosineOfPoiningAngle>0.983 && DCAV0ToPrimVertex>0.09) Xi_Topology_Stand_flag = kTRUE;

				/*
				Bool_t Xi_Topology_Stand_flag = kFALSE; if(V0DecayLength>2.67 && CascDecayLength>0.38 && DCABachToPrimVertex>0.0204
						&& DCANegToPrimVertex>0.073 && DCAPosToPrimVertex>0.073 && V0CosineOfPoiningAngleXi>0.983 && XiCosineOfPoiningAngle>0.983 &&DCAV0ToPrimVertex>0.03) Xi_Topology_Stand_flag = kTRUE;
				Bool_t Xi_Topology_VLoose_flag = kFALSE; if(V0DecayLength>0.02 && CascDecayLength>0.02 && DCABachToPrimVertex>0.01
						&& DCANegToPrimVertex>0.05 && DCAPosToPrimVertex>0.05 && V0CosineOfPoiningAngleXi>0.98 && XiCosineOfPoiningAngle>0.98 &&DCAV0ToPrimVertex>0.01) Xi_Topology_VLoose_flag = kTRUE;
				Bool_t Xi_Topology_Loose_flag = kFALSE; if(V0DecayLength>1.55 && CascDecayLength>0.29 && DCABachToPrimVertex>0.0146
						&& DCANegToPrimVertex>0.061 && DCAPosToPrimVertex>0.061 && V0CosineOfPoiningAngleXi>0.981 && XiCosineOfPoiningAngle>0.981 &&DCAV0ToPrimVertex>0.02) Xi_Topology_Loose_flag = kTRUE;
				Bool_t Xi_Topology_Tight_flag = kFALSE; if(V0DecayLength>3.6 && CascDecayLength>0.53 && DCABachToPrimVertex>0.037
						&& DCANegToPrimVertex>0.088 && DCAPosToPrimVertex>0.088 && V0CosineOfPoiningAngleXi>0.9839 && XiCosineOfPoiningAngle>0.9839 &&DCAV0ToPrimVertex>0.04) Xi_Topology_Tight_flag = kTRUE;
				Bool_t Xi_Topology_VTight_flag = kFALSE; if(V0DecayLength>4.39 && CascDecayLength>0.72 && DCABachToPrimVertex>0.037
						&& DCANegToPrimVertex>0.102 && DCAPosToPrimVertex>0.102 && V0CosineOfPoiningAngleXi>0.985 && XiCosineOfPoiningAngle>0.985 &&DCAV0ToPrimVertex>0.06) Xi_Topology_VTight_flag = kTRUE;

				Bool_t Xi_Topology_Stand_V0 = kFALSE; if(V0DecayLength>0.02 && CascDecayLength>0.38 && DCABachToPrimVertex>0.0204
						&& DCANegToPrimVertex>0.073 && DCAPosToPrimVertex>0.073 && V0CosineOfPoiningAngleXi>0.983 && XiCosineOfPoiningAngle>0.983 &&DCAV0ToPrimVertex>0.03) Xi_Topology_Stand_V0 = kTRUE;
				Bool_t Xi_Topology_Stand_Xi = kFALSE; if(V0DecayLength>2.67 && CascDecayLength>0.02 && DCABachToPrimVertex>0.0204
						&& DCANegToPrimVertex>0.073 && DCAPosToPrimVertex>0.073 && V0CosineOfPoiningAngleXi>0.983 && XiCosineOfPoiningAngle>0.983 &&DCAV0ToPrimVertex>0.03) Xi_Topology_Stand_Xi = kTRUE;
				Bool_t Xi_Topology_Stand_DCAb = kFALSE; if(V0DecayLength>2.67 && CascDecayLength>0.38 && DCABachToPrimVertex>0.01
						&& DCANegToPrimVertex>0.073 && DCAPosToPrimVertex>0.073 && V0CosineOfPoiningAngleXi>0.983 && XiCosineOfPoiningAngle>0.983 &&DCAV0ToPrimVertex>0.03) Xi_Topology_Stand_DCAb = kTRUE;
						*/

				Bool_t Xi_Recon_VTight_flag = kFALSE;
				Bool_t Xi_Recon_Tight_flag = kFALSE;
				Bool_t Xi_Recon_Stand_flag = kFALSE;
				Bool_t Xi_Recon_Loose_flag = kFALSE;
				Bool_t Xi_Recon_VLoose_flag = kFALSE;

				if(pion_crossedratio/pion_findable>0.77 && proton_crossedratio/proton_findable>0.77 && bpion_crossedratio/bpion_findable>0.77 && pion_crossedratio>70 && proton_crossedratio>70 && bpion_crossedratio>70) Xi_Recon_Stand_flag = kTRUE;

				/*
				if(pion_crossedratio/pion_findable>0.81 && proton_crossedratio/proton_findable>0.81 && bpion_crossedratio/bpion_findable>0.81 && pion_crossedratio>80 && proton_crossedratio>80 && bpion_crossedratio>80) Xi_Recon_VTight_flag = kTRUE;
				if(pion_crossedratio/pion_findable>0.79 && proton_crossedratio/proton_findable>0.79 && bpion_crossedratio/bpion_findable>0.79 && pion_crossedratio>75 && proton_crossedratio>75 && bpion_crossedratio>75) Xi_Recon_Tight_flag = kTRUE;
				if(pion_crossedratio/pion_findable>0.77 && proton_crossedratio/proton_findable>0.77 && bpion_crossedratio/bpion_findable>0.77 && pion_crossedratio>70 && proton_crossedratio>70 && bpion_crossedratio>70) Xi_Recon_Stand_flag = kTRUE;
				if(pion_crossedratio/pion_findable>0.75 && proton_crossedratio/proton_findable>0.75 && bpion_crossedratio/bpion_findable>0.75 && pion_crossedratio>65 && proton_crossedratio>65 && bpion_crossedratio>65) Xi_Recon_Loose_flag = kTRUE;
				if(pion_crossedratio/pion_findable>0.70 && proton_crossedratio/proton_findable>0.70 && bpion_crossedratio/bpion_findable>0.70 && pion_crossedratio>65 && proton_crossedratio>65 && bpion_crossedratio>65) Xi_Recon_VLoose_flag = kTRUE;
				*/


				Bool_t e_Recon_Stand_flag = kFALSE; if(e_crossedratio>70 && TPCPIDCluster>50 && e_crossedratio/e_findable>0.8 && ITSCluster>=3) e_Recon_Stand_flag = kTRUE;
				/*
				Bool_t e_Recon_VTight_flag = kFALSE;  if(e_crossedratio>85 && TPCPIDCluster>60 && e_crossedratio/e_findable>0.9 && ITSCluster>=3) e_Recon_VTight_flag = kTRUE;
				Bool_t e_Recon_Tight_flag = kFALSE;  if(e_crossedratio>75 && TPCPIDCluster>55 && e_crossedratio/e_findable>0.85 && ITSCluster>=3) e_Recon_Tight_flag = kTRUE;
				Bool_t e_Recon_Stand_flag = kFALSE; if(e_crossedratio>70 && TPCPIDCluster>50 && e_crossedratio/e_findable>0.8 && ITSCluster>=3) e_Recon_Stand_flag = kTRUE;
				Bool_t e_Recon_Loose_flag = kFALSE;  if(e_crossedratio>65 && TPCPIDCluster>45 && e_crossedratio/e_findable>0.8 && ITSCluster>=3) e_Recon_Loose_flag = kTRUE;
				Bool_t e_Recon_VLoose_flag = kFALSE;  if(e_crossedratio>65 && TPCPIDCluster>40 && e_crossedratio/e_findable>0.75 && ITSCluster>=3) e_Recon_VLoose_flag = kTRUE;
				*/

				Bool_t e_PID_Stand_flag = kFALSE; if(fabs(nSigmaTOF)<=3 && nSigmaTPC>=S_e_nsigma_cut && nSigmaTPC<=3) e_PID_Stand_flag = kTRUE;

				/*
				Bool_t e_PID_VTight_flag = kFALSE; if(fabs(nSigmaTOF)<=2 && nSigmaTPC>=VT_e_nsigma_cut && nSigmaTPC<=3) e_PID_VTight_flag = kTRUE;
				Bool_t e_PID_Tight_flag = kFALSE; if(fabs(nSigmaTOF)<=3 && nSigmaTPC>=T_e_nsigma_cut && nSigmaTPC<=3) e_PID_Tight_flag = kTRUE;
				Bool_t e_PID_Stand_flag = kFALSE; if(fabs(nSigmaTOF)<=3 && nSigmaTPC>=S_e_nsigma_cut && nSigmaTPC<=3) e_PID_Stand_flag = kTRUE;
				Bool_t e_PID_Loose_flag = kFALSE; if(fabs(nSigmaTOF)<=3 && nSigmaTPC>=L_e_nsigma_cut && nSigmaTPC<=3) e_PID_Loose_flag = kTRUE;
				Bool_t e_PID_VLoose_flag = kFALSE; if(fabs(nSigmaTOF)<=3 && nSigmaTPC>=VL_e_nsigma_cut && nSigmaTPC<=3) e_PID_VLoose_flag = kTRUE;
				*/

				Bool_t OPAngle_Tight_flag = kFALSE; if(cosoa>cos(70*(3.141592/180))) OPAngle_Tight_flag = kTRUE;
				Bool_t OPAngle_Stand_flag = kFALSE; if(cosoa>cos(90*(3.141592/180))) OPAngle_Stand_flag = kTRUE;
				Bool_t OPAngle_Loose_flag = kFALSE; if(cosoa>cos(110*(3.141592/180)))OPAngle_Loose_flag = kTRUE;

				Bool_t PairMass_Tight_flag = kFALSE; if(In_Mass<2.3) PairMass_Tight_flag = kTRUE;
				Bool_t PairMass_Stand_flag = kFALSE; if(In_Mass<2.5) PairMass_Stand_flag = kTRUE;
				Bool_t PairMass_Loose_flag = kFALSE; if(In_Mass<2.7) PairMass_Loose_flag = kTRUE;

				//float wt = fmb->Eval(mcptXic0)/fmc->GetBinContent(fmc->FindBin(mcptXic0));
				float wt = 1./fmc->GetBinContent(fmc->FindBin(mcptXic0));

				if ( mcptXic0>0 ){

					hXic0pT_V0DL->Fill(mcptXic0, V0DecayLength);
					hXic0pT_XiDL->Fill(mcptXic0, CascDecayLength);

					if ( PairMass_Stand_flag && OPAngle_Stand_flag ){
						h2dRecoSPDVCut1[iset]->Fill(mcptXic0, t2_SPDValue);
						if ( e_Recon_Stand_flag ){
							h2dRecoSPDVCut2[iset]->Fill(mcptXic0, t2_SPDValue);
							if ( e_PID_Stand_flag ){
								h2dRecoSPDVCut3[iset]->Fill(mcptXic0, t2_SPDValue);
								if ( Xi_Recon_Stand_flag ){
									h2dRecoSPDVCut4[iset]->Fill(mcptXic0, t2_SPDValue);
									if ( Xi_Topology_Stand_flag ){
										h2dRecoSPDVCut5[iset]->Fill(mcptXic0, t2_SPDValue);
									}
								}
							}
						}
					}

					if ( 
							e_Recon_Stand_flag 
							&& Xi_Recon_Stand_flag 
							&& e_PID_Stand_flag 
							&& Xi_Topology_Stand_flag
							&& OPAngle_Stand_flag 
							&& PairMass_Stand_flag
						 ){

						h2dReco[iset]->Fill(mcptXic0, t2_V0MCentrality);

						if ( mcc_flag>0.5 ){
							h2dRecoSPDVCharm[iset]->Fill(mcptXic0, t2_SPDValue);
							h2dRecoV0MVCharm[iset]->Fill(mcptXic0, t2_V0MValue);
							h2dRecoV0MPCharm[iset]->Fill(mcptXic0, t2_V0MCentrality);
							h2dRecoVtxZCharm[iset]->Fill(mcptXic0, t2_VtxZ);
							for (int ii=0; ii<5; ii++){
								float wt = fW_zvtx[ii]->Eval(t2_VtxZ);
								h2dRecoVtxZWTCharm[iset][ii]->Fill(mcptXic0, t2_VtxZ, wt);
							}
						}else if ( mcb_flag>0.5 ){
							h2dRecoSPDVBottom[iset]->Fill(mcptXic0, t2_SPDValue);
							h2dRecoV0MVBottom[iset]->Fill(mcptXic0, t2_V0MValue);
							h2dRecoV0MPBottom[iset]->Fill(mcptXic0, t2_V0MCentrality);
							h2dRecoVtxZBottom[iset]->Fill(mcptXic0, t2_VtxZ);
							for (int ii=0; ii<5; ii++){
								float wt = fW_zvtx[ii]->Eval(t2_VtxZ);
								h2dRecoVtxZWTBottom[iset][ii]->Fill(mcptXic0, t2_VtxZ, wt);
							}
						}
					}//cuts with pair mass

					if ( 
							e_Recon_Stand_flag 
							&& Xi_Recon_Stand_flag 
							&& e_PID_Stand_flag 
							&& Xi_Topology_Stand_flag
							&& OPAngle_Stand_flag 
						 ){

						if ( mcc_flag>0.5 ){
							h2dRecoMassCharm[iset]->Fill(mcptXic0, In_Mass);
							h2dRecoPtMassCharm[iset]->Fill(Pt, In_Mass);
							h2dRecoPtMassWTCharm[iset]->Fill(Pt, In_Mass, wt);
						}else if ( mcb_flag>0.5 ){
							h2dRecoMassBottom[iset]->Fill(mcptXic0, In_Mass);
							h2dRecoPtMassBottom[iset]->Fill(Pt, In_Mass);
							h2dRecoPtMassWTBottom[iset]->Fill(Pt, In_Mass, wt);
						}

						if ( In_Mass<2.5 ){
							hRPM_un[iset][0]->Fill(Pt, mcptXic0);
						}
						if ( In_Mass<2.4 ){
							hRPM_un[iset][1]->Fill(Pt, mcptXic0);
						}
						if ( In_Mass<2.3 ){
							hRPM_un[iset][2]->Fill(Pt, mcptXic0);
						}
						if ( In_Mass<2.2 ){
							hRPM_un[iset][3]->Fill(Pt, mcptXic0);
						}
						if ( In_Mass<2.1 ){
							hRPM_un[iset][4]->Fill(Pt, mcptXic0);
						}
						if ( In_Mass<2.0 ){
							hRPM_un[iset][5]->Fill(Pt, mcptXic0);
						}
						if ( In_Mass<1.9 ){
							hRPM_un[iset][6]->Fill(Pt, mcptXic0);
						}
						if ( In_Mass<1.8 ){
							hRPM_un[iset][7]->Fill(Pt, mcptXic0);
						}

						if ( In_Mass>1.3 ){
							hRPM_un[iset][8]->Fill(Pt, mcptXic0);
						}
						if ( In_Mass>1.4 ){
							hRPM_un[iset][9]->Fill(Pt, mcptXic0);
						}
						if ( In_Mass>1.5 ){
							hRPM_un[iset][10]->Fill(Pt, mcptXic0);
						}
						if ( In_Mass>1.6 ){
							hRPM_un[iset][11]->Fill(Pt, mcptXic0);
						}
						if ( In_Mass>1.7 ){
							hRPM_un[iset][12]->Fill(Pt, mcptXic0);
						}
						if ( In_Mass>1.8 ){
							hRPM_un[iset][13]->Fill(Pt, mcptXic0);
						}
						if ( In_Mass>1.9 ){
							hRPM_un[iset][14]->Fill(Pt, mcptXic0);
						}
						if ( In_Mass>2.0 ){
							hRPM_un[iset][15]->Fill(Pt, mcptXic0);
						}

					}//cuts without pair mass

				}//mcptXic0

				/*
				if ( 
						e_Recon_Stand_flag 
						&& Xi_Recon_Stand_flag 
						&& e_PID_Stand_flag 
						&& Xi_Topology_Stand_flag
						&& OPAngle_Stand_flag 
						&& PairMass_Stand_flag
					 ){
					if ( echarge>0 && vcharge>0 ){
						h2dReco_mass_pt[iset][1]->Fill(In_Mass, Pt);
					}else if ( echarge<0 && vcharge<0 ){
						h2dReco_mass_pt[iset][2]->Fill(In_Mass, Pt);
					}else{
						h2dReco_mass_pt[iset][0]->Fill(In_Mass, Pt);
					}

					if ( mcptXic0>0 ){
						h2dReco_truth_mass_pt[iset]->Fill(In_Mass, Pt);
					}
				}
				*/

			}//ien

			infile->Close();
			delete infile;

		}//while

		flist.close();

		//cout << "Number Generated Xic0 " << NumGen << " " << h2dGen[iset]->Integral() << endl;  

	}//iset

	//return;

	TFile *outfile = new TFile("outfile_efficiencyV5_LHC18MCpythia6.root","recreate");
	//TFile *outfile = new TFile("outfile_efficiencyV5_AOD235_LHC20MCj7all.root","recreate");
	for (int iset=0; iset<nset; iset++){

		h1dGenCharm[iset]->Write();
		h1dGenBottom[iset]->Write();

		h1dGenWTCharm[iset]->Write();
		h1dGenWTBottom[iset]->Write();

		h2dRecoSPDVCut1[iset]->Write();
		h2dRecoSPDVCut2[iset]->Write();
		h2dRecoSPDVCut3[iset]->Write();
		h2dRecoSPDVCut4[iset]->Write();
		h2dRecoSPDVCut5[iset]->Write();

		h2dRecoSPDVCharm[iset]->Write();
		h2dRecoSPDVBottom[iset]->Write();

		h2dRecoV0MVCharm[iset]->Write();
		h2dRecoV0MVBottom[iset]->Write();
		h2dRecoV0MPCharm[iset]->Write();
		h2dRecoV0MPBottom[iset]->Write();

		h2dRecoVtxZCharm[iset]->Write();
		h2dRecoVtxZBottom[iset]->Write();
		for (int ii=0; ii<5; ii++){
			h2dRecoVtxZWTCharm[iset][ii]->Write();
			h2dRecoVtxZWTBottom[iset][ii]->Write();
		}

		h2dRecoMassCharm[iset]->Write();
		h2dRecoMassBottom[iset]->Write();

		h2dRecoPtMassCharm[iset]->Write();
		h2dRecoPtMassBottom[iset]->Write();

		h2dRecoPtMassWTCharm[iset]->Write();
		h2dRecoPtMassWTBottom[iset]->Write();

		h2dGenSPDVCharm[iset]->Write();
		h2dGenSPDVBottom[iset]->Write();

		h2dGenV0MVCharm[iset]->Write();
		h2dGenV0MVBottom[iset]->Write();
		h2dGenV0MPCharm[iset]->Write();
		h2dGenV0MPBottom[iset]->Write();

		h2dGenVtxZCharm[iset]->Write();
		h2dGenVtxZBottom[iset]->Write();
		for (int ii=0; ii<5; ii++){
			h2dGenVtxZWTCharm[iset][ii]->Write();
			h2dGenVtxZWTBottom[iset][ii]->Write();
		}

		for (int ii=0; ii<16; ii++){
			hRPM_un[iset][ii]->Write();
		}

		/*
		h3dGenSPDV[iset]->Write();
		for (int ichg=0; ichg<3; ichg++){
			h2dReco_mass_pt[iset][ichg]->Write();
		}
		h2dReco_truth_mass_pt[iset]->Write();
		*/
	}

	hXic0pT_V0DL->Write();
	hXic0pT_XiDL->Write();

}
