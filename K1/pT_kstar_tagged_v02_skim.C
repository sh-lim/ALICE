//iHadron rescattering effect
//kstar, Kstar, phi

#include <iostream>
#include <fstream>
#include <vector>

#include <TFile.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TTree.h>
#include <TProfile.h>
#include <TVector3.h>
#include <TLorentzVector.h>

using namespace std;

void pT_kstar_tagged_v02_skim(const char *fname="file_kstar.lst"){

	gRandom = new TRandom3(0);

	const int nmult = 1;
	float multcut[nmult+1] = {0, 1000};

	const int npt = 1;
	const float ptbin[npt+1] = {1.0, 2.0};

	const float const_pi = TMath::Pi();

	const float mass_pion = 0.13957; 
	const float mass_kaon = 0.49367; 

	float _rho_pt;
	float _rho_pi1_pt;
	float _rho_pi2_pt;
	float _rho_pt_asym;
	float _rho_oa;
	float _k1_ka_pt;
	float _k1_pt_asym;
	float _k1_oa;
	float _k1_pt;
	float _k1_rap;
	float _k1_mass;
	float _k1_tag;

	TFile *outfile = new TFile("outfile_tree_pT_tagged.root","recreate");

	TTree *Tsig = new TTree("Tsig","Tsig");
	Tsig->Branch("rho_pt",&_rho_pt,"rho_pt/F");
	Tsig->Branch("rho_pi1_pt",&_rho_pi1_pt,"rho_pi1_pt/F");
	Tsig->Branch("rho_pi2_pt",&_rho_pi2_pt,"rho_pi2_pt/F");
	Tsig->Branch("rho_pt_asym",&_rho_pt_asym,"rho_pt_asym/F");
	Tsig->Branch("rho_oa",&_rho_oa,"rho_oa/F");
	Tsig->Branch("k1_ka_pt",&_k1_ka_pt,"k1_ka_pt/F");
	Tsig->Branch("k1_pt_asym",&_k1_pt_asym,"k1_pt_asym/F");
	Tsig->Branch("k1_oa",&_k1_oa,"k1_oa/F");
	Tsig->Branch("k1_pt",&_k1_pt,"k1_pt/F");
	Tsig->Branch("k1_rap",&_k1_rap,"k1_rap/F");
	Tsig->Branch("k1_mass",&_k1_mass,"k1_mass/F");
	Tsig->Branch("k1_tag",&_k1_tag,"k1_tag/F");

	TTree *Tbkg = new TTree("Tbkg","Tbkg");
	Tbkg->Branch("rho_pt",&_rho_pt,"rho_pt/F");
	Tbkg->Branch("rho_pi1_pt",&_rho_pi1_pt,"rho_pi1_pt/F");
	Tbkg->Branch("rho_pi2_pt",&_rho_pi2_pt,"rho_pi2_pt/F");
	Tbkg->Branch("rho_pt_asym",&_rho_pt_asym,"rho_pt_asym/F");
	Tbkg->Branch("rho_oa",&_rho_oa,"rho_oa/F");
	Tbkg->Branch("k1_ka_pt",&_k1_ka_pt,"k1_ka_pt/F");
	Tbkg->Branch("k1_pt_asym",&_k1_pt_asym,"k1_pt_asym/F");
	Tbkg->Branch("k1_oa",&_k1_oa,"k1_oa/F");
	Tbkg->Branch("k1_pt",&_k1_pt,"k1_pt/F");
	Tbkg->Branch("k1_rap",&_k1_rap,"k1_rap/F");
	Tbkg->Branch("k1_mass",&_k1_mass,"k1_mass/F");
	Tbkg->Branch("k1_tag",&_k1_tag,"k1_tag/F");

	TTree *Tdata = new TTree("Tdata","Tdata");
	Tdata->Branch("rho_pt",&_rho_pt,"rho_pt/F");
	Tdata->Branch("rho_pi1_pt",&_rho_pi1_pt,"rho_pi1_pt/F");
	Tdata->Branch("rho_pi2_pt",&_rho_pi2_pt,"rho_pi2_pt/F");
	Tdata->Branch("rho_pt_asym",&_rho_pt_asym,"rho_pt_asym/F");
	Tdata->Branch("rho_oa",&_rho_oa,"rho_oa/F");
	Tdata->Branch("k1_ka_pt",&_k1_ka_pt,"k1_ka_pt/F");
	Tdata->Branch("k1_pt_asym",&_k1_pt_asym,"k1_pt_asym/F");
	Tdata->Branch("k1_oa",&_k1_oa,"k1_oa/F");
	Tdata->Branch("k1_pt",&_k1_pt,"k1_pt/F");
	Tdata->Branch("k1_rap",&_k1_rap,"k1_rap/F");
	Tdata->Branch("k1_mass",&_k1_mass,"k1_mass/F");
	Tdata->Branch("k1_tag",&_k1_tag,"k1_tag/F");

	ifstream flist;
	flist.open(fname);

	char ffname[300];

	int i_np;
	int i_p_id[2000];
	float f_p_pt[2000], f_p_eta[2000], f_p_phi[2000];
	float f_p_vt[2000];

	int i_p_mom_id[2000];
	int i_p_mom_index[2000];
	float f_p_mom_mass[2000];

	int i_p_grandma_id[2000];
	int i_p_grandma_index[2000];
	float f_p_grandma_mass[2000];

	while ( flist >> ffname ){

		cout << "OPEN: " << ffname << endl;

		TFile *infile = new TFile(ffname,"read");

		TTree *T = (TTree*)infile->Get("T");
		if (!T){
			infile->Close();
			delete infile;
			continue;
		}

		T->SetBranchAddress("np",&i_np);
		T->SetBranchAddress("p_id",i_p_id);
		T->SetBranchAddress("p_eta",f_p_eta);
		T->SetBranchAddress("p_phi",f_p_phi);
		T->SetBranchAddress("p_pt",f_p_pt);
		T->SetBranchAddress("p_vt",f_p_vt);

		T->SetBranchAddress("p_mom_id",i_p_mom_id);
		T->SetBranchAddress("p_mom_index",i_p_mom_index);
		T->SetBranchAddress("p_mom_mass",f_p_mom_mass);
		T->SetBranchAddress("p_grandma_id",i_p_grandma_id);
		T->SetBranchAddress("p_grandma_index",i_p_grandma_index);
		T->SetBranchAddress("p_grandma_mass",f_p_grandma_mass);

		int nentries = T->GetEntries();

		for (int ien=0; ien<nentries; ien++){
			T->GetEntry(ien);

			if (i_np>=2000) continue;

			TLorentzVector lvec_ii, lvec_jj, lvec_kk;

			//rho-K mode
			for (int ip=0; ip<i_np; ip++){
				if ( f_p_pt[ip]<0.15 || abs(f_p_eta[ip])>=0.8 ) continue;
				//if ( f_p_vt[ip]>0.1 ) continue;
				if ( abs(i_p_id[ip])!=211 ) continue; 

				lvec_ii.SetPtEtaPhiM(f_p_pt[ip], f_p_eta[ip], f_p_phi[ip], mass_pion);

				for (int jp=0; jp<i_np; jp++){
					if ( ip==jp ) continue;
					if ( f_p_pt[jp]<0.15 || abs(f_p_eta[jp])>=0.8 ) continue;
					//if ( f_p_vt[jp]>0.1 ) continue;
					if ( abs(i_p_id[jp])!=211 ) continue; 

					lvec_jj.SetPtEtaPhiM(f_p_pt[jp], f_p_eta[jp], f_p_phi[jp], mass_pion);
					TLorentzVector lvec_rho = lvec_ii + lvec_jj;

					//rho0 mass selection
					if ( lvec_rho.M()<0.62 || lvec_rho.M()>0.92 ) continue;

					float pt_asym_rho = fabs(lvec_ii.Pt() - lvec_jj.Pt()) / (lvec_ii.Pt() + lvec_jj.Pt());
					float cosoa_rho = (lvec_ii.Px()*lvec_jj.Px() + lvec_ii.Py()*lvec_jj.Py() + lvec_ii.Pz()*lvec_jj.Pz())/lvec_ii.P()/lvec_jj.P();
					float oa_rho = acos(cosoa_rho)*180.0/const_pi;

					for (int kp=0; kp<i_np; kp++){
						if ( kp==ip || kp==jp ) continue;
						if ( f_p_pt[kp]<0.15 || abs(f_p_eta[kp])>=0.8 ) continue;
						//if ( f_p_vt[kp]>0.1 ) continue;
						if ( abs(i_p_id[kp])!=321 ) continue;

						lvec_kk.SetPtEtaPhiM(f_p_pt[kp], f_p_eta[kp], f_p_phi[kp], mass_kaon);
						TLorentzVector lvec_k1 = lvec_rho + lvec_kk;

						if ( fabs(lvec_k1.Rapidity())>0.5 ) continue;
						if ( lvec_k1.M()>2.0 ) continue;
						if ( lvec_k1.Pt()<3 || lvec_k1.Pt()>5 ) continue;

						float pt_asym = (lvec_rho.Pt() - lvec_kk.Pt()) / (lvec_rho.Pt() + lvec_kk.Pt());
						float cosoa = (lvec_rho.Px()*lvec_kk.Px() + lvec_rho.Py()*lvec_kk.Py() + lvec_rho.Pz()*lvec_kk.Pz())/lvec_rho.P()/lvec_kk.P();
						float oa = acos(cosoa)*180.0/const_pi;

						_rho_pt = lvec_rho.Pt();
						_rho_pi1_pt = lvec_ii.Pt();
						_rho_pi2_pt = lvec_jj.Pt();
						_rho_pt_asym = pt_asym_rho;
						_rho_oa = oa_rho; 
						_k1_ka_pt = lvec_kk.Pt();
						_k1_pt_asym = pt_asym;
						_k1_oa = oa;
						_k1_pt = lvec_k1.Pt();
						_k1_rap = lvec_k1.Rapidity();
						_k1_mass = lvec_k1.M();

						if ( i_p_id[ip]*i_p_id[jp]<0 ){
							//signal
							if ( abs(i_p_grandma_id[ip])==10323 
									&& abs(i_p_grandma_id[jp])==10323 
									&& abs(i_p_mom_id[ip])==113
									&& abs(i_p_mom_id[jp])==113
									&& abs(i_p_mom_id[kp])==10323 
									&& i_p_mom_index[kp]==i_p_grandma_index[ip]
									&& i_p_mom_index[kp]==i_p_grandma_index[jp]
								 ){
								_k1_mass = lvec_k1.M() + gRandom->Gaus(0, 0.12);
								_k1_tag = 1.0;
								Tsig->Fill();
								Tdata->Fill();
							}else{
								if ( gRandom->Rndm()<0.5 ){
									_k1_tag = 0.0;
									Tdata->Fill();
								}
							}
						}else{
							if ( gRandom->Rndm()<0.005 ){
								_k1_tag = 0.0;
								Tbkg->Fill();
							}
						}

					}//kp
				}//jp
			}//ip

		}//ien

		delete infile;

	}//

	outfile->cd();
	Tsig->Write();
	Tbkg->Write();
	Tdata->Write();

	outfile->Close();

	return;

}
