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

void pT_kstar_tagged_v01_skim(const char *fname="file_kstar.lst"){

	gRandom = new TRandom3(0);

	const int nmult = 1;
	float multcut[nmult+1] = {0, 1000};

	const int npt = 1;
	const float ptbin[npt+1] = {1.0, 2.0};

	const float const_pi = TMath::Pi();

	const float mass_pion = 0.13957; 
	const float mass_kaon = 0.49367; 

	float _kstar_pt;
	float _kstar_pi_pt;
	float _kstar_ka_pt;
	float _kstar_pt_asym;
	float _kstar_oa;
	float _k1_ka_pt;
	float _k1_pt_asym;
	float _k1_oa;
	float _k1_kpi_mass;
	float _k1_pipi_mass;
	float _k1_pt;
	float _k1_rap;
	float _k1_mass;
	float _k1_tag;

	TFile *outfile = new TFile("outfile_tree_pT_tagged.root","recreate");

	TTree *Tsig = new TTree("Tsig","Tsig");
	Tsig->Branch("kstar_pt",&_kstar_pt,"kstar_pt/F");
	Tsig->Branch("kstar_pi_pt",&_kstar_pi_pt,"kstar_pi_pt/F");
	Tsig->Branch("kstar_ka_pt",&_kstar_ka_pt,"kstar_ka_pt/F");
	Tsig->Branch("kstar_pt_asym",&_kstar_pt_asym,"kstar_pt_asym/F");
	Tsig->Branch("kstar_oa",&_kstar_oa,"kstar_oa/F");
	Tsig->Branch("k1_ka_pt",&_k1_ka_pt,"k1_ka_pt/F");
	Tsig->Branch("k1_pt_asym",&_k1_pt_asym,"k1_pt_asym/F");
	Tsig->Branch("k1_oa",&_k1_oa,"k1_oa/F");
	Tsig->Branch("k1_pipi_mass",&_k1_pipi_mass,"k1_pipi_mass/F");
	Tsig->Branch("k1_kpi_mass",&_k1_kpi_mass,"k1_kpi_masss/F");
	Tsig->Branch("k1_pt",&_k1_pt,"k1_pt/F");
	Tsig->Branch("k1_rap",&_k1_rap,"k1_rap/F");
	Tsig->Branch("k1_mass",&_k1_mass,"k1_mass/F");
	Tsig->Branch("k1_tag",&_k1_tag,"k1_tag/F");

	TTree *Tbkg = new TTree("Tbkg","Tbkg");
	Tbkg->Branch("kstar_pt",&_kstar_pt,"kstar_pt/F");
	Tbkg->Branch("kstar_pi_pt",&_kstar_pi_pt,"kstar_pi_pt/F");
	Tbkg->Branch("kstar_ka_pt",&_kstar_ka_pt,"kstar_ka_pt/F");
	Tbkg->Branch("kstar_pt_asym",&_kstar_pt_asym,"kstar_pt_asym/F");
	Tbkg->Branch("kstar_oa",&_kstar_oa,"kstar_oa/F");
	Tbkg->Branch("k1_ka_pt",&_k1_ka_pt,"k1_ka_pt/F");
	Tbkg->Branch("k1_pt_asym",&_k1_pt_asym,"k1_pt_asym/F");
	Tbkg->Branch("k1_oa",&_k1_oa,"k1_oa/F");
	Tbkg->Branch("k1_pipi_mass",&_k1_pipi_mass,"k1_pipi_mass/F");
	Tbkg->Branch("k1_kpi_mass",&_k1_kpi_mass,"k1_kpi_masss/F");
	Tbkg->Branch("k1_pt",&_k1_pt,"k1_pt/F");
	Tbkg->Branch("k1_rap",&_k1_rap,"k1_rap/F");
	Tbkg->Branch("k1_mass",&_k1_mass,"k1_mass/F");
	Tbkg->Branch("k1_tag",&_k1_tag,"k1_tag/F");

	TTree *Tdata = new TTree("Tdata","Tdata");
	Tdata->Branch("kstar_pt",&_kstar_pt,"kstar_pt/F");
	Tdata->Branch("kstar_pi_pt",&_kstar_pi_pt,"kstar_pi_pt/F");
	Tdata->Branch("kstar_ka_pt",&_kstar_ka_pt,"kstar_ka_pt/F");
	Tdata->Branch("kstar_pt_asym",&_kstar_pt_asym,"kstar_pt_asym/F");
	Tdata->Branch("kstar_oa",&_kstar_oa,"kstar_oa/F");
	Tdata->Branch("k1_ka_pt",&_k1_ka_pt,"k1_ka_pt/F");
	Tdata->Branch("k1_pt_asym",&_k1_pt_asym,"k1_pt_asym/F");
	Tdata->Branch("k1_oa",&_k1_oa,"k1_oa/F");
	Tdata->Branch("k1_pipi_mass",&_k1_pipi_mass,"k1_pipi_mass/F");
	Tdata->Branch("k1_kpi_mass",&_k1_kpi_mass,"k1_kpi_masss/F");
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

			//kstar-pi mode
			for (int ip=0; ip<i_np; ip++){
				if ( f_p_pt[ip]<0.15 || abs(f_p_eta[ip])>=0.8 ) continue;
				//if ( f_p_vt[ip]>0.1 ) continue;
				if ( abs(i_p_id[ip])!=211 ) continue; 

				lvec_ii.SetPtEtaPhiM(f_p_pt[ip], f_p_eta[ip], f_p_phi[ip], mass_pion);

				for (int jp=0; jp<i_np; jp++){
					if ( ip==jp ) continue;
					if ( f_p_pt[jp]<0.15 || abs(f_p_eta[jp])>=0.8 ) continue;
					//if ( f_p_vt[jp]>0.1 ) continue;
					if ( abs(i_p_id[jp])!=321 ) continue; 

					lvec_jj.SetPtEtaPhiM(f_p_pt[jp], f_p_eta[jp], f_p_phi[jp], mass_kaon);
					TLorentzVector lvec_kstar = lvec_ii + lvec_jj;

					//K*0 mass selection
					if ( lvec_kstar.M()<0.84 || lvec_kstar.M()>0.94 ) continue;

					float pt_asym_kstar = fabs(lvec_ii.Pt() - lvec_jj.Pt()) / (lvec_ii.Pt() + lvec_jj.Pt());
					float cosoa_kstar = (lvec_ii.Px()*lvec_jj.Px() + lvec_ii.Py()*lvec_jj.Py() + lvec_ii.Pz()*lvec_jj.Pz())/lvec_ii.P()/lvec_jj.P();
					float oa_kstar = acos(cosoa_kstar)*180.0/const_pi;

					for (int kp=0; kp<i_np; kp++){
						if ( kp==ip || kp==jp ) continue;
						if ( f_p_pt[kp]<0.15 || abs(f_p_eta[kp])>=0.8 ) continue;
						//if ( f_p_vt[kp]>0.1 ) continue;
						if ( abs(i_p_id[kp])!=211 ) continue;

						lvec_kk.SetPtEtaPhiM(f_p_pt[kp], f_p_eta[kp], f_p_phi[kp], mass_pion);
						TLorentzVector lvec_k1 = lvec_kstar + lvec_kk;
						TLorentzVector lvec_pipi = lvec_ii + lvec_kk;
						TLorentzVector lvec_kpi = lvec_jj + lvec_kk;

						if ( fabs(lvec_k1.Rapidity())>0.5 ) continue;
						if ( lvec_k1.M()>2.0 ) continue;
						if ( lvec_k1.Pt()<3 || lvec_k1.Pt()>5 ) continue;

						float pt_asym = (lvec_kstar.Pt() - lvec_kk.Pt()) / (lvec_kstar.Pt() + lvec_kk.Pt());
						float cosoa = (lvec_kstar.Px()*lvec_kk.Px() + lvec_kstar.Py()*lvec_kk.Py() + lvec_kstar.Pz()*lvec_kk.Pz())/lvec_kstar.P()/lvec_kk.P();
						float oa = acos(cosoa)*180.0/const_pi;

						_kstar_pt = lvec_kstar.Pt();
						_kstar_pi_pt = lvec_ii.Pt();
						_kstar_ka_pt = lvec_jj.Pt();
						_kstar_pt_asym = pt_asym_kstar;
						_kstar_oa = oa_kstar; 
						_k1_ka_pt = lvec_kk.Pt();
						_k1_pt_asym = pt_asym;
						_k1_oa = oa;
						_k1_pipi_mass = lvec_pipi.M();
						_k1_kpi_mass = lvec_kpi.M();
						_k1_pt = lvec_k1.Pt();
						_k1_rap = lvec_k1.Rapidity();
						_k1_mass = lvec_k1.M();

						if ( i_p_id[ip]*i_p_id[jp]<0 && i_p_id[ip]*i_p_id[kp]<0 ){
							//signal
							if ( abs(i_p_grandma_id[ip])==10323 
									&& abs(i_p_grandma_id[jp])==10323 
									&& abs(i_p_mom_id[ip])==313
									&& abs(i_p_mom_id[jp])==313
									&& abs(i_p_mom_id[kp])==10323 
									&& i_p_mom_index[kp]==i_p_grandma_index[ip]
									&& i_p_mom_index[kp]==i_p_grandma_index[jp]
								 ){
								//_k1_mass = lvec_k1.M() + gRandom->Gaus(0, 0.12);
								_k1_tag = 1.0;
								Tsig->Fill();
								Tdata->Fill();
							}else{
								if ( gRandom->Rndm()<0.99 ){
									_k1_tag = 0.0;
									Tdata->Fill();
								}
							}
						}else if ( i_p_id[ip]*i_p_id[jp]>0 && i_p_id[ip]*i_p_id[kp]<0 ){
							if ( gRandom->Rndm()<0.02 ){
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

}
