
// kstar pT cut not applied

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

void make_hist_rho(const char *fname="file.lst"){

	const float const_pi = TMath::Pi();

	ifstream flist;
	flist.open(fname);

	char ffname[300];

//	bool _IsMB;
	int i_np;
	float f_vtxz;
	int i_p_ch[2000];//, i_p_id[2000];
	float f_TOF_pion[2000], f_TOF_kaon[2000];
	float f_TPC_pion[2000], f_TPC_kaon[2000];
	float f_p_pt[2000], f_p_eta[2000], f_p_phi[2000], f_p_dca_xy[2000], f_p_dca_z[2000];

//	const float ptcut = 2.0; //min
	const float oacut = 30.0; //max
	const float pipiptasymcut = 0.7; //max
	const float rhoKptasymcut_min = -0.2; //max
	const float rhoKptasymcut_max = 0.55; //max
	const float mpipicut_min= 0.6; //min
	const float mpipicut_max = 0.84; //max
	const float mKpicut= 1.1; //max

	const float mass_pion = 0.13957; 
	const float mass_kaon = 0.49367; 
	const float mass_rho = 0.77526;
	const float mass_k1 = 1.270;

	const float eta_max = 0.8;
	const float pt_min = 0.2;
	/*
	const int npt = 6;
	const float ptbin[npt+1] = {1, 2, 4, 6, 8, 12, 24};
*/

	const int npt = 1; 
	const float ptbin[npt+1] = {0, 100};
	const int nvtxz = 10;
	const int nbuffer = 1;

	TH1D *hmult = new TH1D("hmult","",500,0,500);
	TH2D *hmrho_mm = new TH2D("hmrho_mm","",1500,0,1.5,150,0,15);
	TH2D *hmrho_pp = new TH2D("hmrho_pp","",1500,0,1.5,150,0,15);
	TH2D *hmrho_pm = new TH2D("hmrho_pm","",1500,0,1.5,150,0,15);

	TH2D *hmK1_unlike_pmp = new TH2D("hmK1_unlike_pmp","",200,1.0,2.0,150,0,15);
	TH2D *hmK1_unlike_pmm = new TH2D("hmK1_unlike_pmm","",200,1.0,2.0,150,0,15);
	TH2D *hmK1_like_ppp = new TH2D("hmK1_like_ppp","",200,1.0,2.0,150,0,15);
	TH2D *hmK1_like_ppm = new TH2D("hmK1_like_ppm","",200,1.0,2.0,150,0,15);
	TH2D *hmK1_mixed_pmp = new TH2D("hmK1_mixed_pmp","",200,1.0,2.0,150,0,15);
	TH2D *hmK1_mixed_pmm = new TH2D("hmK1_mixed_pmm","",200,1.0,2.0,150,0,15);
	TH2D *hmK1_mixed_ppp = new TH2D("hmK1_mixed_ppp","",200,1.0,2.0,150,0,15);
	TH2D *hmK1_mixed_ppm = new TH2D("hmK1_mixed_ppm","",200,1.0,2.0,150,0,15);

	TH1D *hevent = new TH1D("hevent","",1,0,2);
	
	TH3D *hpt = new TH3D("hpt","pT_hist_corr",500,0,100,500,0,100,500,0,100); 
	TH3D *heta = new TH3D("heta","eta_hist_corr",500,0,100,500,0,100,500,0,100); 

	TH2D *hpt_K1_rhoK1[2];
	TH2D *hpt_K1_kaonK1[2];
	TH2D *hpt_K1_kaonrho[2];
	TH2D *hpt_K1_OA[2];
	TH2D *hpt_K1_ptasym[2];
	TH2D *hpt_K1_pi1Kmass[2];
	TH2D *hpt_K1_pi2Kmass[2];


	for (int ichg=0; ichg<2; ichg++){
		hpt_K1_rhoK1[ichg] = new TH2D(Form("hpt_K1_rhoK1_chg%d",ichg),"",10,0,10,100,0,2);
		hpt_K1_kaonK1[ichg] = new TH2D(Form("hpt_K1_kaonK1_chg%d",ichg),"",10,0,10,100,0,2);
		hpt_K1_kaonrho[ichg] = new TH2D(Form("hpt_K1_kaonrho_chg%d",ichg),"",10,0,10,100,0,2);
		hpt_K1_OA[ichg] = new TH2D(Form("hpt_K1_OA_chg%d",ichg),"",10,0,10,100,0,180);
		hpt_K1_ptasym[ichg] = new TH2D(Form("hpt_K1_ptasym_chg%d",ichg),"",10,0,10,100,-1,1);
		hpt_K1_pi1Kmass[ichg] = new TH2D(Form("hpt_K1_pi1Kmass_chg%d",ichg),"",10,0,10,100,0,2);
		hpt_K1_pi2Kmass[ichg] = new TH2D(Form("hpt_K1_pi2Kmass_chg%d",ichg),"",10,0,10,100,0,2);
	}


	vector<float> vec_pt[nbuffer][nvtxz];
	vector<float> vec_eta[nbuffer][nvtxz];
	vector<float> vec_phi[nbuffer][nvtxz];
	vector<float> vec_ch[nbuffer][nvtxz];

	int n_filled_buffer[nvtxz] = {0};
	int n_filled2_buffer[nvtxz] = {0};


	while ( flist >> ffname ){
	
		cout << "OPEN : " << ffname << endl;

		TFile *infile = new TFile(ffname,"read");

		TTree *T = (TTree*)infile->Get("tree");
		if (!T){
			infile->Close();
			delete infile;

			continue;
		}

//		T->SetBranchAddress("IsMB",&_IsMB);
		T->SetBranchAddress("VertexZ",&f_vtxz);
		T->SetBranchAddress("nTrack",&i_np);
		T->SetBranchAddress("TrackCharge",i_p_ch);
		T->SetBranchAddress("TrackNosTofPion",f_TOF_pion);
		T->SetBranchAddress("TrackNosTofKaon",f_TOF_kaon);
		T->SetBranchAddress("TrackNosTpcPion",f_TPC_pion);
		T->SetBranchAddress("TrackNosTpcKaon",f_TPC_kaon);
		T->SetBranchAddress("TrackEta",f_p_eta);
		T->SetBranchAddress("TrackPhi",f_p_phi);
		T->SetBranchAddress("TrackPt",f_p_pt);
		T->SetBranchAddress("TrackDCAXY",f_p_dca_xy);
		T->SetBranchAddress("TrackDCAZ",f_p_dca_z);

		int nentries = T->GetEntries();

		for (int ien=0; ien<nentries; ien++){
			T->GetEntry(ien);

//			if ( _IsMB==0 ) continue;
			if ( fabs(f_vtxz)>10.0 ) continue;

			int zbin = int((f_vtxz+10.0)/2.0);
			if ( zbin<0 || zbin>=10 ) continue;

			int mult_count = 0;
			for (int ip=0;ip<i_np; ip++){
				
				if ( (f_p_eta[ip]>2.8 && f_p_eta[ip]<5.1) || (f_p_eta[ip]>-3.7 && f_p_eta[ip]<-1.7)  ){
					mult_count++;
				}
			}

			hmult->Fill(mult_count);

			int nfilled = n_filled_buffer[zbin];

			for (int ip=0; ip<i_np; ip++){

				TLorentzVector lvec_pion;

				float nos_sum_pion = sqrt(pow(f_TOF_pion[ip],2)+pow(f_TPC_pion[ip],2));
				float nos_sum_kaon = sqrt(pow(f_TOF_kaon[ip],2)+pow(f_TPC_kaon[ip],2));

				float nos_tpc_pion = fabs(f_TPC_pion[ip]);
				float nos_tof_pion = fabs(f_TOF_pion[ip]);
				float nos_tpc_kaon = fabs(f_TPC_kaon[ip]);
				float nos_tof_kaon = fabs(f_TOF_kaon[ip]);

				if ( ( nos_tof_pion!=-999 ) && (( nos_sum_pion >= 3.0 ) || ( nos_sum_pion > nos_sum_kaon )) ) continue;
				else if ( ( nos_tof_pion==-999 ) && (( nos_tpc_pion >= 2.0 ) || ( nos_tpc_pion > nos_tpc_kaon )) ) continue;

				if ( f_p_pt[ip]<=pt_min || fabs(f_p_eta[ip])>=eta_max ) continue;
				/*
				   if ( fabs(f_p_dca_z[ip])>0.1 ) continue;
				   if ( fabs(f_p_dca_xy[ip])>0.1 ) continue;
				   */
				if ( fabs(f_p_dca_z[ip])>2.0 ) continue;
				if ( fabs(f_p_dca_xy[ip])>7.0*(0.0026+0.0050/f_p_pt[ip]) ) continue;

				lvec_pion.SetPtEtaPhiM(f_p_pt[ip],f_p_eta[ip],f_p_phi[ip],mass_pion);


				for (int jp=0; jp<i_np; jp++){
					if (ip == jp) continue;

					nos_sum_pion = sqrt(pow(f_TOF_pion[jp],2)+pow(f_TPC_pion[jp],2));
					nos_sum_kaon = sqrt(pow(f_TOF_kaon[jp],2)+pow(f_TPC_kaon[jp],2));

					nos_tpc_pion = fabs(f_TPC_pion[jp]);
					nos_tof_pion = fabs(f_TOF_pion[jp]);
					nos_tpc_kaon = fabs(f_TPC_kaon[jp]);
					nos_tof_kaon = fabs(f_TOF_kaon[jp]);

					if ( ( nos_tof_pion!=-999 ) && (( nos_sum_pion >= 3.0 ) || ( nos_sum_pion > nos_sum_kaon )) ) continue;
					else if ( ( nos_tof_pion==-999 ) && (( nos_tpc_pion >= 2.0 ) || ( nos_tpc_pion > nos_tpc_kaon )) ) continue;


					//if ( i_p_ch[ip]==i_p_ch[jp] ) continue;
					if ( f_p_pt[jp]<=pt_min || fabs(f_p_eta[jp])>=eta_max ) continue;
					/*
					   if ( fabs(f_p_dca_z[jp])>0.1 ) continue;
					   if ( fabs(f_p_dca_xy[jp])>0.1 ) continue;
					*/
					if ( fabs(f_p_dca_z[jp])>2.0 ) continue;
					if ( fabs(f_p_dca_xy[jp])>7.0*(0.0026+0.0050/f_p_pt[jp]) ) continue;

					TLorentzVector lvec_pion2;
					lvec_pion2.SetPtEtaPhiM(f_p_pt[jp],f_p_eta[jp],f_p_phi[jp],mass_pion);

					TLorentzVector lvec_rho = lvec_pion + lvec_pion2;

					//		if ( lvec_rho.Pt()<ptcut ) continue;
					
					float pt_asym_pipi = (lvec_pion.Pt() - lvec_pion2.Pt()) / (lvec_pion.Pt() + lvec_pion2.Pt());
					if ( pt_asym_pipi>pipiptasymcut ) continue;

					if ( i_p_ch[ip]>0 && i_p_ch[jp]>0 ){
						hmrho_pp->Fill(lvec_rho.M(), lvec_rho.Pt());
					}else if ( i_p_ch[ip]<0 && i_p_ch[jp]<0 ){
						hmrho_mm->Fill(lvec_rho.M(), lvec_rho.Pt());
					}else{
						hmrho_pm->Fill(lvec_rho.M(), lvec_rho.Pt());
					}

					if ( lvec_rho.M()<0.62 || lvec_rho.M()>0.92 ) continue;

					
					for (int kp=0; kp<i_np; kp++){

						if ( kp == ip || kp == jp ) continue;

						nos_sum_pion = sqrt(pow(f_TOF_pion[kp],2)+pow(f_TPC_pion[kp],2));
						nos_sum_kaon = sqrt(pow(f_TOF_kaon[kp],2)+pow(f_TPC_kaon[kp],2));

						nos_tpc_pion = fabs(f_TPC_pion[kp]);
						nos_tof_pion = fabs(f_TOF_pion[kp]);
						nos_tpc_kaon = fabs(f_TPC_kaon[kp]);
						nos_tof_kaon = fabs(f_TOF_kaon[kp]);

						if ( ( nos_tof_kaon!=-999 ) && (( nos_sum_kaon >= 3.0 ) || ( nos_sum_kaon > nos_sum_pion )) ) continue;
						else if ( ( nos_tof_kaon==-999 ) && (( nos_tpc_kaon >= 2.0 ) || ( nos_tpc_kaon > nos_tpc_pion )) ) continue;

						if ( f_p_pt[kp]<=pt_min || fabs(f_p_eta[kp])>=eta_max ) continue;
						/*
						   if ( fabs(f_p_dca_z[kp])>0.1 ) continue;
						   if ( fabs(f_p_dca_xy[kp])>0.1 ) continue;
						*/
						if ( fabs(f_p_dca_z[kp])>2.0 ) continue;
						if ( fabs(f_p_dca_xy[kp])>7.0*(0.0026+0.0050/f_p_pt[kp]) ) continue;

						TLorentzVector lvec_kaon;
						lvec_kaon.SetPtEtaPhiM(f_p_pt[kp],f_p_eta[kp],f_p_phi[kp],mass_kaon);

						TLorentzVector lvec_K1 = lvec_rho + lvec_kaon;

						float y = lvec_K1.Rapidity();
						float pt = lvec_K1.Pt();

						if ( fabs(y)>=0.5 ) continue;

						int ind_chg = 0;

						TLorentzVector lvec_pi1K = lvec_pion + lvec_kaon; 
						TLorentzVector lvec_pi2K = lvec_pion2 + lvec_kaon; 

						float cosoa = (lvec_rho.Px()*lvec_kaon.Px() + lvec_rho.Py()*lvec_kaon.Py() + lvec_rho.Pz()*lvec_kaon.Pz())/lvec_rho.P()/lvec_kaon.P();
						float oa = acos(cosoa)*180.0/const_pi;
						float pt_asym_rhoK = (lvec_rho.Pt() - lvec_kaon.Pt()) / (lvec_rho.Pt() + lvec_kaon.Pt());

						hpt_K1_rhoK1[ind_chg]->Fill(lvec_K1.Pt(), lvec_rho.Pt()/lvec_K1.Pt());
						hpt_K1_kaonK1[ind_chg]->Fill(lvec_K1.Pt(), lvec_kaon.Pt()/lvec_K1.Pt());
						hpt_K1_kaonrho[ind_chg]->Fill(lvec_K1.Pt(), lvec_kaon.Pt()/lvec_rho.Pt());
						hpt_K1_OA[ind_chg]->Fill(lvec_K1.Pt(), oa);
						hpt_K1_ptasym[ind_chg]->Fill(lvec_K1.Pt(), pt_asym_rhoK);
						hpt_K1_pi1Kmass[ind_chg]->Fill(lvec_K1.Pt(), lvec_pi1K.M());
						hpt_K1_pi2Kmass[ind_chg]->Fill(lvec_K1.Pt(), lvec_pi2K.M());

						if ( pt_asym_rhoK<rhoKptasymcut_min || pt_asym_rhoK>rhoKptasymcut_max ) continue;
						if ( oa>oacut ) continue;
//						if ( lvec_pi1K.M()>mKpicut ) continue;
//						if ( lvec_pi2K.M()>mKpicut ) continue;


						if ( i_p_ch[ip]==i_p_ch[jp] ){
							if ( i_p_ch[ip]==i_p_ch[kp] ){
								hmK1_like_ppp->Fill(lvec_K1.M(), lvec_K1.Pt());
							}else if ( i_p_ch[ip]!=i_p_ch[kp] ){
								hmK1_like_ppm->Fill(lvec_K1.M(), lvec_K1.Pt());
							}
						}else if ( i_p_ch[ip]!=i_p_ch[jp] ){
							if ( i_p_ch[ip]==i_p_ch[kp] ){
								hmK1_unlike_pmp->Fill(lvec_K1.M(), lvec_K1.Pt());
							}else if ( i_p_ch[ip]!=i_p_ch[kp] ){
								hmK1_unlike_pmm->Fill(lvec_K1.M(), lvec_K1.Pt());
							}
						} 

					}//kp

					for (int ibff=0; ibff<nfilled; ibff++){

						for (unsigned int kk=0; kk<vec_pt[ibff][zbin].size(); kk++){

					//		if ( i_p_ch[ip]==vec_ch[ibff][zbin][kk] ) continue;

							TLorentzVector lvec_kaon;
							lvec_kaon.SetPtEtaPhiM(vec_pt[ibff][zbin][kk],vec_eta[ibff][zbin][kk],vec_phi[ibff][zbin][kk],mass_kaon);

							TLorentzVector lvec_K1 = lvec_rho + lvec_kaon;

							float y = lvec_K1.Rapidity();
							float pt = lvec_K1.Pt();

							if ( fabs(y)>=0.5 ) continue;

							float cosoa = (lvec_rho.Px()*lvec_kaon.Px() + lvec_rho.Py()*lvec_kaon.Py() + lvec_rho.Pz()*lvec_kaon.Pz())/lvec_rho.P()/lvec_kaon.P();
							float oa = acos(cosoa)*180.0/const_pi;
							float pt_asym = (lvec_rho.Pt() - lvec_kaon.Pt()) / (lvec_rho.Pt() + lvec_kaon.Pt());

							TLorentzVector lvec_pi1K = lvec_pion + lvec_kaon; 
							TLorentzVector lvec_pi2K = lvec_pion2 + lvec_kaon; 

//							if ( pt_asym<ptasymcut ) continue;
//							if ( oa>oacut ) continue;
					//		if ( lvec_pi1K.M()>mKpicut ) continue;
//							if ( lvec_pi2K.M()>mKpicut ) continue;


							if ( i_p_ch[ip]==i_p_ch[jp] ){
								if ( i_p_ch[ip]==vec_ch[ibff][zbin][kk] ){
									hmK1_mixed_ppp->Fill(lvec_K1.M(), lvec_K1.Pt());
								}else if ( i_p_ch[ip]!=vec_ch[ibff][zbin][kk] ){
									hmK1_mixed_ppm->Fill(lvec_K1.M(), lvec_K1.Pt());
								}
							}else if ( i_p_ch[ip]!=i_p_ch[jp] ){
								if ( i_p_ch[ip]==vec_ch[ibff][zbin][kk] ){
									hmK1_mixed_pmp->Fill(lvec_K1.M(), lvec_K1.Pt());
								}else if ( i_p_ch[ip]!=vec_ch[ibff][zbin][kk] ){
									hmK1_mixed_pmm->Fill(lvec_K1.M(), lvec_K1.Pt());
								}
							} 
						}//kk
					}//ibff
				}//jp
			}//ip


			if ( nfilled<nbuffer ){

				for (int ip=0; ip<i_np; ip++){
					//if ( i_p_id[ip]!=0 ) continue;

					float nos_sum_pion = sqrt(pow(f_TOF_pion[ip],2)+pow(f_TPC_pion[ip],2));
					float nos_sum_kaon = sqrt(pow(f_TOF_kaon[ip],2)+pow(f_TPC_kaon[ip],2));
                          
					float nos_tpc_pion = fabs(f_TPC_pion[ip]);
					float nos_tof_pion = fabs(f_TOF_pion[ip]);
					float nos_tpc_kaon = fabs(f_TPC_kaon[ip]);
					float nos_tof_kaon = fabs(f_TOF_kaon[ip]);

					if ( ( nos_tof_kaon!=-999 ) && (( nos_sum_kaon >= 3.0 ) || ( nos_sum_kaon > nos_sum_pion )) ) continue;
					else if ( ( nos_tof_kaon==-999 ) && (( nos_tpc_kaon >= 2.0 ) || ( nos_tpc_kaon > nos_tpc_pion )) ) continue;


					if ( f_p_pt[ip]<=pt_min || fabs(f_p_eta[ip])>=eta_max ) continue;
					/*
					   if ( fabs(f_p_dca_z[ip])>0.1 ) continue;
					   if ( fabs(f_p_dca_xy[ip])>0.1 ) continue;
					 */
					if ( fabs(f_p_dca_z[ip])>2.0 ) continue;
					if ( fabs(f_p_dca_xy[ip])>7.0*(0.0026+0.0050/f_p_pt[ip]) ) continue;


					vec_pt[nfilled][zbin].push_back(f_p_pt[ip]);
					vec_eta[nfilled][zbin].push_back(f_p_eta[ip]);
					vec_phi[nfilled][zbin].push_back(f_p_phi[ip]);
					vec_ch[nfilled][zbin].push_back(i_p_ch[ip]);
				}

				n_filled_buffer[zbin]++;

			}else{

				int nfilled2 = n_filled2_buffer[zbin];

				vec_pt[nfilled2][zbin].clear();
				vec_eta[nfilled2][zbin].clear();
				vec_phi[nfilled2][zbin].clear();
				vec_ch[nfilled2][zbin].clear();

				for (int ip=0; ip<i_np; ip++){
					//if ( i_p_id[ip]!=0 ) continue;

					float nos_sum_pion = sqrt(pow(f_TOF_pion[ip],2)+pow(f_TPC_pion[ip],2));
					float nos_sum_kaon = sqrt(pow(f_TOF_kaon[ip],2)+pow(f_TPC_kaon[ip],2));
                          
					float nos_tpc_pion = fabs(f_TPC_pion[ip]);
					float nos_tof_pion = fabs(f_TOF_pion[ip]);
					float nos_tpc_kaon = fabs(f_TPC_kaon[ip]);
					float nos_tof_kaon = fabs(f_TOF_kaon[ip]);

					if ( ( nos_tof_kaon!=-999 ) && (( nos_sum_kaon >= 3.0 ) || ( nos_sum_kaon > nos_sum_pion )) ) continue;
					else if ( ( nos_tof_kaon==-999 ) && (( nos_tpc_kaon >= 2.0 ) || ( nos_tpc_kaon > nos_tpc_pion )) ) continue;

					if ( f_p_pt[ip]<=pt_min || fabs(f_p_eta[ip])>=eta_max ) continue;
					/*
					   if ( fabs(f_p_dca_z[ip])>0.1 ) continue;
					   if ( fabs(f_p_dca_xy[ip])>0.1 ) continue;
					 */
					if ( fabs(f_p_dca_z[ip])>2.0 ) continue;
					if ( fabs(f_p_dca_xy[ip])>7.0*(0.0026+0.0050/f_p_pt[ip]) ) continue;

					vec_pt[nfilled2][zbin].push_back(f_p_pt[ip]);
					vec_eta[nfilled2][zbin].push_back(f_p_eta[ip]);
					vec_phi[nfilled2][zbin].push_back(f_p_phi[ip]);
					vec_ch[nfilled2][zbin].push_back(i_p_ch[ip]);

				}

				n_filled2_buffer[zbin]++;

				if (n_filled2_buffer[zbin] == nbuffer) n_filled2_buffer[zbin]=0;
			}

		}//ien

		delete infile;

	}//while

	TFile *outfile = new TFile("outfile_hist.root","recreate");

	hmult->Write();
	hmrho_mm->Write();
	hmrho_pp->Write();
	hmrho_pm->Write();

	hmK1_unlike_pmp->Write();
	hmK1_unlike_pmm->Write();
	hmK1_like_ppp->Write();
	hmK1_like_ppm->Write();

	hmK1_mixed_pmp->Write();
	hmK1_mixed_pmm->Write();
	hmK1_mixed_ppp->Write();
	hmK1_mixed_ppm->Write();

	for (int ichg=0; ichg<2; ichg++){
		hpt_K1_rhoK1[ichg]->Write();
		hpt_K1_kaonK1[ichg]->Write();
		hpt_K1_kaonrho[ichg]->Write();
		hpt_K1_OA[ichg]->Write();
		hpt_K1_ptasym[ichg]->Write();
		hpt_K1_pi1Kmass[ichg]->Write();
		hpt_K1_pi2Kmass[ichg]->Write();
	}
	
	outfile->Close();

}//void
