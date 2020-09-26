#include "TH2F.h"
#include "TProfile.h"
#include "TFile.h"
#include "TMatrixTSym.h"
#include "TMath.h"
#include "TVector.h"
#include "TLorentzVector.h"
#include "TStyle.h"
#include "TBenchmark.h"
#include "TClonesArray.h"
#include "TSystem.h"
#include "TSQLServer.h"
#include "TSQLResult.h"
#include "TSQLRow.h"

#include <iostream>
#include <string>

using namespace std;

#include "CmdReader.h"
#include "CmdBGOGeometry.h"
#include "CmdRawHit.h"
#include "CmdDCTrack.h"
#include "CmdDCVertex.h"
#include "CmdZCOuterCluster.h"
#include "CmdZCOuterSector.h"
#include "CmdZCSectorHit.h"
#include "CmdZCStripHit.h"
#include "CmdLXeTrack.h"
#include "CmdLXeStripHit.h"
#include "CmdLXeStripCluster.h"
#include "CmdLXeStripClusterExt.h"
#include "CmdLXeCluster.h"
#include "CmdLXeClusterExt.h"
#include "CmdLXeTrack.h"
#include "CmdLXeGeometry.h"
#include "CmdLXeRecoTools.h"
#include "CmdLXeUtils.h"

#include "CmdFullCluster.h"
#include "CmdCsICluster.h"
#include "CmdBGOCluster.h"
#include "CmdMuClbrHit.h"
#include "CmdAnTClbrHit.h"
#include "CmdTrgRefpointRawHit.h"
#include "CmdTrgDCSCRawHit.h"
#include "CmdTrgClusterRawHit.h"
#include "CmdTrgDecision.h"
#include "fitvert.c"

#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IteratorRange.h"


using namespace cmd3;

namespace tr_ph_nm{

  TLorentzVector convert2ROOT(const HepMC::FourVector& in){ 
    TLorentzVector T(in.x(),in.y(),in.z(),in.t());
    return T;
  }

  //
  // Main storage of event data
  //
  const unsigned int MAX_VERTEX = 5;
  const unsigned int MAX_TRACK = 10;
  const unsigned int MAX_LXETRACK = 10;
  const unsigned int MAX_CLUSTER = 50;
  const unsigned int MAX_ZCS = 20;
  const unsigned int MAX_ZCC = 20;
  const unsigned int MAX_ANT = 175;
  const unsigned int MAX_MU = 48;
  const unsigned int MAX_KS = 20;
  const unsigned int MAX_SIM = 100;
  const unsigned int MAX_GOLD = 3000;
  const unsigned int MAX_CORR = 3000;

  typedef struct {

    // Main flags
    Float_t time;
    Float_t ebeam;
    Float_t emeas;
    Float_t demeas;
    Float_t emeas0;
    Float_t demeas0;
    Float_t xbeam;
    Float_t ybeam;
    int runnum;
    int finalstate_id;
    int evnum;
    int trigbits;
    int trigmchs;
    Float_t trigtime;
    Float_t dcfittime;
    Float_t anttime;
    Float_t mutime;
    int is_coll;
    int is_bhabha;
    Float_t ecaltot;
    Float_t ecalneu;
    Float_t z0;
    Float_t psumch; //vector sum of charged track momenta 
    Float_t psumnu;//vector sum of photon momenta
    Float_t lumoff; // offline luminosity from db
    Float_t lumofferr; //stat error offline luminosity from db

    // Vertex data
    int nv_total;
    int nv;
    int vtrk[MAX_VERTEX];// number of tracks in vertex
    int vind[MAX_VERTEX][MAX_TRACK];//indeces of tracks from vertex
    Float_t vchi[MAX_VERTEX];
    Float_t vxyz[MAX_VERTEX][3];

    // KS vertex data
    int nks_total;
    int nks;
    int ksvind[MAX_KS][MAX_TRACK];//indeces of tracks from vertex
    int kstype[MAX_KS];
    int ksfstatus[MAX_KS]; // Status of MINUIT fit     //new1
    Float_t ksvchi[MAX_KS];
    Float_t ksvxyz[MAX_KS][3];
    Float_t ksminv[MAX_KS];
    Float_t ksalign[MAX_KS];
    Float_t ksdpsi[MAX_KS];
    Float_t kstlen[MAX_KS];
    Float_t kslen[MAX_KS];
    Float_t ksz0[MAX_KS];
    Float_t ksphi[MAX_KS];
    Float_t ksth[MAX_KS];
    Float_t ksptot[MAX_KS];
    Float_t kspiphi[MAX_KS][2];
    Float_t kspith[MAX_KS][2];
    Float_t kspipt[MAX_KS][2];
    //Float_t kserr[MAX_KS][3][3];

    // Tracks data
    int nt_total;
    int nt;
    Int_t it[2];
    Int_t   tnhit[MAX_TRACK];
    Float_t tlength[MAX_TRACK];
    Float_t tphi[MAX_TRACK];
    Float_t tth[MAX_TRACK];
    Float_t tptot[MAX_TRACK];
    Float_t tphiv[MAX_TRACK];
    Float_t tthv[MAX_TRACK];
    Float_t tptotv[MAX_TRACK];
    Float_t trho[MAX_TRACK];
    Float_t tdedx[MAX_TRACK];
    Float_t tz[MAX_TRACK];
    Float_t tt0[MAX_TRACK];
    Float_t tant[MAX_TRACK];
    Float_t tchi2r[MAX_TRACK];
    Float_t tchi2z[MAX_TRACK];
    Float_t tchi2ndf[MAX_TRACK];
    Int_t   tcharge[MAX_TRACK];
    Float_t ten[MAX_TRACK];
    Float_t tenlxe[MAX_TRACK];
    Float_t tlengthlxe[MAX_TRACK];
    //Float_t tenslxe[MAX_TRACK]; // by strips == sun of  tenslxe_layers
    Float_t tenslxe_layers[MAX_TRACK][14]; // by strips by layers
    
    Float_t tencsi[MAX_TRACK];
    Float_t tenbgo[MAX_TRACK];
    Float_t tclth[MAX_TRACK];//theta and phi of connected cluster
    Float_t tclphi[MAX_TRACK];
    Float_t terr[MAX_TRACK][3][3];
    Float_t terr0[MAX_TRACK][5][5];

    Int_t tindlxe[MAX_TRACK]; // reference to LXe track number, == -1 if LXe track not found

    Float_t tfc[MAX_TRACK];// size (number of crystals) of the clusters connected to track
    Float_t tzcc[MAX_TRACK][2];  // track point in zc
    Float_t txyzatcl[MAX_TRACK][3]; // track point at cluster conversion point
    Float_t txyzatlxe[MAX_TRACK][3]; // point of LXE cluster
    Int_t   tenconv[MAX_TRACK];// 1 - if LXe has strip info, 0 - if no LXe strip info

    // LXe tracks
    Int_t ntlxe_total; // full number of tracks
    Int_t ntlxe; // num of lxe tracks under data below
    Int_t tlxenhit[MAX_LXETRACK]; // num of points belong to track
    Float_t tlxelength[MAX_LXETRACK]; // len in cm
    Float_t tlxededx[MAX_LXETRACK]; // effective DeDx value of track in LXe, MeV/cm
    Int_t ntlxelayers[MAX_LXETRACK]; // number of LXe layers passed by track

    Float_t tlxeir[MAX_LXETRACK];     // lxe track: start point. R coordinate
    Float_t tlxeitheta[MAX_LXETRACK]; // lxe track: start point. theta coordinate
    Float_t tlxeiphi[MAX_LXETRACK];   // lxe track: start point. Phi coordinate

    Float_t tlxevtheta[MAX_LXETRACK];   // lxe track: svector Phi coordinate
    Float_t tlxevphi[MAX_LXETRACK];   // lxe track: svector Phi coordinate
    Float_t tlxechi2[MAX_LXETRACK];   // lxe track: chi2 value of track approximation
    Float_t tlxesen[MAX_LXETRACK];   // Full strips energy of lxe track
    Float_t tlxesen_layers[MAX_LXETRACK][14]; // by strips by layers
    
    //  Int_t itlxe[MAX_LXETRACK]; // reference to DC track number, == -1 if DC track not found

    // LXE Calorimeter data
    Int_t nlxe_total;
    Int_t nlxe;
    Float_t lxeen[MAX_CLUSTER];
    Float_t lxeth[MAX_CLUSTER];
    Float_t lxephi[MAX_CLUSTER];
    Float_t lxeentot[MAX_CLUSTER];
    Int_t   lxenzcs[MAX_CLUSTER];
    Int_t   lxeflag[MAX_CLUSTER];


    // CSI Calorimeter data
    int ncsi_total;
    int ncsi;
    Float_t csien[MAX_CLUSTER];
    Float_t csith[MAX_CLUSTER];
    Float_t csiphi[MAX_CLUSTER];
    Int_t   csiflag[MAX_CLUSTER];

    // BGO Calorimeter data
    int nbgo_total;
    int nbgo;
    Float_t bgoen[MAX_CLUSTER];
    Float_t bgoth[MAX_CLUSTER];
    Float_t bgophi[MAX_CLUSTER];
    Int_t   bgoflag[MAX_CLUSTER];
	
    // Full clusters data
    int nfc_total;
    int nfc;
    Float_t fcen[MAX_CLUSTER];
    Float_t fcth[MAX_CLUSTER];
    Float_t fcphi[MAX_CLUSTER];
    Int_t   fcconv[MAX_CLUSTER];
    Float_t fclxe[MAX_CLUSTER];
    Float_t fccsi[MAX_CLUSTER];
    Float_t fcbgo[MAX_CLUSTER];
    Int_t   fctype[MAX_CLUSTER];
    Int_t   fcflag[MAX_CLUSTER];
  
    // Photons data
    int nph_total;
    int nph;
    // photon parameters (corrected)
    Float_t phen[MAX_CLUSTER];
    Float_t phth[MAX_CLUSTER];
    Float_t phphi[MAX_CLUSTER];
    Float_t phrho[MAX_CLUSTER];
    // resolution
    Float_t pherr[MAX_CLUSTER][3];
    // measured values (before correction)
    Float_t phen0[MAX_CLUSTER]; 
    Float_t phth0[MAX_CLUSTER]; 
    Float_t phphi0[MAX_CLUSTER];
    Float_t phrho0[MAX_CLUSTER];
    // energy deposition by subsistems
    Float_t phlxe[MAX_CLUSTER];
    Float_t phcsi[MAX_CLUSTER];
    Float_t phbgo[MAX_CLUSTER];
    // clusters type (1 - LXe,  2 - CsI, 3 - BGO based cluster)
    Int_t phflag[MAX_CLUSTER];
    Int_t   phnzcs[MAX_CLUSTER];//number of Z  chamber sectors
    // Conversion point in LXe calorimeter
    Int_t phconv[MAX_CLUSTER];// 1 - if LXe has strip info, 0 - if no LXe strip info
    //Float_t phslxe[MAX_CLUSTER]; // by strips
    Float_t phslxe_layers[MAX_CLUSTER][14]; // by strips by layers
    Int_t phfc[MAX_CLUSTER];//number of calorimeter elements in cluster (cluster size).
    
    // ZC sector
    int nzcs_total;
    int nzcs;
    Int_t zcsch[MAX_ZCS];
    Int_t zcsstat[MAX_ZCS];
    Float_t zcsphi[MAX_ZCS];
    Float_t zcsamp[MAX_ZCS];
    Float_t zcstime[MAX_ZCS];

    // ZC strips
    int nzcc_total;
    int nzcc;
    Int_t zccl[MAX_ZCC];
    Int_t zccns[MAX_ZCC];
    Float_t zccz[MAX_ZCC];
    Float_t zccamp[MAX_ZCC];
    Int_t zcct[MAX_ZCC];
    Int_t zccvalid[MAX_ZCC];

    // Ant conters
    int nant;
    int antch[MAX_ANT];
    Float_t antt0[MAX_ANT];
    Float_t antt1[MAX_ANT];
    Float_t anta0[MAX_ANT];
    Float_t anta1[MAX_ANT];
    int antst[MAX_ANT];

    //Mu system
    int nmu;
    int much[MAX_MU];
    Float_t mut0[MAX_MU];
    Float_t mut1[MAX_MU];
    Float_t mut2[MAX_MU];
    Float_t mut3[MAX_MU];
    Float_t mua0[MAX_MU];
    Float_t mua1[MAX_MU];
    Float_t mua2[MAX_MU];
    Float_t mua3[MAX_MU];
    int must[MAX_MU];

    //Simulation data

    int nsim;
    int simtype[MAX_SIM];// particle ID from GEANT
    int simorig[MAX_SIM];// Origin of particle. 0 - from e+e- or particle ID from GEANT  
    Float_t simmom[MAX_SIM];
    Float_t simphi[MAX_SIM];
    Float_t simtheta[MAX_SIM];
    Float_t simvtx[MAX_SIM];// vertex X for origin
    Float_t simvty[MAX_SIM];//vertex Y
    Float_t simvtz[MAX_SIM];//vertex Z
    
    //Correction data
    int ncorr;
    int idcorr[MAX_CORR];
    int bitcorr[MAX_CORR];
    
    //Corrupted banks data
    int nbadbank;
    int nbadbankg;
    int nbadbanks[MAX_GOLD];
    int nlostbanks;
    int ncorruptedbanks;

  } tree_data;

  //
  // Some global parameters
  //
  vector<double> xyzbeam;
  double Ebeam;
  vector<double> Emeas;
  vector<double> Emeas0;
  double Bfield;
  double SigmaBeam = 0.007;
  double mKs = 497.67;
  double SigmaKsVert = 100000.0;
  double kswindow = 80.;
  double eeangle = 1.;
  double kschi2cut = 50.;
  double ChiTrackMax = 100.;
  double RhoTrackMax = 6.0;// in cm, skip far tracks from tree

  bool tracktrue[MAX_TRACK];
  std::map<cmd3::CmdDCTrack*,int> track2id;
  std::map<int,cmd3::CmdDCTrack*> id2track;

  //
  //
  // Some control histograms
  //
  TH2D *h_dphi_lxe;
  TH2D *h_dphi_lxe0;
  TH2D *h_dphi_csi;
  TH2D *h_dphi_bgo;

  TH2D *h_dphi_lxe1;
  TH2D *h_dphi_lxe01;
  TH2D *h_dphi_csi1;
  TH2D *h_dphi_bgo1;

  TH2D *h_dphi_tlxe;
  TH2D *h_dr_tlxe;

  TH1D *h_dphi_cal;
  TH1D *h_dphi_ant;
  TH1D *h_ksminv;
  TH2D *h_lxe_csi_ph;

  //
  // Clean event data storage
  //
  void clear_event(tree_data &data) {
    data.is_coll=0;
    data.is_bhabha=0;
    data.trigbits=0;
    data.trigmchs=0;
    data.trigtime=0;
    data.time = 0.;
    data.dcfittime=0;
    data.anttime=0;
    data.mutime=0;
    data.evnum=0;
    data.finalstate_id=0;
    data.nt_total=0;
    data.nt=0;
    data.nv_total=0;
    data.nv=0;
    data.nks_total=0;
    data.nks=0;
    data.ecaltot=0;
    data.ecalneu=0;
    data.it[0]=data.it[1]=0;
    data.nlxe_total=0;
    data.ncsi_total=0;
    data.nbgo_total=0;
    data.nfc_total=0;
    data.nlxe=0;
    data.ncsi=0;
    data.nbgo=0;
    data.nfc=0;
    data.nph=0;
    data.nph_total=0;
    data.nzcs_total=0;
    data.nzcs=0;
    data.nzcc_total=0;
    data.nzcc=0;
    data.nant=0;
    data.nmu=0;
    data.psumch=0;
    data.psumnu=0;

    data.ntlxe_total = 0;
    data.ntlxe = 0;

    data.nsim = 0;
    
    data.ncorr = 0;
    data.nbadbank = 0;
    data.nbadbankg = 0;
    data.ncorruptedbanks = 0;
    data.nlostbanks = 0;

  }
  //
  // Fill DC data
  //
  void fill_tracks_data(tree_data &data, const cmd3::CmdEvent* event) {

    const double vlight = 2.99792458;

    track2id.clear();
    id2track.clear();
    //
    //call all needed collections
    //
    const CmdList<CmdDCTrack*> *DCtracks = event->Get<CmdList<CmdDCTrack*> >("dc_tracks");
    const CmdVector<CmdBGOCluster*> *clustersBGO = event->Get<const CmdVector<CmdBGOCluster*> >("bgo_clusters");
    const CmdVector<CmdLXeCluster> * lxe_towers_clusters = event->Get<CmdVector<CmdLXeCluster> >("lxe_clusters");

    const CmdVector<CmdVector<CmdVector<CmdLXeStripCluster> > > * lxe_strips_clusters = event->Get<CmdVector<CmdVector<CmdVector<CmdLXeStripCluster> > > >("lxe_strips_clusters");


    const CmdVector<CmdVector<CmdVector<CmdLXeStripHit> > > * lxe_strips_hits = event->Get<CmdVector<CmdVector<CmdVector<CmdLXeStripHit> > > >("lxe_strips_hits");

    const CmdVector<CmdCsICluster> *clusters = event->Get<const CmdVector<CmdCsICluster> >("csi_clusters");
	
    const cmd3::CmdVector<cmd3::CmdFullCluster> *full_clusters = event->Get<const cmd3::CmdVector<cmd3::CmdFullCluster> >("full_clusters");

    //vector pointers to tracks
    vector<CmdDCTrack*> trks;

    //
    // zero cluster flags
    //
    for (unsigned int t = 0; t < MAX_CLUSTER; ++t) {
      data.bgoflag[t] = 0;
      data.csiflag[t] = 0;
      data.lxeflag[t] = 0;
      data.fcflag[t] = 0;
      data.phflag[t] = 0;
      data.phfc[t] = 0;
    }

    if ( !DCtracks ) return;

    data.nt_total = DCtracks->GetSize();

    if (DCtracks->GetSize() > (int) MAX_TRACK)
      return;

    double zav = 0;
    int itrack=0;
    int itr0=0;

    for(CmdList<CmdDCTrack*>::const_iterator j1=DCtracks->begin();j1!=DCtracks->end(); j1++) {
      CmdDCTrack *trk = (*j1) ;

      //clear pointer vector and fill with last track
      trks.clear();

      //trk->GetInfo();
      //cout<<trk->GetNHits()<<endl;
      if(trk->GetNHits()<5)continue;
      itr0++;
      tracktrue[itr0-1]=true;

      // Get track parameters
      double xpca,ypca,zpca, k, ctg, rho ,t0, phi;
      trk->TrkParams(xyzbeam[0],xyzbeam[1],k, phi, rho, ctg,xpca,ypca,zpca, t0);

      double g = sqrt(1 + ctg*ctg);
      double theta = (ctg>0)?(TMath::ASin(1/TMath::Sqrt(1+ctg*ctg))):
	(M_PI - TMath::ASin(1/TMath::Sqrt(1+ctg*ctg)));
      double p_trk = fabs(vlight*Bfield/k*g);
      double p_perp = fabs(vlight*Bfield/k);

      double chi2r = trk->GetRPhiChi2()/(trk->GetNHits()-3);
      double chi2z = trk->GetZChi2()/(trk->GetNHits()-2);
      double chi2ndf = trk->GetChi2NDF();
 
      // skip bad reconstructed tracks

      //   if(chi2r>ChiTrackMax || chi2z>ChiTrackMax )continue;

      //TVector3 PCA1(xpca-xyzbeam[0],ypca-xyzbeam[1],zpca);
      //rho = PCA1.Perp();

      //remember "good" tracks ID

      track2id[trk]=itrack;
      id2track[itrack]=trk;

      data.tnhit[itrack] = trk->GetNHits();
      data.tlength[itrack] = trk->GetLength();
                
      // covariant matrix
      const TMatrixDSym Dtmp =trk->GetCovMatrix(xyzbeam[0],xyzbeam[1]);

      data.terr[itrack][0][0]=Dtmp[0][0]*(-p_perp/k)*(-p_perp/k); // Move to MeV units
      data.terr[itrack][0][1]=Dtmp[0][1]*(-p_perp/k)*1.;
      data.terr[itrack][1][0]=data.terr[itrack][0][1];
      data.terr[itrack][0][2]=Dtmp[0][3]*(-p_perp/k)/(-g*g);
      data.terr[itrack][2][0]=data.terr[itrack][0][2];
      data.terr[itrack][1][1]=Dtmp[1][1];
      data.terr[itrack][1][2]=Dtmp[1][3]/(-g*g);
      data.terr[itrack][2][1]=data.terr[itrack][1][2];
      data.terr[itrack][2][2]=Dtmp[3][3]/(g*g*g*g);// move to theta from ctg

      //error matrix curveture moved to MeV
      double to_mev_tmp=1;
      for(int i0=0;i0<5;++i0){
	for(int j0=0;j0<5;++j0){
	  to_mev_tmp=1;
	  if(i0==0) to_mev_tmp=to_mev_tmp*(-p_perp/k);
	  if(j0==0) to_mev_tmp=to_mev_tmp*(-p_perp/k);
	  data.terr0[itrack][i0][j0] = Dtmp[i0][j0]*to_mev_tmp;
	}
      }
    

      // scip bad reconstructed tracks
      //if(chi2r>100. || chi2z>100. )continue;

      // reference to track for fitvert
      trks.push_back(trk);

      //----------------------------------------------------------------
      // fit with beam point

      TVector3 vertfit;
      double chi2fit;
      vector<TLorentzVector> Pmfit;
      int fsta = -1;
      double t0fit;
      //		cout<<"fit trks="<<trks.size()<<endl;
      chi2fit=fitvertex(&trks,fsta,xyzbeam[0],xyzbeam[1],SigmaBeam,SigmaBeam,
		        vertfit,Pmfit,t0fit,
                        Bfield,139.57018);
      //		cout<<"fittetd "<<chi2fit<<endl;

      // cout<<"x="<<xyzbeam[0]<<" y="<<xyzbeam[1]<<" p="<<p_trk<<" pfit="<<Pmfit[0].P()<<" Sigma="<<SigmaBeam<<endl;

      data.tphiv[itrack]  = Pmfit[0].Phi();
      if(Pmfit[0].Phi()<0)data.tphiv[itrack]  = Pmfit[0].Phi()+2*M_PI;
      data.tthv[itrack]   = Pmfit[0].Theta();
      data.tptotv[itrack] = Pmfit[0].P();
      //--------------------------------------------------------------

      zav += zpca;

      data.tphi[itrack] = phi;
      if(phi<0)data.tphi[itrack] = phi+2*M_PI;
      data.tth[itrack] = theta;
      //    data.trho[itrack] = sqrt((xpca-xyzbeam[0])*(xpca-xyzbeam[0])-(ypca-xyzbeam[1])*(ypca-xyzbeam[1]));
      data.trho[itrack] = rho;
      data.tptot[itrack] = p_trk;
      data.tdedx[itrack] = trk->GetdEdX();
      data.tz[itrack] = zpca;
      data.tt0[itrack] = t0;
      if( k>0 ) data.tcharge[itrack]=1; else data.tcharge[itrack]=-1;
      data.tchi2r[itrack] = chi2r;
      data.tchi2z[itrack] = chi2z;
      data.tchi2ndf[itrack] = chi2ndf;

      TLorentzVector vec(p_trk*TMath::Cos(phi)*TMath::Sin(theta),
			 p_trk*TMath::Sin(phi)*TMath::Sin(theta),
			 p_trk*TMath::Cos(theta),
			 TMath::Hypot(p_trk,139.57));

      data.tenlxe[itrack] = 0.0;
      //data.tenslxe[itrack] = 0.0;
      data.tencsi[itrack] = 0.0;
      data.tenbgo[itrack] = 0.0;
      data.ten[itrack] = 0.0;
      data.tfc[itrack] = 0.0;
      data.tlengthlxe[itrack] = 0;
      data.txyzatcl[itrack][0] = 150.0;
      data.txyzatcl[itrack][1] = 150.0;
      data.txyzatcl[itrack][2] = 150.0;
      data.txyzatlxe[itrack][0] = 150.0;
      data.txyzatlxe[itrack][1] = 150.0;
      data.txyzatlxe[itrack][2] = 150.0;
      data.tenconv[itrack] = 0;
      
      //------zc radius z--------------- 
      std::pair<std::pair<TVector3,TVector3>,std::pair<TVector3,TVector3> > PointZ_I;
      TVector3 Position_0;
      //data.tzcc[itrack][0]=100.0;
      //data.tzcc[itrack][1]=100.0;

      PointZ_I=trk->GetPointDirectionAtR(0.5*62.07);
      Position_0=(vec.Angle(PointZ_I.first.first)<vec.Angle(PointZ_I.second.first))?PointZ_I.first.first:PointZ_I.second.first;
      data.tzcc[itrack][0]=Position_0.Z();
      //data.tpos_zcin[itrack]=Position_0.Z();

      PointZ_I=trk->GetPointDirectionAtR(0.5*63.87);
      Position_0=(vec.Angle(PointZ_I.first.first)<vec.Angle(PointZ_I.second.first))?PointZ_I.first.first:PointZ_I.second.first;
      data.tzcc[itrack][1]=Position_0.Z();
      //data.tpos_zcout[itrack]=Position_0.Z();
      //------zc radius z----end----------

      // --- LXe clusters for LXe Z-scale studies ---- ---- -----------------
      const cmd3::CmdVector<cmd3::CmdLXeCluster> * lxe_clusters1 = event->Get<CmdVector<CmdLXeCluster> >("lxe_clusters");
      if (lxe_clusters1)
      for (unsigned int t = 0; t < lxe_clusters1->GetSize(); ++t) {
	const cmd3::CmdLXeCluster & lclust = (*lxe_clusters1)[t];
	if(lclust.GetEnergy()>3*data.ebeam)continue;
	//if (lclust.GetEnergy()>50)
	double dphi_min2 = 10.; //for Track Point in LXe
	double dth_min2 = 10.;  //for Track Point in LXe
	double phic = lclust.GetPosition().c3;
	double thc  = lclust.GetPosition().c2;
	//double rhoc = lclust.GetPosition().c1 * sin(thc);

	double x11 = CmdGeomObject::GetXYZfromRThetaPhi(lclust.GetPosition()).c1;
  	double y11 = CmdGeomObject::GetXYZfromRThetaPhi(lclust.GetPosition()).c2;
	double rhoc = sqrt(x11*x11+y11*y11); 

	//get track positions at LXe cluster radius	
	std::pair<std::pair<TVector3,TVector3>,std::pair<TVector3,TVector3> 
	    > ptt=trk->GetPointDirectionAtR(rhoc);
	TVector3 postr=(vec.Angle(ptt.first.first)<vec.Angle(ptt.second.first))?ptt.first.first:ptt.second.first;
	
	double dthc  = postr.Theta()-thc; 
	double dphic = postr.Phi()-phic;
	if( dphic > (TMath::Pi()) ) dphic-=(TMath::Pi()*2);
	if( dphic < -(TMath::Pi()) ) dphic+=(TMath::Pi()*2);

	if(fabs(dphic)<0.3 && fabs(dthc)<0.3){
	    if( (fabs(dphic)<dphi_min2)&&(fabs(dthc)<dth_min2) ){
	      dphi_min2 = fabs(dphic);
	      dth_min2 = fabs(dthc); 
		if (lclust.HasConvPoint() ){
			data.txyzatlxe[itrack][0] = (CmdGeomObject::GetXYZfromRThetaPhi(lclust.GetPosition()).c1)/10;
	      		data.txyzatlxe[itrack][1] = (CmdGeomObject::GetXYZfromRThetaPhi(lclust.GetPosition()).c2)/10;
	      		data.txyzatlxe[itrack][2] = (CmdGeomObject::GetXYZfromRThetaPhi(lclust.GetPosition()).c3)/10;
		}
		else {data.txyzatlxe[itrack][0] = 222;data.txyzatlxe[itrack][1] = 222;data.txyzatlxe[itrack][2] = 222;}
	    }
	}  
      }
      // --- END of LXe Z-scale studies paramrters --------------------------
      
      //////////////////////////////// calculation of tlengthlxe
      double smin,smax,rho1,z1,rho2,z2;
      trk->GetLengths(smin,smax);
      TVector3 ip,fp,p1,p2;
      trk->GetPosition(smin,p1);trk->GetPosition(smax,p2);
      rho1=abs(p1.Perp());z1=p1.Z();
      rho2=abs(p2.Perp());z2=p2.Z();
      if (rho1<=rho2){ip=p1;fp=p2;}else{ip=p2;fp=p1;}
      /////////////////////////////////
      TVector3 p_lxe_in,p_lxe_out;
      std::pair<std::pair<TVector3,TVector3>,std::pair<TVector3,TVector3> > ptt=trk->GetPointDirectionAtR(37.85-1.02);
      TVector3 postr,postr1 = ptt.first.first,postr2 = ptt.second.first;      
      if (abs(trk->GetPosAlongTrack(postr1)-trk->GetPosAlongTrack(fp))<abs(trk->GetPosAlongTrack(postr2)-trk->GetPosAlongTrack(fp))){postr=postr1;}else{postr=postr2;}
      p_lxe_in=postr;
      
      ptt=trk->GetPointDirectionAtR(50.15+1.02);
      postr1 = ptt.first.first;
      postr2 = ptt.second.first;
      if (abs(trk->GetPosAlongTrack(postr1)-trk->GetPosAlongTrack(fp))<abs(trk->GetPosAlongTrack(postr2)-trk->GetPosAlongTrack(fp))){postr=postr1;}else{postr=postr2;}
      p_lxe_out=postr;      
      
      data.tlengthlxe[itrack] = (p_lxe_in-p_lxe_out).Mag();
      /////////////////////////////////
      
      for (unsigned int i=0; i<14; i++) {data.tenslxe_layers[itrack][i] = 0;}
      
      if ( full_clusters ) {
	int ngfc = 0;

	double phi_min = 10.;
	double th_min = 10.;
	double dphi_min = 10.;
	double dth_min = 10.;

	double dphi_min2 = 10.; //for Track Point in LXe
	double dth_min2 = 10.;  //for Track Point in LXe

	for (unsigned int t = 0; t < full_clusters->GetSize(); ++t) {
	  const cmd3::CmdFullCluster & tclust = (*full_clusters)[t];
	
	  if(tclust.GetEnergy()>3*data.ebeam)continue;
	
	  //double phic = tclust.GetPosition().c3;
	  //double thc  = tclust.GetPosition().c2;

	  double phic = tclust.getGammaPhi();
	  double thc  = tclust.getGammaTheta();

	  //double rhoc = tclust.GetPosition().c1 * sin(thc);
	  double rhoc = tclust.getRho(); 
	
	  //get track positions at LXe cluster radius	
	  std::pair<std::pair<TVector3,TVector3>,std::pair<TVector3,TVector3> 
	    > ptt=trk->GetPointDirectionAtR(rhoc);
	
	  TVector3 postr=(vec.Angle(ptt.first.first)<vec.Angle(ptt.second.first))?ptt.first.first:ptt.second.first;
	  //if(fabs(postr.Pt())<30) continue;
	
	  double dthc  = postr.Theta()-thc; 
	  double dphic = postr.Phi()-phic;
	
	  if( dphic > (TMath::Pi()) ) dphic-=(TMath::Pi()*2);
	  if( dphic < -(TMath::Pi()) ) dphic+=(TMath::Pi()*2);
	
	
	  //double cenlxe = tclust.getLXeEnergy();
	  //double cencsi = tclust.getCsIEnergy();
	  //double cenbgo = tclust.getBGOEnergy();
	
	  int ctype = tclust.getCType();
	  //chech where is more deposit energy and change ctype
	  if(ctype == 3 && tclust.getLXeEnergy()>tclust.getCsIEnergy()) ctype =1;
	  if(ctype == 1 && tclust.getLXeEnergy()<tclust.getCsIEnergy()) ctype =3;


	  int cconv = tclust.haveConvPoint();

	  if(ctype == 1 && cconv == 1)h_dphi_lxe->Fill(dphic,dthc);
	  if(ctype == 2)h_dphi_csi->Fill(dphic,dthc);
	  if(ctype == 3)h_dphi_bgo->Fill(dphic,dthc);
	  if(ctype == 1 && cconv == 0)h_dphi_lxe0->Fill(dphic,dthc);

	  // use errors to select connected cluster

	  //double ClThEr = tclust.getThetaError();
	  //double ClPhiEr = tclust.getPhiError();

	  //double deltaTheta = 10.*ClThEr+0.2;
	  //double deltaPhi   = 10.*ClPhiEr+0.2;

	  //if(deltaTheta>0.3)deltaTheta=0.3;
	  //if(deltaPhi>0.3)deltaPhi=0.3;

	  double deltaTheta = 0.2;
	  double deltaPhi   = 0.2;

	  if(ctype == 2){
	    deltaTheta=0.5;
	    deltaPhi=0.5;
	  }

	  if(ctype == 3){
	    deltaTheta=0.15;
	    deltaPhi=0.15;
	  }

	  if(ctype == 1 &&  cconv == 0){
	    deltaTheta=0.3;
	    deltaPhi=0.3;
	  }

	  //	if( fabs(dphic) < deltaPhi && fabs(dthc) < deltaTheta ){


	  //	  if(dphic*dphic/deltaPhi/deltaPhi+dthc*dthc/deltaTheta/deltaTheta<1){
	  if(fabs(dphic)<deltaPhi && fabs(dthc)<deltaTheta){

	    //------Track Point in LXe ------------------------------
	    //-- lxe cluster(Full cluster) rho=rhoc--
	    PointZ_I=trk->GetPointDirectionAtR(rhoc);
	    Position_0=(vec.Angle(PointZ_I.first.first)<vec.Angle(PointZ_I.second.first))?PointZ_I.first.first:PointZ_I.second.first;
	    if( (fabs(dphic)<dphi_min2)&&(fabs(dthc)<dth_min2) ){
	      dphi_min2 = fabs(dphic);
	      dth_min2 = fabs(dthc); 
		if (tclust.haveConvPoint() ){
			data.txyzatcl[itrack][0] = Position_0.X();
	      		data.txyzatcl[itrack][1] = Position_0.Y();
	      		data.txyzatcl[itrack][2] = Position_0.Z();
		}
		else {data.txyzatcl[itrack][0] = 333;data.txyzatcl[itrack][1] = 333;data.txyzatcl[itrack][2] = 333;}
	    }		
	    //------------------------------------------------------

	    if(ctype == 1 && cconv == 1)h_dphi_lxe1->Fill(dphic,dthc);
	    if(ctype == 2)h_dphi_csi1->Fill(dphic,dthc);
	    if(ctype == 3)h_dphi_bgo1->Fill(dphic,dthc);
	    if(ctype == 1 && cconv == 0)h_dphi_lxe01->Fill(dphic,dthc);

	    // remember phi and theta of closest cluster if more than one

	    if(fabs(dphic) < dphi_min){phi_min = phic; dphi_min = fabs(dphic);}
	    if(fabs(dthc) < dth_min){th_min = thc; dth_min = fabs(dthc);} 

	    data.tfc[itrack] += tclust.getNLXe();
	    data.tfc[itrack] += tclust.getNCsI();
	    data.tfc[itrack] += tclust.getNBGO();

	    data.ten[itrack] += tclust.getEnergy();
	    data.tenlxe[itrack] += tclust.getLXeEnergy();
	    data.tencsi[itrack] += tclust.getCsIEnergy();
	    data.tenbgo[itrack] += tclust.getBGOEnergy();
 	    data.tenconv[itrack] = tclust.haveConvPoint();// 1 if LXe strips present, 0 - if no

	    if(ngfc<MAX_CLUSTER)
	      data.fcflag[ngfc] = 1;
	    if(data.tchi2ndf[itrack]>ChiTrackMax || fabs(data.trho[itrack])>RhoTrackMax ) 
		data.fcflag[ngfc] = 0;

	    data.tclphi[itrack] = phi_min;
	    data.tclth[itrack] = th_min;
	    
	    int lxe_index = tclust.getLXeClusterIndex();
	    if (lxe_strips_clusters && lxe_strips_hits && lxe_index > -1) {
	      const cmd3::CmdLXeClusterExt & tclust_lxe = (*lxe_towers_clusters)[lxe_index];
	      const std::set<int> & ind_sclust = tclust_lxe.GetStripClusters();	
	      
	      for (std::set<int>::const_iterator k = ind_sclust.begin(); k != ind_sclust.end(); k++){	    
		const CmdLXeStripClusterExt & sclust = CmdLXeRecoTools::GetStripElem(*k,*lxe_strips_clusters);
		const CmdLXeUtils::PACKED_ID pid = CmdLXeUtils::UnPackValues(sclust.GetIndex());    
		int layer_num=2*(pid.layer)+pid.direction;	    
		double senergy = sclust.GetTotalAmplitude(*lxe_strips_hits);
		data.tenslxe_layers[itrack][layer_num] += senergy;	    
	      }
	    } 
	  }
	  ngfc ++;
	}
      }
      
      // find out proper LXe tracks
      for (int t=0; t<data.ntlxe; t++) {
	double Rotlxe = data.tlxeir[t]*sin(data.tlxeitheta[t]);
	std::pair<std::pair<TVector3,TVector3>,std::pair<TVector3,TVector3> > ptt = trk->GetPointDirectionAtR(Rotlxe);
	TVector3 postr = vec.Angle(ptt.first.first)<vec.Angle(ptt.second.first) ? ptt.first.first : ptt.second.first;

	double dphi = data.tlxeiphi[t] - postr.Phi();
	double dtheta = data.tlxeitheta[t] - postr.Theta();

	if (dphi > TMath::Pi()) dphi  -= TMath::Pi()*2;
	if (dphi < -TMath::Pi()) dphi += TMath::Pi()*2;

	h_dphi_tlxe->Fill(dphi, dtheta);
	CmdPointData plxe = CmdGeomObject::GetXYZfromRThetaPhi(CmdPointData(data.tlxeir[t], data.tlxeitheta[t], data.tlxeiphi[t]));
	double xdist = plxe.GetDistance(CmdPointData(postr.x(), postr.y(), postr.z()));

	h_dr_tlxe->Fill(Rotlxe, xdist);

	//if (fabs(dphi)<0.3 && fabs(dtheta)<0.3) {
	//  data.itlxe[t] = itrack;
	//}

	data.tindlxe[itrack] = -1; // if no LXe track found

	if (xdist < 2.5) {
	  //        data.itlxe[t] = itrack;  // FIX ME LATER to proper LXe-DC linking: find out the closest
	  data.tindlxe[itrack] = t; // DC track reference to LXe track
	}

      }


      // Add ANT data

      const cmd3::CmdVector<cmd3::CmdAnTClbrHit>* AnTHit =
	event->Get<cmd3::CmdVector<cmd3::CmdAnTClbrHit> >("ant_clbr_hits");

      if (AnTHit){

        data.tant[itrack] = 0;
        double tavant = 0;
        int ngood = 0;

	for (unsigned int j=0; j < AnTHit->GetSize(); ++j )
	  {
	    cmd3::CmdAnTClbrHit antHit = (cmd3::CmdAnTClbrHit) AnTHit->Get(j);
	    if ( antHit.GetChannelNumber() != 0 )
	      {
		int ach = antHit.GetChannelNumber();
		double at0 = antHit.GetTime0();
		double at1 = antHit.GetTime1();

		double achphi = ach/16.*6.28;

		double dphi = data.tphi[itrack]-achphi-0.25;

		//cout<<"ach "<<ach<<"achphi "<<achphi<<endl;

		if( dphi >  (TMath::Pi()) ) dphi-=(TMath::Pi()*2);
		if( dphi < -(TMath::Pi()) ) dphi+=(TMath::Pi()*2);
		h_dphi_ant->Fill(dphi);
		if( fabs(dphi)<1.5 ) {
		  if(fabs(at0-at1)*0.14<20){
		    tavant += (at0+at1)/2;
		    ngood++;
		  }
		}
	      }
	  }

	if(ngood>0)data.tant[itrack] = tavant/ngood*0.14;
      }

      itrack++;
    }

    data.nt = itrack;
    data.z0 = zav/itrack;
  }

  //
  // Fill LXE calorimeter data
  //
  void fill_lxe_data(tree_data &data, const cmd3::CmdEvent* event) {

    const cmd3::CmdVector<cmd3::CmdLXeCluster> * lxe_towers_clusters = event->Get<CmdVector<CmdLXeCluster> >("lxe_clusters");
    const CmdVector<CmdVector<CmdVector<CmdLXeStripCluster> > > * lxe_strips_clusters = event->Get<CmdVector<CmdVector<CmdVector<CmdLXeStripCluster> > > >("lxe_strips_clusters");
    const CmdVector<CmdVector<CmdVector<CmdLXeStripHit> > > * lxe_strips_hits = event->Get<CmdVector<CmdVector<CmdVector<CmdLXeStripHit> > > >("lxe_strips_hits");
    const CmdVector<CmdLXeTrack> * lxe_tracks = event->Get<CmdVector<CmdLXeTrack> >("lxe_tracks");

    //const CmdVector<CmdCsICluster> *clusters = event->Get<const CmdVector<CmdCsICluster> >("csi_clusters");
    const CmdVector <CmdZCOuterSector>* Sector = event->Get<CmdVector<CmdZCOuterSector> >("zc_sectors");

    if ( ! lxe_towers_clusters ) return;


    // fill LXe tracks data
    if (lxe_tracks && lxe_strips_hits && lxe_strips_clusters) {
      CmdLXeGeometry *geom = CmdLXeGeometry::GetInstance();
      data.ntlxe_total = lxe_tracks->GetSize();
      data.ntlxe = 0;
      if (data.ntlxe_total <= (int) MAX_LXETRACK) {
	for (unsigned int i=0; i<lxe_tracks->size(); i++) {
	  const CmdLXeTrack & track = (*lxe_tracks)[i];
	  data.tlxenhit[data.ntlxe] = track.GetPoints().size();
	  const CmdLine3D & line = track.GetLine();
	  
	  CmdPointData sp_rthetaphi = CmdGeomObject::GetRThetaPhifromXYZ(line.GetStartPoint());
	  CmdPointData v_rtheta_phi = CmdGeomObject::GetRThetaPhifromXYZ(line.GetSVector());
	  
	  // calculate track length: TODO move to the CmdLXeTrack later
	  // FIX me later
	  CmdPointData cross_romin_xyz, cross_romax_xyz;
	  //std::cout<<"track="<<track<<std::endl;
	  CmdPointData pmin, pmax;
	  double romin = 1e5;
	  double romax = 0;
	  CmdPointData rovals = geom->GetTowerRoGeometry(); // should be fixed to proper track length calulation
	  for (unsigned int ip=0; ip<track.GetPoints().size(); ip++) {
	    CmdPointData pt = track.GetPoints()[ip];
	    double ro = CmdGeomObject::GetRoPhiZfromXYZ(pt).c1;
	    if ( ro <= romin) {
	      romin = ro;
	      pmin = pt;
	    }
	    
	    if (ro>= romax) {
	      romax = ro;
	      pmax = pt;
	    }
	    
	  }
	  
	  data.tlxelength[data.ntlxe] = pmin.GetDistance(pmax) * 0.1; // in cm
	  
	  CmdPointData pmin_rthetaphi = CmdGeomObject::GetRThetaPhifromXYZ(pmin);
	  CmdPointData pmax_rthetaphi = CmdGeomObject::GetRThetaPhifromXYZ(pmax);
	  
	  double dro = fabs((pmin_rthetaphi.c1)*sin(pmin_rthetaphi.c2) - (pmax_rthetaphi.c1)*sin(pmax_rthetaphi.c2));
	  data.ntlxelayers[data.ntlxe] = int(dro/20.4+0.1)+1;
	  
	  data.tlxeir[data.ntlxe] = pmin_rthetaphi.c1 * 0.1; // in cm     // lxe track: start point. R coordinate
	  data.tlxeitheta[data.ntlxe] = pmin_rthetaphi.c2; // lxe track: start point. theta coordinate
	  data.tlxeiphi[data.ntlxe] = pmin_rthetaphi.c3;   // lxe track: start point. Phi coordinate
	  
	  data.tlxevtheta[data.ntlxe] = v_rtheta_phi.c2;   // lxe track: svector theta coordinate
	  data.tlxevphi[data.ntlxe] = v_rtheta_phi.c3;   // lxe track: svector Phi coordinate
	  
	  data.tlxechi2[data.ntlxe] = track.GetChi2();
	  
	  // calc strips energy of track: TODO move to the CmdLXeTrack layer
	  data.tlxesen[data.ntlxe] = 0;
	  for (unsigned int k=0; k<14; k++){data.tlxesen_layers[data.ntlxe][k] = 0;}
	  
	  if (lxe_strips_clusters) {
	    const std::vector<int> & track_clusters = track.GetClusters();
	    for (unsigned int c=0; c<track_clusters.size(); c++) {
	      //data.tlxesen[data.ntlxe] += clust.GetEnergy(); // FIX me later: do upgrade LXeReco
	      const CmdLXeStripClusterExt & clust = CmdLXeRecoTools::GetStripElem(track_clusters[c], *lxe_strips_clusters);
	      const CmdLXeUtils::PACKED_ID pid = CmdLXeUtils::UnPackValues(clust.GetIndex()); 	      
	      int layer_num=2*(pid.layer)+pid.direction;
	      
	      if (lxe_strips_hits) {
		double senergy = clust.GetTotalAmplitude(*lxe_strips_hits);
		data.tlxesen_layers[data.ntlxe][layer_num] = senergy;  
		data.tlxesen[data.ntlxe] += senergy;
	      }
	    }
	  }
	  
	  data.tlxededx[data.ntlxe] = data.tlxesen[data.ntlxe]/data.tlxelength[data.ntlxe];
	  //data.itlxe[data.ntlxe] = -1; // no ref by default
	  data.ntlxe ++;
	}
      }
    }

    data.nlxe_total =  lxe_towers_clusters->GetSize();

    if( data.nlxe_total> (int) MAX_CLUSTER ) return;

    data.nlxe = 0;
    
    for (unsigned int t = 0; t < lxe_towers_clusters->GetSize(); ++t) {

      const CmdLXeClusterExt & tclust = (*lxe_towers_clusters)[t];

      if(tclust.GetEnergy()>3*data.ebeam)continue;

      data.lxeen[data.nlxe]  = tclust.GetEnergy();
      data.lxephi[data.nlxe] = tclust.GetPosition().c3;
      data.lxeth[data.nlxe]  = tclust.GetPosition().c2;
      data.lxeentot[data.nlxe] = tclust.GetEnergy();

      if(data.lxeflag[data.nlxe]==0 && data.nph < (int) MAX_CLUSTER){
	//			data.ecalneu += tclust.GetEnergy();

	//
	// here photon data were filled - removed by V.Shebalin
	// photons are filled in the full_clusters section
	//
      }

      // Add data from ZC sectors
      data.lxenzcs[data.nlxe] = 0;
      if( Sector ) {
	for (unsigned int j=0; j < Sector->GetSize(); ++j ) {
	  cmd3::CmdZCOuterSector SectHit = (cmd3::CmdZCOuterSector) Sector->Get(j);
	  double dphi = data.lxephi[data.nlxe]-SectHit.Get_Phi()+0.2;
	  if( dphi >  (TMath::Pi()) ) dphi-=(TMath::Pi()*2);
	  if( dphi < -(TMath::Pi()) ) dphi+=(TMath::Pi()*2);
	  if( fabs(dphi)<0.5 ) {
	    data.lxenzcs[data.nlxe] += 1;
	  }
	}
      }
      data.nlxe ++;
    }
  }



  //
  // Fill CSI calorimeter data
  //
  void fill_csi_data(tree_data &data, const cmd3::CmdEvent* event) {

    const cmd3::CmdVector<cmd3::CmdCsICluster> *clusters =
      event->Get<const cmd3::CmdVector<cmd3::CmdCsICluster> >("csi_clusters");

    if( ! clusters ) return;

    data.ncsi_total =  clusters->GetSize();

    if( data.ncsi_total> (int) MAX_CLUSTER ) return;

    data.ncsi = 0;

    for (unsigned int t = 0; t < clusters->GetSize(); ++t) {

      const cmd3::CmdCsICluster& cluster = (*clusters)[t];

      if(cluster.GetEnergy()>3*data.ebeam)continue;

      data.csien[data.ncsi]  = cluster.GetEnergy();
      data.csiphi[data.ncsi] = cluster.getPhi();
      data.csith[data.ncsi]  = cluster.getTheta();


      // count neutral energy and form photons from
      // remaining csi clusters
      //

      if(data.csiflag[data.ncsi]==0 && data.nph< (int) MAX_CLUSTER){
	//      data.ecalneu += cluster.GetEnergy();

	//
	// here photon data were filled - removed by V.Shebalin
	// photons are filled in the full_clusters section
	//
      }
      data.ncsi ++;
    }

  }

  //
  // Fill BGO calorimeter data
  //
  void fill_bgo_data(tree_data &data, const cmd3::CmdEvent* event) {

    const cmd3::CmdVector<cmd3::CmdBGOCluster*> *clusters =
      event->Get<const cmd3::CmdVector<cmd3::CmdBGOCluster*> >("bgo_clusters");

    //  const cmd3::CmdVector<cmd3::CmdLXeCluster> * lxe_towers_clusters =
    //    event->Get<cmd3::CmdVector<cmd3::CmdLXeCluster> >("lxe_clusters");

    if( ! clusters ) return;

    data.nbgo_total =  clusters->GetSize();

    if( data.nbgo_total> (int) MAX_CLUSTER ) return;

    data.nbgo = 0;

    for (unsigned int t = 0; t < clusters->GetSize(); ++t) {

      const cmd3::CmdBGOCluster& cluster = *((*clusters)[t]);

      if(cluster.GetEnergy()>3*data.ebeam)continue;

      double phic,thc,rhoc;
      cluster.GetRhoThetaPhi(rhoc,thc,phic);//rhoc in mm!

      data.bgoen[data.nbgo]  = cluster.GetEnergy();
      data.bgophi[data.nbgo] = phic;
      data.bgoth[data.nbgo]  = thc;
      //
      // form photons data
      //
      if(data.bgoflag[data.nbgo]==0 && data.nph< (int) MAX_CLUSTER){
	//      data.ecalneu += cluster.GetEnergy();
	//
	// here photon data were filled - removed by V.Shebalin
	// photons are filled in the full_clusters section
	//
      }
      data.nbgo ++;
    }
  }





  //Fill full_clusters data

  void fill_fullcluster_data(tree_data &data, const cmd3::CmdEvent* event) {
    const cmd3::CmdVector<cmd3::CmdLXeCluster> * lxe_towers_clusters = event->Get<CmdVector<CmdLXeCluster> >("lxe_clusters");
    const CmdVector<CmdVector<CmdVector<CmdLXeStripCluster> > > * lxe_strips_clusters = event->Get<CmdVector<CmdVector<CmdVector<CmdLXeStripCluster> > > >("lxe_strips_clusters");
    const CmdVector<CmdVector<CmdVector<CmdLXeStripHit> > > * lxe_strips_hits = event->Get<CmdVector<CmdVector<CmdVector<CmdLXeStripHit> > > >("lxe_strips_hits");
    //const CmdVector<CmdLXeTrack> * lxe_tracks = event->Get<CmdVector<CmdLXeTrack> >("lxe_tracks");

    const cmd3::CmdVector<cmd3::CmdFullCluster> *clusters = 
      event->Get<const cmd3::CmdVector<cmd3::CmdFullCluster> >("full_clusters");
  
    const cmd3::CmdVector <cmd3::CmdZCOuterSector>* Sector  = event->Get<cmd3::CmdVector<cmd3::CmdZCOuterSector> >("zc_sectors");

    if ( ! clusters ) return;

    data.nfc_total =  clusters->GetSize();

    if( data.nfc_total> (int) MAX_CLUSTER ) return;

    data.nfc = 0;
    data.nph = 0;

    for (unsigned int t = 0; t < MAX_CLUSTER; ++t) {
      data.phflag[t] = 0;
      data.phfc[t] = 0;
    }

    for (unsigned int t = 0; t < clusters->GetSize(); ++t) {

      const cmd3::CmdFullCluster & tclust = (*clusters)[t];

      if (tclust.GetEnergy()>3*data.ebeam) continue;

      data.fcen[data.nfc]  = tclust.GetEnergy();
      data.fclxe[data.nfc]  = tclust.getLXeEnergy();
      data.fccsi[data.nfc]  = tclust.getCsIEnergy();
      data.fcbgo[data.nfc]  = tclust.getBGOEnergy();
      data.fclxe[data.nfc]  = tclust.getLXeEnergy();
      data.fcphi[data.nfc] = tclust.GetPosition().c3;
      data.fcth[data.nfc]  = tclust.GetPosition().c2;
      data.ecaltot += tclust.GetEnergy();

      //
      // form photons data
      //

      if(data.fcflag[data.nfc]==0 && data.nph < (int) MAX_CLUSTER){

	//calculate total energy of free clusters
	data.ecalneu += tclust.GetEnergy();
      
	double thc0 = tclust.getGammaTheta();
	double rhoc = tclust.GetPosition().c1*sin(thc0);// in CM - spherical system - radius from 0,0
      
	// correct cluster theta to track Z (IS IT CORRECT?)
	double clz = rhoc/tan(thc0);
	double thc = atan2(rhoc,clz-data.z0);
      
	//corrected photon parameters      
	data.phphi[data.nph]  = tclust.getGammaPhi();
	data.phen[data.nph] = tclust.getGammaEnergy();
	data.phrho[data.nph] = rhoc;
	if(data.nt>0)data.phth[data.nph]  = thc; else data.phth[data.nph] = thc0;


	//not corrected for photon parameters
	data.phen0[data.nph] = tclust.GetEnergy();
	data.phrho0[data.nph]  = tclust.GetPosition().c1;
	data.phth0[data.nph]   = tclust.GetPosition().c2;
	data.phphi0[data.nph]  = tclust.GetPosition().c3;

	data.phlxe[data.nph] = tclust.getLXeEnergy();
	data.phcsi[data.nph] = tclust.getCsIEnergy(); // zero csi and bgo with the same index
	data.phbgo[data.nph] = tclust.getBGOEnergy();
	data.phflag[data.nph] = tclust.getCType();// 1-  LXe based, 2- CsI based, 3-BGO based cluster
	data.phconv[data.nph] = tclust.haveConvPoint();// 1 if LXe strips present, 0 - if no

	// size of the cluster
      
	data.phfc[data.nph] += tclust.getNLXe();
	data.phfc[data.nph] += tclust.getNCsI();
	data.phfc[data.nph] += tclust.getNBGO();

	//covariant matrix
	data.pherr[data.nph][0] = tclust.getEnergyError();
	data.pherr[data.nph][1] = tclust.getThetaError();
	data.pherr[data.nph][2] = tclust.getPhiError();
      
      
	//data.phslxe[data.nph] = 0;
	for (unsigned int ilayer=0; ilayer<14; ilayer++) data.phslxe_layers[data.nph][ilayer] = 0;
	
	int lxe_index = tclust.getLXeClusterIndex();
	if (lxe_strips_clusters && lxe_strips_hits && lxe_index > -1) {
	  
	  const cmd3::CmdLXeClusterExt & tclust_lxe = (*lxe_towers_clusters)[lxe_index];	  
	  const std::set<int> & ind_sclust = tclust_lxe.GetStripClusters();	  
	  
	  for (std::set<int>::const_iterator k = ind_sclust.begin(); k != ind_sclust.end(); k++){	    
	    const CmdLXeStripClusterExt & sclust = CmdLXeRecoTools::GetStripElem(*k,*lxe_strips_clusters);
	    const CmdLXeUtils::PACKED_ID pid = CmdLXeUtils::UnPackValues(sclust.GetIndex());    
	    int layer_num=2*(pid.layer)+pid.direction;	    
	    double senergy = sclust.GetTotalAmplitude(*lxe_strips_hits);
	    data.phslxe_layers[data.nph][layer_num] += senergy;	    
	  }
	}
	
	data.fcflag[data.nfc]=2;// change flag for photon
	// Add data from ZC sectors
	data.phnzcs[data.nph] = 0;
	if( Sector ){
	  for (unsigned int j=0; j < Sector->GetSize(); ++j ) {
	    cmd3::CmdZCOuterSector SectHit = (cmd3::CmdZCOuterSector) Sector->Get(j);
	    double dphi = data.phphi[data.nph]-SectHit.Get_Phi()+0.2;
	    if( dphi >  (TMath::Pi()) ) dphi-=(TMath::Pi()*2);
	    if( dphi < -(TMath::Pi()) ) dphi+=(TMath::Pi()*2);	
	    if( fabs(dphi)<0.5 ) {
	      data.phnzcs[data.nph] += 1;//number of outer Z sectors connected to photon
	    }
	  }
	}

      }
    
      if(data.fcflag[data.nfc]==2)data.nph++;
    
      data.nfc ++;
    }
  }


  //
  // Fill ZC sectors data
  //
  void fill_zcs_data(tree_data &data, const cmd3::CmdEvent* event) {

    const cmd3::CmdVector <cmd3::CmdZCSectorHit>* Sector  = event->Get<cmd3::CmdVector<cmd3::CmdZCSectorHit> >("zcsector_clbr_hits");
    if( ! Sector ) return;

    data.nzcs_total = Sector->GetSize();

    if( Sector->GetSize()> (int) MAX_ZCS ) return;

    data.nzcs = Sector->GetSize();

    for (unsigned int j=0; j < Sector->GetSize(); ++j )
      {
	cmd3::CmdZCSectorHit SectHit = (cmd3::CmdZCSectorHit) Sector->Get(j);

	data.zcsch[j] = SectHit.GetChannel();
	data.zcsstat[j] = SectHit.GetStatus();
	data.zcsphi[j] = SectHit.GetPhi();
	data.zcsamp[j] = SectHit.GetAmplitude();
	data.zcstime[j] = SectHit.GetTime();
      }
    const cmd3::CmdVector<cmd3::CmdZCOuterCluster>* ZCCCluster = event->Get<cmd3::CmdVector<cmd3::CmdZCOuterCluster> > ("zc_clusters");
    if( ! ZCCCluster ) return;
    data.nzcc_total=ZCCCluster->GetSize();
    if( ZCCCluster->GetSize()>MAX_ZCC ) return;
    data.nzcc=ZCCCluster->GetSize();

    for ( int j1 = 0; j1 < ZCCCluster->GetSize(); j1++ )
      {
	cmd3::CmdZCOuterCluster ZCluster1 = ZCCCluster->Get(j1);

	data.zccamp[j1] = ZCluster1.Get_Amplitude();
	data.zcct[j1] = ZCluster1.GetType();
	data.zccz[j1] = ZCluster1.GetZCoordinate();
	data.zccl[j1] = ZCluster1.Get_NoLayer();
	data.zccns[j1]= ZCluster1.Get_NStrips();
	data.zccvalid[j1] = ZCluster1.GetIsValid();
      }


  }
  //
  // Fill ks candidates data
  //

  void fill_ks_data(tree_data &data, const cmd3::CmdEvent* event) {

    //
    //call all needed collections
    //
    const cmd3::CmdList<cmd3::CmdDCTrack*> *DCtracks =
      event->Get<cmd3::CmdList<cmd3::CmdDCTrack*> >("dc_tracks");

    //vector pointers to tracks
    vector<cmd3::CmdDCTrack*> trks;

    //
    // consider only events with 2 tracks
    //
    if( (DCtracks->GetSize())<2 || (DCtracks->GetSize())> (int) MAX_TRACK ) return;

    int iks=0;
    int iks_total = 0;
    //  int itr0 = 0;
    //  int itr1 = 0;
    int ch0,ch1;
    //
    // cycle ove pair of tracks
    //



    for(cmd3::CmdList<cmd3::CmdDCTrack*>::const_iterator j0=DCtracks->begin();j0!=DCtracks->end(); j0++) {
      cmd3::CmdDCTrack *trk0 = (*j0) ;

      if(trk0->GetNHits()<5)continue;
      std::map<cmd3::CmdDCTrack*,int>::iterator it0=track2id.find(trk0);
      if(it0==track2id.end()||(*it0).second<0) continue;

      // Get track parameters
      double xpca0,ypca0,zpca0, k0, ctg0, rho0 ,t00, phi0;

      /*
	double xfp, yfp, zfp;
	trk0->GetFPT(xfp, yfp, zfp);
	bool inv = (xfp*xfp+yfp*yfp)>(15*15); // invert track if the first point is far away from IP
	trk0->TrkParams(xyzbeam[0],xyzbeam[1],inv,k, phi, rho, ctg,xpca,ypca,zpca, t0);
      */

      trk0->TrkParams(xyzbeam[0],xyzbeam[1],k0, phi0, rho0, ctg0,xpca0,ypca0,zpca0, t00);

      if(k0>0)ch0=1;else ch0=-1;

      cmd3::CmdList<cmd3::CmdDCTrack*>::const_iterator j1=j0;
      j1++;
      for(;j1!=DCtracks->end(); j1++) {
	cmd3::CmdDCTrack *trk1 = (*j1) ;

	if(trk1->GetNHits()<5)continue;
	std::map<cmd3::CmdDCTrack*,int>::iterator it1=track2id.find(trk1);
	if(it1==track2id.end()||(*it1).second<0) continue;

	// Get track parameters
	double xpca1,ypca1,zpca1, k1, ctg1, rho1 ,t01, phi1;
	//    double xfp, yfp, zfp;
	//      trk1->GetFPT(xfp, yfp, zfp);
	//      bool inv = (xfp*xfp+yfp*yfp)>(15*15); // invert track if the first point is far away from IP
	//      trk1->TrkParams(xyzbeam[0],xyzbeam[1],inv,k, phi, rho, ctg,xpca,ypca,zpca, t0);

	trk1->TrkParams(xyzbeam[0],xyzbeam[1],k1, phi1, rho1, ctg1,xpca1,ypca1,zpca1, t01);

	if(k1>0)ch1=1;else ch1=-1;
	//
	// check correct charge
	//
	if( ch0+ch1 != 0)continue;

	// mark events with small and large relative angle

	double dphi = phi0 - phi1;
	if(dphi > M_PI)dphi = dphi - 2.*M_PI;
	if(dphi < -M_PI)dphi = dphi + 2.*M_PI;

	double xfit0 = xyzbeam[0];
	double yfit0 = xyzbeam[1];

	double SigmaVee = SigmaKsVert;

	//	double rtube = 2.5;
	double rtube = 1.75;

	int type = 0;
	double m = 139.57018;
	if(fabs(dphi)<1.0) {
	  m = 0.511;
	  type = 1;

	  rtube = 0.85;  //new1        
	  xfit0 = rtube*cos(phi0-dphi/2.);
	  yfit0 = rtube*sin(phi0-dphi/2.);

	  //SigmaVee = 0.2;
	  SigmaVee = 2.0; //new1

	  /*
	    if(data.nv>0){
	    xfit0 = data.vxyz[0][0];
	    yfit0 = data.vxyz[0][1];
	    }
	  */

	}
	//skip ks finding if below threshold

	if(data.ebeam < 490 && type == 0)continue;

	//clear pointer vector and fill with pair of candidate tracks track
	trks.clear();

	//trk->GetInfo();
	//cout<<trk->GetNHits()<<endl;

	// fill list of tracks for fitvert
	trks.push_back(trk0);
	trks.push_back(trk1);

	//----------------------------------------------------------------
	// fit with beam point as first approximation with sigma = SigmaKsVert, get new vertex

	TVector3 vertfit;
	double chi2ks;
	vector<TLorentzVector> Pmfit;
	int fsta = -1;
        double t0fit;
	//		cout<<"fit trks="<<trks.size()<<endl;
	//
	//chi2ks=fitvertex(&trks,xyzbeam[0],xyzbeam[1],SigmaKsVert,SigmaKsVert,
	//            chi2ks=fitvertex(&trks,xfit0,yfit0,SigmaKsVert,SigmaKsVert,
	chi2ks=fitvertex(&trks,fsta,xfit0,yfit0,SigmaVee,SigmaVee,
			 vertfit,Pmfit,t0fit,
			 Bfield,m);
	//		cout<<"fittetd "<<chi2fit<<endl;

	// cout<<"x="<<xyzbeam[0]<<" y="<<xyzbeam[1]<<" p="<<p_trk<<" pfit="<<Pmfit[0].P()<<" Sigma="<<SigmaBeam<<endl;

	// check if chi2 is good and mass is in KS window

	if(chi2ks>kschi2cut)continue;

	TLorentzVector sum2pi = Pmfit[0] + Pmfit[1];

	//calculate open angle

	double dpsi = acos((Pmfit[0].Px()*Pmfit[1].Px()+Pmfit[0].Py()*Pmfit[1].Py()+Pmfit[0].Pz()*Pmfit[1].Pz())/(Pmfit[0].P()*Pmfit[1].P()));

	//
	// check if pair is in ks window or angle is small
	//
	if( !((fabs(sum2pi.M()-mKs)<kswindow && type == 0 && sum2pi.P() > 20.&& sum2pi.P()<1.5*data.ebeam )
	      || (dpsi<eeangle && type == 1 && sum2pi.M()<100 && data.ecaltot>data.ebeam&&sum2pi.P()<1.5*data.ebeam))
	    )continue;

	iks_total++;

	if(iks>=MAX_KS) continue;

	//fill info for ks or ee candidate

	data.kstype[iks] = type;
	data.ksfstatus[iks] = fsta;
	data.ksminv[iks] = sum2pi.M();
	data.ksphi[iks] =  sum2pi.Phi();
	if(sum2pi.Phi()<0) data.ksphi[iks] =  sum2pi.Phi()+2*M_PI;
	data.ksth[iks] =  sum2pi.Theta();
	data.ksptot[iks] =  sum2pi.P();

	//vertex parameters

	data.ksvchi[iks] = chi2ks;
	data.ksvxyz[iks][0] = vertfit[0];
	data.ksvxyz[iks][1] = vertfit[1];
	data.ksvxyz[iks][2] = vertfit[2];

	data.ksdpsi[iks] = dpsi;

	TVector3 TLen(vertfit[0]-xyzbeam[0],vertfit[1]-xyzbeam[1],vertfit[2]-xyzbeam[2]);

	data.kstlen[iks] = TLen.Perp();
	data.ksalign[iks] = (TLen.X()*sum2pi.Px()+TLen.Y()*sum2pi.Py())/TLen.Perp()/sum2pi.Perp();

	// point of origin for ks

	double z0ks =  vertfit[2] - TLen.Perp()/tan(sum2pi.Theta());// sign "-" verified by A.Semenov

	TVector3 Len(vertfit[0]-xyzbeam[0],vertfit[1]-xyzbeam[1],vertfit[2]-z0ks);
	data.kslen[iks] = Len.Mag();

	data.ksz0[iks] = z0ks;

	//fill info for pions

	data.kspiphi[iks][0]  = Pmfit[0].Phi();
	if(Pmfit[0].Phi()<0) data.kspiphi[iks][0]  = Pmfit[0].Phi()+2*M_PI;
	data.kspith[iks][0]   = Pmfit[0].Theta();
	data.kspipt[iks][0]   = Pmfit[0].P();
	data.kspiphi[iks][1]  = Pmfit[1].Phi();
	if(Pmfit[1].Phi()<0) data.kspiphi[iks][1]  = Pmfit[1].Phi()+2*M_PI;
	data.kspith[iks][1]   = Pmfit[1].Theta();
	data.kspipt[iks][1]   = Pmfit[1].P();



	data.ksvind[iks][0]=track2id[trk0];
	data.ksvind[iks][1]=track2id[trk1];

	// data.kserr[iks][0][0] = ????

	//--------------------------------------------------------------
	iks++;

	h_ksminv->Fill(sum2pi.M());

      }
    }
    data.nks_total = iks_total;
    data.nks = iks;
  }


  //
  // Fill Verteces sectors data
  //
  void fill_vtx_data(tree_data &data, const cmd3::CmdEvent* event) {

    const cmd3::CmdList<cmd3::CmdDCVertex*>* DCVertices = event->Get<cmd3::CmdList<cmd3::CmdDCVertex*> >("dc_verticies");
    if(! DCVertices) return;

    data.nv_total = DCVertices->GetSize();
    if( DCVertices -> GetSize()> (int) MAX_VERTEX ) return;
    DCVertices->Begin();

    //const cmd3::CmdList<cmd3::CmdDCTrack*> *DCtracks = event->Get<cmd3::CmdList<cmd3::CmdDCTrack*> >("dc_tracks");
    //const cmd3::CmdList<cmd3::CmdDCTrack*>* DCTracks;
    //cmd3::CmdDCTrack* DCTrack;
    for(int j=0;j< (int) MAX_VERTEX;j++)
      for(int i=0;i< (int) MAX_TRACK;i++) data.vind[j][i]=-1;

    for (unsigned int j=0; j < DCVertices->GetSize(); ++j )
      {
	cmd3::CmdDCVertex* DCVertex = DCVertices->Next();


	TVector3 VertexCoordinates = DCVertex->GetXYZVertex();
	data.vxyz[data.nv][0] = VertexCoordinates[0];
	data.vxyz[data.nv][1] = VertexCoordinates[1];
	data.vxyz[data.nv][2] = VertexCoordinates[2];
	data.vtrk[data.nv] = 0;
	data.vchi[data.nv] = DCVertex->GetChi2();

	const cmd3::CmdList<cmd3::CmdDCTrack*> *dctracks2 =DCVertex->GetTracks();
	dctracks2->Begin();
	for(unsigned int itr=0;itr<dctracks2->GetSize();itr++){
	  cmd3::CmdDCTrack *dctrack = dctracks2->Next();
	  std::map<cmd3::CmdDCTrack*,int>::iterator it=track2id.find(dctrack);
	  if(it!=track2id.end()&&(*it).second>=0){
	    int idx=(*it).second;
	    data.vind[data.nv][data.vtrk[data.nv]]=idx;
	    data.vtrk[data.nv]++;
	  }else{
	    //cout<<"can not find track from vertex"<<endl;
	  }
	  if(data.vtrk[data.nv]>=(int)MAX_TRACK) break;
	}


	if(data.vtrk[data.nv]>1)data.nv++;
      }

  }


  //
  // Fill Ant data
  //
  void fill_ant_data(tree_data &data, const cmd3::CmdEvent* event) {

    const cmd3::CmdVector<cmd3::CmdAnTClbrHit>* AnTHit =
      event->Get<cmd3::CmdVector<cmd3::CmdAnTClbrHit> >("ant_clbr_hits");

    if (! AnTHit) return;

    data.nant = AnTHit->GetSize();
    if(data.nant  > MAX_ANT ) data.nant=MAX_ANT;
    data.anttime = 0;
    double tavant = 0;
    int ngood = 0;

    for (unsigned int j=0; j < data.nant; ++j )
      {
	cmd3::CmdAnTClbrHit antHit = (cmd3::CmdAnTClbrHit) AnTHit->Get(j);

	data.antch[j] = antHit.GetChannelNumber();
	data.antt0[j] = antHit.GetTime0();
	data.antt1[j] = antHit.GetTime1();
	data.anta0[j] = antHit.GetAmplitude0();
	data.anta1[j] = antHit.GetAmplitude1();
	data.antst[j] = antHit.GetStatus();

	if(fabs(data.antt0[j]-data.antt1[j])*0.14<20){
	  tavant += (data.antt0[j]+data.antt1[j])/2;
	  ngood++;
	}
                    
      }

    if(ngood>0)data.anttime = tavant/ngood*0.14;

  }

  //
  // Fill Mu data
  //
  void fill_mu_data(tree_data &data, const cmd3::CmdEvent* event) {

    const cmd3::CmdVector<cmd3::CmdMuClbrHit>* MuHit =
      event->Get<cmd3::CmdVector<cmd3::CmdMuClbrHit> >("mu_clbr_hits");

    if (! MuHit) return;

    data.nmu = MuHit->GetSize();
    if(data.nmu  > MAX_MU ) data.nmu=MAX_MU;
    data.mutime = 0;
    double tavmu = 0;
    int ngood = 0;

    for (unsigned int j=0; j < data.nmu; ++j )
      {
	cmd3::CmdMuClbrHit muHit = (cmd3::CmdMuClbrHit) MuHit->Get(j);
	data.much[j] = muHit.GetChannelNumber();
	data.mut0[j] = muHit.GetTime0();
	data.mut1[j] = muHit.GetTime1();
        data.mut2[j] = muHit.GetTime2();
        data.mut3[j] = muHit.GetTime3();
	data.mua0[j] = muHit.GetAmplitude0();
	data.mua1[j] = muHit.GetAmplitude1();
	data.mua2[j] = muHit.GetAmplitude2();
	data.mua3[j] = muHit.GetAmplitude3();
	data.must[j] = muHit.GetStatus();
		    
	if(fabs(data.mut0[j]-data.mut1[j])*0.14<40){
	  tavmu += (data.mut0[j]+data.mut1[j])/2;
	  ngood++;
	}
      }

    if(ngood>0)data.mutime = tavmu/ngood*0.14;

  }
  //
  // Fill MCTRUTH information for the simulation data
  //
  void fill_mctruth_data(tree_data &data, const cmd3::CmdEvent* event)
  {

    const cmd3::CmdScalar<HepMC::GenEvent*>* CmdScalGenEvent = NULL;

    if (event->Has ("event_gen")) {
      CmdScalGenEvent = event->Get<cmd3::CmdScalar<HepMC::GenEvent* > > ("event_gen");
    }
    HepMC::GenEvent* _GenEvent = CmdScalGenEvent->Get();
    data.finalstate_id=_GenEvent->signal_process_id();
    if(event->Has ("event_mctruth")) {
      CmdScalGenEvent = event->Get<cmd3::CmdScalar<HepMC::GenEvent* > > ("event_mctruth");
    }

    if(!CmdScalGenEvent)return;

    //  CmdScalGenEvent = dynamic_cast<cmd3::CmdScalar<HepMC::GenEvent*>*>event->Get("event_mctruth");
    _GenEvent = CmdScalGenEvent->Get();
    //  _GenEvent->print();
    data.nsim = 0;
    for ( HepMC::GenEvent::particle_iterator p = _GenEvent->particles_begin(); 
	  p != _GenEvent->particles_end(); ++p ) { // Cycle over all particles in a event
      if ( (*p)->production_vertex()->particles_in_size() == 0 ) { // beam particle is found
	//      cout << "Inside cycle over beam particles " << endl;
	HepMC::GenParticle *pt=(*p);
	data.simtype[data.nsim] = pt->pdg_id();
	data.simorig[data.nsim] = 0;
	TLorentzVector Mom = convert2ROOT(pt->momentum());
	data.simmom[data.nsim] = Mom.Rho();
	data.simphi[data.nsim] = (Mom.Phi() > 0)? Mom.Phi() : Mom.Phi() + 2.0*TMath::Pi();
	data.simtheta[data.nsim] = Mom.Theta();
	HepMC::FourVector vertex = pt->production_vertex()->position();
	data.simvtx[data.nsim] = vertex.x()/10.0;
	data.simvty[data.nsim] = vertex.y()/10.0;
	data.simvtz[data.nsim] = vertex.z()/10.0;
	//      cout << data.simvtx[data.nsim] << " " << data.simvty[data.nsim] <<" " << data.simvtz[data.nsim] << endl;
	data.nsim++;
	if(data.nsim>=(int)MAX_SIM) break;
	//      cout << "Before searching for decayed vertex " << endl;
	HepMC::GenVertex *v1 = _GenEvent->barcode_to_vertex(-((*p)->barcode()));
	//      cout << "After searching for decayed vertex " << endl;
	//      if ( ( pt->pdg_id() == 111 || pt->pdg_id() == 310 || pt->pdg_id() == 130 || pt->pdg_id() == 221 || pt->pdg_id() == 331 ) && v1 != NULL )  { // We save decay particles from pi0/KL/KS/Eta/Eta'
	if ( ( pt->pdg_id() == 111 || pt->pdg_id() == 310 || pt->pdg_id() == 221 || pt->pdg_id() == 331 ) && v1 != NULL )  { // We save decay particles from pi0/KL/KS/Eta/Eta'
	  for(HepMC::GenVertex::particle_iterator p1=v1->particles_begin(HepMC::children);p1 != v1->particles_end(HepMC::children); ++p1 ) {
	    //	  cout << "Inside cycle over decayed pi0 or KS " << endl;
	    HepMC::GenParticle *pt1=(*p1);
	    data.simtype[data.nsim] = pt1->pdg_id();          
	    data.simorig[data.nsim] = pt->pdg_id();
	    TLorentzVector Mom1 = convert2ROOT(pt1->momentum());
	    data.simmom[data.nsim] = Mom1.Rho();
	    data.simphi[data.nsim] = (Mom1.Phi() > 0)? Mom1.Phi() : Mom1.Phi() + 2.0*TMath::Pi();
	    data.simtheta[data.nsim] = Mom1.Theta();
	    HepMC::FourVector vertex1 = v1->position();
	    data.simvtx[data.nsim] = vertex1.x()/10.0;
	    data.simvty[data.nsim] = vertex1.y()/10.0;
	    data.simvtz[data.nsim] = vertex1.z()/10.0;
	    data.nsim++;
	    if(data.nsim>=(int)MAX_SIM) break;
	    HepMC::GenVertex *v2 = _GenEvent->barcode_to_vertex(-(pt1->barcode()));
	    //	  if ( ( pt1->pdg_id() == 111 || pt1->pdg_id() == 310 || pt1->pdg_id() == 130 || pt1->pdg_id() == 221 || pt1->pdg_id() == 331 ) && v2 != NULL )  { // We save decay particles from pi0/KL/KS/Eta/Eta'
	    if ( ( pt1->pdg_id() == 111 || pt1->pdg_id() == 310 || pt1->pdg_id() == 221 || pt1->pdg_id() == 331 ) && v2 != NULL )  { // We save decay particles from pi0/KL/KS/Eta/Eta'
	      for (  HepMC::GenVertex::particle_iterator p2 = v2->particles_begin(HepMC::children); 
		     p2 != v2->particles_end(HepMC::children); ++p2 ) {
		//	      cout << "Inside cycle over decayed pi0" << endl;
		HepMC::GenParticle *pt2=(*p2);
		data.simtype[data.nsim] = pt2->pdg_id();
		data.simorig[data.nsim] =  pt1->pdg_id();
		TLorentzVector Mom2 = convert2ROOT(pt2->momentum());
		data.simmom[data.nsim] = Mom2.Rho();
		data.simphi[data.nsim] = (Mom2.Phi() > 0)? Mom2.Phi() : Mom2.Phi() + 2.0*TMath::Pi();
		data.simtheta[data.nsim] = Mom2.Theta();
		HepMC::FourVector vertex2 = v2->position();
		data.simvtx[data.nsim] = vertex2.x()/10.0;
		data.simvty[data.nsim] = vertex2.y()/10.0;
		data.simvtz[data.nsim] = vertex2.z()/10.0;
		data.nsim++;
		if(data.nsim>=(int)MAX_SIM) break;
		HepMC::GenVertex *v3 = _GenEvent->barcode_to_vertex(-(pt2->barcode()));
		//	      if ( ( pt2->pdg_id() == 111 || pt2->pdg_id() == 310 || pt2->pdg_id() == 130 || pt2->pdg_id() == 221 || pt2->pdg_id() == 331 ) && v3 != NULL )  { // We save decay particles from pi0/KL/KS/Eta/Eta'
		if ( ( pt2->pdg_id() == 111 || pt2->pdg_id() == 310 || pt2->pdg_id() == 221 || pt2->pdg_id() == 331 ) && v3 != NULL )  { // We save decay particles from pi0/KL/KS/Eta/Eta'
		  for ( HepMC::GenVertex::particle_iterator p3 = v3->particles_begin(HepMC::children);
			p3 != v3->particles_end(HepMC::children); ++p3 ) {
		    HepMC::GenParticle *pt3=(*p3);
		    data.simtype[data.nsim] = pt3->pdg_id();
		    data.simorig[data.nsim] =  pt2->pdg_id();
		    TLorentzVector Mom3 = convert2ROOT(pt3->momentum());
		    data.simmom[data.nsim] = Mom3.Rho();
		    data.simphi[data.nsim] = (Mom3.Phi() > 0)? Mom3.Phi() : Mom3.Phi() + 2.0*TMath::Pi();
		    data.simtheta[data.nsim] = Mom3.Theta();
		    HepMC::FourVector vertex3 = v3->position();
		    data.simvtx[data.nsim] = vertex3.x()/10.0;
		    data.simvty[data.nsim] = vertex3.y()/10.0;
		    data.simvtz[data.nsim] = vertex3.z()/10.0;
		    data.nsim++;
		    if(data.nsim>=(int)MAX_SIM) break;
		  }
		}
		if(data.nsim>=(int)MAX_SIM) break;
	      }
	    }
	    if(data.nsim>=(int)MAX_SIM) break;
	  }
	}
	if(data.nsim>=(int)MAX_SIM) break;
      }
    }
  }
  
  void fill_correction_data(tree_data &data, const cmd3::CmdEvent* event) {
    if (event->Has("bits_correction") && event->Has("banks_correction")) {
      const CmdVector<int>* bits_cor = event->Get<CmdVector<int> >("bits_correction");
      const CmdVector<int>* banks_cor = event->Get<CmdVector<int> >("banks_correction");
      data.ncorr = bits_cor->size();
      for ( int i = 0; i < data.ncorr; i++ ) {
	data.idcorr[i] = banks_cor->at(i);
	data.bitcorr[i] = bits_cor->at(i);
      }
    }
  }

  void fill_corrupted_data(tree_data &data, const cmd3::CmdEvent* event) {
    if (event->Has("corrupted_banks") && event->Has("corrupted_status")) {
      const CmdVector<int>* corrupted_b = event->Get<CmdVector<int> >("corrupted_banks");
      const CmdVector<int>* corrupted_bs = event->Get<CmdVector<int> >("corrupted_status");
      data.nbadbank = corrupted_b->size();
      data.nbadbankg = 0;
      for (int i = 0; i < corrupted_b->size(); i++) {
	int linkid = corrupted_b->at(i);
	int link_number = -1;
	if ( linkid == 81 || linkid == 84 || linkid == 95 || linkid == 166 )
	  link_number = 0; //Trig
	else if ( linkid >= 1 && linkid <= 80 )
	  link_number = 1; //DC
	else if ( (linkid >= 201 && linkid <= 209) || (linkid >= 401 && linkid <= 479) )
	  link_number = 2; //LXe
	else if ( linkid >= 121 && linkid <= 159 )
	  link_number = 3; //CSi
	else if ( linkid >= 301 && linkid <= 329 )
	  link_number = 4; //BGO
	else if ( linkid >= 251 && linkid <= 294 )
	  link_number = 5; //TOF
	else if ( linkid >= 221 && linkid <= 232 )
	  link_number = 6; //Mu
	else if ( (linkid >= 111 && linkid <= 114) || (linkid >= 361 && linkid <= 376) )
	  link_number = 7; //ZC
	if (link_number != -1) {
	  data.nbadbanks[data.nbadbankg] = link_number;
	  data.nbadbankg++;
	}
	if (corrupted_bs->at(i) == 1) data.nlostbanks++;
	else if (corrupted_bs->at(i) == 2) data.ncorruptedbanks++;
      }
    }
  }
  
  //
  // Main function: create tree and histograms, select events and fill tree
  //

  TFile *top;
  TTree *tout;
  tree_data data;
  bool if_sim_mode=false;
  bool fillclusters = false;
    
  double cutdTh  = 0.25;  // maximum acollinearity polar angle
  double cutdPhi = 0.15;  // maximum acollinearity asimuth angle
  double cutMinNhit = 5;    // number of points on the track
  double  cutMinEph = 15.;   // cut on minimum photon energy
  //double cutEtot = 200.;         // cut on total energy in calorimeter
  double cutEtot = Ebeam;         // cut on total energy in calorimeter

  void tr_ph_init(const cmd3::CmdRunHeader* run_header){
	
    std::cout<<"init"<<std::endl;
	
    top = new TFile(Form("tr_ph_run%05ld.root",run_header->GetRunID()), "recreate");
    
    h_dphi_lxe = new TH2D("h_dphi_lxe","DPhivsDTheta between LXE cluster and track with strips",500,-3.2,3.2,500,-1.8,1.8);
    h_dphi_lxe0 = new TH2D("h_dphi_lxe0","DPhivsDTheta between LXE cluster and track if no strips",500,-3.2,3.2,500,-1.8,1.8);
    h_dphi_csi = new TH2D("h_dphi_csi","DPhivsDTheta between CSI cluster and track",500,-3.2,3.2,500,-1.8,1.8);
    h_dphi_bgo = new TH2D("h_dphi_bgo","DPhivsDTheta between BGO cluster and track",500,-3.2,3.2,500,-2.4,2.3);

    h_dphi_lxe1 = new TH2D("h_dphi_lxe1","DPhivsDTheta between LXE cluster and track with strips",500,-3.2,3.2,500,-1.8,1.8);
    h_dphi_lxe01 = new TH2D("h_dphi_lxe01","DPhivsDTheta between LXE cluster and track if no strips",500,-3.2,3.2,500,-1.8,1.8);
    h_dphi_csi1 = new TH2D("h_dphi_csi1","DPhivsDTheta between CSI cluster and track",500,-3.2,3.2,500,-1.8,1.8);
    h_dphi_bgo1 = new TH2D("h_dphi_bgo1","DPhivsDTheta between BGO cluster and track",500,-3.2,3.2,500,-2.4,2.3);

    h_dphi_tlxe = new TH2D("h_dphi_tlxe","DPhi vs DTheta between LXE track and DC track; dphi; dtheta",500, -3.2, 3.2, 500,-2, 2);
    h_dr_tlxe = new TH2D("h_dr_tlxe","Dr vs Ro between LXE track and DC track at LXe; ro, cm; dr, cm",100, 36, 52, 350, 0, 50);

    h_lxe_csi_ph = new TH2D("h_lxe_csi_ph","DPhivsDTheta between LXe photon and CsI photon",500,-3.2,3.2,500,-2.4,2.3);
    h_dphi_cal = new TH1D("h_dphi_cal","DPhi between two LXE clusters",100,-1,1);
    h_dphi_ant = new TH1D("h_dphi_ant","DPhi between track and ANT counter",100,-6.2,6.2);
    h_ksminv = new TH1D("h_ksminv","KS candidate invariant mass",200,400.,600.);
    

    tout= new TTree("tr_ph","Tree with the non-collinear events");
    tout->Branch("ebeam",     &data.ebeam,     "ebeam/F");
    tout->Branch("emeas",     &data.emeas,     "emeas/F");
    tout->Branch("demeas",    &data.demeas,    "demeas/F");
    tout->Branch("emeas0",     &data.emeas0,     "emeas0/F");
    tout->Branch("demeas0",    &data.demeas0,    "demeas0/F");
    tout->Branch("xbeam",     &data.xbeam,     "xbeam/F");
    tout->Branch("ybeam",     &data.ybeam,     "ybeam/F");
    tout->Branch("runnum",    &data.runnum,    "runnum/I");
    tout->Branch("finalstate_id", &data.finalstate_id, "finalstate_id/I");
    tout->Branch("evnum",     &data.evnum,     "evnum/I");
    tout->Branch("trigbits",  &data.trigbits,  "trigbits/I");
    tout->Branch("trigmchs",  &data.trigmchs,  "trigmchs/I");
    tout->Branch("trigtime",  &data.trigtime,  "trigtime/F");
    tout->Branch("time",      &data.time,      "time/F");
    tout->Branch("dcfittime", &data.dcfittime, "dcfittime/F");
    tout->Branch("anttime",   &data.anttime,   "anttime/F");
    tout->Branch("mutime",    &data.mutime,    "mutime/F");
    tout->Branch("is_coll",   &data.is_coll,   "is_coll/I");
    tout->Branch("is_bhabha", &data.is_bhabha, "is_bhabha/I");
    tout->Branch("nt_total",  &data.nt_total,  "nt_total/I");
    tout->Branch("ecaltot",   &data.ecaltot,   "ecaltot/F");
    tout->Branch("ecalneu",   &data.ecalneu,   "ecalneu/F");
    tout->Branch("z0",        &data.z0,        "z0/F");
    tout->Branch("psumch",    &data.psumch,    "psumch/F");
    tout->Branch("psumnu",    &data.psumnu,    "psumnu/F");
    tout->Branch("lumoff",    &data.lumoff,    "lumoff/F");
    tout->Branch("lumofferr",    &data.lumofferr,    "lumofferr/F");

    tout->Branch("nv_total", &data.nv_total,   "nv_total/I");
    tout->Branch("nv",       &data.nv,         "nv/I");
    tout->Branch("vtrk",      data.vtrk,       "vtrk[nv]/I");
    tout->Branch("vind",      data.vind,       "vind[nv][10]/I");
    tout->Branch("vchi",      data.vchi,       "vchi[nv]/F");
    tout->Branch("vxyz",      data.vxyz,       "vxyz[nv][3]/F");

    tout->Branch("nt",        &data.nt,        "nt/I");
    tout->Branch("it",         data.it,        "it[2]/I");
    tout->Branch("tnhit",      data.tnhit,     "tnhit[nt]/I");
    tout->Branch("tlength",    data.tlength,      "tlength[nt]/F");
    tout->Branch("tphi",       data.tphi,      "tphi[nt]/F");
    tout->Branch("tth",        data.tth,       "tth[nt]/F");
    tout->Branch("tptot",      data.tptot,     "tptot[nt]/F");
    tout->Branch("tphiv",      data.tphiv,     "tphiv[nt]/F");
    tout->Branch("tthv",       data.tthv,      "tthv[nt]/F");
    tout->Branch("tptotv",     data.tptotv,    "tptotv[nt]/F");
    tout->Branch("trho",       data.trho,      "trho[nt]/F");
    tout->Branch("tdedx",      data.tdedx,     "tdedx[nt]/F");
    tout->Branch("tz",         data.tz,        "tz[nt]/F");
    tout->Branch("tt0",        data.tt0,       "tt0[nt]/F");
    tout->Branch("tant",       data.tant,      "tant[nt]/F");
    tout->Branch("tchi2r",     data.tchi2r,    "tchi2r[nt]/F");
    tout->Branch("tchi2z",     data.tchi2z,    "tchi2z[nt]/F");
    tout->Branch("tchi2ndf",   data.tchi2ndf,  "tchi2ndf[nt]/F");
    tout->Branch("tcharge",    data.tcharge,   "tcharge[nt]/I");
    tout->Branch("ten",        data.ten,       "ten[nt]/F");
    tout->Branch("tfc",        data.tfc,       "tfc[nt]/F");
    tout->Branch("tenlxe",     data.tenlxe,    "tenlxe[nt]/F");
    tout->Branch("tlengthlxe",data.tlengthlxe,"tlengthlxe[nt]/F");
    //tout->Branch("tenslxe",     data.tenslxe,    "tenslxe[nt]/F");
    tout->Branch("tenslxe_layers",data.tenslxe_layers,  "tenslxe_layers[nt][14]/F");
    
    tout->Branch("tencsi",     data.tencsi,    "tencsi[nt]/F");
    tout->Branch("tenbgo",     data.tenbgo,    "tenbgo[nt]/F");
    tout->Branch("tclth",      data.tclth,     "tclth[nt]/F");
    tout->Branch("tclphi",     data.tclphi,    "tclphi[nt]/F");
    tout->Branch("terr",       data.terr,      "terr[nt][3][3]/F");
    tout->Branch("terr0",      data.terr0,     "terr0[nt][5][5]/F");

    tout->Branch("tindlxe",    data.tindlxe,   "tindlxe[nt]/I");
    tout->Branch("tzcc",       data.tzcc,      "tzcc[nt][2]/F");
    tout->Branch("txyzatcl",   data.txyzatcl,  "txyzatcl[nt][3]/F");
    tout->Branch("txyzatlxe",  data.txyzatlxe, "txyzatlxe[nt][3]/F");
    tout->Branch("tenconv",    data.tenconv,   "tenconv[nt]/I");

    tout->Branch("nks_total", &data.nks_total, "nks_total/I");
    tout->Branch("nks",       &data.nks,       "nks/I");
    tout->Branch("ksvind",     data.ksvind,    "ksvind[nks][20]/I");
    tout->Branch("kstype",     data.kstype,    "kstype[nks]/I");
    tout->Branch("ksfstatus",     data.ksfstatus,    "ksfstatus[nks]/I");
    tout->Branch("ksvchi",     data.ksvchi,    "ksvchi[nks]/F");
    tout->Branch("ksvxyz",     data.ksvxyz,    "ksvxyz[nks][3]/F");
    tout->Branch("ksminv",     data.ksminv,    "ksminv[nks]/F");
    tout->Branch("ksalign",    data.ksalign,   "ksalign[nks]/F");
    tout->Branch("kstlen",     data.kstlen,    "kstlen[nks]/F");
    tout->Branch("ksdpsi",     data.ksdpsi,    "ksdpsi[nks]/F");
    tout->Branch("kslen",      data.kslen,     "kslen[nks]/F");
    tout->Branch("ksz0",       data.ksz0,      "ksz0[nks]/F");
    tout->Branch("ksphi",      data.ksphi,     "ksphi[nks]/F");
    tout->Branch("ksth",       data.ksth,      "ksth[nks]/F");
    tout->Branch("ksptot",     data.ksptot,    "ksptot[nks]/F");
    tout->Branch("kspiphi",    data.kspiphi,   "kspiphi[nks][2]/F");
    tout->Branch("kspith",     data.kspith,    "kspith[nks][2]/F");
    tout->Branch("kspipt",     data.kspipt,    "kspipt[nks][2]/F");
    //tout->Branch("kserr",      data.kserr,     "kserr[nks][3][3]/F");

    // LXe treacks
    tout->Branch("ntlxe_total", &data.ntlxe_total, "ntlxe_total/I");
    tout->Branch("ntlxe",       &data.ntlxe,       "ntlxe/I");
    tout->Branch("ntlxelayers",       data.ntlxelayers,       "ntlxelayers[ntlxe]/I");
    tout->Branch("tlxenhit",    data.tlxenhit,     "tlxenhit[ntlxe]/I");
    tout->Branch("tlxelength",    data.tlxelength,      "tlxelength[ntlxe]/F");
    tout->Branch("tlxededx",    data.tlxededx,      "tlxededx[ntlxe]/F");

    tout->Branch("tlxeir",    data.tlxeir,      "tlxeir[ntlxe]/F");
    tout->Branch("tlxeitheta",    data.tlxeitheta,      "tlxeitheta[ntlxe]/F");
    tout->Branch("tlxeiphi",    data.tlxeiphi,      "tlxeiphi[ntlxe]/F");
    tout->Branch("tlxevtheta",    data.tlxevtheta,      "tlxevtheta[ntlxe]/F");
    tout->Branch("tlxevphi",    data.tlxevphi,      "tlxevphi[ntlxe]/F");
    tout->Branch("tlxechi2",    data.tlxechi2,      "tlxechi2[ntlxe]/F");
    tout->Branch("tlxesen",    data.tlxesen,      "tlxesen[ntlxe]/F");
    tout->Branch("tlxesen_layers",data.tlxesen_layers,  "tlxesen_layers[ntlxe][14]/F");
    //tout->Branch("itlxe",    data.itlxe,     "itlxe[ntlxe]/I");
    
    if(fillclusters){// make true is need detaoled data from calorimeters
      
      tout->Branch("nfc_total", &data.nfc_total, "nfc_total/I");
      tout->Branch("nfc",       &data.nfc,       "nfc/I");
      tout->Branch("fcen",       data.fcen,      "fcen[nfc]/F");
      tout->Branch("fclxe",       data.fclxe,      "fclxe[nfc]/F");
      tout->Branch("fccsi",       data.fccsi,      "fccsi[nfc]/F");
      tout->Branch("fcbgo",       data.fcbgo,      "fcbgo[nfc]/F");
      tout->Branch("fcth",       data.fcth,      "fcth[nfc]/F");
      tout->Branch("fcphi",      data.fcphi,     "fcphi[nfc]/F");
      //    tout->Branch("fcentot",    data.fcentot,   "fcentot[nfc]/F");
      tout->Branch("fcflag",     data.fcflag,    "fcflag[nfc]/I");
      
      tout->Branch("nlxe_total", &data.nlxe_total, "nlxe_total/I");
      tout->Branch("nlxe",       &data.nlxe,       "nlxe/I");
      tout->Branch("lxeen",       data.lxeen,      "lxeen[nlxe]/F");
      tout->Branch("lxeth",       data.lxeth,      "lxeth[nlxe]/F");
      tout->Branch("lxephi",      data.lxephi,     "lxephi[nlxe]/F");
      tout->Branch("lxeentot",    data.lxeentot,   "lxeentot[nlxe]/F");
      tout->Branch("lxenzcs",     data.lxenzcs,    "lxenzcs[nlxe]/I");
      tout->Branch("lxeflag",     data.lxeflag,    "lxeflag[nlxe]/I");

      tout->Branch("ncsi_total", &data.ncsi_total, "ncsi_total/I");
      tout->Branch("ncsi",       &data.ncsi,       "ncsi/I");
      tout->Branch("csien",       data.csien,      "csien[ncsi]/F");
      tout->Branch("csith",       data.csith,      "csith[ncsi]/F");
      tout->Branch("csiphi",      data.csiphi,     "csiphi[ncsi]/F");
      tout->Branch("csiflag",     data.csiflag,    "csiflag[ncsi]/I");

      tout->Branch("nbgo_total", &data.nbgo_total, "nbgo_total/I");
      tout->Branch("nbgo",       &data.nbgo,       "nbgo/I");
      tout->Branch("bgoen",       data.bgoen,      "bgoen[nbgo]/F");
      tout->Branch("bgoth",       data.bgoth,      "bgoth[nbgo]/F");
      tout->Branch("bgophi",      data.bgophi,     "ngophi[nbgo]/F");
      tout->Branch("bgoflag",     data.bgoflag,    "bgoflag[nbgo]/I");

    }

    tout->Branch("nph_total",  &data.nph_total,"nph_total/I");
    tout->Branch("nph",        &data.nph,      "nph/I");
    tout->Branch("phen",       data.phen,      "phen[nph]/F");
    tout->Branch("phth",       data.phth,      "phth[nph]/F");
    tout->Branch("phphi",      data.phphi,     "phphi[nph]/F");
    tout->Branch("phrho",      data.phrho,     "phrho[nph]/F");
    tout->Branch("phen0",       data.phen0,      "phen0[nph]/F");
    tout->Branch("phth0",       data.phth0,      "phth0[nph]/F");
    tout->Branch("phphi0",      data.phphi0,     "phphi0[nph]/F");
    tout->Branch("phlxe",      data.phlxe,     "phlxe[nph]/F");
    //tout->Branch("phslxe",      data.phslxe,     "phslxe[nph]/F");
    tout->Branch("phslxe_layers",  data.phslxe_layers,     "phslxe_layers[nph][14]/F");
    
    tout->Branch("pherr",      data.pherr,     "pherr[nph][3]/F");
    tout->Branch("phcsi",      data.phcsi,     "phcsi[nph]/F");
    tout->Branch("phbgo",      data.phbgo,     "phbgo[nph]/F");
    tout->Branch("phflag",     data.phflag,    "phflag[nph]/I");
    tout->Branch("phconv",     data.phconv,    "phconv[nph]/I");
    tout->Branch("phfc",       data.phfc,      "phfc[nph]/I");


    tout->Branch("nzcs_total", &data.nzcs_total, "nzcs_total/I");
    tout->Branch("nzcs",       &data.nzcs,       "nzcs/I");
    tout->Branch("zcsch",       data.zcsch,      "zcsch[nzcs]/I");
    tout->Branch("zcsstat",     data.zcsstat,    "zcsstat[nzcs]/I");
    tout->Branch("zcsamp",      data.zcsamp,     "zcsamp[nzcs]/F");
    tout->Branch("zcstime",     data.zcstime,    "zcstime[nzcs]/F");
    tout->Branch("zcsphi",      data.zcsphi,     "zcsphi[nzcs]/F");

    tout->Branch("nzcc_total", &data.nzcc_total, "nzcc_total/I");
    tout->Branch("nzcc",       &data.nzcc,       "nzcc/I");
    tout->Branch("zccl",       data.zccl,        "zccl[nzcc]/I");
    tout->Branch("zccns",      data.zccns,       "zccns[nzcc]/I");
    tout->Branch("zccamp",     data.zccamp,      "zccamp[nzcc]/F");
    tout->Branch("zcct",       data.zcct,        "zcct[nzcc]/I");
    tout->Branch("zccz",       data.zccz,        "zccz[nzcc]/F");
    tout->Branch("zccvalid",   data.zccvalid,    "zccvalid[nzcc]/I");

    tout->Branch("nant",       &data.nant,       "nant/I");
    tout->Branch("antch",      data.antch,     "antch[nant]/I");
    tout->Branch("antt0",      data.antt0,     "antt0[nant]/F");
    tout->Branch("antt1",      data.antt1,     "antt1[nant]/F");
    tout->Branch("anta0",      data.anta0,     "anta0[nant]/F");
    tout->Branch("anta1",      data.anta1,     "anta1[nant]/F");
    tout->Branch("antst",      data.antst,     "antst[nant]/I");

    tout->Branch("nmu",       &data.nmu,       "nmu/I");
    tout->Branch("much",      data.much,     "much[nmu]/I");
    tout->Branch("mut0",      data.mut0,     "mut0[nmu]/F");
    tout->Branch("mut1",      data.mut1,     "mut1[nmu]/F");
    tout->Branch("mut2",      data.mut2,     "mut2[nmu]/F");
    tout->Branch("mut3",      data.mut3,     "mut3[nmu]/F");
    tout->Branch("mua0",      data.mua0,     "mua0[nmu]/F");
    tout->Branch("mua1",      data.mua1,     "mua1[nmu]/F");
    tout->Branch("mua2",      data.mua2,     "mua2[nmu]/F");
    tout->Branch("mua3",      data.mua3,     "mua3[nmu]/F");
    tout->Branch("must",      data.must,     "must[nmu]/I");

    tout->Branch("nsim",&data.nsim,"nsim/I");
    tout->Branch("simtype",data.simtype,"simtype[nsim]/I");
    tout->Branch("simorig",data.simorig,"simorig[nsim]/I");
    tout->Branch("simmom",data.simmom,"simmom[nsim]/F");
    tout->Branch("simphi",data.simphi,"simphi[nsim]/F");
    tout->Branch("simtheta",data.simtheta,"simtheta[nsim]/F");
    tout->Branch("simvtx",data.simvtx,"simvtx[nsim]/F");
    tout->Branch("simvty",data.simvty,"simvty[nsim]/F");
    tout->Branch("simvtz",data.simvtz,"simvtz[nsim]/F");
    
    tout->Branch("ncorr", &data.ncorr, "ncorr/I");
    tout->Branch("idcorr", data.idcorr, "idcorr[ncorr]/I");
    tout->Branch("bitcorr", data.bitcorr, "bitcorr[ncorr]/I");

    tout->Branch("nbadbank", &data.nbadbank, "nbadbank/I");
    tout->Branch("nbadbankg", &data.nbadbankg, "nbadbankg/I");
    tout->Branch("nbadbanks", data.nbadbanks, "nbadbanks[nbadbankg]/I");
    tout->Branch("nlostbanks", &data.nlostbanks, "nlostbanks/I");
    tout->Branch("ncorruptedbanks", &data.ncorruptedbanks, "ncorruptedbanks/I");

    data.runnum = run_header->GetRunID();
    cout<<"runnum "<<data.runnum<<endl;

    const cmd3::CmdParam& param =run_header->GetParam();
    xyzbeam=param.GetDoubleArray("BeamPosition");
    if(xyzbeam.size()<3) xyzbeam.resize(3,0);
    Ebeam = param.GetDouble("BeamEnergy");

    std::vector<double> energyy;
    energyy.push_back(Ebeam);
    if (! param.GetInt("SimMode")) {
      Emeas0 = param.IsSet("Compton1", CMDPARAM_DOUBLEARRAY_TYPE) ? param.GetDoubleArray("Compton1") : energyy;
      Emeas = param.IsSet("Compton", CMDPARAM_DOUBLEARRAY_TYPE) ? param.GetDoubleArray("Compton") : energyy;
      cout << " param.IsSet(Compton1, CMDPARAM_DOUBLEARRAY_TYPE) = " << param.IsSet("Compton1", CMDPARAM_DOUBLEARRAY_TYPE) << endl;
    } else {
      Emeas.push_back(Ebeam);
      Emeas0.push_back(Ebeam);
    }
    if(Emeas.size()<2) Emeas.resize(2,0);
    if(Emeas0.size()<2) Emeas0.resize(2,0);   
    cout<<"Ebeam "<<Ebeam<< " Emeasured "<< Emeas0[0]  <<" dEmeas0 "<< Emeas0[1] << " Emeasured averaged"<< Emeas[0] <<" dEmeas "<< Emeas[1] << endl;

    cutEtot = Ebeam;         // cut on total energy in calorimeter
    //cutEtot = 0; // no cut on total energy for latest experimental  tr_ph_fc*.root 
    
    data.ebeam = Ebeam;
    data.emeas = Emeas[0];
    data.demeas = Emeas[1];
    data.emeas0 = Emeas0[0];
    data.demeas0 = Emeas0[1];
    data.xbeam = xyzbeam[0];
    data.ybeam = xyzbeam[1];
    
    Bfield = param.GetDouble("MagneticField")/10000;
    cout<<"BField "<<Bfield<<endl;
    
    // --- LUMoff ---
    char comand[4096];
    double  a1 = 0.0, a2 = 0.0;
    const char *ch=0;
    
    // open connection to MySQL server on localhost
    /*
    TSQLServer *db = TSQLServer::Connect("mysql://sl10cmd.inp.nsk.su/online", "daqread","daqread");
    TSQLResult *res;
    TSQLRow *row; 
    cout<<"online db connect code "<< db-> IsConnected()<<" error code "<<db->IsError()<<endl;
    
    sprintf(comand,"select value, error from RunparOff where run=%d and param=0 and name=\"lumoff\"",data.runnum);
    cout<<"Comand: "<<comand<<endl;
    res = db->Query(comand);
    while ((row = res->Next())) {
      cout<<"lumdb ok "<<endl;
      a1=atof(row->GetField(0));
      a2=atof(row->GetField(1));
      delete row;}   
    delete res;
    delete db;
    */

    data.lumoff = a1; 
    data.lumofferr = a2;
    cout<<"LUMoff "<<data.lumoff<<" +/- "<<data.lumofferr<<endl;
    // --- LUMoff end ---

    /// Selection of collinear tracks-----------------------------------------------

    // Parameters for collinear selection

    if(param.GetInt("SimMode")==1)
      if_sim_mode=true;
    else
      if_sim_mode=false;

  }

  void tr_ph_event(const cmd3::CmdEvent* event){
    clear_event(data);
    int EvNum = event->GetID();//event number
    data.evnum = EvNum;
    //
    // Get trigger bits
    //
    const cmd3::CmdScalar<cmd3::CmdTrgDecision> *decisions_collection = 
      event->Get< cmd3::CmdScalar<cmd3::CmdTrgDecision> >("trg_decision");

    // check if collection exists
    try{
      const cmd3::CmdScalar<double> *trg_guessed =
	event->Get< cmd3::CmdScalar<double> >("trg_guessed");

      if( trg_guessed ) data.dcfittime = trg_guessed->Get();

      const cmd3::CmdScalar<double> *trg_time =
	event->Get< cmd3::CmdScalar<double> >("trg_time");

      const cmd3::CmdScalar<double> *iptz_trg_time =
	event->Get< cmd3::CmdScalar<double> >("iptz_trg_time");

      if( trg_time ) data.trigtime = trg_time->Get();
      if(iptz_trg_time && data.trigtime == 0)data.trigtime = iptz_trg_time->Get();
	      
      const cmd3::CmdScalar<long> *time_fromadis = event->Get< cmd3::CmdScalar<long> >("mchs_vepp_time");
      data.time = time_fromadis->Get();

    }catch(...){};


    if( decisions_collection ){
      cmd3::CmdTrgDecision decision = decisions_collection->Get();
      data.trigbits = decision.GetBitMask();
      data.trigmchs = decision.GetFirstTrg();
    }

    fill_lxe_data(data,event);

    fill_tracks_data(data,event);
    fill_zcs_data(data,event);
    fill_ant_data(data,event);
    fill_mu_data(data,event);
    fill_csi_data(data,event);
    fill_bgo_data(data,event);
    fill_fullcluster_data(data,event);
    fill_correction_data(data, event);
    fill_corrupted_data(data, event);
    if(if_sim_mode) fill_mctruth_data(data,event);

    // select only close tracks to be stored
    int nt_new =0;
    double z0_new = 0;

    for(int t =0; t<data.nt; ++t ){

      //	if(data.tchi2r[t]>100 || data.tchi2z[t]>100 || fabs(data.trho[t])>6. )continue;
      //	if(data.tchi2r[t]>ChiTrackMax || data.tchi2z[t]>ChiTrackMax || fabs(data.trho[t])>RhoTrackMax ){
      if(data.tchi2ndf[t]>ChiTrackMax || fabs(data.trho[t])>RhoTrackMax ){
	//!! to save cluster, the fcflag[] is zeroed on the line 830, with the same selections!!!
	track2id[id2track[t]]=-1;
	tracktrue[t]=false;
	continue;
      }

      data.tlength[nt_new]  = data.tlength[t];
      data.tphiv[nt_new]  = data.tphiv[t];
      data.tthv[nt_new]   = data.tthv[t];
      data.tptotv[nt_new] = data.tptotv[t];
      data.tphi[nt_new] = data.tphi[t];
      data.tth[nt_new] = data.tth[t];
      data.trho[nt_new] = data.trho[t];
      data.tptot[nt_new] =data.tptot[t];
      data.tdedx[nt_new] = data.tdedx[t];
      data.tz[nt_new] = data.tz[t] ;
      z0_new += data.tz[nt_new];
      data.tt0[nt_new] =data.tt0[t] ;
      data.tchi2r[nt_new] = data.tchi2r[t];
      data.tchi2z[nt_new] = data.tchi2z[t] ;
      data.tchi2ndf[nt_new] = data.tchi2ndf[t] ;
      data.tenlxe[nt_new] = data.tenlxe[t];
      data.tlengthlxe[nt_new] = data.tlengthlxe[t];
      //data.tenslxe[nt_new] = data.tenslxe[t];
      for (int i=0; i<14; i++) data.tenslxe_layers[nt_new][i] = data.tenslxe_layers[t][i];
      data.tencsi[nt_new] = data.tencsi[t];
      data.tenbgo[nt_new] = data.tenbgo[t];
      data.ten[nt_new] = data.ten[t];
      data.tclth[nt_new] = data.tclth[t];
      data.tclphi[nt_new] = data.tclphi[t];
      data.tfc[nt_new] = data.tfc[t];
      data.tant[nt_new] =data.tant[t];
      data.tnhit[nt_new] = data.tnhit[t];
      data.tcharge[nt_new] = data.tcharge[t];
      data.tindlxe[nt_new] = data.tindlxe[t];
      data.tzcc[nt_new][0] = data.tzcc[t][0];
      data.tzcc[nt_new][1] = data.tzcc[t][1];
      data.txyzatcl[nt_new][0] = data.txyzatcl[t][0];
      data.txyzatcl[nt_new][1] = data.txyzatcl[t][1];
      data.txyzatcl[nt_new][2] = data.txyzatcl[t][2];
      data.txyzatlxe[nt_new][0] = data.txyzatlxe[t][0];
      data.txyzatlxe[nt_new][1] = data.txyzatlxe[t][1];
      data.txyzatlxe[nt_new][2] = data.txyzatlxe[t][2];
      data.tenconv[nt_new] = data.tenconv[t];

      //error matrix
      for(int i=0;i<3;++i){
	for(int j=0;j<3;++j){
	  data.terr[nt_new][i][j] = data.terr[t][i][j];
	}
      }
      //full error matrix from fit

      for(int i0=0;i0<5;++i0){
	for(int j0=0;j0<5;++j0){
	  data.terr0[nt_new][i0][j0] = data.terr0[t][i0][j0];
	}
      }

      if(nt_new!=t){
	id2track[nt_new]=id2track[t];
	track2id[id2track[nt_new]]=nt_new;
	id2track[t]=0;
      }

      nt_new ++;
    }

    data.nt = nt_new;

    if(nt_new>0){
      data.z0 = z0_new/nt_new;

      fill_vtx_data(data,event);

      fill_ks_data(data,event);
    }
    
    TLorentzVector psum;
    
    for(int i = 0; i< data.nt; ++i){
      TLorentzVector vec(data.tptot[i]*TMath::Cos(data.tphi[i])*TMath::Sin(data.tth[i]),
			 data.tptot[i]*TMath::Sin(data.tphi[i])*TMath::Sin(data.tth[i]),
			 data.tptot[i]*TMath::Cos(data.tth[i]),
			 TMath::Hypot(data.tptot[i],139.67));

      psum += vec;
    }
    data.psumch = psum.P();
    
    
    //single track filter
    bool is_one_track = false;
    //      if(data.nt==1 && fabs(data.z0)<12.&&data.ecalneu>0) is_one_track = true;
    //      if(data.nt==1 && fabs(data.z0)<12.&&fabs(data.trho[0])<0.5&&data.ten[0]>cutMinEph) is_one_track = true;
    if(data.nt==1 && fabs(data.z0)<12.&&fabs(data.trho[0])<0.5) is_one_track = true;
    //
    // set thresholt to the photon energies and exclude them from the list
    // always increase photon threshold if no tracks or only one track
    //
    bool is_no_track = false;
    if(data.nt==0 && data.ecaltot>cutEtot) is_no_track = true;
    
    //generator induced trigger for background study
    
    bool is_bkg_trigger = false;
    if(data.trigmchs==16)is_bkg_trigger = true;
    
    
    //remove flags == 0 and set photon threshold
    data.nph_total = data.nph;
    int nph_new = 0;

    for (int t = 0; t < data.nph; ++t) {
      if(data.phflag[t] > 0 ){
	if(data.phen0[t] > cutMinEph){
	  data.phen[nph_new] = data.phen[t];
	  data.phen0[nph_new] = data.phen0[t];
	  data.phphi[nph_new] = data.phphi[t];
	  data.pherr[nph_new][0] = data.pherr[t][0];
	  data.pherr[nph_new][1] = data.pherr[t][1];
	  data.pherr[nph_new][2] = data.pherr[t][2];
	  data.phphi0[nph_new] = data.phphi0[t];
	  data.phth[nph_new] = data.phth[t];
	  data.phth0[nph_new] = data.phth0[t];
	  data.phrho[nph_new] = data.phrho[t];
	  data.phrho0[nph_new] = data.phrho0[t];
	  data.phconv[nph_new] = data.phconv[t];
	  data.phlxe[nph_new] = data.phlxe[t];
	  //data.phslxe[nph_new] = data.phslxe[t];
	  for (unsigned int ilayer=0; ilayer<14; ilayer++) data.phslxe_layers[nph_new][ilayer] = data.phslxe_layers[t][ilayer];
	  data.phcsi[nph_new] = data.phcsi[t];
	  data.phbgo[nph_new] = data.phbgo[t];
	  data.phflag[nph_new] = data.phflag[t];
	  data.phfc[nph_new] = data.phfc[t];
	  nph_new++;
	}
      }
    }
    data.nph = nph_new;
    
    TLorentzVector psum0;
    
    for(int i = 0; i< data.nph; ++i){
      TLorentzVector vec(data.phen[i]*TMath::Cos(data.phphi[i])*TMath::Sin(data.phth[i]),
			 data.phen[i]*TMath::Sin(data.phphi[i])*TMath::Sin(data.phth[i]),
			 data.phen[i]*TMath::Cos(data.phth[i]),
			 TMath::Hypot(data.phen[i],0));

      psum0 += vec;
    }
    data.psumnu = psum0.P();

    // Select collinear events

    while( true ) {

      data.is_coll=-1;
      if( data.nt<2 ) break;

      data.is_coll=-2;
      if( data.nt>10 ) break;// was >4

      // Loop over all track pairs and try to find the collinear one
      for( int i=0; i<(data.nt-1); i++ ){
	for( int j=(i+1); j<data.nt; j++ ) {

	  // --- cut on number of points on track
	  if( data.is_coll > -3 ) data.is_coll = -3;
	  if( data.tnhit[i]<cutMinNhit || data.tnhit[j]<cutMinNhit ) continue;

	  // --- cut on total charge
	  if( data.is_coll > -4 ) data.is_coll = -4;
	  if( (data.tcharge[i]+data.tcharge[j])!=0 ) continue;

	  // --- cut on delta_phi
	  if( data.is_coll > -5 ) data.is_coll = -5;
	  double dphi = fabs(data.tphi[i]-data.tphi[j])-TMath::Pi();
	  if( fabs(dphi)>cutdPhi ) continue;

	  // --- cut on delta_theta
	  if( data.is_coll > -6 ) data.is_coll = -6;
	  double dth = data.tth[i]+data.tth[j]-TMath::Pi();
	  if( fabs(dth)>cutdTh ) continue;

	  data.is_coll=1;
	  data.it[0]=i;
	  data.it[1]=j;
	  if(data.tcharge[i]<0){
	    data.it[0]=j;
	    data.it[1]=i;
	  }
	  break;
	}
	if(data.is_coll==1)break;
      }

      break;
    }

    // Select bhabha events using LXe calorimeter data only

    while( true ) {

      if( data.nlxe!=2 ) break;

      double thavr = (data.lxeth[0]-data.lxeth[1]+TMath::Pi())/2.0;
      //if( thavr<1.0 || thavr>(TMath::Pi()-1.0) ) break;

      if( data.ecaltot < Ebeam/2 ) break;
      //	if( data.fcen[1]<(Ebeam/2) ) break;

      if( data.lxenzcs[0]<=0 ) break;
      if( data.lxenzcs[1]<=0 ) break;

      // Calculate rotation to the calorimeter
      double pt = Ebeam*sin(thavr)*1.0e-3;
      double sinphi = 0.40*0.3*Bfield/(2*pt);
      double phi0 = sinphi<1.0?TMath::ASin(sinphi):0.0;

      double dphi = fabs(data.lxephi[0]-data.lxephi[1])-TMath::Pi();
      if( dphi > (TMath::Pi()) ) dphi-=(TMath::Pi()*2);
      if( dphi < -(TMath::Pi()) ) dphi+=(TMath::Pi()*2);

      dphi = fabs(dphi)-2*phi0;

      h_dphi_cal->Fill(dphi);

      if( fabs(dphi)>0.2 ) break;

      double dth = data.lxeth[0]+data.lxeth[1]-TMath::Pi();
      if( fabs(dth)>0.3 ) break;

      data.is_bhabha=1;
      break;

    }


    //      if( data.is_coll>0 || data.is_bhabha>0 || data.nt>1 || is_one_track ) tout->Fill();
    //            if( data.nt>1 || is_one_track || (is_no_track&&data.nph>0)) tout->Fill();
    if( is_bkg_trigger || data.nt>1 || is_one_track || (is_no_track&&data.nph>0) || if_sim_mode) 
      {
	int nwrite=tout->Fill();
	if(nwrite<0) throw std::ios_base::failure("tr_ph_event: IO Error to fill tree");
      }
    //if((is_no_track&&data.nph>0)) tout->Fill();
  }

  void tr_ph_stop(){
    // save histogram hierarchy in the file
    top->Write();
    std::cout<<"Outputs data saved to file="<<top->GetName()<<std::endl;
    std::cout<<gSystem->pwd()<<std::endl;

    top->Close();

    //delete top;
    std::cout<<"Hist production succeffuly done!"<<std::endl;
  }
  void tr_ph_fc(const std::string& filename){
    gBenchmark->Start("tr_phdo");

    cmd3::CmdReader player;
    player.Open(filename.c_str(), "READ");
    player.GetRunHeader()->GetInfo();

    tr_ph_init(player.GetRunHeader());

    int size = player.GetSize();
    for(int i = 0; i != size; ++i)
      {
        if((i%1000)==0) cout<<"Event "<<i<<endl;
	const cmd3::CmdEvent* event = player.GetEvent(i);
	tr_ph_event(event);
      }
    tr_ph_stop();

    player.Close();

    gBenchmark->Stop("tr_phdo");
    gBenchmark->Print("tr_phdo");

  }
}
