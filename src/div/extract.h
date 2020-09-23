#ifndef EXTRACT_H
#define EXTRACT_H

namespace tr_ph_fc {

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
		Float_t psumch; // vector sum of charged track momenta 
		Float_t psumnu; // vector sum of photon momenta
		Float_t lumoff; // offline luminosity from db
		Float_t lumofferr; // stat error offline luminosity from db

		// Vertex data
		int nv_total;
		int nv;
		int vtrk[MAX_VERTEX]; // number of tracks in vertex
		int vind[MAX_VERTEX][MAX_TRACK]; // indices of tracks from vertex
		Float_t vchi[MAX_VERTEX];
		Float_t vxyz[MAX_VERTEX][3];

		// KS vertex data
		int nks_total;
		int nks;
		int ksvind[MAX_KS][MAX_TRACK]; // indices of tracks from vertex
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
		//Float_t tenslxe[MAX_TRACK]; // by strips == sun of tenslxe_layers
		Float_t tenslxe_layers[MAX_TRACK][14]; // by strips by layers

		Float_t tencsi[MAX_TRACK];
		Float_t tenbgo[MAX_TRACK];
		Float_t tclth[MAX_TRACK]; // theta and phi of connected cluster
		Float_t tclphi[MAX_TRACK];
		Float_t terr[MAX_TRACK][3][3];
		Float_t terr0[MAX_TRACK][5][5];

		Int_t tindlxe[MAX_TRACK]; // reference to LXe track number, == -1 if LXe track not found

		Float_t tfc[MAX_TRACK]; // size (number of crystals) of the clusters connected to track
		Float_t tzcc[MAX_TRACK][2]; // track point in zc
		Float_t txyzatcl[MAX_TRACK][3]; // track point at cluster conversion point
		Float_t txyzatlxe[MAX_TRACK][3]; // point of LXE cluster
		Int_t   tenconv[MAX_TRACK]; // 1 - if LXe has strip info, 0 - if no LXe strip info

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
		Int_t   phnzcs[MAX_CLUSTER]; // number of Z  chamber sectors
		// Conversion point in LXe calorimeter
		Int_t phconv[MAX_CLUSTER]; // 1 - if LXe has strip info, 0 - if no LXe strip info
		//Float_t phslxe[MAX_CLUSTER]; // by strips
		Float_t phslxe_layers[MAX_CLUSTER][14]; // by strips by layers
		Int_t phfc[MAX_CLUSTER]; // number of calorimeter elements in cluster (cluster size).

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
		int simtype[MAX_SIM]; // particle ID from GEANT
		int simorig[MAX_SIM]; // origin of particle. 0 - from e+e- or particle ID from GEANT  
		Float_t simmom[MAX_SIM];
		Float_t simphi[MAX_SIM];
		Float_t simtheta[MAX_SIM];
		Float_t simvtx[MAX_SIM]; // vertex X for origin
		Float_t simvty[MAX_SIM]; // vertex Y
		Float_t simvtz[MAX_SIM]; // vertex Z

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

	std::vector< double > xyzbeam;
	double Ebeam;
	std::vector< double > Emeas;
	std::vector< double > Emeas0;
	double Bfield;
	double SigmaBeam = 0.007;
	double mKs = 497.67;
	double SigmaKsVert = 100000.0;
	double kswindow = 80.;
	double eeangle = 1.;
	double kschi2cut = 50.;
	double ChiTrackMax = 100.;
	double RhoTrackMax = 6.0; // in cm, skip far tracks from tree
	bool tracktrue[ MAX_TRACK ];
	std::map< cmd3::CmdDCTrack*, int > track2id;
	std::map< int, cmd3::CmdDCTrack* > id2track;

}

#endif // EXTRACT_H
