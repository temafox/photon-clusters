/*
  gerasim example
  "gera" - must be unic name of your gerasim
  it is neccessary to define
  namespace gera_nm{
    void gera_init(const cmd3::CmdRunHeader* run_header){}
    void gera_event(const cmd3::CmdEvent* event){}
    void gera_stop(){}
  }

  void gerasim(const std::string& filename){} can be used for standalone running of macro
*/

#include "TMath.h"
#include "TLorentzVector.h"
#include "TBenchmark.h"
#include "TH1D.h"

#include "CmdReader.h"

#include "CmdLXeStripHit.h"
#include "CmdLXeStripHitExt.h"
#include "CmdLXeStripCluster.h"
#include "CmdLXeStripClusterExt.h"
#include "CmdSpiral.h"
#include "CmdPointData.h"

#include "HepMC/GenParticle.h"
#include "HepMC/GenVertex.h"
#include "HepMC/GenEvent.h"
#include "HepMC/IteratorRange.h"

namespace gera_nm{
  TFile *fout;

  TH1D *hist;

  const unsigned int MAX_SIM = 100;

  // common data about strips
  struct strip_data{
    std::vector<int> packedID;
    std::vector<int> innerID;
    std::vector<int> layer;
    std::vector<int> direction;
    std::vector<int> cluster_id;
    std::vector<double> amp;
    std::vector<double> sigma;
  };

  // input two packed index to get a position of two hits
  struct cross_data{
    std::vector<int> id1;
    std::vector<int> id2; // id1 is less than id2, that is collection order
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> z;
  };

  typedef struct {
    int nsim;
    int finalstate_id;
    int simtype[MAX_SIM];// particle ID from GEANT
    int simorig[MAX_SIM];// Origin of particle. 0 - from e+e- or particle ID from GEANT  
    Float_t simmom[MAX_SIM];
    Float_t simphi[MAX_SIM];
    Float_t simtheta[MAX_SIM];
    Float_t simvtx[MAX_SIM];// vertex X for origin
    Float_t simvty[MAX_SIM];//vertex Y
    Float_t simvtz[MAX_SIM];//vertex Z
   
   } tree_data;

  vector<double> xyzbeam;
  double ebeam;

  void fill_strips(strip_data &data, cross_data &cross, const cmd3::CmdEvent *event){

    const cmd3::CmdVector<cmd3::CmdVector<cmd3::CmdVector<cmd3::CmdLXeStripHit> > > * strip_hits
      =event->Get<const cmd3:: CmdVector<cmd3:: CmdVector<cmd3::CmdVector<cmd3::CmdLXeStripHit> > > >("lxe_strips_hits");
    const cmd3::CmdVector<cmd3::CmdVector<cmd3::CmdVector<cmd3::CmdLXeStripCluster> > > * strip_clusters
      =event->Get< cmd3::CmdVector<cmd3::CmdVector<cmd3::CmdVector<cmd3::CmdLXeStripCluster> > > >("lxe_strips_clusters");

    //clear all vectors

    data.packedID.clear();
    data.innerID.clear();
    data.layer.clear();
    data.direction.clear();
    data.cluster_id.clear();
    data.amp.clear();
    data.sigma.clear();

    cross.id1.clear();
    cross.id2.clear();
    cross.x.clear();
    cross.y.clear();
    cross.z.clear();

    std::vector<cmd3::CmdSpiral> spiral_massive;
    
    // fill strips data and also fill their spiral positions
    for(unsigned layer=0; layer < strip_hits->size(); layer++){
      for(unsigned drc=0; drc < (*strip_hits)[layer].size(); drc++){
	cmd3::CmdVector<cmd3::CmdLXeStripHit> strips = (*strip_hits)[layer][drc];
	std::sort(strips.begin(),strips.end(),cmd3::CmdLXeStripHitExt::cmp);
	for (unsigned int id = 0; id <strips.size(); id++){                                                                                                                                                          
          const cmd3::CmdLXeStripHitExt & hit = strips[id];                                                         
          double ampl = hit.GetAmplitude();
	  if(ampl>0){

	    data.packedID.push_back(hit.GetID());
	    data.innerID.push_back(id);
	    data.layer.push_back(layer);
	    data.direction.push_back(drc);
	    data.amp.push_back(hit.GetAmplitude());
	    data.sigma.push_back(hit.GetSigma());
	    data.cluster_id.push_back(-1);

	    spiral_massive.push_back(hit.GetPosition());

	    for(unsigned clayer=0; clayer < strip_clusters->size(); clayer++){
	      for(unsigned cdrc=0; cdrc < (*strip_clusters)[clayer].size(); cdrc++){
	        cmd3::CmdVector<cmd3::CmdLXeStripCluster> cstrips = (*strip_clusters)[clayer][cdrc];
		for(unsigned int cid = 0; cid < cstrips.size(); cid++){
		  const cmd3::CmdLXeStripClusterExt & cluster = cstrips[cid];
		  for(unsigned int i=0; i<cluster.GetHits().size(); i++){
		    if( hit.GetIndex() == cluster.GetHits().at(i)){
			data.cluster_id.back() = cluster.GetIndex();
		    }
		  }
		}
	      }
	    }


	  }
	}
      }
    }

    // find the crosses of the spirals
    if(spiral_massive.size() > 0){
      for(unsigned int i=0; i<spiral_massive.size()-1; i++){
        for(unsigned int j=i+1; j<spiral_massive.size(); j++){
          std::vector<cmd3::CmdPointData> coordinate = spiral_massive[i].CalculateCrossPoints(spiral_massive[j], -fabs(spiral_massive[i].GetStartPoint().c3), fabs(spiral_massive[i].GetStartPoint().c3),1);
	  if(coordinate.size() > 0){ // Does it have a size more than 1?
            const cmd3::CmdPointData cross_point = coordinate[0].GetXYZfromRoPhiZ(coordinate[0]);

	    cross.id1.push_back(data.packedID.at(i));
	    cross.id2.push_back(data.packedID.at(j));
  	    cross.x.push_back(cross_point.c1);
	    cross.y.push_back(cross_point.c2);
	    cross.z.push_back(cross_point.c3);
	  }
        }
      }
    }
  }

  TLorentzVector convert2ROOT(const HepMC::FourVector& in){ 
    TLorentzVector T(in.x(),in.y(),in.z(),in.t());
    return T;         
  }                   


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
            //    cout << "Inside cycle over decayed pi0 or KS " << endl;
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
            //    if ( ( pt1->pdg_id() == 111 || pt1->pdg_id() == 310 || pt1->pdg_id() == 130 || pt1->pdg_id() == 221 || pt1->pdg_id() == 331 ) && v2 != NULL )  { // We save decay particles from pi0/KL/KS/Eta/Eta'
            if ( ( pt1->pdg_id() == 111 || pt1->pdg_id() == 310 || pt1->pdg_id() == 221 || pt1->pdg_id() == 331 ) && v2 != NULL )  { // We save decay particles from pi0/KL/KS/Eta/Eta'
              for (  HepMC::GenVertex::particle_iterator p2 = v2->particles_begin(HepMC::children); 
                     p2 != v2->particles_end(HepMC::children); ++p2 ) {
                //            cout << "Inside cycle over decayed pi0" << endl;
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
                //            if ( ( pt2->pdg_id() == 111 || pt2->pdg_id() == 310 || pt2->pdg_id() == 130 || pt2->pdg_id() == 221 || pt2->pdg_id() == 331 ) && v3 != NULL )  { // We save decay particles from pi0/KL/KS/Eta/Eta'
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

  tree_data data;
  strip_data strips;
  cross_data cross;
  TTree *tout;

  void gera_init(const cmd3::CmdRunHeader* run_header){
    fout=new TFile(Form("strips_run%05ld.root",run_header->GetRunID()),"recreate");
    hist=new TH1D("hist","histogram;event number",1000,0,300000);
    tout = new TTree("tr_lxe","Tree with info for the lxe strips");
    tout->Branch("ebeam",&ebeam,"ebeam/F");
    tout->Branch("xyzbam",&xyzbeam);
    tout->Branch("nsim",&data.nsim,"nsim/I");
    tout->Branch("simtype",data.simtype,"simtype[nsim]/I");
    tout->Branch("simorig",data.simorig,"simorig[nsim]/I");
    tout->Branch("simmom",data.simmom,"simmom[nsim]/F");
    tout->Branch("simphi",data.simphi,"simphi[nsim]/F");
    tout->Branch("simtheta",data.simtheta,"simtheta[nsim]/F");
    tout->Branch("simvtx",data.simvtx,"simvtx[nsim]/F");
    tout->Branch("simvty",data.simvty,"simvty[nsim]/F");
    tout->Branch("simvtz",data.simvtz,"simvtz[nsim]/F");
    tout->Branch("strips",&strips);
    tout->Branch("cross_pos",&cross);

    const cmd3::CmdParam& param =run_header->GetParam();
    xyzbeam=param.GetDoubleArray("BeamPosition");
    if(xyzbeam.size()<3) xyzbeam.resize(3,0);
    ebeam = param.GetDouble("BeamEnergy");
  }

  void gera_event(const cmd3::CmdEvent* event){
    fill_mctruth_data(data, event);
    fill_strips(strips, cross, event);

    int nwrite = hist->Fill(event->GetID());
    if(nwrite<0) throw std::ios_base::failure("gerasim_event: IO Error to fill tree");
    nwrite = tout->Fill();                                                  
    if(nwrite<0) throw std::ios_base::failure("tr_lxe_event: IO Error to fill tree");

  }

  void gera_stop(){
    fout->Write();
  }

  void gerasim(const std::string& filename){
    gBenchmark->Start("gerado");

    cmd3::CmdReader player;
    player.Open(filename.c_str(), "READ");
    player.GetRunHeader()->GetInfo();

    gera_init(player.GetRunHeader());

    int size = player.GetSize();
    for(int i = 0; i !=size; ++i)
    {
        if((i%1000)==0) cout<<"Event "<<i<<endl;
	const cmd3::CmdEvent* event = player.GetEvent(i);
	gera_event(event);
    }
    gera_stop();

    player.Close();

    gBenchmark->Stop("gerado");
    gBenchmark->Print("gerado");

  }
}
