#include "ClusterShapeAnalysis/ClusterShapeHistProc.h"

#include "ClusterShapeAnalysis/ClusterHists.h"
#include "ClusterShapeAnalysis/TrackerHitResoHists.h"

#include <set>

#include <EVENT/LCCollection.h>
#include <EVENT/LCEvent.h>
#include <EVENT/LCRelation.h>
#include <EVENT/MCParticle.h>
#include <EVENT/Track.h>
#include <EVENT/TrackerHit.h>
#include <EVENT/SimTrackerHit.h>
#include <EVENT/TrackerHitPlane.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <UTIL/CellIDDecoder.h>
#include <UTIL/LCTrackerConf.h>

#include <marlin/AIDAProcessor.h>

#include <AIDA/ITree.h>

ClusterShapeHistProc  aClusterShapeHistProc ;

ClusterShapeHistProc::ClusterShapeHistProc()
  : Processor("ClusterShapeHistProc")
{  
  // modify processor description
  _description = "ClusterShapeHistProc creates a series of output histograms for tracking measurements (clusters) performance studies when realistic digitization is used." ;

  // register steering parameters: name, description, class-variable, default value
  // --- List of tracker cluster collections
  registerInputCollection( LCIO::TRACKERHIT,
			   "VBTrackerHitsCollection" , 
			   "Name of vertex barrel tracker hits collection",
			   _vbtrkhitColName,
			   ""
			   );     

  registerInputCollection( LCIO::TRACKERHIT,
			   "IBTrackerHitsCollection" , 
			   "Name of inner barrel tracker hits collection",
			   _ibtrkhitColName,
			   ""
			   );  

  registerInputCollection( LCIO::TRACKERHIT,
			   "OBTrackerHitsCollection" , 
			   "Name of outer barrel tracker hits collection",
			   _obtrkhitColName,
			   ""
			   );  

  registerInputCollection( LCIO::TRACKERHIT,
			   "VETrackerHitsCollection" , 
			   "Name of vertex endcap tracker hits collection",
			   _vetrkhitColName,
			   ""
			   );  

  registerInputCollection( LCIO::TRACKERHIT,
			   "IETrackerHitsCollection" , 
			   "Name of inner endcap tracker hits collection",
			   _ietrkhitColName,
			   ""
			   );     

  registerInputCollection( LCIO::TRACKERHIT,
			   "OETrackerHitsCollection" , 
			   "Name of outer endcap tracker hits collection",
			   _oetrkhitColName,
			   ""
			   );
  // --- Collection of MC particles
  registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticleCollection" , 
			   "Name of the MCParticle collection"  ,
			   _mcpColName,
			   ""
			   );  

  // --- Relationship reco-truth for tracker clusters
  registerInputCollection( LCIO::LCRELATION,
		     "VBRelationCollection" ,
			   "Name of the input vertex barrel relation collection",
			   _VBRelationCollection,
			   ""
		 	    );

    registerInputCollection( LCIO::LCRELATION,
		     "IBRelationCollection" ,
			   "Name of the input inner tracker barrel relation collection",
			   _IBRelationCollection,
			     ""
		 	    );

    registerInputCollection( LCIO::LCRELATION,
		     "OBRelationCollection" ,
			   "Name of the input outer tracker barrel relation collection",
			   _OBRelationCollection,
			     ""
		 	    );

    registerInputCollection( LCIO::LCRELATION,
		     "VERelationCollection" ,
			   "Name of the input vertex endcap relation collection",
			   _VERelationCollection,
			     ""
		 	    );

    registerInputCollection( LCIO::LCRELATION,
		     "IERelationCollection" ,
			   "Name of the input inner tracker endcap relation collection",
			   _IERelationCollection,
			     ""
		 	    );

    registerInputCollection( LCIO::LCRELATION,
		     "OERelationCollection" ,
			   "Name of the input outer tracker endcap relation collection",
			   _OERelationCollection,
			     ""
		 	    );
    
    registerProcessorParameter("muColDet",
                               "Muon collider detector geometry",
                               _muDet,
                               (int)0);

}


void ClusterShapeHistProc::init()
{
  // Print the initial parameters
  printParameters() ;

  // Create histograms
  AIDA::ITree* tree=marlin::AIDAProcessor::tree(this);
  marlin::AIDAProcessor::histogramFactory(this);

  // Create ROOT histograms, with location setup by the above factory
  tree->mkdir("clusters_vb" ); tree->cd("clusters_vb" );
  _resolution_vb=std::make_shared<TrackerHitResoHists>();
  _clusters_vb=std::make_shared<ClusterHists>();

  tree->cd("../");
  tree->mkdir("clusters_ve" ); tree->cd("clusters_ve" );
  _resolution_ve=std::make_shared<TrackerHitResoHists>();
  _clusters_ve=std::make_shared<ClusterHists>();

  tree->cd("../");
  tree->mkdir("clusters_ib" ); tree->cd("clusters_ib" );
  _resolution_ib=std::make_shared<TrackerHitResoHists>();
  _clusters_ib=std::make_shared<ClusterHists>();

  tree->cd("../");
  tree->mkdir("clusters_ie" ); tree->cd("clusters_ie" );
  _resolution_ie=std::make_shared<TrackerHitResoHists>();
  _clusters_ie=std::make_shared<ClusterHists>();

  tree->cd("../");
  tree->mkdir("clusters_ob" ); tree->cd("clusters_ob" );
  _resolution_ob=std::make_shared<TrackerHitResoHists>();
  _clusters_ob=std::make_shared<ClusterHists>();

  tree->cd("../");
  tree->mkdir("clusters_oe" ); tree->cd("clusters_oe" );
  _resolution_oe=std::make_shared<TrackerHitResoHists>();
  _clusters_oe=std::make_shared<ClusterHists>();

  tree->cd("../");
  h_trackerhit_timing = new TH1F("hit_timing", "Time of arrival of hits [ns]", 110, -10, 100);
  h_clusters_by_bLayer   = new TH1F("numClusters_by_bLayer"      , ";Barrel Layer Index; Number of Clusters",30,0,30);
  h_hits_by_bLayer   = new TH1F("numhits_by_bLayer"      , ";Barrel Layer Index; Number of Hits",30,0,30);
  h_clusters_by_bLayer_BX   = new TH1F("numClusters_by_bLayer_BX"      , ";Barrel Layer Index; Number of Clusters / 1 BX",30,0,30);
  h_hits_by_bLayer_BX   = new TH1F("numhits_by_bLayer_BX"      , ";Barrel Layer Index; Number of Hits / 1 BX",30,0,30);
  h_clusters_by_eLayer   = new TH1F("numClusters_by_eLayer"      , ";Endcap Layer Index; Number of Clusters",60,-30,30);
  h_hits_by_eLayer   = new TH1F("numhits_by_eLayer"      , ";Endcap Layer Index; Number of Hits",60,-30,30);
  h_clusters_by_eLayer_BX   = new TH1F("numClusters_by_eLayer_BX"      , ";Endcap Layer Index; Number of Clusters / 1 BX",60,-30,30);
  h_hits_by_eLayer_BX   = new TH1F("numhits_by_eLayer_BX"      , ";Endcap Layer Index; Number of Hits / 1 BX",60,-30,30);
  h_clusterDensity_bLayer = new TH1F("ClusterDensity_bLayer", ";Barrel Layer Index; Number of Clusters / 1 BX / cm^2",30,0,30);
  h_clusterDensity_eLayer = new TH1F("ClusterDensity_eLayer", ";Endcap Layer Index; Number of Clusters / 1 BX / cm^2",60,-30,30);
  h_hitDensity_bLayer = new TH1F("HitDensity_bLayer", ";Barrel Layer Index; Number of Hits / 1 BX / cm^2",30,0,30);
  h_hitDensity_eLayer = new TH1F("HitDensity_eLayer", ";Endcap Layer Index; Number of Hits / 1 BX / cm^2",60,-30,30);
  h_inPixPU     = new TH1F("Hits_inPixPU"          , ";#Hits in same pixel; Number of pixels" ,50,0,50);
  nEvtTotal = 0;

  if(_muDet==0 || _muDet==1){ //MuCol_v1 or MAIA_v0 geometry
    vxb_area[0] = 270.4; vxb_area[1] = 270.4; vxb_area[2] = 448.5; vxb_area[3] = 448.5; vxb_area[4] = 655.2; vxb_area[5] = 655.2; vxb_area[6] = 904.8; vxb_area[7] = 904.8;
    vxe_area[0] = 389.0; vxe_area[1] = 389.0; vxe_area[2] = 378.96; vxe_area[3] = 378.96; vxe_area[4] = 364.36; vxe_area[5] = 364.36; vxe_area[6] = 312.48; vxe_area[7] = 312.48;
    itb_area[0] = 8117.85; itb_area[1] = 22034.16; itb_area[2] = 51678.81;
    ite_area[0] = 6639.65; ite_area[1] = 10611.59; ite_area[2] = 10078.04; ite_area[3] = 9900.19; ite_area[4] = 9307.37; ite_area[5] = 8595.98; ite_area[6] = 8299.56;
    otb_area[0] = 140032.91; otb_area[1] = 194828.39; otb_area[2] = 249623.88;
    ote_area[0] = 69545.45; ote_area[1] = 69545.45; ote_area[2] = 69545.45; ote_area[3] = 69545.45;
  }
  if(_muDet==2){ //MuSIC_v2 geometry
    vxb_area[0] = 540.8; vxb_area[2] = 702.0; vxb_area[4] = 897.0; vxb_area[6] = 1310.4; vxb_area[8] = 1809.6;
    vxe_area[0] = 365.73; vxe_area[2] = 343.28; vxe_area[4] = 304.72; vxe_area[6] = 257.0;
    itb_area[0] = 10437.24; itb_area[1] = 22034.16; itb_area[2] = 51678.81;
    ite_area[0] = 6639.65; ite_area[1] = 10611.59; ite_area[2] = 10078.04; ite_area[3] = 9900.19; ite_area[4] = 9307.37; ite_area[5] = 8595.98; ite_area[6] = 8299.56;
    otb_area[0] = 140032.91; otb_area[1] = 194828.39; otb_area[2] = 286154.2;
    ote_area[0] = 69545.45; ote_area[1] = 69545.45; ote_area[2] = 69545.45; ote_area[3] = 69545.45;
  }
}

void ClusterShapeHistProc::processRunHeader( LCRunHeader* /*run*/)
{ } 


void ClusterShapeHistProc::processEvent( LCEvent * evt )
{
  //
  // Get object required collections and create lists
  // to keep track of unsaved objects.

  nEvtTotal++;

  // --- MCParticles

  LCCollection* mcpCol  =evt->getCollection(_mcpColName);

  if( mcpCol->getTypeName() != lcio::LCIO::MCPARTICLE )
    { throw EVENT::Exception( "Invalid collection type: " + mcpCol->getTypeName() ) ; }


  LCCollection* vbtrkhitCol = nullptr;
  if (_vbtrkhitColName != "") {
    vbtrkhitCol = evt->getCollection(_vbtrkhitColName);
  }

  LCCollection* ibtrkhitCol = nullptr;
  if (_ibtrkhitColName != "") {
    ibtrkhitCol = evt->getCollection(_ibtrkhitColName);
  }
  
  LCCollection* obtrkhitCol = nullptr;
  if (_obtrkhitColName != "") {
    obtrkhitCol = evt->getCollection(_obtrkhitColName);
  }
  
  LCCollection* vetrkhitCol = nullptr;
  if (_vetrkhitColName != "") {
    vetrkhitCol = evt->getCollection(_vetrkhitColName);
  }
  
  LCCollection* ietrkhitCol = nullptr;
  if (_ietrkhitColName != "") {
    ietrkhitCol = evt->getCollection(_ietrkhitColName);
  }
  
  LCCollection* oetrkhitCol = nullptr;
  if (_oetrkhitColName != "") {
    oetrkhitCol = evt->getCollection(_oetrkhitColName);
  }
  
  LCCollection* VBRelationCollection = nullptr;
  if (_VBRelationCollection != "") {
    VBRelationCollection = evt->getCollection(_VBRelationCollection);
  }
	
  
  LCCollection* IBRelationCollection = nullptr;
  if (_IBRelationCollection != "") {
    IBRelationCollection = evt->getCollection(_IBRelationCollection);
  }
	
  
  LCCollection* OBRelationCollection = nullptr;
  if (_OBRelationCollection != "") {
    OBRelationCollection = evt->getCollection(_OBRelationCollection);
  }
	
  
  LCCollection* VERelationCollection = nullptr;
  if (_VERelationCollection != "") {
    VERelationCollection = evt->getCollection(_VERelationCollection);
  }
	
  
  LCCollection* IERelationCollection = nullptr;
  if (_IERelationCollection != "") {
    IERelationCollection = evt->getCollection(_IERelationCollection);
  }
	
  
  LCCollection* OERelationCollection = nullptr;
  if (_OERelationCollection != "") {
    OERelationCollection = evt->getCollection(_OERelationCollection);
  }

  // --- Tracker hits multiplicity
  
  // vertex barrel tracker hits
  streamlog_out(DEBUG3) << "Num Events in VB Hit Collection: " << vbtrkhitCol->getNumberOfElements() << std::endl;
  int maxTrkHits=0;
  allPixels.clear();
  if (vbtrkhitCol) maxTrkHits = vbtrkhitCol->getNumberOfElements();
  for (int i=0; i<maxTrkHits; ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(vbtrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());
      streamlog_out(DEBUG9) << "Filling VB clusters with VB track hits..." << std::endl;
      _clusters_vb->fill(trkhit);
      LayerInfo(trkhit, 0);
    }
  for (const auto& pair : allPixels) {
    uint64_t key = pair.first;
    const std::vector<lcio::SimTrackerHit*>& vec = pair.second;
    h_inPixPU->Fill(vec.size());
  }
  _clusters_vb->h_cluster_edep_BX->Scale(1/_clusters_vb->h_cluster_edep_BX->Integral());

  // vertex endcap tracker hits
  streamlog_out(DEBUG3) << "Num Events in VE Hit Collection: " << vetrkhitCol->getNumberOfElements() << std::endl;
  maxTrkHits=0;
  if (vetrkhitCol) maxTrkHits = vetrkhitCol->getNumberOfElements();
  for (int i=0; i<maxTrkHits; ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(vetrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());
      streamlog_out(DEBUG9) << "Filling VE clusters with VE track hits..." << std::endl;
      _clusters_ve->fill(trkhit);
      LayerInfo(trkhit, 0);
    }

  // inner tracker barrel
  maxTrkHits=0;
  if (ibtrkhitCol) maxTrkHits = ibtrkhitCol->getNumberOfElements();
  for (int i=0; i<maxTrkHits; ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(ibtrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());
      streamlog_out(DEBUG9) << "Filling IB clusters with IB track hits..." << std::endl;
      _clusters_ib->fill(trkhit);
      LayerInfo(trkhit, 10);
    }

  // inner tracker endcap
  maxTrkHits=0;
  if (ietrkhitCol) maxTrkHits = ietrkhitCol->getNumberOfElements();
  for (int i=0; i<maxTrkHits; ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(ietrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());
      streamlog_out(DEBUG9) << "Filling IE clusters with IE track hits..." << std::endl;
      _clusters_ie->fill(trkhit);
      LayerInfo(trkhit, 10);
    }
  
  // outer tracker barrel
  maxTrkHits=0;
  if (obtrkhitCol) maxTrkHits = obtrkhitCol->getNumberOfElements();
  for (int i=0; i<maxTrkHits; ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(obtrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());
      streamlog_out(DEBUG9) << "Filling OB clusters with OB track hits..." << std::endl;
      _clusters_ob->fill(trkhit);
      LayerInfo(trkhit, 20);
    }      

  // outer tracker endcap  
  maxTrkHits=0;
  if (oetrkhitCol) maxTrkHits = oetrkhitCol->getNumberOfElements();
  for (int i=0; i<maxTrkHits; ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(oetrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());
      streamlog_out(DEBUG9) << "Filling OE clusters with OE track hits..." << std::endl;
      _clusters_oe->fill(trkhit);
      LayerInfo(trkhit, 20);
    }

   
  // --- tracker hit resolution histograms
  // vertex barrel resolution
  maxTrkHits=0;
  if (VBRelationCollection) maxTrkHits = VBRelationCollection->getNumberOfElements();
  for (int i=0; i<maxTrkHits; ++i)
    {
      streamlog_out(DEBUG3) << "Events in VB Relation Collection: " << VBRelationCollection->getNumberOfElements() << std::endl;
      EVENT::LCRelation *rel=static_cast<EVENT::LCRelation*>(VBRelationCollection->getElementAt(i));
      EVENT::TrackerHit *trkhit=dynamic_cast<EVENT::TrackerHit*>(rel->getFrom());
      EVENT::SimTrackerHit *simtrkhit=dynamic_cast<EVENT::SimTrackerHit*>(rel->getTo());
      IMPL::TrackerHitPlaneImpl *trkhitplane=dynamic_cast<IMPL::TrackerHitPlaneImpl*>(trkhit);
      if(trkhit==nullptr or simtrkhit==nullptr or trkhitplane==nullptr){
        std::cout << "Warning: Failed to dynamic cast to planar sensor" << std::endl;
        std::cout << "- Trackhit: " << trkhit << std::endl;
        std::cout << "- Simtrackhit: " << simtrkhit << std::endl;
        std::cout << "- Trackhitplane: " << trkhitplane << std::endl;
        continue;
      }
      _resolution_vb->fill(trkhit,simtrkhit,trkhitplane);
    }

  // vertex endcap resolution
  maxTrkHits=0;
  if (VERelationCollection) maxTrkHits = VERelationCollection->getNumberOfElements();
  for (int i=0; i<maxTrkHits; ++i)
    {
      streamlog_out(DEBUG3) << "Events in VE Relation Collection: " << VERelationCollection->getNumberOfElements() << std::endl;
      EVENT::LCRelation *rel=static_cast<EVENT::LCRelation*>(VERelationCollection->getElementAt(i));
      EVENT::TrackerHit *trkhit=dynamic_cast<EVENT::TrackerHit*>(rel->getFrom());
      EVENT::SimTrackerHit *simtrkhit=dynamic_cast<EVENT::SimTrackerHit*>(rel->getTo());
      IMPL::TrackerHitPlaneImpl *trkhitplane=dynamic_cast<IMPL::TrackerHitPlaneImpl*>(trkhit);
      if(trkhit==nullptr or simtrkhit==nullptr or trkhitplane==nullptr){
        std::cout << "Warning: Failed to dynamic cast to planar sensor" << std::endl;
        std::cout << "- Trackhit: " << trkhit << std::endl;
        std::cout << "- Simtrackhit: " << simtrkhit << std::endl;
        std::cout << "- Trackhitplane: " << trkhitplane << std::endl;
        continue;
      }
      _resolution_ve->fill(trkhit,simtrkhit,trkhitplane);
    }

  // inner tracker barrel hits resolution
  maxTrkHits=0;
  if (IBRelationCollection) maxTrkHits = IBRelationCollection->getNumberOfElements();
  for (int i=0; i<maxTrkHits; ++i)
    {
      streamlog_out(DEBUG3) << "Events in IB Relation Collection: " << IBRelationCollection->getNumberOfElements() << std::endl;
      EVENT::LCRelation *rel=static_cast<EVENT::LCRelation*>(IBRelationCollection->getElementAt(i));
      EVENT::TrackerHit *trkhit=dynamic_cast<EVENT::TrackerHit*>(rel->getFrom());
      EVENT::SimTrackerHit *simtrkhit=dynamic_cast<EVENT::SimTrackerHit*>(rel->getTo());
      IMPL::TrackerHitPlaneImpl *trkhitplane=dynamic_cast<IMPL::TrackerHitPlaneImpl*>(trkhit);
      if(trkhit==nullptr or simtrkhit==nullptr or trkhitplane==nullptr){
        std::cout << "Warning: Failed to dynamic cast to planar sensor" << std::endl;
        std::cout << "- Trackhit: " << trkhit << std::endl;
        std::cout << "- Simtrackhit: " << simtrkhit << std::endl;
        std::cout << "- Trackhitplane: " << trkhitplane << std::endl;
        continue;
      }
      _resolution_ib->fill(trkhit,simtrkhit,trkhitplane);
    }
  

  // outer tracker barrel hits resolution
  maxTrkHits=0;
  if (OBRelationCollection) maxTrkHits = OBRelationCollection->getNumberOfElements();
  for (int i=0; i<maxTrkHits; ++i)
    {
      streamlog_out(DEBUG3) << "Events in OB Relation Collection: " << OBRelationCollection->getNumberOfElements() << std::endl;
      EVENT::LCRelation *rel=static_cast<EVENT::LCRelation*>(OBRelationCollection->getElementAt(i));
      EVENT::TrackerHit *trkhit=dynamic_cast<EVENT::TrackerHit*>(rel->getFrom());
      EVENT::SimTrackerHit *simtrkhit=dynamic_cast<EVENT::SimTrackerHit*>(rel->getTo());
      IMPL::TrackerHitPlaneImpl *trkhitplane=dynamic_cast<IMPL::TrackerHitPlaneImpl*>(trkhit);
      if(trkhit==nullptr or simtrkhit==nullptr or trkhitplane==nullptr){
        std::cout << "Warning: Failed to dynamic cast to planar sensor" << std::endl;
        std::cout << "- Trackhit: " << trkhit << std::endl;
        std::cout << "- Simtrackhit: " << simtrkhit << std::endl;
        std::cout << "- Trackhitplane: " << trkhitplane << std::endl;
        continue;
      }
      _resolution_ob->fill(trkhit,simtrkhit,trkhitplane);
    }


  // inner tracker endcap hits resolution
  maxTrkHits=0;
  if (IERelationCollection) maxTrkHits = IERelationCollection->getNumberOfElements();
  for (int i=0; i<maxTrkHits; ++i)
    {
      streamlog_out(DEBUG3) << "Events in IE Relation Collection: " << IERelationCollection->getNumberOfElements() << std::endl;
      EVENT::LCRelation *rel=static_cast<EVENT::LCRelation*>(IERelationCollection->getElementAt(i));
      EVENT::TrackerHit *trkhit=dynamic_cast<EVENT::TrackerHit*>(rel->getFrom());
      EVENT::SimTrackerHit *simtrkhit=dynamic_cast<EVENT::SimTrackerHit*>(rel->getTo());
      IMPL::TrackerHitPlaneImpl *trkhitplane=dynamic_cast<IMPL::TrackerHitPlaneImpl*>(trkhit);
      if(trkhit==nullptr or simtrkhit==nullptr or trkhitplane==nullptr){
        std::cout << "Warning: Failed to dynamic cast to planar sensor" << std::endl;
        std::cout << "- Trackhit: " << trkhit << std::endl;
        std::cout << "- Simtrackhit: " << simtrkhit << std::endl;
        std::cout << "- Trackhitplane: " << trkhitplane << std::endl;
        continue;
      }
      _resolution_ie->fill(trkhit,simtrkhit,trkhitplane);
    }


  // outer tracker endcap hits resolution
  maxTrkHits=0;
  if (OERelationCollection) maxTrkHits = OERelationCollection->getNumberOfElements();
  for (int i=0; i<maxTrkHits; ++i)
    {
      streamlog_out(DEBUG3) << "Events in OE Relation Collection: " << OERelationCollection->getNumberOfElements() << std::endl;
      EVENT::LCRelation *rel=static_cast<EVENT::LCRelation*>(OERelationCollection->getElementAt(i));
      EVENT::TrackerHit *trkhit=dynamic_cast<EVENT::TrackerHit*>(rel->getFrom());
      EVENT::SimTrackerHit *simtrkhit=dynamic_cast<EVENT::SimTrackerHit*>(rel->getTo());
      IMPL::TrackerHitPlaneImpl *trkhitplane=dynamic_cast<IMPL::TrackerHitPlaneImpl*>(trkhit);
      if(trkhit==nullptr or simtrkhit==nullptr or trkhitplane==nullptr){
        std::cout << "Warning: Failed to dynamic cast to planar sensor" << std::endl;
        std::cout << "- Trackhit: " << trkhit << std::endl;
        std::cout << "- Simtrackhit: " << simtrkhit << std::endl;
        std::cout << "- Trackhitplane: " << trkhitplane << std::endl;
        continue;
      }
      _resolution_oe->fill(trkhit,simtrkhit,trkhitplane);
    }
}

void ClusterShapeHistProc::LayerInfo(const EVENT::TrackerHit* trkhit, int offset)
{
  const lcio::LCObjectVec &rawHits = trkhit->getRawHits();
  float z = trkhit->getPosition()[2];
  float loopsize = rawHits.size();

  //Get hit layer
  std::string _encoderString = lcio::LCTrackerCellID::encoding_string();
  UTIL::CellIDDecoder<lcio::TrackerHit> decoder(_encoderString);
  uint32_t systemID = decoder(trkhit)["system"];
  uint32_t layerID = decoder(trkhit)["layer"];
  int z_dir = z > 0. ? 1 : -1;

  if(systemID%2==0){ //Endcap layers
    h_clusters_by_eLayer->Fill(z_dir*(layerID+1) + z_dir*offset);
    h_clusters_by_eLayer_BX->Fill(z_dir*(layerID+1) + z_dir*offset);
    h_hits_by_eLayer->Fill(z_dir*(layerID+1) + z_dir*offset, std::min((int)loopsize,30));
    h_hits_by_eLayer_BX->Fill(z_dir*(layerID+1) + z_dir*offset, std::min((int)loopsize, 30));
    h_clusterDensity_eLayer->Fill(z_dir*(layerID+1) + z_dir*offset);
    h_hitDensity_eLayer->Fill(z_dir*(layerID+1) + z_dir*offset, std::min((int)loopsize, 30));
  }
  else{ //barrel layers
    h_clusters_by_bLayer->Fill(layerID + offset);
    h_clusters_by_bLayer_BX->Fill(layerID + offset);
    h_hits_by_bLayer->Fill(layerID + offset, std::min((int)loopsize,30));
    h_hits_by_bLayer_BX->Fill(layerID + offset, std::min((int)loopsize, 30));
    h_clusterDensity_bLayer->Fill(layerID + offset);
    h_hitDensity_bLayer->Fill(layerID + offset, std::min((int)loopsize, 30));
  }

  if(layerID==0 && systemID==1){ //only first layer of vertex barrel
    for (size_t j=0; j<loopsize; ++j) {
      lcio::SimTrackerHit *hitConstituent = dynamic_cast<lcio::SimTrackerHit*>( rawHits[j] );
      const double *localPos = hitConstituent->getPosition();
      int16_t x_local = localPos[0]>0. ? static_cast<int16_t>(localPos[0]) : static_cast<int16_t>(localPos[0]+10000);
      int16_t y_local = localPos[1]>0. ? static_cast<int16_t>(localPos[1]) : static_cast<int16_t>(localPos[1]+10000);
      UTIL::CellIDDecoder<lcio::SimTrackerHit> simdecoder(_encoderString);
      int16_t ladderID = simdecoder(hitConstituent)["module"];
      int16_t modID = 0;
      if(z>=-130. && z<-104.)
	modID = 0;
      else if(z>=-104. && z<-78.)
	modID = 1;
      else if(z>=-78. && z<-52.)
	modID = 2;
      else if(z>=-52. && z<-26.)
	modID = 3;
      else if(z>=-26. && z<0)
	modID = 4;
      else if(z>=0. && z<26.)
	modID = 5;
      else if(z>=26. && z<52.)
	modID = 6;
      else if(z>=52. && z<78.)
	modID = 7;
      else if(z>=78. && z<104.)
	modID = 8;
      else
	modID = 9;

      uint64_t pixel_hash=0;
      pixel_hash |= ((uint64_t)x_local)<<48;
      pixel_hash |= ((uint64_t)y_local)<<32;
      pixel_hash |= ((uint64_t)modID)<<16;
      pixel_hash |= (uint64_t)ladderID;

      if(allPixels.find(pixel_hash) == allPixels.end())
	allPixels[pixel_hash] = {};

      allPixels[pixel_hash].push_back(hitConstituent);
    }
  }


  return;
}
void ClusterShapeHistProc::check( LCEvent * /*evt*/ )
{ }

void ClusterShapeHistProc::end()
{
  std::cout<<"Total events = "<<nEvtTotal<<std::endl;
  h_clusters_by_bLayer_BX->Scale(1.0/nEvtTotal);
  h_hits_by_bLayer_BX->Scale(1.0/nEvtTotal);
  h_clusters_by_eLayer_BX->Scale(1.0/nEvtTotal);
  h_hits_by_eLayer_BX->Scale(1.0/nEvtTotal);
  h_clusterDensity_bLayer->Scale(1.0/nEvtTotal);
  h_clusterDensity_eLayer->Scale(1.0/nEvtTotal);
  h_hitDensity_bLayer->Scale(1.0/nEvtTotal);
  h_hitDensity_eLayer->Scale(1.0/nEvtTotal);

  float x_val;
  int z_dir, layerID;

  //barrel density
  for(int ibin=1; ibin<=h_clusters_by_bLayer_BX->GetNbinsX(); ibin++){
    x_val = h_clusters_by_bLayer_BX->GetBinLowEdge(ibin);
    if(x_val<10 && x_val>0 && h_clusters_by_bLayer_BX->GetBinContent(ibin)>0){
      layerID = x_val;
      h_clusterDensity_bLayer->SetBinContent(ibin, h_clusterDensity_bLayer->GetBinContent(ibin)/vxb_area.at(layerID));
      h_clusterDensity_bLayer->SetBinError(ibin, h_clusterDensity_bLayer->GetBinError(ibin)/vxb_area.at(layerID));
      h_hitDensity_bLayer->SetBinContent(ibin, h_hitDensity_bLayer->GetBinContent(ibin)/vxb_area.at(layerID));
      h_hitDensity_bLayer->SetBinError(ibin, h_hitDensity_bLayer->GetBinError(ibin)/vxb_area.at(layerID));
    }
    else if(x_val<20 && x_val>10 && h_clusters_by_bLayer_BX->GetBinContent(ibin)>0){
      layerID = x_val - 10;
      h_clusterDensity_bLayer->SetBinContent(ibin, h_clusterDensity_bLayer->GetBinContent(ibin)/itb_area.at(layerID));
      h_clusterDensity_bLayer->SetBinError(ibin, h_clusterDensity_bLayer->GetBinError(ibin)/itb_area.at(layerID));
      h_hitDensity_bLayer->SetBinContent(ibin, h_hitDensity_bLayer->GetBinContent(ibin)/itb_area.at(layerID));
      h_hitDensity_bLayer->SetBinError(ibin, h_hitDensity_bLayer->GetBinError(ibin)/itb_area.at(layerID));
    }
    else if(x_val<30 && x_val>20 && h_clusters_by_bLayer_BX->GetBinContent(ibin)>0){
      layerID = x_val - 20;
      h_clusterDensity_bLayer->SetBinContent(ibin, h_clusterDensity_bLayer->GetBinContent(ibin)/otb_area.at(layerID));
      h_clusterDensity_bLayer->SetBinError(ibin, h_clusterDensity_bLayer->GetBinError(ibin)/otb_area.at(layerID));
      h_hitDensity_bLayer->SetBinContent(ibin, h_hitDensity_bLayer->GetBinContent(ibin)/otb_area.at(layerID));
      h_hitDensity_bLayer->SetBinError(ibin, h_hitDensity_bLayer->GetBinError(ibin)/otb_area.at(layerID));
    }
    else{}
  }
    
  //endcap density
  for(int ibin=1; ibin<=h_clusters_by_eLayer_BX->GetNbinsX(); ibin++){
    x_val = h_clusters_by_eLayer_BX->GetBinCenter(ibin);
    z_dir = x_val > 0 ? 1 : -1;
    if(((x_val>0 && x_val<10) || (x_val<0 && x_val>-10)) && h_clusters_by_eLayer_BX->GetBinContent(ibin)>0){
      layerID = int(x_val/z_dir) - 1;
      h_clusterDensity_eLayer->SetBinContent(ibin, h_clusterDensity_eLayer->GetBinContent(ibin)/vxe_area.at(layerID));
      h_clusterDensity_eLayer->SetBinError(ibin, h_clusterDensity_eLayer->GetBinError(ibin)/vxe_area.at(layerID));
      h_hitDensity_eLayer->SetBinContent(ibin, h_hitDensity_eLayer->GetBinContent(ibin)/vxe_area.at(layerID));
      h_hitDensity_eLayer->SetBinError(ibin, h_hitDensity_eLayer->GetBinError(ibin)/vxe_area.at(layerID));
    }
    else if(((x_val>10 && x_val<20) || (x_val<-10 && x_val>-20)) && h_clusters_by_eLayer_BX->GetBinContent(ibin)>0){
      layerID = int((x_val - z_dir*10)/z_dir) - 1;
      h_clusterDensity_eLayer->SetBinContent(ibin, h_clusterDensity_eLayer->GetBinContent(ibin)/ite_area.at(layerID));
      h_clusterDensity_eLayer->SetBinError(ibin, h_clusterDensity_eLayer->GetBinError(ibin)/ite_area.at(layerID));
      h_hitDensity_eLayer->SetBinContent(ibin, h_hitDensity_eLayer->GetBinContent(ibin)/ite_area.at(layerID));
      h_hitDensity_eLayer->SetBinError(ibin, h_hitDensity_eLayer->GetBinError(ibin)/ite_area.at(layerID));
    }
    else if(((x_val>20 && x_val<30) || (x_val<-20 && x_val>-30)) && h_clusters_by_eLayer_BX->GetBinContent(ibin)>0){
      layerID = int((x_val - z_dir*20)/z_dir) - 1;
      h_clusterDensity_eLayer->SetBinContent(ibin, h_clusterDensity_eLayer->GetBinContent(ibin)/ote_area.at(layerID));
      h_clusterDensity_eLayer->SetBinError(ibin, h_clusterDensity_eLayer->GetBinError(ibin)/ote_area.at(layerID));
      h_hitDensity_eLayer->SetBinContent(ibin, h_hitDensity_eLayer->GetBinContent(ibin)/ote_area.at(layerID));
      h_hitDensity_eLayer->SetBinError(ibin, h_hitDensity_eLayer->GetBinError(ibin)/ote_area.at(layerID));
    }
    else{}
  }

}
