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
			   _vbtrkhitColName
			   );     

  registerInputCollection( LCIO::TRACKERHIT,
			   "IBTrackerHitsCollection" , 
			   "Name of inner barrel tracker hits collection",
			   _ibtrkhitColName,
			   _ibtrkhitColName
			   );  

  registerInputCollection( LCIO::TRACKERHIT,
			   "OBTrackerHitsCollection" , 
			   "Name of outer barrel tracker hits collection",
			   _obtrkhitColName,
			   _obtrkhitColName
			   );  

  registerInputCollection( LCIO::TRACKERHIT,
			   "VETrackerHitsCollection" , 
			   "Name of vertex endcap tracker hits collection",
			   _vetrkhitColName,
			   _vetrkhitColName
			   );  

  registerInputCollection( LCIO::TRACKERHIT,
			   "IETrackerHitsCollection" , 
			   "Name of inner endcap tracker hits collection",
			   _ietrkhitColName,
			   _ietrkhitColName
			   );     

  registerInputCollection( LCIO::TRACKERHIT,
			   "OETrackerHitsCollection" , 
			   "Name of outer endcap tracker hits collection",
			   _oetrkhitColName,
			   _oetrkhitColName
			   );
  // --- Collection of MC particles
  registerInputCollection( LCIO::MCPARTICLE,
			   "MCParticleCollection" , 
			   "Name of the MCParticle collection"  ,
			   _mcpColName,
			   _mcpColName
			   );  

  // --- Relationship reco-truth for tracker clusters
  registerInputCollection( LCIO::LCRELATION,
		     "VBRelationCollection" ,
			   "Name of the input vertex barrel relation collection",
			   _VBRelationCollection,
		     _VBRelationCollection
		 	    );

    registerInputCollection( LCIO::LCRELATION,
		     "IBRelationCollection" ,
			   "Name of the input inner tracker barrel relation collection",
			   _IBRelationCollection,
		     _IBRelationCollection
		 	    );

    registerInputCollection( LCIO::LCRELATION,
		     "OBRelationCollection" ,
			   "Name of the input outer tracker barrel relation collection",
			   _OBRelationCollection,
		     _OBRelationCollection
		 	    );

    registerInputCollection( LCIO::LCRELATION,
		     "VERelationCollection" ,
			   "Name of the input vertex endcap relation collection",
			   _VERelationCollection,
		     _VERelationCollection
		 	    );

    registerInputCollection( LCIO::LCRELATION,
		     "IERelationCollection" ,
			   "Name of the input inner tracker endcap relation collection",
			   _IERelationCollection,
		     _IERelationCollection
		 	    );

    registerInputCollection( LCIO::LCRELATION,
		     "OERelationCollection" ,
			   "Name of the input outer tracker endcap relation collection",
			   _OERelationCollection,
		     _OERelationCollection
		 	    );
    
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

}

void ClusterShapeHistProc::processRunHeader( LCRunHeader* /*run*/)
{ } 


void ClusterShapeHistProc::processEvent( LCEvent * evt )
{
  //
  // Get object required collections and create lists
  // to keep track of unsaved objects.

  // --- MCParticles

  LCCollection* mcpCol  =evt->getCollection(_mcpColName);

  if( mcpCol->getTypeName() != lcio::LCIO::MCPARTICLE )
    { throw EVENT::Exception( "Invalid collection type: " + mcpCol->getTypeName() ) ; }

  LCCollection* vbtrkhitCol  =evt->getCollection(_vbtrkhitColName);
  LCCollection* ibtrkhitCol  =evt->getCollection(_ibtrkhitColName);
  LCCollection* obtrkhitCol  =evt->getCollection(_obtrkhitColName);
  LCCollection* vetrkhitCol  =evt->getCollection(_vetrkhitColName);
  LCCollection* ietrkhitCol  =evt->getCollection(_ietrkhitColName);
  LCCollection* oetrkhitCol  =evt->getCollection(_oetrkhitColName);
  LCCollection* VBRelationCollection =evt->getCollection(_VBRelationCollection);
  LCCollection* IBRelationCollection =evt->getCollection(_IBRelationCollection);
  LCCollection* OBRelationCollection =evt->getCollection(_OBRelationCollection);
  LCCollection* VERelationCollection =evt->getCollection(_VERelationCollection);
  LCCollection* IERelationCollection =evt->getCollection(_IERelationCollection);
  LCCollection* OERelationCollection =evt->getCollection(_OERelationCollection);

  // --- Tracker hits multiplicity
  
  // vertex barrel tracker hits
  streamlog_out(DEBUG3) << "Num Events in VB Hit Collection: " << vbtrkhitCol->getNumberOfElements() << std::endl;
  for(int i=0; i<vbtrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(vbtrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());
      streamlog_out(DEBUG9) << "Filling VB clusters with VB track hits..." << std::endl;
      _clusters_vb->fill(trkhit);}


  // vertex endcap tracker hits
  streamlog_out(DEBUG3) << "Num Events in VE Hit Collection: " << vetrkhitCol->getNumberOfElements() << std::endl;
  for(int i=0; i<vetrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(vetrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());
      streamlog_out(DEBUG9) << "Filling VE clusters with VE track hits..." << std::endl;
      _clusters_ve->fill(trkhit);}

  // inner tracker barrel
  for(int i=0; i<ibtrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(ibtrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());
      streamlog_out(DEBUG9) << "Filling IB clusters with IB track hits..." << std::endl;
      _clusters_ib->fill(trkhit);}

  // inner tracker endcap
  for(int i=0; i<ietrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(ietrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());
      streamlog_out(DEBUG9) << "Filling IE clusters with IE track hits..." << std::endl;
      _clusters_ie->fill(trkhit);}
  
  // outer tracker barrel
  for(int i=0; i<obtrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(obtrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());
      streamlog_out(DEBUG9) << "Filling OB clusters with OB track hits..." << std::endl;
      _clusters_ob->fill(trkhit);}      

  // outer tracker endcap      
  for(int i=0; i<oetrkhitCol->getNumberOfElements(); ++i)
    {
      const EVENT::TrackerHit *trkhit=static_cast<const EVENT::TrackerHit*>(oetrkhitCol->getElementAt(i));
      h_trackerhit_timing -> Fill(trkhit->getTime());
      streamlog_out(DEBUG9) << "Filling OE clusters with OE track hits..." << std::endl;
      _clusters_oe->fill(trkhit);
    }

   
  // --- tracker hit resolution histograms

  // vertex barrel resolution
  for(int i=0; i<VBRelationCollection->getNumberOfElements(); ++i)
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
  for(int i=0; i<VERelationCollection->getNumberOfElements(); ++i)
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
  for(int i=0; i<IBRelationCollection->getNumberOfElements(); ++i)
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
  for(int i=0; i<OBRelationCollection->getNumberOfElements(); ++i)
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
  for(int i=0; i<IERelationCollection->getNumberOfElements(); ++i)
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
  for(int i=0; i<OERelationCollection->getNumberOfElements(); ++i)
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

void ClusterShapeHistProc::check( LCEvent * /*evt*/ )
{ }

void ClusterShapeHistProc::end()
{ }
