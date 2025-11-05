#pragma once

#include <marlin/Processor.h>
#include <EVENT/SimTrackerHit.h>
#include <UTIL/LCTrackerConf.h>

#include <TH1.h>

namespace EVENT
{
  class TrackerHit;
}

class ClusterHists;
class TrackerHitResoHists;

//! Creates a simple column wise ntuple in a HistProc from LCIO collections.
class ClusterShapeHistProc : public marlin::Processor
{
public:
  virtual Processor*  newProcessor() { return new ClusterShapeHistProc ; }

  ClusterShapeHistProc(const ClusterShapeHistProc &) = delete ;
  ClusterShapeHistProc& operator =(const ClusterShapeHistProc &) = delete ;
  ClusterShapeHistProc() ;

  /** Called at the begin of the job before anything is read.
   * Use to initialize the processor, e.g. book histograms.
   */
  virtual void init() ;

  /** Called for every run.
   */
  virtual void processRunHeader( LCRunHeader* run ) ;

  /** Called for every event - the working horse.
   */
  virtual void processEvent( LCEvent * evt ) ; 

  virtual void check( LCEvent * evt ) ; 

  /** Called after data processing for clean up.
   */
  virtual void end() ;

  void LayerInfo(const EVENT::TrackerHit* trkhit, int offset);
  
private:
  //! Tracker hit collections
  std::string _vbtrkhitColName {};
  std::string _ibtrkhitColName {};
  std::string _obtrkhitColName {};
  std::string _vetrkhitColName {};
  std::string _ietrkhitColName {};
  std::string _oetrkhitColName {};

  //! MC Particle Collection
  std::string _mcpColName {};

  //! Tracker hit relation collection
  std::string _VBRelationCollection {};
  std::string _IBRelationCollection {};
  std::string _OBRelationCollection {};
  std::string _VERelationCollection {};
  std::string _IERelationCollection {};
  std::string _OERelationCollection {};  

  //! Histograms
  std::shared_ptr<ClusterHists> _clusters_vb;
  std::shared_ptr<ClusterHists> _clusters_ve;
  std::shared_ptr<ClusterHists> _clusters_ib;
  std::shared_ptr<ClusterHists> _clusters_ie;
  std::shared_ptr<ClusterHists> _clusters_ob;
  std::shared_ptr<ClusterHists> _clusters_oe;
  std::shared_ptr<TrackerHitResoHists> _resolution_vb;
  std::shared_ptr<TrackerHitResoHists> _resolution_ve;
  std::shared_ptr<TrackerHitResoHists> _resolution_ib;
  std::shared_ptr<TrackerHitResoHists> _resolution_ie;
  std::shared_ptr<TrackerHitResoHists> _resolution_ob;
  std::shared_ptr<TrackerHitResoHists> _resolution_oe;

  TH1* h_trackerhit_timing;  
  TH1* h_clusters_by_bLayer;
  TH1* h_hits_by_bLayer;
  TH1* h_clusters_by_bLayer_BX;
  TH1* h_hits_by_bLayer_BX;
  TH1* h_clusters_by_eLayer;
  TH1* h_hits_by_eLayer;
  TH1* h_clusters_by_eLayer_BX;
  TH1* h_hits_by_eLayer_BX;
  TH1* h_clusterDensity_bLayer;
  TH1* h_clusterDensity_eLayer;
  TH1* h_hitDensity_bLayer;
  TH1* h_hitDensity_eLayer;
  TH1* h_inPixPU[9];
  TH1* h_inPixPUTimeDiff[9];
  
  int nEvtTotal;
  int _muDet;

  //2D map of layer ID and active sensor area
  std::map<int, double> vxb_area;
  std::map<int, double> vxe_area;
  std::map<int, double> itb_area;
  std::map<int, double> ite_area;
  std::map<int, double> otb_area;
  std::map<int, double> ote_area;

  struct PixelData {
    int layer;
    std::vector<lcio::SimTrackerHit*> hits;
  };

  std::map<uint64_t, PixelData> allPixels;
};
