#pragma once

#include <TH2.h>
#include <TH1.h>

//#include <ACTSTracking/GeometryIdMappingTool.hxx>

namespace EVENT
{
  class TrackerHit;
}

//! Histograms for cluster analysis
class ClusterHists
{
public:
  ClusterHists(const ClusterHists &) = delete ;
  ClusterHists& operator =(const ClusterHists &) = delete ;

  //! Initialize empty histograms
  ClusterHists() ;

  // Fill histograms with a single track hit
  void fill(const EVENT::TrackerHit* trkhit);

  TH1* h_cluster_edep_BX;

private:
  TH2* h_size_theta_y;
  TH2* h_size_theta_x;
  TH2* h_size_theta_tot;
  TH2* h_size_theta_tot_layer[9];
  TH2* h_size_theta_x_layer[9];
  TH2* h_size_theta_y_layer[9];
  TH2* h_size_r_tot;
  TH2* h_size_r_tot_0;
  TH2* h_size_r_tot_1;
  TH2* h_size_r_tot_2;
  TH2* h_size_r_tot_3;
  TH2* h_cluster_pos;
  TH2* h_cluster_pos_0;
  TH2* h_cluster_pos_1;
  TH2* h_cluster_pos_2;
  TH2* h_cluster_pos_3;
  TH1* h_clusters_by_layer;
  TH1* h_hits_by_layer;
  TH1* h_theta;
  TH1* h_cluster_edep;
  TH1* h_hit_edep;
  TH2* h_edep_r;
  TH2* h_edep_cluster;
  /*TH1* h_edep_0deg;
  TH1* h_edep_90deg; */
  // positions
  TH1* h_x;
  TH1* h_y;
  TH1* h_z;
  TH1* h_r;
  TH2* h_z_r;
  TH2* h_x_y;
  TH2* h_z_r_hits;
  TH2* h_x_y_hits;
  TH1* h_z_layer[9];
  TH1* h_r_layer[9];
  TH1* h_avgHits;
  TH1* h_sysID;
  TH1* h_cluster_timing;
  TH1* h_hit_timing;
  TH1* h_hits_layer[9];
  TH1* h_cluster_layer[9];

  // vertex
  TH1* h_x_vx;  
  TH1* h_y_vx;  
  TH1* h_z_vx;  
  TH1* h_r_vx;
  TH2* h_z_r_vx;
  TH2* h_x_y_vx;

};
