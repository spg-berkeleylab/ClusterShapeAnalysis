#pragma once

#include <TH3.h>
#include <TH2.h>
#include <TH1.h>

//#include <ACTSTracking/GeometryIdMappingTool.hxx>

namespace EVENT
{
  class TrackerHit;
  class SimTrackerHit;
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
  void fill(const EVENT::SimTrackerHit* simtrkhit);

  TH1* h_cluster_edep_BX;

private:
  TH2* h_size_theta_y;
  TH2* h_size_theta_x;
  TH2* h_size_theta_tot;
  TH2* h_size_theta_tot_0;
  TH2* h_size_theta_tot_1;
  TH2* h_size_theta_tot_2;
  TH2* h_size_theta_tot_3;
  TH2* h_size_theta_tot_4;
  TH2* h_size_theta_tot_5;
  TH2* h_size_theta_tot_6;
  TH2* h_size_theta_tot_7;
  TH2* h_size_theta_tot_8;
  TH2* h_size_theta_x_0;
  TH2* h_size_theta_x_1;
  TH2* h_size_theta_x_2;
  TH2* h_size_theta_x_3;
  TH2* h_size_theta_x_4;
  TH2* h_size_theta_x_5;
  TH2* h_size_theta_x_6;
  TH2* h_size_theta_x_7;
  TH2* h_size_theta_x_8;
  TH2* h_size_theta_y_0;
  TH2* h_size_theta_y_1;
  TH2* h_size_theta_y_2;
  TH2* h_size_theta_y_3;
  TH2* h_size_theta_y_4;
  TH2* h_size_theta_y_5;
  TH2* h_size_theta_y_6;
  TH2* h_size_theta_y_7;
  TH2* h_size_theta_y_8;
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
  TH1* h_z_layer0;
  TH1* h_z_layer1;
  TH1* h_z_layer2;
  TH1* h_z_layer3;
  TH1* h_z_layer4;
  TH1* h_z_layer5;
  TH1* h_z_layer6;
  TH1* h_z_layer7;
  TH1* h_z_layer8;
  TH1* h_r_layer0;
  TH1* h_r_layer1;
  TH1* h_r_layer2;
  TH1* h_r_layer3;
  TH1* h_r_layer4;
  TH1* h_r_layer5;
  TH1* h_r_layer6;
  TH1* h_r_layer7;
  TH1* h_r_layer8;
  TH1* h_avgHits;
  TH1* h_sysID;
  TH1* h_cluster_timing;
  TH1* h_hit_timing;
  
  // vertex
  TH1* h_x_vx;  
  TH1* h_y_vx;  
  TH1* h_z_vx;  
  TH1* h_r_vx;
  TH2* h_z_r_vx;
  TH2* h_x_y_vx;

  //time and edep cut for all layers -- juliet
  TH1* h_cluster_edep_Tcut;
  TH1* h_cluster_edep_Tcut_layer0;
  TH1* h_cluster_edep_Tcut_layer1;
  TH1* h_cluster_edep_Tcut_layer2;   
  TH1* h_cluster_edep_Tcut_layer3;
  TH1* h_cluster_edep_Tcut_layer4;
  TH1* h_cluster_edep_Tcut_layer5;
  TH1* h_cluster_edep_Tcut_layer6;
  TH1* h_cluster_edep_Tcut_layer7;
  TH1* h_cluster_edep_Tcut_layer8;

  TH1* h_trackerhit_time_Tcut;
  TH1* h_trackerhit_time_Tcut_layer0;
  TH1* h_trackerhit_time_Tcut_layer1;
  TH1* h_trackerhit_time_Tcut_layer2;
  TH1* h_trackerhit_time_Tcut_layer3;
  TH1* h_trackerhit_time_Tcut_layer4;
  TH1* h_trackerhit_time_Tcut_layer5;
  TH1* h_trackerhit_time_Tcut_layer6;
  TH1* h_trackerhit_time_Tcut_layer7;

  //time and edep for 1st layer -- juliet
  TH1* h_cluster_edep_layer0;
  TH1* h_cluster_edep_layer1;
  TH1* h_cluster_edep_layer2;   
  TH1* h_cluster_edep_layer3;
  TH1* h_cluster_edep_layer4;
  TH1* h_cluster_edep_layer5;
  TH1* h_cluster_edep_layer6;
  TH1* h_cluster_edep_layer7;
  TH1* h_cluster_edep_layer8;

//hit energy depo per layer
  TH1* h_hit_edep_layer0;
  TH1* h_hit_edep_layer1;
  TH1* h_hit_edep_layer2;   
  TH1* h_hit_edep_layer3;
  TH1* h_hit_edep_layer4;
  TH1* h_hit_edep_layer5;
  TH1* h_hit_edep_layer6;
  TH1* h_hit_edep_layer7;
  TH1* h_hit_edep_layer8;
//time of arriver per layer
  TH1* h_trackerhit_time;
  TH1* h_trackerhit_time_layer0;
  TH1* h_trackerhit_time_layer1;
  TH1* h_trackerhit_time_layer2;
  TH1* h_trackerhit_time_layer3;
  TH1* h_trackerhit_time_layer4;
  TH1* h_trackerhit_time_layer5;
  TH1* h_trackerhit_time_layer6;
  TH1* h_trackerhit_time_layer7;
  TH1* h_trackerhit_time_layer8;
  //hits per cluster per layer: 
  TH1* h_thclen;
  TH1* h_thclen_layer0;
  TH1* h_thclen_layer1;
  TH1* h_thclen_layer2;
  TH1* h_thclen_layer3;
  TH1* h_thclen_layer4;
  TH1* h_thclen_layer5;
  TH1* h_thclen_layer6;
  TH1* h_thclen_layer7;
  TH1* h_thclen_layer8;
  //hits per cluster per layer cut: 
  TH1* h_thclen_cut;
  TH1* h_thclen_layer0_cut;
  TH1* h_thclen_layer1_cut;
  TH1* h_thclen_layer2_cut;
  TH1* h_thclen_layer3_cut;
  TH1* h_thclen_layer4_cut;
  TH1* h_thclen_layer5_cut;
  TH1* h_thclen_layer6_cut;
  TH1* h_thclen_layer7_cut;
  TH1* h_thclen_layer8_cut;
  //2D
  TH2* h_toa_edepCluster;
  TH2* h_cluster_edep_thlen;
  TH2* h_cluster_edep_thlen_cut;
  //3D Histos: 
  TH3* h_3DPosition_digi;
  TH3* h_3DPosition_cdigi;  
  
};
