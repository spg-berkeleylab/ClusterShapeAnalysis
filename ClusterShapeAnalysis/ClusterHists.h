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

private:
  TH2* h_size_theta_y;
  TH2* h_size_theta_x;
  TH2* h_size_theta_tot;
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

  // vertex
  TH1* h_x_vx;  
  TH1* h_y_vx;  
  TH1* h_z_vx;  
  TH1* h_r_vx;
  TH2* h_z_r_vx;
  TH2* h_x_y_vx;

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

//time and edep for 1st layer -- juliet
  TH1* h_truth_cluster_edep_layer0;
  TH1* h_truth_cluster_edep_layer1;
  TH1* h_truth_cluster_edep_layer2;   
  TH1* h_truth_cluster_edep_layer3;
  TH1* h_truth_cluster_edep_layer4;
  TH1* h_truth_cluster_edep_layer5;
  TH1* h_truth_cluster_edep_layer6;
  TH1* h_truth_cluster_edep_layer7;
  TH1* h_truth_cluster_edep_layer8;

  //normalized cluster edp:
  TH1* h_cluster_edep_norm_layer0;
  TH1* h_cluster_edep_norm_layer1;
  TH1* h_cluster_edep_norm_layer2;   
  TH1* h_cluster_edep_norm_layer3;
  TH1* h_cluster_edep_norm_layer4;
  TH1* h_cluster_edep_norm_layer5;
  TH1* h_cluster_edep_norm_layer6;
  TH1* h_cluster_edep_norm_layer7;
  TH1* h_cluster_edep_norm_layer8;

  TH1* h_hit_edep_layer0;
  TH1* h_hit_edep_layer1;
  TH1* h_hit_edep_layer2;   
  TH1* h_hit_edep_layer3;
  TH1* h_hit_edep_layer4;
  TH1* h_hit_edep_layer5;
  TH1* h_hit_edep_layer6;
  TH1* h_hit_edep_layer7;
  TH1* h_hit_edep_layer8;

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
  TH1* h_thclen_layer0;
  TH1* h_thclen_layer1;
  TH1* h_thclen_layer2;
  TH1* h_thclen_layer3;
  TH1* h_thclen_layer4;
  TH1* h_thclen_layer5;
  TH1* h_thclen_layer6;
  TH1* h_thclen_layer7;
  TH1* h_thclen_layer8;
  //New stuff
  TH1* h_thclen;
  TH1* h_cluster_1hits;
  TH1* h_cluster_2hits; 
  TH1* h_cluster_3hits; 
  TH1* h_cluster_4hits;
  TH1* h_cluster_5hits; 
  TH1* h_cluster_6hits; 
  TH1* h_cluster_7hits;
  TH1* h_cluster_8hits; 
  TH1* h_cluster_9hits; 
  TH1* h_cluster_norm;
 
  //diff histos in cluster and hit edep: 
  TH1* h_diffHitCluster_edep_layer0;
  TH1* h_diffHitCluster_edep_layer1;
  TH1* h_diffHitCluster_edep_layer2;
  TH1* h_diffHitCluster_edep_layer3;
  TH1* h_diffHitCluster_edep_layer4;
  TH1* h_diffHitCluster_edep_layer5;
  TH1* h_diffHitCluster_edep_layer6;
  TH1* h_diffHitCluster_edep_layer7;
  TH1* h_diffHitCluster_edep_layer8;

  
  //2D
  TH2* h_toa_edepCluster;
  TH2* h_3hitEDEP_vs_clusterEDEP_1;
  TH2* h_3hitEDEP_vs_clusterEDEP_2;
  TH2* h_3hitEDEP_vs_clusterEDEP_3;

  TH2* h_edepVhits;
  TH2* h_edepVhits_layer0;
  TH2* h_edepVhits_layer1;
  TH2* h_edepVhits_layer2;
  TH2* h_edepVhits_layer3;
  TH2* h_edepVhits_layer4;
  TH2* h_edepVhits_layer5;
  TH2* h_edepVhits_layer6;
  TH2* h_edepVhits_layer7;

  //3D Histos: 
  TH3* h_3DPosition_digi;
  TH3* h_3DPosition_20digi;
  TH3* h_3DPosition_cdigi;
  TH3* h_3DPosition_20cdigi;
  TH3* h_3DPosition_r_z_hit;
  TH2* h_2D_r_hitNum;
  TH2* h_2D_z_hitNum;
  TH3* h_3DPosition_r_z_20hit;
  TH3* h_3DPosition_theta_z_20hit;
  TH3* h_3DPosition_theta_z_hit;
  TH3* h_3DPosition_theta_r_20hit;
  TH3* h_3DPosition_theta_r_hit;
  
  TH2* h_2D_r_20hitNum;
  TH2* h_2D_z_20hitNum;
  TH2* h_2D_theta_hitNum;
  TH2* h_2D_theta_20hitNum;
  
  //theta, r, and z histos for larger than 20 hit clusters: 
  TH1* h_theta_20hit;
  TH1* h_theta_20hit_layer0;
  TH1* h_theta_20hit_layer1;
  TH1* h_theta_20hit_layer2;
  TH1* h_theta_20hit_layer3;
  TH1* h_theta_20hit_layer4;
  TH1* h_theta_20hit_layer5;
  TH1* h_theta_20hit_layer6;
  TH1* h_theta_20hit_layer7;

  TH1* h_r_20hit;
  TH1* h_r_20hit_layer0;
  TH1* h_r_20hit_layer1;
  TH1* h_r_20hit_layer2;
  TH1* h_r_20hit_layer3;
  TH1* h_r_20hit_layer4;
  TH1* h_r_20hit_layer5;
  TH1* h_r_20hit_layer6;
  TH1* h_r_20hit_layer7;

  TH1* h_z_20hit;
  TH1* h_z_20hit_layer0;
  TH1* h_z_20hit_layer1;
  TH1* h_z_20hit_layer2;
  TH1* h_z_20hit_layer3;
  TH1* h_z_20hit_layer4;
  TH1* h_z_20hit_layer5;
  TH1* h_z_20hit_layer6;
  TH1* h_z_20hit_layer7;

  
};
