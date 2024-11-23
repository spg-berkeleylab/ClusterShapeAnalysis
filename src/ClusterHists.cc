#include "ClusterShapeAnalysis/ClusterHists.h"
#include "marlin/VerbosityLevels.h"
#include <EVENT/LCCollection.h>

#include <EVENT/TrackerHit.h>
#include <EVENT/SimTrackerHit.h>
#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <UTIL/LCTrackerConf.h>

#include <TCanvas.h>
#include <TPaletteAxis.h>


ClusterHists::ClusterHists()
{
  h_size_theta_y    = new TH2F("cluster_size_vs_theta_y" , ";Cluster #theta; Cluster size" , 100, 0,  3.14,  11,  -0.5,  10.5  );
  h_size_theta_x    = new TH2F("cluster_size_vs_theta_x" , ";Cluster #theta; Cluster size" , 100, 0,  3.14,  11,  -0.5,  10.5  );
  h_size_theta_tot  = new TH2F("cluster_size_vs_theta_tot" , ";Cluster #theta; Cluster size" , 100, 0,  3.14,  11,  -0.5,  10.5  );
  h_size_r_tot      = new TH2F("cluster_size_vs_R_tot" , ";Cluster R (x^2+y^2)^(1/2) (mm); Cluster size" , 100, 20,  120,  11,  -0.5,  10.5  );
  h_size_r_tot_0     = new TH2F("cluster_size_vs_R_tot_0" , ";Cluster R (x^2+y^2)^(1/2) (mm); Cluster size" , 100, 20,  120,  11,  -0.5,  10.5  );
  h_size_r_tot_1     = new TH2F("cluster_size_vs_R_tot_1" , ";Cluster R (x^2+y^2)^(1/2) (mm); Cluster size" , 100, 20,  120,  11,  -0.5,  10.5  );
  h_size_r_tot_2     = new TH2F("cluster_size_vs_R_tot_2" , ";Cluster R (x^2+y^2)^(1/2) (mm); Cluster size" , 100, 20,  120,  11,  -0.5,  10.5  );
  h_size_r_tot_3     = new TH2F("cluster_size_vs_R_tot_3" , ";Cluster R (x^2+y^2)^(1/2) (mm); Cluster size" , 100, 20,  120,  11,  -0.5,  10.5  );
  h_cluster_pos   = new TH2F("cluster_position"      , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_0 = new TH2F("cluster_position_0"    , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_1 = new TH2F("cluster_position_1"    , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_2 = new TH2F("cluster_position_2"    , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_3 = new TH2F("cluster_position_3"    , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_clusters_by_layer   = new TH1F("numClusters_by_layer"      , ";Layer Index; Number of Clusters",8,0,8);
  h_hits_by_layer   = new TH1F("numhits_by_layer"      , ";Layer Index; Number of Hits",8,0,8);
  h_theta         = new TH1F("theta"                 , ";Theta;Number of Clusters"       ,100,0,3.15);
  h_cluster_edep     = new TH1F("Clusters_edep"          , ";Energy Deposited (GeV);Clusters" ,100,0,0.0005);
  h_hit_edep     = new TH1F("Hits_edep"          , ";Deposited charge (electrons);Hits" ,5000, 0, 50000);//5000, 0, 50000 Change JULIET
  h_edep_r   = new TH2F("edep_vs_r" , ";Cluster R (x^2+y^2)^(1/2) (mm); Energy Deposited (GeV)" , 100, 20,  120,  100,  0,  0.002 );
  h_edep_cluster   = new TH2F("edep_vs_cluster_size" , "; Energy Deposited (GeV); Total Cluster Size" ,100,  0,  0.002, 100, -0.5, 99.5 );
  h_toa_edepCluster = new TH2F("toa_vs_edepCluster", "; Time of Arrival [ns];Energy Deposited [GeV]", 100, 0, 10, 100, 0, 0.0005); //--JULIET
  h_3hitEDEP_vs_clusterEDEP_1 = new TH2F("3hitEDEP_vs_clusterEDEP1", "; Hit EDEP [GeV];Cluster EDEP [GeV]", 100, 0, 0.00005, 100, 0, 0.00005); //--JULIET
  h_3hitEDEP_vs_clusterEDEP_2 = new TH2F("3hitEDEP_vs_clusterEDEP2", "; Hit EDEP [GeV];Cluster EDEP [GeV]", 100, 0, 0.0001, 100, 0.000045, 0.0001);
  h_3hitEDEP_vs_clusterEDEP_3 = new TH2F("3hitEDEP_vs_clusterEDEP3", "; Hit EDEP [GeV];Cluster EDEP [GeV]", 100, 0, 0.0001, 100, 0.00009, 0.0003);

  //2D plot of each layer cluster EDEP vs # of Hits:
  h_edepVhits   = new TH2F("edepVhits" , "; Energy Deposited (GeV); Number of Hits" ,100,  0,  0.0005, 100, 0, 10);
  h_edepVhits_layer0   = new TH2F("edepVhits_layer0" , "; Energy Deposited (GeV); Number of Hits" ,100,  0,  0.0005, 100, 0, 10);
  h_edepVhits_layer1   = new TH2F("edepVhits_layer1" , "; Energy Deposited (GeV); Number of Hits" ,100,  0,  0.0005, 100, 0, 10);
  h_edepVhits_layer2   = new TH2F("edepVhits_layer2" , "; Energy Deposited (GeV); Number of Hits" ,100,  0,  0.0005, 100, 0, 10);
  h_edepVhits_layer3   = new TH2F("edepVhits_layer3" , "; Energy Deposited (GeV); Number of Hits" ,100,  0,  0.0005, 100, 0, 10);
  h_edepVhits_layer4   = new TH2F("edepVhits_layer4" , "; Energy Deposited (GeV); Number of Hits" ,100,  0,  0.0005, 100, 0, 10);
  h_edepVhits_layer5   = new TH2F("edepVhits_layer5" , "; Energy Deposited (GeV); Number of Hits" ,100,  0,  0.0005, 100, 0, 10);
  h_edepVhits_layer6   = new TH2F("edepVhits_layer6" , "; Energy Deposited (GeV); Number of Hits" ,100,  0,  0.0005, 100, 0, 10);
  h_edepVhits_layer7   = new TH2F("edepVhits_layer7" , "; Energy Deposited (GeV); Number of Hits" ,100,  0,  0.0005, 100, 0, 10);
  //h_edepVhits_layer8   = new TH2F("edepVhits_layer8" , "; Energy Deposited (GeV); Number of Hits" ,100,  0,  0.1e-3, 100, 0, 5);


  // Create position histograms for tracker hits
  int numbins_all = 1000;
  int rmin_all = 0;
  int rmax_all = 1600;
  int zmin_all = -2500;
  int zmax_all = 2500;
  h_x   = new TH1F("x  " , ";x   ; Num Hits" , numbins_all, -rmax_all, rmax_all);
  h_y   = new TH1F("y  " , ";y   ; Num Hits" , numbins_all, -rmax_all, rmax_all);
  h_z   = new TH1F("z  " , ";z   ; Num Hits" , numbins_all,  zmin_all, zmax_all);
  h_r   = new TH1F("r  " , ";r   ; Num Hits" , numbins_all,  rmin_all, rmax_all);
  h_z_r = new TH2F("z_r" , ";z_r ; r"        , numbins_all/5,  zmin_all, zmax_all, numbins_all/10, rmin_all, rmax_all);
  h_x_y = new TH2F("x_y" , ";x_y ; r"        , numbins_all, -rmax_all, rmax_all, numbins_all, -rmax_all, rmax_all);
  h_z_r_hits = new TH2F("z_r_hits" , ";z_r ; r"        , numbins_all/5,  zmin_all, zmax_all, numbins_all/10, rmin_all, rmax_all);
  h_x_y_hits = new TH2F("x_y_hits" , ";x_y ; r"        , numbins_all, -rmax_all, rmax_all, numbins_all, -rmax_all, rmax_all);
  h_z_layer0   = new TH1F("hits_vs_z_layer0  " , ";z   ; Num Hits" , numbins_all,  zmin_all, zmax_all);
  h_z_layer1   = new TH1F("hits_vs_z_layer1  " , ";z   ; Num Hits" , numbins_all,  zmin_all, zmax_all);
  h_z_layer2   = new TH1F("hits_vs_z_layer2  " , ";z   ; Num Hits" , numbins_all,  zmin_all, zmax_all);
  h_z_layer3   = new TH1F("hits_vs_z_layer3  " , ";z   ; Num Hits" , numbins_all,  zmin_all, zmax_all);
  h_z_layer4   = new TH1F("hits_vs_z_layer4  " , ";z   ; Num Hits" , numbins_all,  zmin_all, zmax_all);
  h_z_layer5   = new TH1F("hits_vs_z_layer5  " , ";z   ; Num Hits" , numbins_all,  zmin_all, zmax_all);
  h_z_layer6   = new TH1F("hits_vs_z_layer6  " , ";z   ; Num Hits" , numbins_all,  zmin_all, zmax_all);
  h_z_layer7   = new TH1F("hits_vs_z_layer7  " , ";z   ; Num Hits" , numbins_all,  zmin_all, zmax_all);
  h_z_layer8   = new TH1F("hits_vs_z_layer8  " , ";z   ; Num Hits" , numbins_all,  zmin_all, zmax_all);
  h_r_layer0   = new TH1F("hits_vs_r_layer0  " , ";r   ; Num Hits" , numbins_all,  rmin_all, rmax_all);
  h_r_layer1   = new TH1F("hits_vs_r_layer1  " , ";r   ; Num Hits" , numbins_all,  rmin_all, rmax_all);
  h_r_layer2   = new TH1F("hits_vs_r_layer2  " , ";r   ; Num Hits" , numbins_all,  rmin_all, rmax_all);
  h_r_layer3   = new TH1F("hits_vs_r_layer3  " , ";r   ; Num Hits" , numbins_all,  rmin_all, rmax_all);
  h_r_layer4   = new TH1F("hits_vs_r_layer4  " , ";r   ; Num Hits" , numbins_all,  rmin_all, rmax_all);
  h_r_layer5   = new TH1F("hits_vs_r_layer5  " , ";r   ; Num Hits" , numbins_all,  rmin_all, rmax_all);
  h_r_layer6   = new TH1F("hits_vs_r_layer6  " , ";r   ; Num Hits" , numbins_all,  rmin_all, rmax_all);
  h_r_layer7   = new TH1F("hits_vs_r_layer7  " , ";r   ; Num Hits" , numbins_all,  rmin_all, rmax_all);
  h_r_layer8   = new TH1F("hits_vs_r_layer8  " , ";r   ; Num Hits" , numbins_all,  rmin_all, rmax_all);


  // histograms with ranges that will just show vertex tracker hits
  int numbins_vx = 500;
  int rmin_vx = 0;
  int rmax_vx = 120;
  int zmin_vx = -300;
  int zmax_vx = 300;
  h_x_vx   = new TH1F("x_vx  " , ";x   ; Num Hits" , numbins_vx, -rmax_vx, rmax_vx);
  h_y_vx   = new TH1F("y_vx  " , ";y   ; Num Hits" , numbins_vx, -rmax_vx, rmax_vx);
  h_z_vx   = new TH1F("z_vx  " , ";z   ; Num Hits" , numbins_vx,  zmin_vx, zmax_vx);
  h_r_vx   = new TH1F("r_vx  " , ";r   ; Num Hits" , numbins_vx,  rmin_vx, rmax_vx);
  h_z_r_vx = new TH2F("z_r_vx" , ";z_r ; r"        , numbins_vx/10,  zmin_vx, zmax_vx, numbins_vx/10, rmin_vx, rmax_vx);
  h_x_y_vx = new TH2F("x_y_vx" , ";x_y ; r"        , numbins_vx, -rmax_vx, rmax_vx, numbins_vx, -rmax_vx, rmax_vx);

  //create histograms for 1st layer in full silicon tracking detector: EDEP & Time
  h_cluster_edep_layer0   = new TH1F("Clusters_edep_layer0", ";Energy Deposited (GeV);Clusters" ,40,0,0.0005);
  h_cluster_edep_layer1   = new TH1F("Clusters_edep_layer1", ";Energy Deposited (GeV);Clusters" ,40,0,0.0005);
  h_cluster_edep_layer2   = new TH1F("Clusters_edep_layer2", ";Energy Deposited (GeV);Clusters" ,40,0,0.0005);
  h_cluster_edep_layer3   = new TH1F("Clusters_edep_layer3", ";Energy Deposited (GeV);Clusters" ,40,0,0.0005);
  h_cluster_edep_layer4   = new TH1F("Clusters_edep_layer4", ";Energy Deposited (GeV);Clusters" ,40,0,0.0005);
  h_cluster_edep_layer5   = new TH1F("Clusters_edep_layer5", ";Energy Deposited (GeV);Clusters" ,40,0,0.0005);
  h_cluster_edep_layer6   = new TH1F("Clusters_edep_layer6", ";Energy Deposited (GeV);Clusters" ,40,0,0.0005);
  h_cluster_edep_layer7   = new TH1F("Clusters_edep_layer7", ";Energy Deposited (GeV);Clusters" ,40,0,0.0005);
  h_cluster_edep_layer8   = new TH1F("Clusters_edep_layer8", ";Energy Deposited (GeV);Clusters" ,40,0,0.0005);

  //normalized cluster energies:
  h_cluster_edep_norm_layer0   = new TH1F("Clusters_edep_norm_layer0", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  h_cluster_edep_norm_layer1   = new TH1F("Clusters_edep_norm_layer1", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  h_cluster_edep_norm_layer2   = new TH1F("Clusters_edep_norm_layer2", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  h_cluster_edep_norm_layer3   = new TH1F("Clusters_edep_norm_layer3", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  h_cluster_edep_norm_layer4   = new TH1F("Clusters_edep_norm_layer4", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  h_cluster_edep_norm_layer5   = new TH1F("Clusters_edep_norm_layer5", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  h_cluster_edep_norm_layer6   = new TH1F("Clusters_edep_norm_layer6", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  h_cluster_edep_norm_layer7   = new TH1F("Clusters_edep_norm_layer7", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  h_cluster_edep_norm_layer8   = new TH1F("Clusters_edep_norm_layer8", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);

  h_truth_cluster_edep_layer0   = new TH1F("Clusters_truth_edep_norm_layer0", ";Truth Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  h_truth_cluster_edep_layer1   = new TH1F("Clusters_truth_edep_norm_layer1", ";Truth Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  h_truth_cluster_edep_layer2   = new TH1F("Clusters_truth_edep_norm_layer2", ";Truth Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  h_truth_cluster_edep_layer3   = new TH1F("Clusters_truth_edep_norm_layer3", ";Truth Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  h_truth_cluster_edep_layer4   = new TH1F("Clusters_truth_edep_norm_layer4", ";Truth Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  h_truth_cluster_edep_layer5   = new TH1F("Clusters_truth_edep_norm_layer5", ";Truth Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  h_truth_cluster_edep_layer6   = new TH1F("Clusters_truth_edep_norm_layer6", ";Truth Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  h_truth_cluster_edep_layer7   = new TH1F("Clusters_truth_edep_norm_layer7", ";Truth Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  h_truth_cluster_edep_layer8   = new TH1F("Clusters_truth_edep_norm_layer8", ";Truth Energy Deposited (GeV);Clusters" ,200,0,0.0005);

  h_hit_edep_layer0   = new TH1F("hit_edep_layer0", ";Deposited charge (electrons);Hits" ,5000, 0, 50000);
  h_hit_edep_layer1   = new TH1F("hit_edep_layer1", ";Deposited charge (electrons);Hits" ,5000, 0, 50000);
  h_hit_edep_layer2   = new TH1F("hit_edep_layer2", ";Deposited charge (electrons);Hits" ,5000, 0, 50000);
  h_hit_edep_layer3   = new TH1F("hit_edep_layer3", ";Deposited charge (electrons);Hits" ,5000, 0, 50000);
  h_hit_edep_layer4   = new TH1F("hit_edep_layer4", ";Deposited charge (electrons);Hits" ,5000, 0, 50000);
  h_hit_edep_layer5   = new TH1F("hit_edep_layer5", ";Deposited charge (electrons);Hits" ,5000, 0, 50000);
  h_hit_edep_layer6   = new TH1F("hit_edep_layer6", ";Deposited charge (electrons);Hits" ,5000, 0, 50000); //(100, 0, 36000)
  h_hit_edep_layer7   = new TH1F("hit_edep_layer7", ";Deposited charge (electrons);Hits" ,5000, 0, 50000);//5000,0,50000
  h_hit_edep_layer8   = new TH1F("hit_edep_layer8", ";Deposited charge (electrons);Hits" ,5000, 0, 50000);//5000,0,50000

  h_trackerhit_time_layer0  = new TH1F("trackerhit_time_layer0", ";Time (ns);Events" ,200,-10,100);
  h_trackerhit_time_layer1  = new TH1F("trackerhit_time_layer1", ";Time (ns);Events" ,200,-10,100);
  h_trackerhit_time_layer2  = new TH1F("trackerhit_time_layer2", ";Time (ns);Events" ,200,-10,100);
  h_trackerhit_time_layer3  = new TH1F("trackerhit_time_layer3", ";Time (ns);Events" ,200,-10,100);
  h_trackerhit_time_layer4  = new TH1F("trackerhit_time_layer4", ";Time (ns);Events" ,200,-10,100);
  h_trackerhit_time_layer5  = new TH1F("trackerhit_time_layer5", ";Time (ns);Events" ,200,-10,100);
  h_trackerhit_time_layer6  = new TH1F("trackerhit_time_layer6", ";Time (ns);Events" ,200,-10,100);
  h_trackerhit_time_layer7  = new TH1F("trackerhit_time_layer7", ";Time (ns);Events" ,200,-10,100);
  h_trackerhit_time_layer8  = new TH1F("trackerhit_time_layer8", ";Time (ns);Events" ,200,-10,100);

  //number of hits per cluster
  h_thclen = new TH1F("thclen", ";Number of Hit Constituents; Events", 50, 0, 2000); //Total number of hit constituents 
  h_cluster_1hits = new TH1F("cluster_1hits", ";Cluster EDP (GeV);Clusters",200, 0, 0.0005);
  h_cluster_2hits = new TH1F("cluster_2hits", ";Cluster EDP (GeV);Clusters",200, 0, 0.0005);
  h_cluster_3hits = new TH1F("cluster_3hits", ";Cluster EDP (GeV);Clusters",200, 0, 0.0005);
  h_cluster_4hits = new TH1F("cluster_4hits", ";Cluster EDP (GeV);Clusters",200, 0, 0.0005);
  h_cluster_5hits = new TH1F("cluster_5hits", ";Cluster EDP (GeV);Clusters",200, 0, 0.0005);
  h_cluster_6hits = new TH1F("cluster_6hits", ";Cluster EDP (GeV);Clusters",200, 0, 0.0005);
  h_cluster_7hits = new TH1F("cluster_7hits", ";Cluster EDP (GeV);Clusters",200, 0, 0.0005);
  h_cluster_8hits = new TH1F("cluster_8hits", ";Cluster EDP (GeV);Clusters",200, 0, 0.0005);
  h_cluster_9hits = new TH1F("cluster_9hits", ";Cluster EDP (GeV);Clusters",200, 0, 0.0005);
  //cluster edep needs to be divided by total number of hits in the cluster and plotted
  h_cluster_norm = new TH1F("cluster_norm", ";Cluster EDP (GeV);Clusters",100, 0, 0.1e-3);

  //number of hits per cluster per layer: 
  h_thclen_layer0 = new TH1F("thclen_layer0", ";Number of Hit Constituents; Events", 50, 0, 1000);
  h_thclen_layer1 = new TH1F("thclen_layer1", ";Number of Hit Constituents; Events", 50, 0, 1000);
  h_thclen_layer2 = new TH1F("thclen_layer2", ";Number of Hit Constituents; Events", 50, 0, 1000);
  h_thclen_layer3 = new TH1F("thclen_layer3", ";Number of Hit Constituents; Events", 50, 0, 1000);
  h_thclen_layer4 = new TH1F("thclen_layer4", ";Number of Hit Constituents; Events", 50, 0, 1000);
  h_thclen_layer5 = new TH1F("thclen_layer5", ";Number of Hit Constituents; Events", 50, 0, 1000);
  h_thclen_layer6 = new TH1F("thclen_layer6", ";Number of Hit Constituents; Events", 50, 0, 1000);
  h_thclen_layer7 = new TH1F("thclen_layer7", ";Number of Hit Constituents; Events", 50, 0, 1000);
  h_thclen_layer8 = new TH1F("thclen_layer8", ";Number of Hit Constituents; Events", 50, 0, 1000);

//Differnece in hit and cluster edep per layer: 
h_diffHitCluster_edep_layer0   = new TH1F("diffHitCluster_edep_layer0", ";Difference in Hit and Cluster Energy Deposited (GeV);Difference" ,50,-1e-6,1e-6);
h_diffHitCluster_edep_layer1   = new TH1F("diffHitCluster_edep_layer1", ";Difference in Hit and Cluster Energy Deposited (GeV);Difference" ,50,-1e-6,1e-6);
h_diffHitCluster_edep_layer2   = new TH1F("diffHitCluster_edep_layer2", ";Difference in Hit and Cluster Energy Deposited (GeV);Difference" ,50,-1e-6,1e-6);
h_diffHitCluster_edep_layer3   = new TH1F("diffHitCluster_edep_layer3", ";Difference in Hit and Cluster Energy Deposited (GeV);Difference" ,50,-1e-6,1e-6);
h_diffHitCluster_edep_layer4   = new TH1F("diffHitCluster_edep_layer4", ";Difference in Hit and Cluster Energy Deposited (GeV);Difference" ,50,-1e-6,1e-6);
h_diffHitCluster_edep_layer5   = new TH1F("diffHitCluster_edep_layer5", ";Difference in Hit and Cluster Energy Deposited (GeV);Difference" ,50,-1e-6,1e-6);
h_diffHitCluster_edep_layer6   = new TH1F("diffHitCluster_edep_layer6", ";Difference in Hit and Cluster Energy Deposited (GeV);Difference" ,50,-1e-6,1e-6);
h_diffHitCluster_edep_layer7   = new TH1F("diffHitCluster_edep_layer7", ";Difference in Hit and Cluster Energy Deposited (GeV);Difference" ,50,-1e-6,1e-6);
h_diffHitCluster_edep_layer8   = new TH1F("diffHitCluster_edep_layer8", ";Difference in Hit and Cluster Energy Deposited (GeV);Difference" ,50,-1e-6,1e-6);

//3D HISTO for X vs Y vs Z position in digitized and truth 
h_3DPosition_digi = new TH3F("3DPosition_digi", "3D Digitized Position;x[mm];y[mm];z[mm]", 100, 0, 1600, 100, 0, 1600, 100,  -2500, 2500);
h_3DPosition_20digi = new TH3D ("3DPosition_20digi", "3D Digitized Position;x[mm];y[mm];z[mm]", numbins_all/10, -rmax_all, rmax_all, numbins_all/10, -rmax_all, rmax_all, numbins_all/10,  zmin_all, zmax_all);
h_3DPosition_cdigi = new TH3D ("3DPosition_cdigi", "3D Digitized Position;#theta;r[mm];z[mm]", numbins_all/10, 0, 3.14, numbins_all/10, -rmax_all, rmax_all, numbins_all/10,  zmin_all, zmax_all);
h_3DPosition_20cdigi = new TH3D ("3DPosition_20cdigi", "3D Digitized Position;#theta;r[mm];z[mm]", numbins_all/10, 0, 3.14, numbins_all/10, -rmax_all, rmax_all, numbins_all/10,  zmin_all, zmax_all);
 
h_3DPosition_r_z_hit = new TH3D ("3DPosition_r_z_hit", "3D Digitized Position;Hit Num;r[mm];z[mm]", 100, 0, 400, numbins_all/10, 0, rmax_all, numbins_all/10,  zmin_all, zmax_all);
h_3DPosition_r_z_20hit = new TH3D ("3DPosition_r_z_20hit", "3D Digitized Position;Hit Num;r[mm];z[mm]", 100, 0, 400, numbins_all/10, 0, rmax_all, numbins_all/10,  zmin_all, zmax_all);

h_3DPosition_theta_z_20hit = new TH3D ("3DPosition_theta_z_20hit", "3D Digitized Position;#theta;z[mm];Hit Num", numbins_all/10, 0, 3.14, numbins_all/10, zmin_all, zmax_all, numbins_all/10,  0, 400);
h_3DPosition_theta_z_hit = new TH3D ("3DPosition_theta_z_hit", "3D Digitized Position;#theta;z[mm];Hit Num", numbins_all/10, 0, 3.14, numbins_all/10, zmin_all, zmax_all, numbins_all/10,  0, 400);
h_3DPosition_theta_r_20hit = new TH3D ("3DPosition_theta_r_20hit", "3D Digitized Position;#theta;r[mm];Hit Num", numbins_all/10, 0, 3.14, numbins_all/10, 0, rmax_all, numbins_all/10,  0, 400);
h_3DPosition_theta_r_hit = new TH3D ("3DPosition_theta_r_hit", "3D Digitized Position;#theta;r[mm];Hit Num", numbins_all/10, 0, 3.14, numbins_all/10, 0, rmax_all, numbins_all/10,  0, 400);

h_2D_r_hitNum = new TH2D ("2D_r_hitNum", "2D Digitized Radius vs Hit Density;r[mm];Hit Number", numbins_all/10, 0, rmax_all, 100, 0, 400);
h_2D_z_hitNum = new TH2D ("2D_z_hitNum", "2D Digitized Z vs Hit Density;z[mm];Hit Number", numbins_all/10, zmin_all, zmax_all, 100, 0, 400);
h_2D_r_20hitNum = new TH2D ("2D_r_20hitNum", "2D Digitized Radius vs Hit Density;r[mm];Hit Number", numbins_all/10, 0, rmax_all, 100, 0, 400);
h_2D_z_20hitNum = new TH2D ("2D_z_20hitNum", "2D Digitized Z vs Hit Density;z[mm];Hit Number", numbins_all/10, zmin_all, zmax_all, 100, 0,400);
h_2D_theta_hitNum = new TH2D ("2D_theta_hitNum", "2D Digitized #theta vs Hit Density;#theta;Hit Number", numbins_all/10, 0, 3.14, 100, 0, 400);
h_2D_theta_20hitNum = new TH2D ("2D_theta_20hitNum", "2D Digitized #theta vs Hit Density;#theta;Hit Number", numbins_all/10, 0, 3.14, 100, 0, 400);

//individual distributions for theta, r, and z above 20 hit cluster histrogram defs: 
h_theta_20hit = new TH1F("theta_20hit", "#theta (> 20 Hit Clusters);#theta;Events",numbins_all/10, 0, 3.14);
h_theta_20hit_layer0 = new TH1F("theta_20hit_layer0", "#theta (> 20 Hit Clusters);#theta;Events",numbins_all/10, 0, 3.14);
h_theta_20hit_layer1 = new TH1F("theta_20hit_layer1", "#theta (> 20 Hit Clusters);#theta;Events",numbins_all/10, 0, 3.14);
h_theta_20hit_layer2 = new TH1F("theta_20hit_layer2", "#theta (> 20 Hit Clusters);#theta;Events",numbins_all/10, 0, 3.14);
h_theta_20hit_layer3 = new TH1F("theta_20hit_layer3", "#theta (> 20 Hit Clusters);#theta;Events",numbins_all/10, 0, 3.14);
h_theta_20hit_layer4 = new TH1F("theta_20hit_layer4", "#theta (> 20 Hit Clusters);#theta;Events",numbins_all/10, 0, 3.14);
h_theta_20hit_layer5 = new TH1F("theta_20hit_layer5", "#theta (> 20 Hit Clusters);#theta;Events",numbins_all/10, 0, 3.14);
h_theta_20hit_layer6 = new TH1F("theta_20hit_layer6", "#theta (> 20 Hit Clusters);#theta;Events",numbins_all/10, 0, 3.14);
h_theta_20hit_layer7 = new TH1F("theta_20hit_layer7", "#theta (> 20 Hit Clusters);#theta;Events",numbins_all/10, 0, 3.14);

h_r_20hit = new TH1F("r_20hit", "R (> 20 Hit Clusters);r[mm];Events",numbins_all/10, 0, rmax_all);
h_r_20hit_layer0 = new TH1F("r_20hit_layer0", "R (> 20 Hit Clusters);r[mm];Events",numbins_all/10, 0, rmax_all);
h_r_20hit_layer1 = new TH1F("r_20hit_layer1", "R (> 20 Hit Clusters);r[mm];Events",numbins_all/10, 0, rmax_all);
h_r_20hit_layer2 = new TH1F("r_20hit_layer2", "R (> 20 Hit Clusters);r[mm];Events",numbins_all/10, 0, rmax_all);
h_r_20hit_layer3 = new TH1F("r_20hit_layer3", "R (> 20 Hit Clusters);r[mm];Events",numbins_all/10, 0, rmax_all);
h_r_20hit_layer4 = new TH1F("r_20hit_layer4", "R (> 20 Hit Clusters);r[mm];Events",numbins_all/10, 0, rmax_all);
h_r_20hit_layer5 = new TH1F("r_20hit_layer5", "R (> 20 Hit Clusters);r[mm];Events",numbins_all/10, 0, rmax_all);
h_r_20hit_layer6 = new TH1F("r_20hit_layer6", "R (> 20 Hit Clusters);r[mm];Events",numbins_all/10, 0, rmax_all);
h_r_20hit_layer7 = new TH1F("r_20hit_layer7", "R (> 20 Hit Clusters);r[mm];Events",numbins_all/10, 0, rmax_all);

h_z_20hit = new TH1F("z_20hit", "Z (> 20 Hit Clusters);z[mm];Events",numbins_all/10, zmin_all, zmax_all);
h_z_20hit_layer0 = new TH1F("z_20hit_layer0", "Z (> 20 Hit Clusters);z[mm];Events",numbins_all/10, zmin_all, zmax_all);
h_z_20hit_layer1 = new TH1F("z_20hit_layer1", "Z (> 20 Hit Clusters);z[mm];Events",numbins_all/10, zmin_all, zmax_all);
h_z_20hit_layer2 = new TH1F("z_20hit_layer2", "Z (> 20 Hit Clusters);z[mm];Events",numbins_all/10, zmin_all, zmax_all);
h_z_20hit_layer3 = new TH1F("z_20hit_layer3", "Z (> 20 Hit Clusters);z[mm];Events",numbins_all/10, zmin_all, zmax_all);
h_z_20hit_layer4 = new TH1F("z_20hit_layer4", "Z (> 20 Hit Clusters);z[mm];Events",numbins_all/10, zmin_all, zmax_all);
h_z_20hit_layer5 = new TH1F("z_20hit_layer5", "Z (> 20 Hit Clusters);z[mm];Events",numbins_all/10, zmin_all, zmax_all);
h_z_20hit_layer6 = new TH1F("z_20hit_layer6", "Z (> 20 Hit Clusters);z[mm];Events",numbins_all/10, zmin_all, zmax_all);
h_z_20hit_layer7 = new TH1F("z_20hit_layer7", "Z (> 20 Hit Clusters);z[mm];Events",numbins_all/10, zmin_all, zmax_all);

}

void ClusterHists::fill(const EVENT::TrackerHit* trkhit)
{
  //Calculate energy deposited
  float EDep = trkhit->getEDep();
  float toa = trkhit->getTime(); //time of arrival -- Juliet
  
  //Calculating theta
  float x = trkhit->getPosition()[0];
  float y = trkhit->getPosition()[1];
  float z = trkhit->getPosition()[2];
  float r = sqrt(pow(x,2)+pow(y,2));
  float incidentTheta = std::atan(r/z);
  
  if(incidentTheta<0)
    incidentTheta += M_PI;
  streamlog_out(DEBUG6) << "Theta: " << incidentTheta << std::endl;


  //Calculating cluster size
  const lcio::LCObjectVec &rawHits = trkhit->getRawHits(); 
  float ymax = -1000000;
  float xmax = -1000000;
  float ymin =  1000000; 
  float xmin =  1000000; 

  float loopsize = rawHits.size();
  streamlog_out(DEBUG3) << "Number of raw hits: " << loopsize << std::endl;

  for (size_t j=0; j<loopsize; ++j) {
    lcio::SimTrackerHit *hitConstituent = dynamic_cast<lcio::SimTrackerHit*>( rawHits[j] );
    //lcio::TrackerHit *hitTracker = dynamic_cast<lcio::TrackerHit*>( rawHits[j] );
    h_hit_edep->Fill(hitConstituent->getEDep());
    const double *localPos = hitConstituent->getPosition();
    float x_local = localPos[0];
    float y_local = localPos[1];
    streamlog_out(DEBUG6) << "Local y: " << y_local << ", local x: " << x_local << std::endl;

    if (y_local < ymin){
      ymin = y_local;
      }
    if (y_local > ymax){
      ymax = y_local;          
      } 

    if (x_local < xmin){
      xmin = x_local;
      }
    if (x_local > xmax){
      xmax = x_local;
      //streamlog_out(DEBUG2) << "Updated ymin: " << ymin << ", xmin: " << xmin << std::endl;
      //streamlog_out(DEBUG2) << "Updated ymax: " << ymax << ", xmax: " << xmax << std::endl;          
     }
  }
  streamlog_out(DEBUG2) << "Min y pos: " << ymin  << ", max y pos: " << ymax << std::endl;
  streamlog_out(DEBUG2) << "Min x pos: " << xmin  << ", max x pos: " << xmax << std::endl;
  float cluster_size_y = (ymax - ymin)+1;
  float cluster_size_x = (xmax - xmin)+1;
  float cluster_size_tot = loopsize;
  streamlog_out(DEBUG6) << "Cluster size, y direction (parallel to beam line for barrel, radial dir for endcap) " << cluster_size_y << std::endl;
  streamlog_out(DEBUG6) << "Cluster size, x direction (parallel to beam line in ladder plane for barrel, phi dir for endcap) " << cluster_size_x << std::endl;
  streamlog_out(DEBUG6) << "Cluster size, TOTAL " << cluster_size_tot << std::endl;
  /* if (cluster_size_y < 1) {
    streamlog_out(WARNING) << "Cluster calculated y size less than 1. Skip cluster." << std::endl;
    std::stringstream err  ; err << " Cluster calculation failed.";
    throw EVENT::Exception ( err.str() );
  } */

  //Get hit subdetector/layer 
  std::string _encoderString = lcio::LCTrackerCellID::encoding_string();
  UTIL::CellIDDecoder<lcio::TrackerHit> decoder(_encoderString);
  uint32_t systemID = decoder(trkhit)["system"];
  uint32_t layerID = decoder(trkhit)["layer"];
  
  streamlog_out(DEBUG9) << "Layer ID is: " << layerID << std::endl;

  // Fill for all hits
  h_size_theta_y->Fill(incidentTheta, cluster_size_y);
  h_size_theta_x->Fill(incidentTheta, cluster_size_x);
  h_size_theta_tot->Fill(incidentTheta, cluster_size_tot);
  h_size_r_tot->Fill(r, cluster_size_tot);

  h_theta->Fill(incidentTheta);
  h_cluster_pos->Fill(z,r);
  h_clusters_by_layer->Fill(layerID);

  // tracker hit hists
  h_x->Fill(x);  
  h_y->Fill(y);  
  h_z->Fill(z);  
  h_r->Fill(r);
  h_z_r->Fill(z,r);
  h_x_y->Fill(x,y);
  for (size_t j=0; j<loopsize; ++j) {
    lcio::SimTrackerHit *hitConstituent = dynamic_cast<lcio::SimTrackerHit*>( rawHits[j] ); // adding this -- Juliet
    //lcio::TrackerHit *hitTracker = dynamic_cast<lcio::TrackerHit*>( rawHits[j] );
    h_hits_by_layer->Fill(layerID);
    h_z_r_hits->Fill(z,r);
    h_x_y_hits->Fill(x,y);
    h_edepVhits->Fill(EDep, rawHits.size());

    if(rawHits.size() > 20){
      h_3DPosition_20digi->Fill(x, y, z);
      h_3DPosition_20cdigi->Fill(incidentTheta, r, z);
      h_2D_r_20hitNum->Fill(r, rawHits.size());
      h_2D_z_20hitNum->Fill(z, rawHits.size());
      h_2D_theta_20hitNum->Fill(incidentTheta, rawHits.size());
      h_3DPosition_r_z_20hit->Fill(rawHits.size(), r, z);
      h_3DPosition_theta_z_20hit->Fill(incidentTheta, z, rawHits.size());
      h_3DPosition_theta_r_20hit->Fill(incidentTheta, r, rawHits.size());
      //theta, r, z
      h_theta_20hit->Fill(incidentTheta);
      h_r_20hit->Fill(r);
      h_z_20hit->Fill(z);

    }
    h_3DPosition_digi->Fill(x, y, z);
    h_3DPosition_cdigi->Fill(incidentTheta, r, z);
    h_2D_r_hitNum->Fill(r, rawHits.size());
    h_2D_z_hitNum->Fill(z, rawHits.size());
    h_2D_theta_hitNum->Fill(incidentTheta, rawHits.size());
    h_3DPosition_r_z_hit->Fill(rawHits.size(), r, z);
    h_3DPosition_theta_z_hit->Fill(incidentTheta, z, rawHits.size());
      h_3DPosition_theta_r_hit->Fill(incidentTheta, r, rawHits.size());

    if(rawHits.size() == 3){
      if(EDep < 0.05e-3){
        h_3hitEDEP_vs_clusterEDEP_1->Fill(hitConstituent->getEDep()*3.6e-9,EDep);
      }
      else if(EDep >= 0.05e-3 && EDep <= 0.1e-3){
        h_3hitEDEP_vs_clusterEDEP_2->Fill(hitConstituent->getEDep()*3.6e-9,EDep);
      }
      else if(EDep > 0.1e-3){
        h_3hitEDEP_vs_clusterEDEP_3->Fill(hitConstituent->getEDep()*3.6e-9,EDep);
      }
    }

    if(layerID==0){
      //remove events with hit size greater than 200: 
      // if(rawHits.size() > 5 || rawHits.size()< 2) continue;
      // if(EDep > 0.075e-3) continue;

      if(rawHits.size() > 20){
        //theta, r, z
        h_theta_20hit_layer0->Fill(incidentTheta);
        h_r_20hit_layer0->Fill(r);
        h_z_20hit_layer0->Fill(z);
      }

      h_edepVhits_layer0->Fill(EDep, rawHits.size());
      
      h_z_layer0->Fill(z);
      h_r_layer0->Fill(r);
      h_cluster_edep_layer0->Fill(EDep); //energy cluster hits in GeV in layer 1 -- Juliet
      h_cluster_edep_norm_layer0->Fill(EDep/rawHits.size());
      h_trackerhit_time_layer0->Fill(toa); //time of arrive in layer 1 -- Juliet
      h_hit_edep_layer0->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
      h_thclen_layer0->Fill(rawHits.size());
      //diff
      h_diffHitCluster_edep_layer0->Fill(abs(EDep/rawHits.size() - hitConstituent->getEDep()*3.7e-6));
    }
    if(layerID==1){
      // if(rawHits.size() > 5 || rawHits.size() < 2) continue;
      // if(EDep > 0.075e-3) continue;

      if(rawHits.size() > 20){
        //theta, r, z
        h_theta_20hit_layer1->Fill(incidentTheta);
        h_r_20hit_layer1->Fill(r);
        h_z_20hit_layer1->Fill(z);
      }

      h_edepVhits_layer1->Fill(EDep, rawHits.size());
      //h_truth_cluster_edep_layer1->Fill(hitTracker->getEDep());
      h_z_layer1->Fill(z);
      h_r_layer1->Fill(r);
      h_cluster_edep_layer1->Fill(EDep); //energy cluster hits in GeV in layer 2 -- Juliet
      h_cluster_edep_norm_layer1->Fill(EDep/rawHits.size());
      h_trackerhit_time_layer1->Fill(toa); //time of arrive in layer 2 -- Juliet
      h_hit_edep_layer1->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
      h_thclen_layer1->Fill(rawHits.size());
      h_diffHitCluster_edep_layer1->Fill(abs(EDep/rawHits.size() - hitConstituent->getEDep()*3.7e-6));
    }
    if(layerID==2){
      // if(rawHits.size() > 5 || rawHits.size()< 2) continue;
      // if(EDep > 0.075e-3) continue;

      if(rawHits.size() > 20){
        //theta, r, z
        h_theta_20hit_layer2->Fill(incidentTheta);
        h_r_20hit_layer2->Fill(r);
        h_z_20hit_layer2->Fill(z);
      }

      h_edepVhits_layer2->Fill(EDep, rawHits.size());
      //h_truth_cluster_edep_layer2->Fill(hitTracker->getEDep());
      h_z_layer2->Fill(z);
      h_r_layer2->Fill(r);
      h_cluster_edep_layer2->Fill(EDep); //energy cluster hits in GeV in layer 3 -- Juliet
      h_cluster_edep_norm_layer2->Fill(EDep/rawHits.size());
      h_trackerhit_time_layer2->Fill(toa); //time of arrive in layer 3 -- Juliet
      h_hit_edep_layer2->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
      h_thclen_layer2->Fill(rawHits.size());
      h_diffHitCluster_edep_layer2->Fill(abs(EDep/rawHits.size() - hitConstituent->getEDep()*3.7e-6));
    }
    if(layerID==3){
      // if(rawHits.size() > 5 || rawHits.size()< 2) continue;
      // if(EDep > 0.075e-3) continue;

      if(rawHits.size() > 20){
        //theta, r, z
        h_theta_20hit_layer3->Fill(incidentTheta);
        h_r_20hit_layer3->Fill(r);
        h_z_20hit_layer3->Fill(z);
      }

      h_edepVhits_layer3->Fill(EDep, rawHits.size());
      //h_truth_cluster_edep_layer3->Fill(hitTracker->getEDep());
      h_z_layer3->Fill(z);
      h_r_layer3->Fill(r);
      h_cluster_edep_layer3->Fill(EDep); //energy cluster hits in GeV in layer 4 -- Juliet
      h_cluster_edep_norm_layer3->Fill(EDep/rawHits.size());
      h_trackerhit_time_layer3->Fill(toa); //time of arrive in layer 4 -- Juliet
      h_hit_edep_layer3->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
      h_thclen_layer3->Fill(rawHits.size());
      h_diffHitCluster_edep_layer3->Fill(abs(EDep/rawHits.size() - hitConstituent->getEDep()*3.7e-6));
    }
    if(layerID==4){
      // if(rawHits.size() > 5 || rawHits.size()< 2) continue;
      // if(EDep > 0.075e-3) continue;

      if(rawHits.size() > 20){
        //theta, r, z
        h_theta_20hit_layer4->Fill(incidentTheta);
        h_r_20hit_layer4->Fill(r);
        h_z_20hit_layer4->Fill(z);
      }

      h_edepVhits_layer4->Fill(EDep, rawHits.size());
      //h_truth_cluster_edep_layer4->Fill(hitTracker->getEDep());
      h_z_layer4->Fill(z);
      h_r_layer4->Fill(r);
      h_cluster_edep_layer4->Fill(EDep); //energy cluster hits in GeV in layer 5 -- Juliet
      h_cluster_edep_norm_layer4->Fill(EDep/rawHits.size());
      h_trackerhit_time_layer4->Fill(toa); //time of arrive in layer 5 -- Juliet
      h_hit_edep_layer4->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
      h_thclen_layer4->Fill(rawHits.size());
      h_diffHitCluster_edep_layer4->Fill(abs(EDep/rawHits.size() - hitConstituent->getEDep()*3.7e-6));
    }
    if(layerID==5){
      // if(rawHits.size() > 5 || rawHits.size()< 2) continue;
      // if(EDep > 0.075e-3) continue;

      if(rawHits.size() > 20){
        //theta, r, z
        h_theta_20hit_layer5->Fill(incidentTheta);
        h_r_20hit_layer5->Fill(r);
        h_z_20hit_layer5->Fill(z);
      }
      
      h_edepVhits_layer5->Fill(EDep, rawHits.size());
      //h_truth_cluster_edep_layer5->Fill(hitTracker->getEDep());
      h_z_layer5->Fill(z);
      h_r_layer5->Fill(r);
      h_cluster_edep_layer5->Fill(EDep); //energy cluster hits in GeV in layer 6 -- Juliet
      h_cluster_edep_norm_layer5->Fill(EDep/rawHits.size());
      h_trackerhit_time_layer5->Fill(toa); //time of arrive in layer 6 -- Juliet
      h_hit_edep_layer5->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
      h_thclen_layer5->Fill(rawHits.size());
      h_diffHitCluster_edep_layer5->Fill(abs(EDep/rawHits.size() - hitConstituent->getEDep()*3.7e-6));
    }
    if(layerID==6){
      // if(rawHits.size() > 5 || rawHits.size()< 2) continue;
      // if(EDep > 0.075e-3) continue;

      if(rawHits.size() > 20){
        //theta, r, z
        h_theta_20hit_layer6->Fill(incidentTheta);
        h_r_20hit_layer6->Fill(r);
        h_z_20hit_layer6->Fill(z);
      }

      h_edepVhits_layer6->Fill(EDep, rawHits.size());
      //h_truth_cluster_edep_layer6->Fill(hitTracker->getEDep());
      h_z_layer6->Fill(z);
      h_r_layer6->Fill(r);
      h_cluster_edep_layer6->Fill(EDep); //energy cluster hits in GeV in layer 7 -- Juliet
      h_cluster_edep_norm_layer6->Fill(EDep/rawHits.size());
      h_trackerhit_time_layer6->Fill(toa); //time of arrive in layer 7 -- Juliet
      h_hit_edep_layer6->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
      h_thclen_layer6->Fill(rawHits.size());
      h_diffHitCluster_edep_layer6->Fill(abs(EDep/rawHits.size() - hitConstituent->getEDep()*3.7e-6));
    }
    if(layerID==7){
      // if(rawHits.size() > 5 || rawHits.size()< 2) continue;
      // if(EDep > 0.075e-3) continue;

      if(rawHits.size() > 20){
        //theta, r, z
        h_theta_20hit_layer7->Fill(incidentTheta);
        h_r_20hit_layer7->Fill(r);
        h_z_20hit_layer7->Fill(z);
      }

      h_edepVhits_layer7->Fill(EDep, rawHits.size());
      //h_truth_cluster_edep_layer7->Fill(hitTracker->getEDep());
      h_z_layer7->Fill(z);
      h_r_layer7->Fill(r);
      h_cluster_edep_layer7->Fill(EDep); //energy cluster hits in GeV in layer 8 -- Juliet
      h_cluster_edep_norm_layer7->Fill(EDep/rawHits.size());
      h_trackerhit_time_layer7->Fill(toa); //time of arrive in layer 8 -- Julie
      h_hit_edep_layer7->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
      h_thclen_layer7->Fill(rawHits.size());
      h_diffHitCluster_edep_layer7->Fill(abs(EDep/rawHits.size() - hitConstituent->getEDep()*3.7e-6));
    }
    if(layerID==8){
      //if(rawHits.size() > 5 ||rawHits.size()< 2) continue;

      h_z_layer8->Fill(z);
      h_r_layer8->Fill(r);
      h_cluster_edep_layer8->Fill(EDep); //energy cluster hits in GeV in layer 8 -- Juliet
      h_cluster_edep_norm_layer8->Fill(EDep/rawHits.size());
      h_trackerhit_time_layer8->Fill(toa); //time of arrive in layer 8 -- Julie
      h_hit_edep_layer8->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
      h_thclen_layer8->Fill(rawHits.size());
      h_diffHitCluster_edep_layer8->Fill(abs(EDep/rawHits.size() - hitConstituent->getEDep()*3.7e-6));
    }
  }
  

  // vertex hists
  h_x_vx->Fill(x);  
  h_y_vx->Fill(y);  
  h_z_vx->Fill(z);  
  h_r_vx->Fill(r);  
  h_z_r_vx->Fill(z,r);
  h_x_y_vx->Fill(x,y);

  h_cluster_edep->Fill(EDep);
  h_edep_r->Fill(r,EDep);
  h_edep_cluster->Fill(EDep,cluster_size_tot);

  h_thclen->Fill(rawHits.size()); // -- JULIET
  if (rawHits.size() != 0){ //this if statement is pointless apparently - ask Simone about it next time
    h_cluster_norm->Fill(EDep/rawHits.size());
  }
  //std::cout << "EDP of " << EDep << " with raw hit num: " << rawHits.size() << " || and here is the EDEP/hitNum: "<< EDep/rawHits.size() << std::endl;
  if(rawHits.size() == 1){
    h_cluster_1hits->Fill(EDep);
  }
  if(rawHits.size() == 2){
    h_cluster_2hits->Fill(EDep);
  }
  if(rawHits.size() == 3){
    h_cluster_3hits->Fill(EDep);
  }
  if(rawHits.size() == 4){
    h_cluster_4hits->Fill(EDep);
  }
  if(rawHits.size() == 5){
    h_cluster_5hits->Fill(EDep);
  }
  if(rawHits.size() == 6){
    h_cluster_6hits->Fill(EDep);
  }
  if(rawHits.size() == 7){
    h_cluster_7hits->Fill(EDep);
  }
  if(rawHits.size() == 8){
    h_cluster_8hits->Fill(EDep);
  }
  if(rawHits.size() == 9){
    h_cluster_9hits->Fill(EDep);
  }

  //2D hist for clusterEnergy vs Time of arrival with color bar for number of hits per cluster 
  for(int i = 0; i < rawHits.size(); i++){
    h_toa_edepCluster->Fill(toa, EDep);
  }

  // Fill based on which double layer region was hit
  if(layerID==0 or layerID==1){
    h_cluster_pos_0->Fill(z,r);
    h_size_r_tot_0->Fill(r, cluster_size_tot);
    }
  if(layerID==2 or layerID==3){
    h_cluster_pos_1->Fill(z,r);
    h_size_r_tot_1->Fill(r, cluster_size_tot);
    }
  if(layerID==4 or layerID==5){
    h_cluster_pos_2->Fill(z,r);
    h_size_r_tot_2->Fill(r, cluster_size_tot);
    }
  if(layerID==6 or layerID==7){
    h_cluster_pos_3->Fill(z,r);
    h_size_r_tot_3->Fill(r, cluster_size_tot);
    }   
}


void ClusterHists::fill(const EVENT::SimTrackerHit* simtrkhit){

//   TH1* h_truth_cluster_edep_layer0;
//   TH1* h_truth_cluster_edep_layer1;
//   TH1* h_truth_cluster_edep_layer2;   
//   TH1* h_truth_cluster_edep_layer3;
//   TH1* h_truth_cluster_edep_layer4;
//   TH1* h_truth_cluster_edep_layer5;
//   TH1* h_truth_cluster_edep_layer6;
//   TH1* h_truth_cluster_edep_layer7;

//   h_truth_cluster_edep_layer0   = new TH1F("Clusters_truth_edep_norm_layer0", ";Truth Energy Deposited (GeV);Clusters" ,200,0,0.0005);
//   h_truth_cluster_edep_layer1   = new TH1F("Clusters_truth_edep_norm_layer1", ";Truth Energy Deposited (GeV);Clusters" ,200,0,0.0005);
//   h_truth_cluster_edep_layer2   = new TH1F("Clusters_truth_edep_norm_layer2", ";Truth Energy Deposited (GeV);Clusters" ,200,0,0.0005);
//   h_truth_cluster_edep_layer3   = new TH1F("Clusters_truth_edep_norm_layer3", ";Truth Energy Deposited (GeV);Clusters" ,200,0,0.0005);
//   h_truth_cluster_edep_layer4   = new TH1F("Clusters_truth_edep_norm_layer4", ";Truth Energy Deposited (GeV);Clusters" ,200,0,0.0005);
//   h_truth_cluster_edep_layer5   = new TH1F("Clusters_truth_edep_norm_layer5", ";Truth Energy Deposited (GeV);Clusters" ,200,0,0.0005);
//   h_truth_cluster_edep_layer6   = new TH1F("Clusters_truth_edep_norm_layer6", ";Truth Energy Deposited (GeV);Clusters" ,200,0,0.0005);
//   h_truth_cluster_edep_layer7   = new TH1F("Clusters_truth_edep_norm_layer7", ";Truth Energy Deposited (GeV);Clusters" ,200,0,0.0005);


  std::string _encoderString = lcio::LCTrackerCellID::encoding_string();
  UTIL::CellIDDecoder<lcio::SimTrackerHit> decoder(_encoderString);
  uint32_t layerID = decoder(simtrkhit)["layer"];
  //const lcio::LCObjectVec &rawHits = simtrkhit->getRawHits();
  //float loopsize = rawHits.size();

//   for (size_t j=0; j<loopsize; ++j) {
//     lcio::SimTrackerHit* truthHit = static_cast<lcio::SimTrackerHit*>( col->getElementAt(j) ) ; //-JULIET
  float truthEDep = simtrkhit->getEDep();//--JUlIET

  if(layerID==0){
    h_truth_cluster_edep_layer0->Fill(truthEDep);
  }
  if(layerID==1){
    h_truth_cluster_edep_layer1->Fill(truthEDep);
  }
  if(layerID==2){
    h_truth_cluster_edep_layer2->Fill(truthEDep);
  }
  if(layerID==3){
    h_truth_cluster_edep_layer3->Fill(truthEDep);
  }
  if(layerID==4){
    h_truth_cluster_edep_layer4->Fill(truthEDep);
  }
  if(layerID==5){
    h_truth_cluster_edep_layer5->Fill(truthEDep);
  }
  if(layerID==6){
    h_truth_cluster_edep_layer6->Fill(truthEDep);
  }
  if(layerID==7){
    h_truth_cluster_edep_layer7->Fill(truthEDep);
  }

//   }


}
