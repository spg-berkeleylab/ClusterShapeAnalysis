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

#include "DD4hep/Detector.h"
#include "DD4hep/DD4hepUnits.h"
#include "/global/cfs/cdirs/atlas/juliet/work/TrkHitsStudiesWorkspace/packages/MarlinTrkProcessors/source/Utils/include/FilterTimeHits.h"

#include <iostream>
#include <cmath>
#include <set>

using namespace lcio;
using namespace marlin;
using namespace std;


ClusterHists::ClusterHists()
{

  //Calculate the MPV scaling
  float mpv = 74e-6; //GeV
  float edp_rangeMax = 10 * mpv; //GeV 
  //float edp_binWidth = 20e-7; //[GeV] = 2000eV, each bin has a width of 2000 eV 
  float edp_binNum = 60; //edp_rangeMax/edp_binWidth; //bin number will be the same for electrons and GeV units
  //3.7 eV = 1 e, ou, 3.7e-9 GeV = 1e
  float edp_rangeMax_e = 2 * mpv / (3.7e-9); //range in electrons
  float edp_binNum_e = 5000; //10*edp_binNum;
  
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
  h_cluster_edep     = new TH1F("Clusters_edep"          , ";Energy Deposited (GeV);Clusters" ,edp_binNum,0,edp_rangeMax);
  h_hit_edep     = new TH1F("Hits_edep"          , ";Deposited charge (electrons);Hits" ,edp_binNum_e,0,edp_rangeMax_e);//5000, 0, 50000 Change JULIET
  h_edep_r   = new TH2F("edep_vs_r" , ";Cluster R (x^2+y^2)^(1/2) (mm); Energy Deposited (GeV)" , 100, 20,  120,  100,  0,  0.002 );
  h_edep_cluster   = new TH2F("edep_vs_cluster_size" , "; Energy Deposited (GeV); Total Cluster Size" ,100,  0,  0.002, 100, -0.5, 99.5 );
  h_toa_edepCluster = new TH2F("toa_vs_edepCluster", "; Time of Arrival [ns];Energy Deposited [GeV]", 100, 0, 10, 100, 0, 0.0005); //--JULIET
  
  // Create position histograms for tracker hits
  int numbins_all = 1000; //change this for 10 mm bins 
  int numbins_r = 160; //for 10 mm bins
  int numbins_z = 500; //for 10 mm bins
  int rmin_all = 0;
  int rmax_all = 1600;
  int zmin_all = -2500;
  int zmax_all = 2500;
  h_x   = new TH1F("x  " , ";x   ; Num Hits" , numbins_r, -rmax_all, rmax_all); //mm
  h_y   = new TH1F("y  " , ";y   ; Num Hits" , numbins_r, -rmax_all, rmax_all);
  h_z   = new TH1F("z  " , ";z   ; Num Hits" , numbins_z,  zmin_all, zmax_all);
  h_r   = new TH1F("r  " , ";r   ; Num Hits" , numbins_r,  rmin_all, rmax_all);
  h_z_r = new TH2F("z_r" , ";z_r ; r"        , numbins_z,  zmin_all, zmax_all, numbins_r, rmin_all, rmax_all);//numbins_all/5,  zmin_all, zmax_all, numbins_all/10, rmin_all, rmax_all);
  h_x_y = new TH2F("x_y" , ";x_y ; r"        , numbins_r, -rmax_all, rmax_all, numbins_r, -rmax_all, rmax_all);
  h_z_r_hits = new TH2F("z_r_hits" , ";z_r ; r"        , numbins_z,  zmin_all, zmax_all, numbins_r, rmin_all, rmax_all);//numbins_all/5,  zmin_all, zmax_all, numbins_all/10, rmin_all, rmax_all);
  h_x_y_hits = new TH2F("x_y_hits" , ";x_y ; r"        , numbins_r, -rmax_all, rmax_all, numbins_r, -rmax_all, rmax_all);
  h_z_layer0   = new TH1F("hits_vs_z_layer0  " , ";z   ; Num Hits" , numbins_z,  zmin_all, zmax_all);
  h_z_layer1   = new TH1F("hits_vs_z_layer1  " , ";z   ; Num Hits" , numbins_z,  zmin_all, zmax_all);
  h_z_layer2   = new TH1F("hits_vs_z_layer2  " , ";z   ; Num Hits" , numbins_z,  zmin_all, zmax_all);
  h_z_layer3   = new TH1F("hits_vs_z_layer3  " , ";z   ; Num Hits" , numbins_z,  zmin_all, zmax_all);
  h_z_layer4   = new TH1F("hits_vs_z_layer4  " , ";z   ; Num Hits" , numbins_z,  zmin_all, zmax_all);
  h_z_layer5   = new TH1F("hits_vs_z_layer5  " , ";z   ; Num Hits" , numbins_z,  zmin_all, zmax_all);
  h_z_layer6   = new TH1F("hits_vs_z_layer6  " , ";z   ; Num Hits" , numbins_z,  zmin_all, zmax_all);
  h_z_layer7   = new TH1F("hits_vs_z_layer7  " , ";z   ; Num Hits" , numbins_z,  zmin_all, zmax_all);
  h_z_layer8   = new TH1F("hits_vs_z_layer8  " , ";z   ; Num Hits" , numbins_z,  zmin_all, zmax_all);
  h_r_layer0   = new TH1F("hits_vs_r_layer0  " , ";r   ; Num Hits" , numbins_r,  rmin_all, rmax_all);
  h_r_layer1   = new TH1F("hits_vs_r_layer1  " , ";r   ; Num Hits" , numbins_r,  rmin_all, rmax_all);
  h_r_layer2   = new TH1F("hits_vs_r_layer2  " , ";r   ; Num Hits" , numbins_r,  rmin_all, rmax_all);
  h_r_layer3   = new TH1F("hits_vs_r_layer3  " , ";r   ; Num Hits" , numbins_r,  rmin_all, rmax_all);
  h_r_layer4   = new TH1F("hits_vs_r_layer4  " , ";r   ; Num Hits" , numbins_r,  rmin_all, rmax_all);
  h_r_layer5   = new TH1F("hits_vs_r_layer5  " , ";r   ; Num Hits" , numbins_r,  rmin_all, rmax_all);
  h_r_layer6   = new TH1F("hits_vs_r_layer6  " , ";r   ; Num Hits" , numbins_r,  rmin_all, rmax_all);
  h_r_layer7   = new TH1F("hits_vs_r_layer7  " , ";r   ; Num Hits" , numbins_r,  rmin_all, rmax_all);
  h_r_layer8   = new TH1F("hits_vs_r_layer8  " , ";r   ; Num Hits" , numbins_r,  rmin_all, rmax_all);


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
  h_cluster_edep_layer0   = new TH1F("Clusters_edep_layer0", ";Energy Deposited (GeV);Clusters" ,edp_binNum,0,edp_rangeMax);
  h_cluster_edep_layer1   = new TH1F("Clusters_edep_layer1", ";Energy Deposited (GeV);Clusters" ,edp_binNum,0,edp_rangeMax);
  h_cluster_edep_layer2   = new TH1F("Clusters_edep_layer2", ";Energy Deposited (GeV);Clusters" ,edp_binNum,0,edp_rangeMax);
  h_cluster_edep_layer3   = new TH1F("Clusters_edep_layer3", ";Energy Deposited (GeV);Clusters" ,edp_binNum,0,edp_rangeMax);
  h_cluster_edep_layer4   = new TH1F("Clusters_edep_layer4", ";Energy Deposited (GeV);Clusters" ,edp_binNum,0,edp_rangeMax);
  h_cluster_edep_layer5   = new TH1F("Clusters_edep_layer5", ";Energy Deposited (GeV);Clusters" ,edp_binNum,0,edp_rangeMax);
  h_cluster_edep_layer6   = new TH1F("Clusters_edep_layer6", ";Energy Deposited (GeV);Clusters" ,edp_binNum,0,edp_rangeMax);
  h_cluster_edep_layer7   = new TH1F("Clusters_edep_layer7", ";Energy Deposited (GeV);Clusters" ,edp_binNum,0,edp_rangeMax);
  h_cluster_edep_layer8   = new TH1F("Clusters_edep_layer8", ";Energy Deposited (GeV);Clusters" ,edp_binNum,0,edp_rangeMax);

  //time-energy cluster cut: 
  // h_cluster_edep_Tcut     = new TH1F("Clusters_edep_Tcut", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  // h_cluster_edep_Tcut_layer0   = new TH1F("Clusters_edep_Tcut_layer0", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  // h_cluster_edep_Tcut_layer1   = new TH1F("Clusters_edep_Tcut_layer1", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  // h_cluster_edep_Tcut_layer2   = new TH1F("Clusters_edep_Tcut_layer2", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  // h_cluster_edep_Tcut_layer3   = new TH1F("Clusters_edep_Tcut_layer3", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  // h_cluster_edep_Tcut_layer4   = new TH1F("Clusters_edep_Tcut_layer4", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  // h_cluster_edep_Tcut_layer5   = new TH1F("Clusters_edep_Tcut_layer5", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  // h_cluster_edep_Tcut_layer6   = new TH1F("Clusters_edep_Tcut_layer6", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);
  // h_cluster_edep_Tcut_layer7   = new TH1F("Clusters_edep_Tcut_layer7", ";Energy Deposited (GeV);Clusters" ,200,0,0.0005);

  // h_trackerhit_time_Tcut  = new TH1F("trackerhit_time_Tcut", ";Time (ns);Events" ,800,-1,20);
  // h_trackerhit_time_Tcut_layer0  = new TH1F("trackerhit_time_Tcut_layer0", ";Time (ns);Events" ,800,-1,20);
  // h_trackerhit_time_Tcut_layer1  = new TH1F("trackerhit_time_Tcut_layer1", ";Time (ns);Events" ,800,-1,20);
  // h_trackerhit_time_Tcut_layer2  = new TH1F("trackerhit_time_Tcut_layer2", ";Time (ns);Events" ,800,-1,20);
  // h_trackerhit_time_Tcut_layer3  = new TH1F("trackerhit_time_Tcut_layer3", ";Time (ns);Events" ,800,-1,20);
  // h_trackerhit_time_Tcut_layer4  = new TH1F("trackerhit_time_Tcut_layer4", ";Time (ns);Events" ,800,-1,20);
  // h_trackerhit_time_Tcut_layer5  = new TH1F("trackerhit_time_Tcut_layer5", ";Time (ns);Events" ,800,-1,20);
  // h_trackerhit_time_Tcut_layer6  = new TH1F("trackerhit_time_Tcut_layer6", ";Time (ns);Events" ,800,-1,20);
  // h_trackerhit_time_Tcut_layer7  = new TH1F("trackerhit_time_Tcut_layer7", ";Time (ns);Events" ,800,-1,20); 

  h_hit_edep_layer0   = new TH1F("hit_edep_layer0", ";Deposited charge (electrons);Hits" ,edp_binNum_e,0,edp_rangeMax_e); //5000, 0, 50000);
  h_hit_edep_layer1   = new TH1F("hit_edep_layer1", ";Deposited charge (electrons);Hits" ,edp_binNum_e,0,edp_rangeMax_e);
  h_hit_edep_layer2   = new TH1F("hit_edep_layer2", ";Deposited charge (electrons);Hits" ,edp_binNum_e,0,edp_rangeMax_e);
  h_hit_edep_layer3   = new TH1F("hit_edep_layer3", ";Deposited charge (electrons);Hits" ,edp_binNum_e,0,edp_rangeMax_e);
  h_hit_edep_layer4   = new TH1F("hit_edep_layer4", ";Deposited charge (electrons);Hits" ,edp_binNum_e,0,edp_rangeMax_e);
  h_hit_edep_layer5   = new TH1F("hit_edep_layer5", ";Deposited charge (electrons);Hits" ,edp_binNum_e,0,edp_rangeMax_e);
  h_hit_edep_layer6   = new TH1F("hit_edep_layer6", ";Deposited charge (electrons);Hits" ,edp_binNum_e,0,edp_rangeMax_e); //(100, 0, 36000)
  h_hit_edep_layer7   = new TH1F("hit_edep_layer7", ";Deposited charge (electrons);Hits" ,edp_binNum_e,0,edp_rangeMax_e);//5000,0,50000
  h_hit_edep_layer8   = new TH1F("hit_edep_layer8", ";Deposited charge (electrons);Hits" ,edp_binNum_e,0,edp_rangeMax_e);//5000,0,50000
  
  h_trackerhit_time  = new TH1F("trackerhit_time", ";Time (ns);Events" ,200,-0.2,0.5); //-0.2 - 0.5
  h_trackerhit_time_layer0  = new TH1F("trackerhit_time_layer0", ";Time (ns);Events" ,200,-0.2,0.5); //,200000,0,20
  h_trackerhit_time_layer1  = new TH1F("trackerhit_time_layer1", ";Time (ns);Events" ,200,-0.2,0.5);
  h_trackerhit_time_layer2  = new TH1F("trackerhit_time_layer2", ";Time (ns);Events" ,200,-0.2,0.5);
  h_trackerhit_time_layer3  = new TH1F("trackerhit_time_layer3", ";Time (ns);Events" ,200,-0.2,0.5);
  h_trackerhit_time_layer4  = new TH1F("trackerhit_time_layer4", ";Time (ns);Events" ,200,-0.2,0.5);
  h_trackerhit_time_layer5  = new TH1F("trackerhit_time_layer5", ";Time (ns);Events" ,200,-0.2,0.5);
  h_trackerhit_time_layer6  = new TH1F("trackerhit_time_layer6", ";Time (ns);Events" ,200,-0.2,0.5);
  h_trackerhit_time_layer7  = new TH1F("trackerhit_time_layer7", ";Time (ns);Events" ,200,-0.2,0.5);
  h_trackerhit_time_layer8  = new TH1F("trackerhit_time_layer8", ";Time (ns);Events" ,200,-0.2,0.5);

  //number of hits per cluster
  h_thclen = new TH1F("thclen", ";Number of Hit Constituents; Events", 50, 0, 50); //Total number of hit constituents 
  //number of hits per cluster per layer: 
  h_thclen_layer0 = new TH1F("thclen_layer0", ";Number of Hit Constituents; Events", 50, 0, 50);
  h_thclen_layer1 = new TH1F("thclen_layer1", ";Number of Hit Constituents; Events", 50, 0, 50);
  h_thclen_layer2 = new TH1F("thclen_layer2", ";Number of Hit Constituents; Events", 50, 0, 50);
  h_thclen_layer3 = new TH1F("thclen_layer3", ";Number of Hit Constituents; Events", 50, 0, 50);
  h_thclen_layer4 = new TH1F("thclen_layer4", ";Number of Hit Constituents; Events", 50, 0, 50);
  h_thclen_layer5 = new TH1F("thclen_layer5", ";Number of Hit Constituents; Events", 50, 0, 50);
  h_thclen_layer6 = new TH1F("thclen_layer6", ";Number of Hit Constituents; Events", 50, 0, 50);
  h_thclen_layer7 = new TH1F("thclen_layer7", ";Number of Hit Constituents; Events", 50, 0, 50);
  h_thclen_layer8 = new TH1F("thclen_layer8", ";Number of Hit Constituents; Events", 50, 0, 50);

  //number of hits per cluster for hit clusters > 15:
  h_thclen_cut = new TH1F("thclen_cut", ";Number of Hit Constituents; Events", 100, 15, 115); //Total number of hit constituents 
  //number of hits per cluster per layer: 
  h_thclen_layer0_cut = new TH1F("thclen_layer0_cut", ";Number of Hit Constituents; Events", 100, 15, 115);
  h_thclen_layer1_cut = new TH1F("thclen_layer1_cut", ";Number of Hit Constituents; Events", 100, 15, 115);
  h_thclen_layer2_cut = new TH1F("thclen_layer2_cut", ";Number of Hit Constituents; Events", 100, 15, 115);
  h_thclen_layer3_cut = new TH1F("thclen_layer3_cut", ";Number of Hit Constituents; Events", 100, 15, 115);
  h_thclen_layer4_cut = new TH1F("thclen_layer4_cut", ";Number of Hit Constituents; Events", 100, 15, 115);
  h_thclen_layer5_cut = new TH1F("thclen_layer5_cut", ";Number of Hit Constituents; Events", 100, 15, 115);
  h_thclen_layer6_cut = new TH1F("thclen_layer6_cut", ";Number of Hit Constituents; Events", 100, 15, 115);
  h_thclen_layer7_cut = new TH1F("thclen_layer7_cut", ";Number of Hit Constituents; Events", 100, 15, 115);
  h_thclen_layer8_cut = new TH1F("thclen_layer8_cut", ";Number of Hit Constituents; Events", 100, 15, 115);

  //2D histogram for all detector layers for Cluster Edep vs Number of Hits per Cluster:
  h_cluster_edep_thlen = new TH2F("edep_vs_hitNum" , ";Hits/Cluster; Energy Deposited (GeV)" , 100, 0, 100, edp_binNum,0,edp_rangeMax);
  h_cluster_edep_thlen_cut = new TH2F("edep_vs_hitNum_cut" , ";Hits/Cluster; Energy Deposited (GeV)" , 100, 15, 115, edp_binNum,0,edp_rangeMax);

//3D HISTO for X vs Y vs Z position in digitized and truth 
h_3DPosition_digi = new TH3F("3DPosition_digi", "3D Digitized Position;x[mm];y[mm];z[mm]", 100, 0, 1600, 100, 0, 1600, 100,  -2500, 2500);
h_3DPosition_cdigi = new TH3D ("3DPosition_cdigi", "3D Digitized Position;#theta;r[mm];z[mm]", numbins_all/10, 0, 3.14, numbins_all/10, -rmax_all, rmax_all, numbins_all/10,  zmin_all, zmax_all);
 
}

void ClusterHists::fill(const EVENT::TrackerHit* trkhit)
{
 //Calculate energy deposited
  float EDep = trkhit->getEDep();
  float toa = trkhit->getTime(); //time of arrival -- Juliet
  // Correcting for the propagation time
  /*dd4hep::rec::Vector3D pos = trkhit->getPosition();
  double hitR = pos.r();
  double m_beta = 1.0;
  double tmin = -0.09;//-90.0; //ns - get min and max from the config file
  double tmax = 0.15; //90.0; //ns
  double dt = hitR / (TMath::C() * m_beta / 1e6);
  toa -= dt;*/
  
  

  //Calculating theta
  float x = trkhit->getPosition()[0];
  float y = trkhit->getPosition()[1];
  float z = trkhit->getPosition()[2];
  float r = sqrt(pow(x,2)+pow(y,2));
  float incidentTheta = std::atan(r/z);
  float theta_min = std::atan(r/z);
  
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
  h_trackerhit_time->Fill(toa);//if (toa > tmin && toa < tmax) h_trackerhit_time->Fill(toa);

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
    
    //skip events out of the time range -> Correcting for the propagation time
    // if (toa < tmin || toa > tmax){
    //   continue;
    // }

    //hits/cluster cut: 
    if(rawHits.size() > 15){
      h_thclen_cut->Fill(rawHits.size());
      h_cluster_edep_thlen_cut->Fill(rawHits.size(), EDep);
    }
    h_cluster_edep_thlen->Fill(rawHits.size(), EDep);
      
    h_hits_by_layer->Fill(layerID);
    h_z_r_hits->Fill(z,r);
    h_x_y_hits->Fill(x,y);
    
    
    h_3DPosition_digi->Fill(x, y, z);
    h_3DPosition_cdigi->Fill(theta_min, r, z);

    if(layerID==0){
      h_z_layer0->Fill(z);
      h_r_layer0->Fill(r);
      //setting cluster and time of arrive with if j == 0 becuase we do not want fill the cluster edep and toa for clusters multiple times
      if(j==0){
        h_cluster_edep_layer0->Fill(EDep);
        h_trackerhit_time_layer0->Fill(toa);
        h_thclen_layer0->Fill(rawHits.size());
        if(rawHits.size() > 15) h_thclen_layer0_cut->Fill(rawHits.size());
        // if(toa < 0.2){
        //   h_cluster_edep_Tcut_layer0->Fill(EDep);
        //   h_trackerhit_time_Tcut_layer0->Fill(toa);
        // }
      }
      h_hit_edep_layer0->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
      

      //std::cout << "tracker hit time: " << toa << ", and edep cluster: "<< EDep <<", and edep hit : " << hitConstituent->getEDep() << ", with hit num: " << rawHits.size() << std::endl;

    }
    if(layerID==1){
      h_z_layer1->Fill(z);
      h_r_layer1->Fill(r);
      if(j == 0){
        h_cluster_edep_layer1->Fill(EDep);
        h_trackerhit_time_layer1->Fill(toa); //time of arrive in layer 2 -- Juliet
        h_thclen_layer1->Fill(rawHits.size());
        if(rawHits.size() > 15) h_thclen_layer1_cut->Fill(rawHits.size());
        // if(toa < 0.2){
        //   h_cluster_edep_Tcut_layer1->Fill(EDep);
        //   h_trackerhit_time_Tcut_layer1->Fill(toa);
        // }
      }   
      h_hit_edep_layer1->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
    }
    if(layerID==2){
      h_z_layer2->Fill(z);
      h_r_layer2->Fill(r);
      if(j==0){
        h_cluster_edep_layer2->Fill(EDep);
        h_trackerhit_time_layer2->Fill(toa);
        h_thclen_layer2->Fill(rawHits.size());
        if(rawHits.size() > 15) h_thclen_layer2_cut->Fill(rawHits.size());
        // if(toa < 0.2){
        //   h_cluster_edep_Tcut_layer2->Fill(EDep);
        //   h_trackerhit_time_Tcut_layer2->Fill(toa);
        // }
      }
      h_hit_edep_layer2->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
     
    }
    if(layerID==3){
      h_z_layer3->Fill(z);
      h_r_layer3->Fill(r);
      if(j==0){
        h_cluster_edep_layer3->Fill(EDep);
        h_trackerhit_time_layer3->Fill(toa);
        h_thclen_layer3->Fill(rawHits.size());
        if(rawHits.size() > 15) h_thclen_layer3_cut->Fill(rawHits.size());
        // if(toa < 0.2){
        //   h_cluster_edep_Tcut_layer3->Fill(EDep);
        //   h_trackerhit_time_Tcut_layer3->Fill(toa);
        // }
      }
      h_hit_edep_layer3->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
    }
    if(layerID==4){
      h_z_layer4->Fill(z);
      h_r_layer4->Fill(r);
      if(j==0){
        h_cluster_edep_layer4->Fill(EDep);
        h_trackerhit_time_layer4->Fill(toa);
        h_thclen_layer4->Fill(rawHits.size());
        if(rawHits.size() > 15) h_thclen_layer4_cut->Fill(rawHits.size());
        // if(toa < 0.2){
        //   h_cluster_edep_Tcut_layer4->Fill(EDep);
        //   h_trackerhit_time_Tcut_layer4->Fill(toa);
        // }
      }
      h_hit_edep_layer4->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
    }
    if(layerID==5){
      h_z_layer5->Fill(z);
      h_r_layer5->Fill(r);
      if(j==0){
        h_cluster_edep_layer5->Fill(EDep);
        h_trackerhit_time_layer5->Fill(toa);
        h_thclen_layer5->Fill(rawHits.size());
        if(rawHits.size() > 15) h_thclen_layer5_cut->Fill(rawHits.size());
        // if(toa < 0.2){
        //   h_cluster_edep_Tcut_layer5->Fill(EDep);
        //   h_trackerhit_time_Tcut_layer5->Fill(toa);
        // }
      }
      h_hit_edep_layer5->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
    }
    if(layerID==6){
      h_z_layer6->Fill(z);
      h_r_layer6->Fill(r);
      if(j==0){
        h_cluster_edep_layer6->Fill(EDep);
        h_trackerhit_time_layer6->Fill(toa);
        h_thclen_layer6->Fill(rawHits.size());
        if(rawHits.size() > 15) h_thclen_layer6_cut->Fill(rawHits.size());
        // if(toa < 0.2){
        //   h_cluster_edep_Tcut_layer6->Fill(EDep);
        //   h_trackerhit_time_Tcut_layer6->Fill(toa);
        // }
      }
      h_hit_edep_layer6->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
    }
    if(layerID==7){
      h_z_layer7->Fill(z);
      h_r_layer7->Fill(r);
      if(j==0){
        h_cluster_edep_layer7->Fill(EDep);
        h_trackerhit_time_layer7->Fill(toa);
        h_thclen_layer7->Fill(rawHits.size());
        if(rawHits.size() > 15) h_thclen_layer7_cut->Fill(rawHits.size());
        // if(toa < 0.2){
        //   h_cluster_edep_Tcut_layer7->Fill(EDep);
        //   h_trackerhit_time_Tcut_layer7->Fill(toa);
        // }
      }
      h_hit_edep_layer7->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
    }
    if(layerID==8){
      h_z_layer8->Fill(z);
      h_r_layer8->Fill(r);
      if(j==0){
        h_cluster_edep_layer8->Fill(EDep);
        h_trackerhit_time_layer8->Fill(toa);
        h_thclen_layer8->Fill(rawHits.size());
        if(rawHits.size() > 15) h_thclen_layer8_cut->Fill(rawHits.size());
      }
      h_hit_edep_layer8->Fill(hitConstituent->getEDep()); //energy in electrons -- Juliet
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
  //NEW 2/19/2025
  // if (toa > tmin && toa < tmax){
  //   if(toa < 0.2e-3){
  //     h_cluster_edep_Tcut->Fill(EDep);
  //     h_trackerhit_time_Tcut->Fill(toa);
  //   }
  // }
  h_edep_r->Fill(r,EDep);
  h_edep_cluster->Fill(EDep,cluster_size_tot);

  h_thclen->Fill(rawHits.size()); // -- JULIET

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
