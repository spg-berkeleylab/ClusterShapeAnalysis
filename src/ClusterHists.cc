#include "ClusterShapeAnalysis/ClusterHists.h"
#include "marlin/VerbosityLevels.h"

#include <EVENT/TrackerHit.h>
#include <EVENT/SimTrackerHit.h>
#include <UTIL/CellIDDecoder.h>
#include <IMPL/TrackerHitPlaneImpl.h>
#include <UTIL/LCTrackerConf.h>


ClusterHists::ClusterHists()
{
  h_size_theta_y    = new TH2F("cluster_size_vs_theta_y" , ";Cluster #theta; Cluster size" , 100, 0,  3.14,  31,  -0.5,  30.5  );
  h_size_theta_x    = new TH2F("cluster_size_vs_theta_x" , ";Cluster #theta; Cluster size" , 100, 0,  3.14,  31,  -0.5,  30.5  );
  h_size_theta_tot  = new TH2F("cluster_size_vs_theta_tot" , ";Cluster #theta; Cluster size" , 100, 0,  3.14,  31,  -0.5,  30.5  );
  for(int i=0; i<=8; i++){
    h_size_theta_tot_layer[i]  = new TH2F(Form("cluster_size_vs_theta_tot_layer%d", i), ";Cluster #theta; Cluster size" , 100, 0,  3.14,  31,  -0.5,  30.5  );
    h_size_theta_x_layer[i]  = new TH2F(Form("cluster_size_vs_theta_x_layer%d", i), ";Cluster #theta; Cluster size" , 100, 0,  3.14,  31,  -0.5,  30.5  );
    h_size_theta_y_layer[i]  = new TH2F(Form("cluster_size_vs_theta_y_layer%d", i), ";Cluster #theta; Cluster size" , 100, 0,  3.14,  31,  -0.5,  30.5  );
    h_cluster_layer[i] = new TH1F(Form("nClusters_layer%d", i), ";Number of Hits / Cluster; Entries",50,0,50);
    h_hits_layer[i] = new TH1F(Form("nHits_layer%d", i), ";Total number of hits; Entries",50,0,50);
  }
  h_size_r_tot      = new TH2F("cluster_size_vs_R_tot" , ";Cluster R (x^2+y^2)^(1/2) (mm); Cluster size" , 100, 20,  120,  31,  -0.5,  30.5  );
  h_size_r_tot_0     = new TH2F("cluster_size_vs_R_tot_0" , ";Cluster R (x^2+y^2)^(1/2) (mm); Cluster size" , 100, 20,  120,  31,  -0.5,  30.5  );
  h_size_r_tot_1     = new TH2F("cluster_size_vs_R_tot_1" , ";Cluster R (x^2+y^2)^(1/2) (mm); Cluster size" , 100, 20,  120,  31,  -0.5,  30.5  );
  h_size_r_tot_2     = new TH2F("cluster_size_vs_R_tot_2" , ";Cluster R (x^2+y^2)^(1/2) (mm); Cluster size" , 100, 20,  120,  31,  -0.5,  30.5  );
  h_size_r_tot_3     = new TH2F("cluster_size_vs_R_tot_3" , ";Cluster R (x^2+y^2)^(1/2) (mm); Cluster size" , 100, 20,  120,  31,  -0.5,  30.5  );
  h_cluster_pos   = new TH2F("cluster_position"      , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_0 = new TH2F("cluster_position_0"    , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_1 = new TH2F("cluster_position_1"    , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_2 = new TH2F("cluster_position_2"    , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_cluster_pos_3 = new TH2F("cluster_position_3"    , ";z; r"                         , 100, -500, 500, 100, 0, 200);
  h_clusters_by_layer   = new TH1F("numClusters_by_layer"      , ";Layer Index; Number of Clusters",8,0,8);
  h_hits_by_layer   = new TH1F("numhits_by_layer"      , ";Layer Index; Number of Hits",8,0,8);
  h_theta         = new TH1F("theta"                 , ";Theta;Number of Clusters"       ,100,0,3.15);
  h_cluster_edep     = new TH1F("Clusters_edep"          , ";Energy Deposited (GeV);Clusters" ,100,0,0.0005);
  h_cluster_edep_BX     = new TH1F("Clusters_edep_BX"          , ";Energy Deposited (GeV);Clusters/BX" ,100,0,0.0005);
  h_hit_edep     = new TH1F("Hits_edep"          , ";Deposited charge (electrons);Hits" ,5000,0,50000);
  h_edep_r   = new TH2F("edep_vs_r" , ";Cluster R (x^2+y^2)^(1/2) (mm); Energy Deposited (GeV)" , 100, 20,  120,  100,  0,  0.002 );
  h_edep_cluster   = new TH2F("edep_vs_cluster_size" , "; Energy Deposited (GeV); Total Cluster Size" ,100,  0,  0.002, 100, -0.5, 99.5 );
  h_avgHits = new TH1F("nHitsperCluster", ";Number of Hits / Cluster; Clusters",20,0,20);
  h_sysID = new TH1F("systemID", ";System ID; Number of Clusters",20,0,20);
  h_cluster_timing = new TH1F("cluster_timing", "Time of arrival of clusters [ns]", 200, -1, 1);
  h_hit_timing = new TH1F("hit_timing", "Time of arrival of hits [ns]", 200, -1, 1);

  // Create position histograms for tracker hits
  int numbins_z = 1000;
  int numbins_r = 800;
  int rmin_all = 0;
  int rmax_all = 1600;
  int zmin_all = -2500;
  int zmax_all = 2500;
  h_x   = new TH1F("x  " , ";x   ; Num Hits" , numbins_r, -rmax_all, rmax_all);
  h_y   = new TH1F("y  " , ";y   ; Num Hits" , numbins_r, -rmax_all, rmax_all);
  h_z   = new TH1F("z  " , ";z   ; Num Hits" , 5000,  zmin_all, zmax_all);
  h_r   = new TH1F("r  " , ";r   ; Num Hits" , numbins_r,  rmin_all, rmax_all);
  h_z_r = new TH2F("z_r" , ";z_r ; r"        , numbins_z,  zmin_all, zmax_all, numbins_r, rmin_all, rmax_all);
  h_x_y = new TH2F("x_y" , ";x_y ; r"        , numbins_r, -rmax_all, rmax_all, numbins_r, -rmax_all, rmax_all);
  h_z_r_hits = new TH2F("z_r_hits" , ";z_r ; r"        , numbins_z,  zmin_all, zmax_all, numbins_z, rmin_all, rmax_all);
  h_x_y_hits = new TH2F("x_y_hits" , ";x_y ; r"        , numbins_r, -rmax_all, rmax_all, numbins_r, -rmax_all, rmax_all);
  for(int i=0; i<=8; i++){
    h_z_layer[i]   = new TH1F(Form("hits_vs_z_layer%d", i), ";z   ; Num Hits" , numbins_z,  zmin_all, zmax_all);
    h_r_layer[i]   = new TH1F(Form("hits_vs_r_layer%d", i), ";r   ; Num Hits" , numbins_r,  rmin_all, rmax_all);
  }
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
}

void ClusterHists::fill(const EVENT::TrackerHit* trkhit)
{
  //Calculate energy deposited
  float EDep = trkhit->getEDep();
  
  //Calculating theta
  float x = trkhit->getPosition()[0];
  float y = trkhit->getPosition()[1];
  float z = trkhit->getPosition()[2];
  float r = sqrt(pow(x,2)+pow(y,2));
  float incidentTheta = std::atan(r/z);
  
  if(incidentTheta<0)
    incidentTheta += M_PI;
  streamlog_out(DEBUG6) << "Theta: " << incidentTheta << std::endl;

  h_cluster_timing -> Fill(trkhit->getTime());
  
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
    h_hit_edep->Fill(hitConstituent->getEDep());
    h_hit_timing -> Fill(hitConstituent->getTime());
    const double *localPos = hitConstituent->getPosition();
    float x_local = localPos[0];
    float y_local = localPos[1];
    streamlog_out(DEBUG6) << "Hit: "<<j<< " Local y: " << y_local << ", local x: " << x_local << std::endl;

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

  h_sysID->Fill(systemID);

  streamlog_out(DEBUG9) << "Layer ID is: " << layerID << std::endl;

  // Fill for all hits
  h_size_theta_y->Fill(incidentTheta, cluster_size_y);
  h_size_theta_x->Fill(incidentTheta, cluster_size_x);
  h_size_theta_tot->Fill(incidentTheta, cluster_size_tot);
  h_size_r_tot->Fill(r, cluster_size_tot);

  h_theta->Fill(incidentTheta);
  h_cluster_pos->Fill(z,r);
  h_clusters_by_layer->Fill(layerID);
  h_avgHits->Fill(cluster_size_tot);

  // tracker hit hists
  h_x->Fill(x);  
  h_y->Fill(y);  
  h_z->Fill(z);  
  h_r->Fill(r);
  h_z_r->Fill(z,r);
  h_x_y->Fill(x,y);
  h_cluster_layer[layerID]->Fill(loopsize);
  
  for (size_t j=0; j<loopsize; ++j) {
    h_hits_by_layer->Fill(layerID);
    h_z_r_hits->Fill(z,r);
    h_x_y_hits->Fill(x,y);
    h_z_layer[layerID]->Fill(z);
    h_r_layer[layerID]->Fill(r);
    h_hits_layer[layerID]->Fill(loopsize);
    h_size_theta_tot_layer[layerID]->Fill(incidentTheta, cluster_size_tot);
    h_size_theta_x_layer[layerID]->Fill(incidentTheta, cluster_size_x);
    h_size_theta_y_layer[layerID]->Fill(incidentTheta, cluster_size_y);
  }  

  // vertex hists
  h_x_vx->Fill(x);  
  h_y_vx->Fill(y);  
  h_z_vx->Fill(z);  
  h_r_vx->Fill(r);  
  h_z_r_vx->Fill(z,r);
  h_x_y_vx->Fill(x,y);

  h_cluster_edep->Fill(EDep);
  h_cluster_edep_BX->Fill(EDep);
  h_edep_r->Fill(r,EDep);
  h_edep_cluster->Fill(EDep,cluster_size_tot);


  // Fill energy deposition histograms based on angle
  /*float theta_deg = incidentTheta * (180/3.1416);
  if(theta_deg < 5 || theta_deg > 175){
    h_edep_0deg->Fill(EDep);
  }
  if(theta_deg > 89 && theta_deg < 91){
    h_edep_90deg->Fill(EDep);
  } */
  
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
