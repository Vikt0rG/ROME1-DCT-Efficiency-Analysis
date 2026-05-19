#include <string.h>

void postDecoder(std::string configuration_filename="settings.conf"){
  gROOT->ProcessLine("#include <vector>");

  // cabling map (-1 means connector not connected)
  std::map<std::string, std::vector<int>> cabling_maps;
  cabling_maps["BIL680"] = {0,-1,1,2,-1,3};
  cabling_maps["BIL520"] = {0,-1,1,2,-1,-1};
  cabling_maps["BIS1"] = {0,1,2,3,4,5};
  cabling_maps["BIS2_6"] = {-1,0,1,2,3,4};

  // Read configuration file
  std::ifstream config_file(configuration_filename);
  std::string line, key, value;
  std::string rpc_type = "";
  int dt_max = 300, dt_min = -300, daq_window = 128; // default values
  while (std::getline(config_file, line)) {
    size_t comment_pos = line.find('#');
    if (comment_pos != std::string::npos) line = line.substr(0, comment_pos);
    size_t eq_pos = line.find('=');
    if (eq_pos != std::string::npos) {
      key = line.substr(0, eq_pos);
      value = line.substr(eq_pos + 1);
      key.erase(key.find_last_not_of(" \t\n\r") + 1);
      key.erase(0, key.find_first_not_of(" \t\n\r"));
      value.erase(value.find_last_not_of(" \t\n\r") + 1);
      value.erase(0, value.find_first_not_of(" \t\n\r"));
      
      if (key == "dt_max") dt_max = std::stoi(value);
      else if (key == "dt_min") dt_min = std::stoi(value);
      else if (key == "RPC_type") rpc_type = value;
      else if (key == "daq_window") daq_window = std::stoi(value);
    }
  }
  config_file.close();
  
  if (cabling_maps.find(rpc_type) == cabling_maps.end()) {
    std::cerr << "Error: RPC type '" << rpc_type << "' not found in cabling maps!" << std::endl;
    return;
  }

  std::cout<<"Inputs taken from config file:"<<std::endl;
  std::cout<<"\tdt_max = "<<dt_max<<std::endl;
  std::cout<<"\tdt_min = "<<dt_min<<std::endl;
  std::cout<<"\tRPC type: "<<rpc_type<<std::endl;

  std::vector<int> cabling_map = cabling_maps[rpc_type];
  int trig_channel = 143;

  TFile *o = TFile::Open("out.root","RECREATE");
  if (!o) { return; }

  int nHits=0;
  std::vector<int> hit_dct, hit_channel, hit_clk, hit_rawbcid, hit_bcid, hit_rawbcout, hit_bcout;
  std::vector<int> hit_time1, hit_time2, hit_rise, hit_layer, hit_strip;
 
  int nProc=0;
  std::vector<int> proc_layer, proc_strip, proc_time1, proc_time2;
  std::vector<int> proc_tot1, proc_tot2, proc_deltat, proc_dt_trig1, proc_dt_trig2;
  int trig_time=-1;

  int nClus_l0=0, nClus_l1=0, nClus_l2=0;
  std::vector<int> clus_strip_l0, clus_strip_l1, clus_strip_l2;
  std::vector<int> clus_time1_l0, clus_time1_l1, clus_time1_l2;
  std::vector<int> clus_time2_l0, clus_time2_l1, clus_time2_l2;
  std::vector<int> clus_tot1_l0, clus_tot1_l1, clus_tot1_l2;
  std::vector<int> clus_tot2_l0, clus_tot2_l1, clus_tot2_l2;
  std::vector<int> clus_dt_l0, clus_dt_l1, clus_dt_l2;
  std::vector<int> clus_nstrips_l0, clus_nstrips_l1, clus_nstrips_l2;
  std::vector<float> clus_center_l0, clus_center_l1, clus_center_l2;
  std::vector<int> clus_mult_l0, clus_mult_l1, clus_mult_l2;
  
  TTree *tree = new TTree("tree","Events");
  int event=0;
  tree->Branch("event",&event);
  tree->Branch("nHits",&nHits);
  tree->Branch("hit_dct",&hit_dct);
  tree->Branch("hit_channel",&hit_channel);
  tree->Branch("hit_clk",&hit_clk);
  tree->Branch("hit_rawbcid",&hit_rawbcid);
  tree->Branch("hit_bcid",&hit_bcid);
  tree->Branch("hit_rawbcout",&hit_rawbcout);
  tree->Branch("hit_bcout",&hit_bcout);
  tree->Branch("hit_time1",&hit_time1);
  tree->Branch("hit_time2",&hit_time2);
  tree->Branch("hit_rise",&hit_rise);
  tree->Branch("hit_layer",&hit_layer);
  tree->Branch("nProc",&nProc);
  tree->Branch("proc_layer",&proc_layer);
  tree->Branch("proc_strip",&proc_strip);
  tree->Branch("proc_time1",&proc_time1);
  tree->Branch("proc_time2",&proc_time2);
  tree->Branch("proc_tot1",&proc_tot1);
  tree->Branch("proc_tot2",&proc_tot2);
  tree->Branch("proc_deltat",&proc_deltat);
  tree->Branch("proc_dt_trig1",&proc_dt_trig1);
  tree->Branch("proc_dt_trig2",&proc_dt_trig2);
  tree->Branch("trig_time",&trig_time);
  tree->Branch("nClus_l0",&nClus_l0);  tree->Branch("nClus_l1",&nClus_l1);  tree->Branch("nClus_l2",&nClus_l2);
  tree->Branch("clus_strip_l0",&clus_strip_l0); tree->Branch("clus_strip_l1",&clus_strip_l1); tree->Branch("clus_strip_l2",&clus_strip_l2);
  tree->Branch("clus_time1_l0",&clus_time1_l0); tree->Branch("clus_time1_l1",&clus_time1_l1); tree->Branch("clus_time1_l2",&clus_time1_l2);
  tree->Branch("clus_time2_l0",&clus_time2_l0); tree->Branch("clus_time2_l1",&clus_time2_l1); tree->Branch("clus_time2_l2",&clus_time2_l2);
  tree->Branch("clus_tot1_l0",&clus_tot1_l0); tree->Branch("clus_tot1_l1",&clus_tot1_l1); tree->Branch("clus_tot1_l2",&clus_tot1_l2);
  tree->Branch("clus_tot2_l0",&clus_tot2_l0); tree->Branch("clus_tot2_l1",&clus_tot2_l1); tree->Branch("clus_tot2_l2",&clus_tot2_l2);
  tree->Branch("clus_dt_l0",&clus_dt_l0); tree->Branch("clus_dt_l1",&clus_dt_l1); tree->Branch("clus_dt_l2",&clus_dt_l2);
  tree->Branch("clus_nstrips_l0",&clus_nstrips_l0); tree->Branch("clus_nstrips_l1",&clus_nstrips_l1); tree->Branch("clus_nstrips_l2",&clus_nstrips_l2);
  tree->Branch("clus_center_l0",&clus_center_l0);  tree->Branch("clus_center_l1",&clus_center_l1);  tree->Branch("clus_center_l2",&clus_center_l2);
  tree->Branch("clus_mult_l0",&clus_mult_l0);  tree->Branch("clus_mult_l1",&clus_mult_l1);  tree->Branch("clus_mult_l2",&clus_mult_l2);

  TH2D * tot_eta1_ly0_2D = new TH2D("tot_eta1_ly0_2D","tot_eta1_ly0_2D",48,0,48,30,0,30);
  TH2D * tot_eta1_ly1_2D = new TH2D("tot_eta1_ly1_2D","tot_eta1_ly1_2D",48,0,48,30,0,30);
  TH2D * tot_eta1_ly2_2D = new TH2D("tot_eta1_ly2_2D","tot_eta1_ly2_2D",48,0,48,30,0,30);
  TH2D * tot_eta2_ly0_2D = new TH2D("tot_eta2_ly0_2D","tot_eta2_ly0_2D",48,0,48,30,0,30);
  TH2D * tot_eta2_ly1_2D = new TH2D("tot_eta2_ly1_2D","tot_eta2_ly1_2D",48,0,48,30,0,30);
  TH2D * tot_eta2_ly2_2D = new TH2D("tot_eta2_ly2_2D","tot_eta2_ly2_2D",48,0,48,30,0,30);
  TH1D * tot_eta1_ly0_1D = new TH1D("tot_eta1_ly0_1D","tot_eta1_ly0_1D",30,0,30);
  TH1D * tot_eta1_ly1_1D = new TH1D("tot_eta1_ly1_1D","tot_eta1_ly1_1D",30,0,30);
  TH1D * tot_eta1_ly2_1D = new TH1D("tot_eta1_ly2_1D","tot_eta1_ly2_1D",30,0,30);
  TH1D * tot_eta2_ly0_1D = new TH1D("tot_eta2_ly0_1D","tot_eta2_ly0_1D",30,0,30);
  TH1D * tot_eta2_ly1_1D = new TH1D("tot_eta2_ly1_1D","tot_eta2_ly1_1D",30,0,30);
  TH1D * tot_eta2_ly2_1D = new TH1D("tot_eta2_ly2_1D","tot_eta2_ly2_1D",30,0,30);
  TH2D * dt_trig_vs_strip_eta1 = new TH2D("dt_trig_vs_strip_eta1","dt_trig_vs_strip_eta1",144,0,144,350,-300,50);
  TH2D * dt_trig_vs_strip_eta2 = new TH2D("dt_trig_vs_strip_eta2","dt_trig_vs_strip_eta2",144,0,144,350,-300,50);
  TH1D * strips_eta1_layer0 = new TH1D("strips_eta1_layer0","strips_eta1_layer0",48,0,48);
  TH1D * strips_eta1_layer1 = new TH1D("strips_eta1_layer1","strips_eta1_layer1",48,0,48);
  TH1D * strips_eta1_layer2 = new TH1D("strips_eta1_layer2","strips_eta1_layer2",48,0,48);
  TH1D * strips_eta2_layer0 = new TH1D("strips_eta2_layer0","strips_eta2_layer0",48,0,48);
  TH1D * strips_eta2_layer1 = new TH1D("strips_eta2_layer1","strips_eta2_layer1",48,0,48);
  TH1D * strips_eta2_layer2 = new TH1D("strips_eta2_layer2","strips_eta2_layer2",48,0,48);

  strips_eta1_layer0->SetLineColor(kRed); strips_eta1_layer1->SetLineColor(kBlue); strips_eta1_layer2->SetLineColor(kGreen);
  strips_eta2_layer0->SetLineColor(kRed); strips_eta2_layer1->SetLineColor(kBlue); strips_eta2_layer2->SetLineColor(kGreen);

  int clk,word,rawbcout,bcout;

  int evts_eta1[3] = {0, 0, 0}; int evts_eta2[3] = {0, 0, 0};
  int evts_or[3] = {0, 0, 0}; int evts_and[3] = {0, 0, 0}; int trig_evts = 0;
  int evts_eta1_rpc_trig[3] = {0, 0, 0}; int evts_eta2_rpc_trig[3] = {0, 0, 0};
  int evts_OR_rpc_trig[3] = {0, 0, 0}; int evts_AND_rpc_trig[3] = {0, 0, 0};
  int evts_rpc_trig_denom[3] = {0,0,0};
  
  TFile *myFile = TFile::Open("raw.root");
  TTreeReader myReader("tree", myFile);
  TTreeReaderValue<std::vector<int>> my_dct(myReader, "hit_dct");
  TTreeReaderValue<std::vector<int>> my_channel(myReader, "hit_channel");
  TTreeReaderValue<std::vector<int>> my_clk(myReader, "hit_clk");
  TTreeReaderValue<std::vector<int>> my_rawbcid(myReader, "hit_rawbcid");
  TTreeReaderValue<std::vector<int>> my_time1(myReader, "hit_time1");
  TTreeReaderValue<std::vector<int>> my_time2(myReader, "hit_time2");
  TTreeReaderValue<std::vector<int>> my_rise(myReader, "hit_rise");
  TTreeReaderValue<std::vector<int>> my_layer(myReader, "hit_layer");
  TTreeReaderValue<std::vector<int>> my_strip(myReader, "hit_strip");
  Long64_t tot_events = myReader.GetEntries();
  
  o->cd();
  while(myReader.Next()){
    event++;

    if (event % 100 == 0) {
      int barWidth = 50;
      float progress = (float)event / tot_events;
      std::cout << "Progress: [";
      int pos = barWidth * progress;
      for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
      }
      std::cout << "] " << (int)(progress * 100.0) << "% (" << event << "/" << tot_events << " events)\r";
      std::cout.flush();
    } else if (event == tot_events - 1) {
        std::cout << "Progress: [" << std::string(50, '=') << "] 100% (" << tot_events << "/" << tot_events << " events)" << std::endl;
        std::cout.flush();
    }

    nHits = my_channel->size();

    // FAST COPY: Block memory assignment
    hit_dct = *my_dct;
    hit_channel = *my_channel;
    hit_clk = *my_clk;
    hit_rawbcid = *my_rawbcid;
    hit_time1 = *my_time1;
    hit_time2 = *my_time2;
    hit_rise = *my_rise;

    hit_layer.resize(nHits);
    hit_strip.resize(nHits);
    hit_bcid.resize(nHits);

    bool eta1_hit_flag[3] = {0, 0, 0};
    bool eta2_hit_flag[3] = {0, 0, 0};

    int bc0 = 0;
    if (nHits > 0) bc0 = hit_rawbcid[0] % 256;

    for (int ih=0; ih<nHits; ih++){
      int bcid = hit_rawbcid[ih] - bc0;
      if (bcid < -128) bcid += 256;
      if (bcid > 128) bcid -= 256;
      
      int connector = hit_channel[ih] / 24;
      int layer = -1;
      int strip = -1;
      if (connector < cabling_map.size() && cabling_map[connector] >= 0){
        layer = (hit_channel[ih] % 24) / 8;
        strip = 8 * cabling_map[connector] + hit_channel[ih] % 8;
      }
      hit_layer[ih] = layer;
      hit_strip[ih] = strip;
      hit_bcid[ih] = bcid;
    }

    std::vector<int> hit_used1(nHits, 0);
    std::vector<int> hit_used2(nHits, 0);

    proc_layer.reserve(nHits);
    proc_strip.reserve(nHits);
    proc_time1.reserve(nHits);
    proc_time2.reserve(nHits);
    proc_tot1.reserve(nHits);
    proc_tot2.reserve(nHits);
    proc_deltat.reserve(nHits);
    proc_dt_trig1.reserve(nHits);
    proc_dt_trig2.reserve(nHits);

    // FAST LOOKUP: Group hits by channel
    std::vector<std::vector<int>> hits_by_ch(256);
    for (int ih = 0; ih < nHits; ih++) {
        if (hit_channel[ih] >= 0 && hit_channel[ih] < 256) {
            hits_by_ch[hit_channel[ih]].push_back(ih);
        }
    }

    for (int ih=0; ih<nHits; ih++){
      if (hit_channel[ih]==trig_channel && hit_rise[ih]==1 && trig_time==-1) {
        if (hit_time2[ih]!=0) {
          trig_time = (hit_bcid[ih]%256)*30+hit_time2[ih]-1;
        } else {
          cout<<"WARNING: trig time empty."<<std::endl;
        }
      }
    }

    for (int ih=0; ih<nHits; ih++){

      if (hit_rise[ih]==1 && hit_channel[ih]!=trig_channel && hit_time1[ih]>0 &&
          (hit_bcid[ih]%256)*30+hit_time1[ih]-1-trig_time<dt_max &&
          (hit_bcid[ih]%256)*30+hit_time1[ih]-1-trig_time>dt_min){
        eta1_hit_flag[hit_layer[ih]] = true;
      }
      if (hit_rise[ih]==1 && hit_channel[ih]!=trig_channel && hit_time2[ih]>0 &&
          (hit_bcid[ih]%256)*30+hit_time2[ih]-1-trig_time<dt_max &&
          (hit_bcid[ih]%256)*30+hit_time2[ih]-1-trig_time>dt_min){
        eta2_hit_flag[hit_layer[ih]] = true;
      }

      int connector = hit_channel[ih]/24;
      if (hit_rise[ih]==1 && hit_channel[ih]!=trig_channel && hit_layer[ih]>=0){
        int strip = hit_strip[ih];
        int layer = hit_layer[ih];
        int time1=-1, time2=-1;
        int tot1=9999, tot2=9999;
        int ihrise1=-1, ihrise2=-1;
        int deltat= 9999, dt_trig1 = 9999, dt_trig2 = 9999;
        int ch = hit_channel[ih];
        
        for (int ii : hits_by_ch[ch]) {
          if (ii < ih) continue;
          if (time1 == -1 && hit_rise[ii] == 1 && hit_used1[ii] == 0 && hit_time1[ii] != 0) {
            time1 = (hit_bcid[ii] % 256) * 30 + (hit_time1[ii] - 1);
            ihrise1 = ii;
            hit_used1[ii] = 1;
            break;
          }
        }

        if (time1 != -1) {
          for (int ii : hits_by_ch[ch]) {
            if (ii <= ihrise1) continue;
            if (tot1 == 9999 && hit_rise[ii] == 0 && hit_used1[ii] == 0 && hit_time1[ii] != 0) {
              int tmp_tot = (hit_bcid[ii] % 256) * 30 + (hit_time1[ii] - 1) - time1;
              if (tmp_tot > 0) {
                tot1 = tmp_tot;
                hit_used1[ii] = 1;
                break;
              }
            }
          }
        }
        
        if (tot1 != 9999) {
          if (tot1 < -256 * 30 / 2) tot1 += 256 * 30;
          if (tot1 > 256 * 30 / 2) tot1 -= 256 * 30;
        }

        for (int ii : hits_by_ch[ch]) {
          if (ii < ih) continue;
          if (time2 == -1 && hit_rise[ii] == 1 && hit_used2[ii] == 0 && hit_time2[ii] != 0) {
            time2 = (hit_bcid[ii] % 256) * 30 + (hit_time2[ii] - 1);
            ihrise2 = ii;
            hit_used2[ii] = 1;
            break;
          }
        }
      
        if (time2 > 0) {
          for (int ii : hits_by_ch[ch]) {
            if (ii <= ihrise2) continue;
            if (tot2 == 9999 && hit_rise[ii] == 0 && hit_used2[ii] == 0 && hit_time2[ii] != 0) {
              int tmp_tot = (hit_bcid[ii] % 256) * 30 + (hit_time2[ii] - 1) - time2;
              if (tmp_tot > 0) {
                tot2 = tmp_tot;
                hit_used2[ii] = 1;
                break;
              }
            }
          }
        }
        
        if (tot2 != 9999) {
          if (tot2 < -256 * 30 / 2) tot2 += 256 * 30;
          if (tot2 > 256 * 30 / 2) tot2 -= 256 * 30;
        }

        if (time1>-1 && time2>-1) deltat = time2 - time1;
        if (trig_time>0 && time1>-1) dt_trig1 = time1 - trig_time;
        if (trig_time>0 && time2>-1) dt_trig2 = time2 - trig_time;
    
        if (time1>0||time2>0){
          proc_layer.push_back(layer);
          proc_strip.push_back(strip);
          proc_time1.push_back(time1);
          proc_time2.push_back(time2);
          proc_tot1.push_back(tot1);
          proc_tot2.push_back(tot2);
          proc_deltat.push_back(deltat);
          proc_dt_trig1.push_back(dt_trig1);
          proc_dt_trig2.push_back(dt_trig2);
          nProc++;
        }
      }
    }

    // FAST ARRAYS: Using memset
    int hit[3][49];
    int hitmult[3][49];
    memset(hit, -1, sizeof(hit));
    memset(hitmult, 0, sizeof(hitmult));
      
    for (int iprc=0; iprc<nProc; iprc++){
      hitmult[proc_layer[iprc]][proc_strip[iprc]]++;
      if (hit[proc_layer[iprc]][proc_strip[iprc]]<0){
        hit[proc_layer[iprc]][proc_strip[iprc]]=iprc;
      }
    }
    
    vector<int> vclus_w_l[3];
    vector<int> vclus_indx_l[3];
    vector<float> vclus_center_l[3];
    vector<int> vclus_mult_l[3];
      
    for (int ly=0; ly<3; ly++){
      int cl_start=-1;
      for (int jj=0; jj<49; jj++){
        if (cl_start==-1){
          if (hit[ly][jj]>=0) cl_start=jj;
        } else {
          if (hit[ly][jj]<0) {
            vclus_w_l[ly].push_back(jj-cl_start);
            vclus_center_l[ly].push_back(0.5*(float)(jj+cl_start-1.));
        
            int strip=(jj+cl_start-1)/2;
            int indx=hit[ly][strip];
            if (!(proc_time1[indx]>0&&proc_time2[indx]>0)) {
              if (strip+1<jj){
                if (proc_time1[hit[ly][strip+1]]>0&&proc_time2[hit[ly][strip+1]]>0){
                  indx=hit[ly][strip+1];
                } else if (strip-1>=cl_start) {
                  if (proc_time1[hit[ly][strip-1]]>0&&proc_time2[hit[ly][strip-1]]>0){
                    indx=hit[ly][strip+1];
                  }
                }
              }
            }
            vclus_indx_l[ly].push_back(indx);
            int mult=0; for (int is=cl_start; is<=jj; is++) mult+=hitmult[ly][is];
            vclus_mult_l[ly].push_back(mult);
            cl_start=-1;
          }
        }
      }
    }

    for (int ic=0; ic<vclus_indx_l[0].size(); ic++){
      clus_nstrips_l0.push_back(vclus_w_l[0][ic]);
      clus_center_l0.push_back(vclus_center_l[0][ic]);
      clus_mult_l0.push_back(vclus_mult_l[0][ic]);
      int indx=vclus_indx_l[0][ic];
      clus_strip_l0.push_back(proc_strip[indx]);
      clus_tot1_l0.push_back(proc_tot1[indx]);
      clus_tot2_l0.push_back(proc_tot2[indx]);
      clus_time1_l0.push_back(proc_time1[indx]);
      clus_time2_l0.push_back(proc_time2[indx]);
      clus_dt_l0.push_back(proc_deltat[indx]);
      nClus_l0++;
    }
    for (int ic=0; ic<vclus_indx_l[1].size(); ic++){
      clus_nstrips_l1.push_back(vclus_w_l[1][ic]);
      clus_center_l1.push_back(vclus_center_l[1][ic]);
      clus_mult_l1.push_back(vclus_mult_l[1][ic]);
      int indx=vclus_indx_l[1][ic];
      clus_strip_l1.push_back(proc_strip[indx]);
      clus_tot1_l1.push_back(proc_tot1[indx]);
      clus_tot2_l1.push_back(proc_tot2[indx]);
      clus_time1_l1.push_back(proc_time1[indx]);
      clus_time2_l1.push_back(proc_time2[indx]);
      clus_dt_l1.push_back(proc_deltat[indx]);
      nClus_l1++;
    }
    for (int ic=0; ic<vclus_indx_l[2].size(); ic++){
      clus_nstrips_l2.push_back(vclus_w_l[2][ic]);
      clus_center_l2.push_back(vclus_center_l[2][ic]);
      clus_mult_l2.push_back(vclus_mult_l[2][ic]);
      int indx=vclus_indx_l[2][ic];
      clus_strip_l2.push_back(proc_strip[indx]);
      clus_tot1_l2.push_back(proc_tot1[indx]);
      clus_tot2_l2.push_back(proc_tot2[indx]);
      clus_time1_l2.push_back(proc_time1[indx]);
      clus_time2_l2.push_back(proc_time2[indx]);
      clus_dt_l2.push_back(proc_deltat[indx]);
      nClus_l2++;
    }

    for (int ip=0; ip<nProc; ip++){
      if (proc_time1[ip]>-1) dt_trig_vs_strip_eta1->Fill(proc_layer[ip]*48+proc_strip[ip],(proc_time1[ip]-trig_time));
      if (proc_time2[ip]>-1) dt_trig_vs_strip_eta2->Fill(proc_layer[ip]*48+proc_strip[ip],(proc_time2[ip]-trig_time));

      if (proc_layer[ip]==0 && proc_time1[ip]>-1) strips_eta1_layer0->Fill(proc_strip[ip]);
      if (proc_layer[ip]==1 && proc_time1[ip]>-1) strips_eta1_layer1->Fill(proc_strip[ip]);
      if (proc_layer[ip]==2 && proc_time1[ip]>-1) strips_eta1_layer2->Fill(proc_strip[ip]);
      if (proc_layer[ip]==0 && proc_time2[ip]>-1) strips_eta2_layer0->Fill(proc_strip[ip]);
      if (proc_layer[ip]==1 && proc_time2[ip]>-1) strips_eta2_layer1->Fill(proc_strip[ip]);
      if (proc_layer[ip]==2 && proc_time2[ip]>-1) strips_eta2_layer2->Fill(proc_strip[ip]);

      if (proc_layer[ip]==0 && proc_time1[ip]>-1 && proc_tot1[ip]!=9999) {
        tot_eta1_ly0_1D->Fill(proc_tot1[ip]);
        tot_eta1_ly0_2D->Fill(proc_strip[ip],proc_tot1[ip]);
      }
      if (proc_layer[ip]==1 && proc_time1[ip]>-1 && proc_tot1[ip]!=9999) {
        tot_eta1_ly1_1D->Fill(proc_tot1[ip]);
        tot_eta1_ly1_2D->Fill(proc_strip[ip],proc_tot1[ip]);
      }
      if (proc_layer[ip]==2 && proc_time1[ip]>-1 && proc_tot1[ip]!=9999) {
        tot_eta1_ly2_1D->Fill(proc_tot1[ip]);
        tot_eta1_ly2_2D->Fill(proc_strip[ip],proc_tot1[ip]);
      }
      if (proc_layer[ip]==0 && proc_time2[ip]>-1 && proc_tot2[ip]!=9999) {
        tot_eta2_ly0_1D->Fill(proc_tot2[ip]);
        tot_eta2_ly0_2D->Fill(proc_strip[ip],proc_tot2[ip]);
      }
      if (proc_layer[ip]==1 && proc_time2[ip]>-1 && proc_tot2[ip]!=9999) {
        tot_eta2_ly1_1D->Fill(proc_tot2[ip]);
        tot_eta2_ly1_2D->Fill(proc_strip[ip],proc_tot2[ip]);
      }
      if (proc_layer[ip]==2 && proc_time2[ip]>-1 && proc_tot2[ip]!=9999) {
        tot_eta2_ly2_1D->Fill(proc_tot2[ip]);
        tot_eta2_ly2_2D->Fill(proc_strip[ip],proc_tot2[ip]);
      }
    }
      
    tree->Fill();
    nHits=0; nProc=0; nClus_l0=0; nClus_l1=0; nClus_l2=0;
    hit_channel.clear(); hit_clk.clear(); hit_rawbcid.clear();
    hit_bcid.clear(); hit_rawbcout.clear(); hit_bcout.clear();
    hit_time1.clear(); hit_time2.clear(); hit_rise.clear();
    hit_layer.clear(); hit_strip.clear(); proc_layer.clear();
    proc_strip.clear(); proc_time1.clear(); proc_time2.clear();
    proc_tot1.clear(); proc_tot2.clear(); proc_deltat.clear();
    proc_dt_trig1.clear(); proc_dt_trig2.clear();
    clus_strip_l0.clear(); clus_strip_l1.clear();  clus_strip_l2.clear();
    clus_time1_l0.clear(); clus_time1_l1.clear();  clus_time1_l2.clear();
    clus_time2_l0.clear(); clus_time2_l1.clear();  clus_time2_l2.clear();
    clus_tot1_l0.clear();  clus_tot1_l1.clear();   clus_tot1_l2.clear();
    clus_tot2_l0.clear();  clus_tot2_l1.clear();   clus_tot2_l2.clear();
    clus_dt_l0.clear();    clus_dt_l1.clear();     clus_dt_l2.clear();
    clus_nstrips_l0.clear();  clus_nstrips_l1.clear(); clus_nstrips_l2.clear();
    clus_center_l0.clear();  clus_center_l1.clear(); clus_center_l2.clear();
    clus_mult_l0.clear(); clus_mult_l1.clear(); clus_mult_l2.clear();

    if (trig_time>0 && eta2_hit_flag[1] && eta2_hit_flag[2]) {
      evts_rpc_trig_denom[0]++;
      if (eta1_hit_flag[0]) evts_eta1_rpc_trig[0]++;
      if (eta2_hit_flag[0]) evts_eta2_rpc_trig[0]++;
      if (eta1_hit_flag[0] || eta2_hit_flag[0]) evts_OR_rpc_trig[0]++;
      if (eta1_hit_flag[0] && eta2_hit_flag[0]) evts_AND_rpc_trig[0]++;
    }

    if (trig_time>0 && eta2_hit_flag[0] && eta2_hit_flag[2]) {
      evts_rpc_trig_denom[1]++;
      if (eta1_hit_flag[1]) evts_eta1_rpc_trig[1]++;
      if (eta2_hit_flag[1]) evts_eta2_rpc_trig[1]++;
      if (eta1_hit_flag[1] || eta2_hit_flag[1]) evts_OR_rpc_trig[1]++;
      if (eta1_hit_flag[1] && eta2_hit_flag[1]) evts_AND_rpc_trig[1]++;
    }

    if (trig_time>0 && eta2_hit_flag[0] && eta2_hit_flag[1]) {
      evts_rpc_trig_denom[2]++;
      if (eta1_hit_flag[2]) evts_eta1_rpc_trig[2]++;
      if (eta2_hit_flag[2]) evts_eta2_rpc_trig[2]++;
      if (eta1_hit_flag[2] || eta2_hit_flag[2]) evts_OR_rpc_trig[2]++;
      if (eta1_hit_flag[2] && eta2_hit_flag[2]) evts_AND_rpc_trig[2]++;
    }

    if (trig_time>0) {
      trig_evts += 1;
      for (int l = 0; l < 3; l++) {
        if (eta1_hit_flag[l]) evts_eta1[l]++;
        if (eta2_hit_flag[l]) evts_eta2[l]++;
        if (eta1_hit_flag[l] || eta2_hit_flag[l]) evts_or[l]++;
        if (eta1_hit_flag[l] && eta2_hit_flag[l]) evts_and[l]++;
      }
    }

    trig_time=-1;
  }

  double eff_eta1[3] = {0,0,0}, eff_eta2[3] = {0,0,0}, eff_or[3] = {0,0,0}, eff_and[3] = {0,0,0};
  double eff_err_eta1[3] = {0,0,0}, eff_err_eta2[3] = {0,0,0}, eff_err_or[3] = {0,0,0}, eff_err_and[3] = {0,0,0};
  if (trig_evts > 0.0) {
    for (int l = 0; l < 3; l++) {
      eff_eta1[l] = static_cast<double>(evts_eta1[l]) / trig_evts;
      eff_eta2[l] = static_cast<double>(evts_eta2[l]) / trig_evts;
      eff_or[l]  = static_cast<double>(evts_or[l])  / trig_evts;
      eff_and[l]  = static_cast<double>(evts_and[l])  / trig_evts;

      eff_err_eta1[l] = sqrt(eff_eta1[l]*(1.0-eff_eta1[l])/trig_evts);
      eff_err_eta2[l] = sqrt(eff_eta2[l]*(1.0-eff_eta2[l])/trig_evts);
      eff_err_or[l] = sqrt(eff_or[l]*(1.0-eff_or[l])/trig_evts);
      eff_err_and[l] = sqrt(eff_and[l]*(1.0-eff_and[l])/trig_evts);
    }
  }

  std::cout<<std::endl;
  std::cout << "Applied cuts on timing window from trigger (ticks): [" << dt_min << ", " << dt_max << "]" << std::endl;
  std::cout << "***************** Efficiency with external trigger only *****************" << std::endl;
  std::cout << "* eta1:\t ly0 = " << std::fixed << std::setprecision(3) << eff_eta1[0] << " +/- " << std::fixed << std::setprecision(3) << eff_err_eta1[0]
  << "\t ly1 = " << std::fixed << std::setprecision(3) << eff_eta1[1] << " +/- " << std::fixed << std::setprecision(3) << eff_err_eta1[1]
  << "\t ly2 = " << std::fixed << std::setprecision(3) << eff_eta1[2] << " +/- " << std::fixed << std::setprecision(3) << eff_err_eta1[2] << std::endl;
  std::cout << "* eta2:\t ly0 = " << std::fixed << std::setprecision(3) << eff_eta2[0] << " +/- " << std::fixed << std::setprecision(3) << eff_err_eta2[0]
  << "\t ly1 = " << std::fixed << std::setprecision(3) << eff_eta2[1] << " +/- " << std::fixed << std::setprecision(3) << eff_err_eta2[1]
  << "\t ly2 = " << std::fixed << std::setprecision(3) << eff_eta2[2] << " +/- " << std::fixed << std::setprecision(3) << eff_err_eta2[2] << std::endl;
  std::cout << "* OR:\t ly0 = " << std::fixed << std::setprecision(3) << eff_or[0] << " +/- " << std::fixed << std::setprecision(3) << eff_err_or[0]
  << "\t ly1 = " << std::fixed << std::setprecision(3) << eff_or[1] << " +/- " << std::fixed << std::setprecision(3) << eff_err_or[1]
  << "\t ly2 = " << std::fixed << std::setprecision(3) << eff_or[2] << " +/- " << std::fixed << std::setprecision(3) << eff_err_or[2] << std::endl;
  std::cout << "* AND:\t ly0 = " << std::fixed << std::setprecision(3) << eff_and[0] << " +/- " << std::fixed << std::setprecision(3) << eff_err_and[0]
  << "\t ly1 = " << std::fixed << std::setprecision(3) << eff_and[1] << " +/- " << std::fixed << std::setprecision(3) << eff_err_and[1]
  << "\t ly2 = " << std::fixed << std::setprecision(3) << eff_and[2] << " +/- " << std::fixed << std::setprecision(3) << eff_err_and[2] << std::endl;
  std::cout << "**********************************************************************" << std::endl;
  std::cout << "Triggered events: " << trig_evts << std::endl;
  std::cout << std::endl;

  for (int l = 0; l < 3; l++) {
    TH1F *hEff1 = new TH1F(Form("efficiency_ext_trigger_layer%d", l), Form("Efficiency layer %d (with external trigger only);type;efficiency", l), 4, 0.5, 4.5);
    hEff1->GetXaxis()->SetBinLabel(1,"eta1");
    hEff1->GetXaxis()->SetBinLabel(2,"eta2");
    hEff1->GetXaxis()->SetBinLabel(3,"OR");
    hEff1->GetXaxis()->SetBinLabel(4,"AND");
    hEff1->SetBinContent(1, eff_eta1[l]);
    hEff1->SetBinError(1, eff_err_eta1[l]);
    hEff1->SetBinContent(2, eff_eta2[l]);
    hEff1->SetBinError(2, eff_err_eta2[l]);
    hEff1->SetBinContent(3, eff_or[l]);
    hEff1->SetBinError(3, eff_err_or[l]);
    hEff1->SetBinContent(4, eff_and[l]);
    hEff1->SetBinError(4, eff_err_and[l]);
  }
  
  o->Write();
  o->Close();
  myFile->Close();
}
