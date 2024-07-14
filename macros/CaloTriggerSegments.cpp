#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <iostream>
#include <vector>
#include <cmath>
#include <map>
#include <bitset>
#include "TriggerDefs.h"
#include <phool/PHCompositeNode.h>
#include <phool/getClass.h>
#include <calobase/TowerInfoContainer.h>
#include <calobase/TowerInfoDefs.h>
#include <calobase/TowerInfo.h>
#include <unordered_map>

R__LOAD_LIBRARY(libcalotrigger.so)
R__LOAD_LIBRARY(libtriggervalid.so)

void CaloTriggerSegments(const char* filename, const char* output_filename, Long64_t nEvents = 0) {
    TFile *file = TFile::Open(filename);
    if (!file || file->IsZombie()) {
        std::cerr << "Error: Cannot open file " << filename << std::endl;
        return;
    }
    std::cout << "ROOT file opened successfully: " << filename << std::endl;

    TTree *tree = (TTree*)file->Get("ttree");
    if (!tree) {
        std::cerr << "Error: Cannot find TTree 'ttree' in the file." << std::endl;
        file->Close();
        return;
    }

    
    ULong64_t b_gl1_scaledvec, gl1_livevec, b_gl1_rawvec;;
    ULong64_t b_gl1_clock;  // Variable to hold the GL1 clock
    // Define pointers to vectors to hold EMCal data
    std::vector<bool>* b_emcal_good= 0;
    std::vector<float>* b_emcal_energy= 0;
    std::vector<float>* b_emcal_time= 0;
    std::vector<int>* b_emcal_etabin = 0;
    std::vector<int>* b_emcal_phibin = 0;
      
    // Define pointers to vectors to hold HCal inner data
    std::vector<short>* b_hcalin_good= nullptr;
    std::vector<float>* b_hcalin_energy= nullptr;
    std::vector<float>* b_hcalin_time= nullptr;
    std::vector<float>* b_hcalin_etabin= nullptr;
    std::vector<float>* b_hcalin_phibin= nullptr;

    // Define pointers to vectors to hold HCal outer data
    std::vector<short>* b_hcalout_good= nullptr;
    std::vector<float>* b_hcalout_energy= nullptr;
    std::vector<float>* b_hcalout_time= nullptr;
    std::vector<float>* b_hcalout_etabin= nullptr;
    std::vector<float>* b_hcalout_phibin= nullptr;

    // Define variables to hold cluster data
    int b_cluster_n;
    std::vector<float>* b_cluster_prob = nullptr;
    std::vector<float>* b_cluster_chi2 = nullptr;
    std::vector<float>* b_cluster_ecore = nullptr;
    std::vector<float>* b_cluster_pt = nullptr;
    std::vector<float>* b_cluster_phi = nullptr;
    std::vector<float>* b_cluster_eta = nullptr;
    std::vector<float>* b_cluster_iso = nullptr;
    Float_t mbd_vertex_z;

    // Define arrays to hold trigger sums for EMCal LL1 and EMCal
    unsigned int b_trigger_sum_emcal_ll1[384];
    unsigned int b_trigger_sumkey_emcal_ll1[384];
    unsigned int b_trigger_sum_emcal[6144];
    unsigned int b_trigger_sumkey_emcal[6144];
    

    tree->SetBranchAddress("gl1_scaled", gl1_scaled);
    tree->SetBranchAddress("gl1_live", gl1_live);
    tree->SetBranchAddress("gl1_livevec", &gl1_livevec); //equivalent to TriggerVector
    
    /*
     Vectors of length 24576 for 96 * 216 towers
     */
    tree->SetBranchAddress("gl1_scaled", gl1_scaled);
    tree->SetBranchAddress("gl1_live", gl1_live);
    tree->SetBranchAddress("gl1_livevec", &gl1_livevec);
    tree->SetBranchAddress("trigger_sum_emcal_ll1", &b_trigger_sum_emcal_ll1);
    tree->SetBranchAddress("trigger_sumkey_emcal_ll1", &b_trigger_sumkey_emcal_ll1);
    tree->SetBranchAddress("trigger_sum_emcal", &b_trigger_sum_emcal);
    tree->SetBranchAddress("trigger_sumkey_emcal", &b_trigger_sumkey_emcal);
    tree->SetBranchAddress("gl1_clock", &b_gl1_clock);
    tree->SetBranchAddress("gl1_scaledvec", &b_gl1_scaledvec);
    tree->SetBranchAddress("gl1_rawvec", &b_gl1_rawvec);

    /*
     In CaloEmulatorTreeMaker.cc:
     
     _towers = findNode::getClass<TowerInfoContainer>(topNode, "TOWERINFO_CALIB_CEMC");
     if (_towers) {
     size = _towers->size(); //online towers should be the same!
     for (int channel = 0; channel < size;channel++) {
         _tower = _towers->get_tower_at_channel(channel);
         float energy = _tower->get_energy();
         float time = _tower->get_time();
         unsigned int towerkey = _towers->encode_key(channel);
         int ieta = _towers->getTowerEtaBin(towerkey);
         int iphi = _towers->getTowerPhiBin(towerkey);
         short good = (_tower->get_isGood() ? 1:0);
         b_emcal_good.push_back(good);
         b_emcal_energy.push_back(energy);
         b_emcal_time.push_back(time);
         b_emcal_etabin.push_back(ieta);
         b_emcal_phibin.push_back(iphi);
       }
     }
     */
    tree->SetBranchAddress("emcal_good", &b_emcal_good);
    tree->SetBranchAddress("emcal_energy", &b_emcal_energy);
    tree->SetBranchAddress("emcal_time", &b_emcal_time);
    tree->SetBranchAddress("emcal_phibin", &b_emcal_phibin);
    tree->SetBranchAddress("emcal_etabin", &b_emcal_etabin);
    tree->SetBranchAddress("hcalin_good", &b_hcalin_good);
    tree->SetBranchAddress("hcalin_energy", &b_hcalin_energy);
    tree->SetBranchAddress("hcalin_time", &b_hcalin_time);
    tree->SetBranchAddress("hcalin_phibin", &b_hcalin_phibin);
    tree->SetBranchAddress("hcalin_etabin", &b_hcalin_etabin);
    tree->SetBranchAddress("hcalout_good", &b_hcalout_good);
    tree->SetBranchAddress("hcalout_energy", &b_hcalout_energy);
    tree->SetBranchAddress("hcalout_time", &b_hcalout_time);
    tree->SetBranchAddress("hcalout_phibin", &b_hcalout_phibin);
    tree->SetBranchAddress("hcalout_etabin", &b_hcalout_etabin);
    tree->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z);
    tree->SetBranchAddress("cluster_prob", &b_cluster_prob);
    tree->SetBranchAddress("cluster_chi2", &b_cluster_chi2);
    tree->SetBranchAddress("cluster_ecore", &b_cluster_ecore);
    tree->SetBranchAddress("cluster_pt", &b_cluster_pt);
    tree->SetBranchAddress("cluster_phi", &b_cluster_phi);
    tree->SetBranchAddress("cluster_eta", &b_cluster_eta);
    tree->SetBranchAddress("cluster_iso", &b_cluster_iso);
    tree->SetBranchAddress("cluster_n", &b_cluster_n);

    // Print out the branch addresses and their statuses
    std::cout << "\033[1mBranch addresses set successfully:\033[0m" << std::endl;

    TH1D *h_emcal_energy_gl1[64];
    TH1D *h_cluster_energy_gl1[64];

    for (int i = 0; i < 64; ++i) {
        h_emcal_energy_gl1[i] = new TH1D(Form("h_emcal_energy_gl1_%d", i), "; Energy;", 40, 0,  20);
        h_cluster_energy_gl1[i] = new TH1D(Form("h_cluster_energy_gl1_%d", i), "; Energy;", 40, 0,  20);
    }
    
    Long64_t nentries = std::min(tree->GetEntries(), (nEvents > 0 ? nEvents : tree->GetEntries()));
    std::cout << "Number of entries to process: " << nentries << std::endl;

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        
        std::cout << "Processing entry " << i << std::endl;
        
        std::bitset<64> bits(b_gl1_scaledvec);
        std::cout << "Processing entry " << i << ", gl1_scaledvec (bits): " << bits.to_string() << std::endl;

        std::vector<int> trig_bits{};

        // Loop over each bit position from 0 to 63
        for (unsigned int bit = 0; bit < 64; bit++) {
            /* 
               Check if the bit at the specified position 'bit' in the variable 'b_gl1_scaledvec' is set to 1.
               This is done by shifting 'b_gl1_scaledvec' right by 'bit' positions, which moves the bit of interest to the
               least significant bit position. 
             
               Then, it is bitwise ANDed with 0x1U to isolate this bit. If the result is 1,
               then the bit at position 'bit' was originally set to 1.
             */
            if (((b_gl1_scaledvec >> bit) & 0x1U) == 0x1U) {
                /* If the bit is set to 1, add its position to the vector 'trig_bits',
                   which stores the positions of all set bits found. */
                trig_bits.push_back(bit);
            }
        }
        std::cout << "\033[31mScaled Trigger Bits: [ ";
        for (auto bit : trig_bits) {
            std::cout << bit << " ";
        }
        std::cout << "]\033[0m" << std::endl;

        
        // Loop over each row in the map (35 rows) -- each sector in phi
        for (int j = 0; j < 35; j++) {
            // Loop over each column in the map (12 columns) -- each ib board across eta
            for (int k = 0; k < 12; k++) {
                if (j < 32) {
                    energymap[k][j] = 0.0;
                }
                energymap_emcal[k][j] = 0.0;
            }
        }
        std::cout << "Energy and good maps initialized for entry " << i << std::endl;
        
        //b_emcal_energy is the size of all towers in the EMCal -- 256*96 = 24576
        for (int ie = 0; ie < b_emcal_energy->size(); ie++) {
            // Calculate the eta bin and phi bin for the current entry
            // Each IB board is 8x8 towers, hence we divide by 8 to get the IB bin index
            int ebin = b_emcal_etabin->at(ie) / 8;
            int pbin = b_emcal_phibin->at(ie) / 8;
            
            if (!b_emcal_good->at(ie)) continue;
            energymap[ebin][pbin] += b_emcal_energy->at(ie);
        }

        float max_energy = 0.0;
        float max_energy_clus = 0.0;
        int ebin = 0;
        int pbin = 0;

        // Loop over all clusters to find the maximum cluster energy (ECore)
        for (int j = 0; j < b_cluster_n; j++) {
            if (b_cluster_ecore->at(j) > max_energy_clus) {
                max_energy_clus = b_cluster_ecore->at(j); // Update the maximum cluster energy if current cluster has higher energy
            }
        }
        std::cout << "Maximum cluster energy (ECore): " << max_energy_clus << " for entry " << i << std::endl;

        //emcal has sectors 0 through 31 on north side and 32 through 63 on south
        //12 IB boards across eta
        for (int j = 0; j < 32; j++) { //go through all sectors in phi
            for (int k = 0; k < 12; k++) { //loop over each IB board through eta
                if (energymap[k][j] > max_energy) {
                    max_energy = energymap[k][j]; // Update maximum EMCal energy
                    ebin = k; // Update EMCal eta bin index
                    pbin = j; // Update EMCal phi bin index
                }
            }
        }
        h_emcal_energy->Fill(max_energy); // Fill EMCal energy histogram
        h_cluster_energy->Fill(max_energy_clus); // Fill cluster energy histogram


        std::cout << "Filled histograms for entry " << i << std::endl;

        // Loop over scaled trigger bits and fill histograms
        for (auto &b : trig_bits) {
            h_cluster_energy_gl1[b]->Fill(max_energy_clus); // Fill GL1 cluster energy histogram for each trigger bit
            h_emcal_energy_gl1[b]->Fill(max_energy); // Fill GL1 EMCal energy histogram for each trigger bit -- I think this is 8x8
        }
    }

    TFile *outfile = new TFile(output_filename, "RECREATE");
    TTree *outputTree = new TTree("SummaryTree", "Summary data including accumulated counts");

    outputTree->Fill();

    for (int i = 0; i < 64; ++i) {
        h_emcal_energy_gl1[i]->Write();
        h_cluster_energy_gl1[i]->Write();
    }
    outputTree->Write();
    outfile->Close();
    file->Close();
    std::cout << "Histograms and accumulated counts written to file and files closed." << std::endl;
}
