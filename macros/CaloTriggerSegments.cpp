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

/*
make list file: ./createListFile /sphenix/u/patsfan753/scratch/analysis/calotriggeremulator/output/44630 44630
 */

// Define a key for the eta-phi bin pair
struct EtaPhiKey {
    int eta;
    int phi;

    bool operator==(const EtaPhiKey &other) const {
        return eta == other.eta && phi == other.phi;
    }
};

// Define a hash function for the EtaPhiKey
namespace std {
    template <>
    struct hash<EtaPhiKey> {
        std::size_t operator()(const EtaPhiKey &k) const {
            return ((std::hash<int>()(k.eta) ^ (std::hash<int>()(k.phi) << 1)) >> 1);
        }
    };
}

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

    
    ULong64_t gl1_scaled[64], gl1_live[64], gl1_livevec, b_gl1_rawvec;;
    ULong64_t b_gl1_scaledvec; // Variable to hold the GL1 scaled vector
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

    
    TH1F *histograms_ecore[64];
    TH1F *histograms_energy[64];
    TH1F *histograms_ecore_liveCounter[64];
    TH1F *histograms_energy_liveCounter[64];
    std::vector<float> accumulated_live(64, 0);
    std::vector<float> accumulated_scaled(64, 0);
    
    
    float energymap[12][32] = {0};
    float energymap_hcalin[12][35] = {0};
    float energymap_hcalout[12][35] = {0};
    float energymap_emcal[12][35] = {0};
    float energymap_jet[9][32] = {0};
    float energymap_jet_emcal[9][32] = {0};
    float energymap_jet_hcalin[9][32] = {0};
    float energymap_jet_hcalout[9][32] = {0};
    float energymap_extend[12][35] = {0};

    int goodmap[12][32] = {0};
    int goodmap_hcalin[12][35] = {0};
    int goodmap_hcalout[12][35] = {0};
    int goodmap_emcal[12][35] = {0};
    int goodmap_jet[9][32] = {0};
    int goodmap_jet_emcal[9][32] = {0};
    int goodmap_jet_hcalin[9][32] = {0};
    int goodmap_jet_hcalout[9][32] = {0};
    int goodmap_extend[12][35] = {0};

    int n_primitives = 384;
    int n_sums = 6144;
    int energy_bins = 100;
    float energy_max = 30.0;
    int turnon_bins = 40;
    float turnon_max = 10.0;
    int lut_bins = 256;

    
    // clock v sum
    TH1D *h_triggers = new TH1D("h_triggers", "", 64, -0.5, 63.5);
    TH1D *h_mbd_triggers = new TH1D("h_mbd_triggers", "", 64, -0.5, 63.5);
    TH1D *h_emcal_energy = new TH1D("h_emcal_energy", "; Energy;", 40, 0,  20);
    TH1D *h_cluster_energy = new TH1D("h_cluster_energy", "; Energy;", 40, 0,  20);
    TH1D *h_hcal_energy = new TH1D("h_hcal_energy", "; Energy;", 40, 0,  20);
    TH1D *h_jet_energy = new TH1D("h_hmcal_energy", "; Energy;", 50, 0,  50);
    TH1D *h_emcal_energy_ref = new TH1D("h_emcal_energy_ref", "; Energy;", 40, 0,20);
    TH1D *h_jet_energy_ref = new TH1D("h_jet_energy_ref", "; Energy;", 50, 0,  50);

    TH1D *h_cluster_phi = new TH1D("h_cluster_phi", "; #phi;", 64, -1*TMath::Pi(), TMath::Pi());
    TH1D *h_cluster_eta = new TH1D("h_cluster_eta", "; #eta;", 64, -1.2, 1.2);
    TH1D *h_cluster_iso = new TH1D("h_cluster_iso", "; Energy;", 100, -50, 50);
    

    TH1D *h_cluster_phi_gl1[64];
    TH1D *h_cluster_eta_gl1[64];
    TH1D *h_cluster_iso_gl1[64];

    TH1D *h_cluster_phi_chi2 = new TH1D("h_cluster_phi_chi2", "; #phi;", 64, -1*TMath::Pi(), TMath::Pi());
    TH1D *h_cluster_eta_chi2 = new TH1D("h_cluster_eta_chi2", "; #eta;", 64, -1.2, 1.2);
    TH1D *h_cluster_iso_chi2 = new TH1D("h_cluster_iso_chi2", "; Energy;", 100, -50, 50);

    TH1D *h_cluster_phi_gl1_chi2[64];
    TH1D *h_cluster_eta_gl1_chi2[64];
    TH1D *h_cluster_iso_gl1_chi2[64];

    TH1D *h_emcal_energy_gl1[64];
    TH1D *h_emcal_energy_raw_gl1[64];
    TH1D *h_cluster_energy_gl1[64];
    TH1D *h_cluster_energy_raw_gl1[64];
    TH1D *h_emcal_energy_scaled_gl1[64];
    TH1D *h_cluster_energy_scaled_gl1[64];
    TH1D *h_hcal_energy_gl1[64];
    TH1D *h_jet_emcal_energy_gl1[64];
    TH1D *h_jet_hcalin_energy_gl1[64];
    TH1D *h_jet_hcalout_energy_gl1[64];
    TProfile *h_jet_emcal_fraction_gl1[64];
    TProfile *h_jet_hcalin_fraction_gl1[64];
    TProfile *h_jet_hcalout_fraction_gl1[64];
    TH2F *h2_good = new TH2F("h2_good","; Good towers; energy", 65, -0.5, 64.5, 200, 0, 20);
    TH2F *h2_good_jet = new TH2F("h2_good_jet","; Good towers; energy", 64*16, -0.5, 64*16 + .5, 200, 0, 20);
    TH2F *h2_good_gl1[64];
    TH2F *h2_good_jet_gl1[64];
    TH1D *h_jet_energy_gl1[64];
    TH1D *h_jet_energy_scaled_gl1[64];
    TH1D *h_jet_energy_raw_gl1[64];
    TH1D *h_hcalout_imbalance[64];
    TH1D *h_emcal_imbalance[64];
    TProfile *h_emcal_spread = new TProfile("h","h", 96*256, -0.5, 96*256 - 0.5);

    for (int i = 0; i < 64; ++i) {
        histograms_ecore[i] = new TH1F(Form("trigger_%d_ecore_mycalculation", i), Form("Trigger %d Maximum Cluster ECore", i), 24, 0, 12);
        histograms_energy[i] = new TH1F(Form("trigger_%d_energy_mycalculation", i), Form("Trigger %d Maximum EMCAL Energy", i), 24, 0, 12);
        histograms_ecore_liveCounter[i] = new TH1F(Form("trigger_%d_ecore_liveCounter_mycalculation", i), Form("Trigger %d Maximum Cluster ECore, Live Counter", i), 24, 0, 12);
        histograms_energy_liveCounter[i] = new TH1F(Form("trigger_%d_energy_liveCounter_mycalculation_", i), Form("Trigger %d Maximum EMCAL Energy, Live Counter", i), 24, 0, 12);
        
        
        h_cluster_phi_gl1[i] = new TH1D(Form("h_cluster_phi_gl1_%d", i), "; #phi;", 64, -1*TMath::Pi(), TMath::Pi());
        h_cluster_eta_gl1[i] = new TH1D(Form("h_cluster_eta_gl1_%d", i), "; #eta;", 64, -1.2, 1.2);
        h_cluster_iso_gl1[i] = new TH1D(Form("h_cluster_iso_gl1_%d", i), "; Energy;", 100, -50, 50);
        
        h2_good_gl1[i] = new TH2F(Form("h2_good_gl1_%d", i),"; Good towers; energy", 65, -0.5, 64.5, 200, 0, 20);
        h2_good_jet_gl1[i] = new TH2F(Form("h2_good_jet_gl1_%d", i),"; Good towers; energy", 64*16, -0.5, 64*16 + .5, 40, 0, 20);
        h_jet_emcal_fraction_gl1[i] = new TProfile(Form("h_jet_emcal_fraction_gl1_%d", i), "; Energy; Fraction", 40, 0, 20);
        h_jet_hcalin_fraction_gl1[i] = new TProfile(Form("h_jet_hcalin_fraction_gl1_%d", i), "; Energy; Fraction", 40, 0, 20);
        h_jet_hcalout_fraction_gl1[i] = new TProfile(Form("h_jet_hcalout_fraction_gl1_%d", i), "; Energy; Fraction", 40, 0, 20);
        h_emcal_energy_gl1[i] = new TH1D(Form("h_emcal_energy_gl1_%d", i), "; Energy;", 40, 0,  20);
        h_emcal_energy_raw_gl1[i] = new TH1D(Form("h_emcal_energy_raw_gl1_%d", i), "; Energy;", 40, 0,  20);
        h_cluster_energy_gl1[i] = new TH1D(Form("h_cluster_energy_gl1_%d", i), "; Energy;", 40, 0,  20);
        h_cluster_energy_raw_gl1[i] = new TH1D(Form("h_cluster_energy_raw_gl1_%d", i), "; Energy;", 40, 0,  20);
        h_hcal_energy_gl1[i] = new TH1D(Form("h_hcal_energy_gl1_%d", i), "; Energy;", 40, 0,  20);
        
        h_hcalout_imbalance[i] = new TH1D(Form("h_hcalout_imbalance_%d", i), "; Energy;", 40, 0,  20);
        h_emcal_imbalance[i] = new TH1D(Form("h_emcal_imbalance_%d", i), "; Energy;", 40, 0,  20);
        h_jet_emcal_energy_gl1[i] = new TH1D(Form("h_jet_emcal_energy_gl1_%d", i), "; Energy;", 40, 0,  20);
        h_jet_hcalin_energy_gl1[i] = new TH1D(Form("h_jet_hcalin_energy_gl1_%d", i), "; Energy;", 40, 0,  20);
        h_jet_hcalout_energy_gl1[i] = new TH1D(Form("h_jet_hcalout_energy_gl1_%d", i), "; Energy;", 40, 0,  20);
        h_jet_energy_gl1[i] = new TH1D(Form("h_jet_energy_gl1_%d", i), "; Energy;", 50, 0,  50);
        h_jet_energy_raw_gl1[i] = new TH1D(Form("h_jet_energy_raw_gl1_%d", i), "; Energy;", 50, 0,  50);
    }

    std::vector<float> maxEcoreDansCode;
    std::vector<float> max8x8DansCalc;
    std::vector<float> maxEcoreMyCalc;
    std::vector<float> max8x8MyCalc;

    float maxEcoreMyCalcVal = 0.0;
    float max8x8MyCalcVal = 0.0;
    
    Long64_t nentries = std::min(tree->GetEntries(), (nEvents > 0 ? nEvents : tree->GetEntries()));
    std::cout << "Number of entries to process: " << nentries << std::endl;

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        
        std::cout << "Processing entry " << i << std::endl;
        
//        // Apply the z vertex cut
//        if (std::abs(mbd_vertex_z) >= 40) {
//            std::cout << "Skipping entry " << i << " due to z vertex cut: " << mbd_vertex_z << std::endl;
//            continue;
//        }
        
        std::bitset<64> bits(b_gl1_scaledvec);
        std::bitset<64> bits_live(gl1_livevec);
        
        std::cout << "Processing entry " << i << ", gl1_scaledvec (bits): " << bits.to_string() << std::endl;
        std::cout << "Processing entry " << i << ", trigger_vector (bits): " << bits_live.to_string() << std::endl;
        

        std::unordered_map<EtaPhiKey, float> energy_map;
        std::map<int, std::map<int, float>> emcal_energies;


        for (size_t k = 0; k < b_emcal_energy->size(); ++k) {
            if ((*b_emcal_good)[k]) {
                EtaPhiKey key = {static_cast<int>((*b_emcal_etabin)[k]), static_cast<int>((*b_emcal_phibin)[k])};
                energy_map[key] += (*b_emcal_energy)[k];
            }
        }
        

        // Loop over all trigger sums
        for (int idx = 0; idx < 6144; ++idx) {
            unsigned int sumk = b_trigger_sumkey_emcal[idx];
            uint16_t sum_phi = TriggerDefs::getSumPhiId(sumk) + 2 * TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
            uint16_t sum_eta = TriggerDefs::getSumEtaId(sumk);

            float energy_sum = 0.0;

            for (int itower = 0; itower < 64; ++itower) {
                int ieta = sum_eta * 8 + itower % 8;
                int iphi = sum_phi * 8 + itower / 8;

                EtaPhiKey key = {ieta, iphi};
                if (energy_map.find(key) != energy_map.end()) {
                    energy_sum += energy_map[key];
                }
            }

            emcal_energies[sum_eta][sum_phi] = energy_sum;
        }

        
        for (int j = 0; j < 64; ++j) {
            if (bits.test(j)) {
                accumulated_live[j] += gl1_live[j];
                accumulated_scaled[j] += gl1_scaled[j];

                if (b_cluster_ecore && !b_cluster_ecore->empty()) {
                    float maxEcore = *std::max_element(b_cluster_ecore->begin(), b_cluster_ecore->end());
                    histograms_ecore[j]->Fill(maxEcore);
                    maxEcoreMyCalc.push_back(maxEcore);
                    std::cout << "\033[1;34mMax ECore: " << maxEcore << " added to histogram.\033[0m" << std::endl;
                }

                float max8x8Energy = 0.0;
                for (const auto& entry : emcal_energies) {
                    for (const auto& sub_entry : entry.second) {
                        // Skip entries with emcal energy < 0.1
                        if (sub_entry.second < 0.1) {
                            continue;
                        }

                        if (sub_entry.second > max8x8Energy) {
                            max8x8Energy = sub_entry.second;
                        }
                    }
                }
                histograms_energy[j]->Fill(max8x8Energy);
                std::cout << "\033[1;31mMax 8x8 Energy: " << max8x8Energy << " added to histogram for trigger " << j << ".\033[0m" << std::endl;
                max8x8MyCalcVal = max8x8Energy;
            }

            if (bits_live.test(j)) {
                if (b_cluster_ecore && !b_cluster_ecore->empty()) {
                    float maxEcore = *std::max_element(b_cluster_ecore->begin(), b_cluster_ecore->end());
                    histograms_ecore_liveCounter[j]->Fill(maxEcore);
                }

                float max8x8Energy = 0.0;
                for (const auto& entry : emcal_energies) {
                    for (const auto& sub_entry : entry.second) {
                        if (sub_entry.second > max8x8Energy) {
                            max8x8Energy = sub_entry.second;
                        }
                    }
                }
                histograms_energy_liveCounter[j]->Fill(max8x8Energy);
            }
        }

        
        
        // Initialize a vector to store scaled trigger bits that are set (i.e., bits that are 1)
        std::vector<int> trig_bits{};

        // Initialize a vector to store raw trigger bits that are set (i.e., bits that are 1)
        std::vector<int> trig_bits_raw{};

        // Initialize a counter for raw MBD (Minimum Bias Detector) triggers
        int raw_mbd = 0;

        // Loop over each bit position from 0 to 63
        for (unsigned int bit = 0; bit < 64; bit++) {
            // Check if the current bit is set in the GL1 scaled vector
            // (i.e., the bit at position 'bit' is 1)
            if (((b_gl1_scaledvec >> bit) & 0x1U) == 0x1U) {
                // If the bit is set, add it to the vector of scaled trigger bits
                trig_bits.push_back(bit);
            }

            // Check if the current bit is set in the GL1 raw vector
            // (i.e., the bit at position 'bit' is 1)
            if (((b_gl1_rawvec >> bit) & 0x1U) == 0x1U) {
                // If the bit is 8 or 9, increment the raw MBD counter
                // These bits correspond to specific MBD triggers
                if (bit == 8 || bit == 9) raw_mbd++;

                // Fill the histogram with the bit value
                // This records the frequency of each trigger bit being set
                h_triggers->Fill(bit);

                // Add the bit to the vector of raw trigger bits
                trig_bits_raw.push_back(bit);
            }
        }
        std::cout << "\033[31mScaled Trigger Bits: [ ";
        for (auto bit : trig_bits) {
            std::cout << bit << " ";
        }
        std::cout << "]\033[0m" << std::endl;

        std::cout << "\033[31mRaw Trigger Bits: [ ";
        for (auto bit : trig_bits_raw) {
            std::cout << bit << " ";
        }
        std::cout << "]\033[0m" << std::endl;

        std::cout << "\033[31mRaw MBD Count: " << raw_mbd << "\033[0m" << std::endl;
        
        

        // Loop over each row in the map (35 rows)
        for (int j = 0; j < 35; j++) {
            // Loop over each column in the map (12 columns)
            for (int k = 0; k < 12; k++) {
                // For the first 32 rows, initialize energymap and goodmap
                if (j < 32) {
                    energymap[k][j] = 0.0; // Initialize the energy map entry to 0.0
                    goodmap[k][j] = 0; // Initialize the good map entry to 0
                }

                // Initialize all entries in the HCAL inner energy map to 0.0
                energymap_hcalin[k][j] = 0.0;

                // Initialize all entries in the HCAL outer energy map to 0.0
                energymap_hcalout[k][j] = 0.0;

                // Initialize all entries in the EMCal energy map to 0.0
                energymap_emcal[k][j] = 0.0;

                // Initialize all entries in the extended energy map to 0.0
                energymap_extend[k][j] = 0.0;

                // Initialize all entries in the extended good map to 0
                goodmap_extend[k][j] = 0;
            }
        }
        std::cout << "Energy and good maps initialized for entry " << i << std::endl;
        
        // Loop over each entry in the EMCal energy vector
        for (int ie = 0; ie < b_emcal_energy->size(); ie++) {
            // Calculate the eta bin and phi bin for the current entry
            int ebin = b_emcal_etabin->at(ie) / 8;
            int pbin = b_emcal_phibin->at(ie) / 8;

            // Update the goodmap and goodmap_extend with the good status of the current entry
            goodmap[ebin][pbin] += b_emcal_good->at(ie);
            goodmap_extend[ebin][pbin] += b_emcal_good->at(ie);
            // If phi bin is less than 3, update the extended region in goodmap_extend
            if (pbin < 3) {
                goodmap_extend[ebin][pbin + 32] += b_emcal_good->at(ie);
            }

            // Skip entries that are not marked as good
            if (!b_emcal_good->at(ie)) continue;

            // Update the energy maps with the energy of the current entry
            energymap[ebin][pbin] += b_emcal_energy->at(ie);
            energymap_emcal[ebin][pbin] += b_emcal_energy->at(ie);
            energymap_extend[ebin][pbin] += b_emcal_energy->at(ie);

            // Fill the EMCal spread histogram with the current entry index and energy
            h_emcal_spread->Fill(ie, b_emcal_energy->at(ie));

            // If phi bin is less than 3, update the extended region in energymap_emcal and energymap_extend
            if (pbin < 3) {
                energymap_emcal[ebin][pbin + 32] += b_emcal_energy->at(ie);
                energymap_extend[ebin][pbin + 32] += b_emcal_energy->at(ie);
            }
        }

        std::cout << "Processed EMCal energy entries for entry " << i << std::endl;

        // Loop over each entry in the HCAL inner energy vector
        for (int ie = 0; ie < b_hcalin_energy->size(); ie++) {
            // Calculate the eta bin and phi bin for the current entry
            int ebin = b_hcalin_etabin->at(ie) / 2;
            int pbin = b_hcalin_phibin->at(ie) / 2;

            // If phi bin is less than 3, update the extended region in goodmap_extend
            if (pbin < 3) {
                goodmap_extend[ebin][pbin + 32] += b_hcalin_good->at(ie);
            }
            // Update the goodmap_extend with the good status of the current entry
            goodmap_extend[ebin][pbin] += b_hcalin_good->at(ie);

            // Skip entries that are not marked as good
            if (!b_hcalin_good->at(ie)) continue;

            // Update the energy maps with the energy of the current entry
            energymap_extend[ebin][pbin] += b_hcalin_energy->at(ie);
            energymap_hcalin[ebin][pbin] += b_hcalin_energy->at(ie);

            // If phi bin is less than 3, update the extended region in energymap_extend and energymap_hcalin
            if (pbin < 3) {
                energymap_extend[ebin][pbin + 32] += b_hcalin_energy->at(ie);
                energymap_hcalin[ebin][pbin + 32] += b_hcalin_energy->at(ie);
            }
        }
        std::cout << "Processed HCAL inner energy entries for entry " << i << std::endl;

        
        // Loop over each entry in the HCAL outer energy vector
        for (int ie = 0; ie < b_hcalout_energy->size(); ie++) {
            // Calculate the eta bin and phi bin for the current entry
            int ebin = b_hcalout_etabin->at(ie) / 2;
            int pbin = b_hcalout_phibin->at(ie) / 2;

            // If phi bin is less than 3, update the extended region in goodmap_extend
            if (pbin < 3) {
                goodmap_extend[ebin][pbin + 32] += b_hcalout_good->at(ie);
            }
            // Update the goodmap_extend with the good status of the current entry
            goodmap_extend[ebin][pbin] += b_hcalout_good->at(ie);

            // Skip entries that are not marked as good
            if (!b_hcalout_good->at(ie)) continue;

            // Update the energy maps with the energy of the current entry
            energymap_extend[ebin][pbin] += b_hcalout_energy->at(ie);
            energymap_hcalout[ebin][pbin] += b_hcalout_energy->at(ie);

            // If phi bin is less than 3, update the extended region in energymap_extend and energymap_hcalout
            if (pbin < 3) {
                energymap_hcalout[ebin][pbin + 32] += b_hcalout_energy->at(ie);
                energymap_extend[ebin][pbin + 32] += b_hcalout_energy->at(ie);
            }
        }
        std::cout << "Processed HCAL outer energy entries for entry " << i << std::endl;


        // Loop over eta bins for jets
        for (int ie = 0; ie < 9; ie++) {
            // Loop over phi bins for jets
            for (int ip = 0; ip < 32; ip++) {
                // Initialize energy and good maps for jets to 0.0 and 0 respectively
                energymap_jet[ie][ip] = 0.0;
                energymap_jet_emcal[ie][ip] = 0.0;
                energymap_jet_hcalin[ie][ip] = 0.0;
                energymap_jet_hcalout[ie][ip] = 0.0;
                goodmap_jet[ie][ip] = 0;

                // Loop over sub-bins within the current eta and phi bins
                for (int is = 0; is < 16; is++) {
                    // Calculate the sub-bin indices for eta and phi
                    int sub_eta = ie + is % 4;
                    int sub_phi = ip + is / 4;

                    // Accumulate good and energy values from extended maps into jet maps
                    goodmap_jet[ie][ip] += goodmap_extend[sub_eta][sub_phi];
                    energymap_jet[ie][ip] += energymap_extend[sub_eta][sub_phi];
                    energymap_jet_emcal[ie][ip] += energymap_emcal[sub_eta][sub_phi];
                    energymap_jet_hcalin[ie][ip] += energymap_hcalin[sub_eta][sub_phi];
                    energymap_jet_hcalout[ie][ip] += energymap_hcalout[sub_eta][sub_phi];
                }
            }
        }
        std::cout << "Processed jet energy and good maps for entry " << i << std::endl;

        
        // Initialize variables to store maximum energy values and corresponding bin indices
        float max_energy = 0.0;
        float max_energy_clus = 0.0;
        float max_energy_hcal = 0.0;
        float max_energy_jet = 0.0;
        int jet_ebin = 0;
        int jet_pbin = 0;
        int ebin = 0;
        int pbin = 0;
        int hcal_ebin = 0;
        int hcal_pbin = 0;

        // Loop over all clusters to find the maximum cluster energy (ECore)
        for (int j = 0; j < b_cluster_n; j++) {
            if (b_cluster_ecore->at(j) > max_energy_clus) {
                max_energy_clus = b_cluster_ecore->at(j); // Update the maximum cluster energy if current cluster has higher energy
            }
        }
        std::cout << "Maximum cluster energy (ECore): " << max_energy_clus << " for entry " << i << std::endl;


        // Loop over eta bins
        for (int j = 0; j < 32; j++) {
            if (nentries <= 10) std::cout << j << " \t |"; // Print the eta bin index if entries are 10 or fewer

            // Loop over phi bins
            for (int k = 0; k < 12; k++) {
                if (k < 9) {
                    // Check if current jet energy is the highest found so far
                    if (energymap_jet[k][j] > max_energy_jet) {
                        max_energy_jet = energymap_jet[k][j]; // Update maximum jet energy
                        jet_ebin = k; // Update jet eta bin index
                        jet_pbin = j; // Update jet phi bin index
                    }
                }
                if (nentries <= 10) std::cout << energymap[k][j] << "\t"; // Print the energy value if entries are 10 or fewer

                // Check if current EMCal energy is the highest found so far
                if (energymap[k][j] > max_energy) {
                    max_energy = energymap[k][j]; // Update maximum EMCal energy
                    ebin = k; // Update EMCal eta bin index
                    pbin = j; // Update EMCal phi bin index
                }

                // Check if current HCal (inner + outer) energy is the highest found so far
                if (energymap_hcalout[k][j] + energymap_hcalin[k][j] > max_energy_hcal) {
                    max_energy_hcal = energymap_hcalin[k][j] + energymap_hcalout[k][j]; // Update maximum HCal energy
                    hcal_ebin = k; // Update HCal eta bin index
                    hcal_pbin = j; // Update HCal phi bin index
                }
            }
            
            if (nentries <= 10) std::cout << " | " << std::endl; // Print a separator if entries are 10 or fewer
        }
        std::cout << "Maximum energies for entry " << i << " - EMCal: " << max_energy << ", HCal: " << max_energy_hcal << ", Jet: " << max_energy_jet << std::endl;


        // Fill histograms with maximum energy values
        h_emcal_energy->Fill(max_energy); // Fill EMCal energy histogram
        h_cluster_energy->Fill(max_energy_clus); // Fill cluster energy histogram
        h_hcal_energy->Fill(max_energy_hcal); // Fill HCal energy histogram
        h_jet_energy->Fill(max_energy_jet); // Fill jet energy histogram
        
        
        // Save calculated values for tabulated output
        maxEcoreDansCode.push_back(max_energy_clus);
        max8x8DansCalc.push_back(max_energy);
        max8x8MyCalc.push_back(max8x8MyCalcVal);

        // Fill 2D histograms with goodmap values and maximum energy
        h2_good->Fill(goodmap[ebin][pbin], max_energy); // Fill 2D histogram for EMCal
        h2_good_jet->Fill(goodmap_jet[jet_ebin][jet_pbin], max_energy_jet); // Fill 2D histogram for jets

        std::cout << "Filled histograms for entry " << i << std::endl;
        // If raw MBD (Minimum Bias Detector) triggers are set
        if (raw_mbd == 2) {
            //Maximum 8x8 Energy Sum (EMCAL) [GeV]
            h_emcal_energy_ref->Fill(max_energy); // Fill reference EMCal energy histogram
            h_jet_energy_ref->Fill(max_energy_jet); // Fill reference jet energy histogram
            
            // Loop over raw trigger bits and fill histograms
            for (auto &b : trig_bits_raw) {
                h_mbd_triggers->Fill(b); // Fill MBD triggers histogram
                h_emcal_energy_raw_gl1[b]->Fill(max_energy); // Fill raw GL1 EMCal energy histogram for each trigger bit
                h_cluster_energy_raw_gl1[b]->Fill(max_energy_clus); // Fill raw GL1 cluster energy histogram for each trigger bit
                h_jet_energy_raw_gl1[b]->Fill(max_energy_jet); // Fill raw GL1 jet energy histogram for each trigger bit
                
            }
            std::cout << "Filled raw GL1 histograms for entry " << i << std::endl;
        }

        // Loop over scaled trigger bits and fill histograms
        for (auto &b : trig_bits) {
            h2_good_gl1[b]->Fill(goodmap[ebin][pbin], max_energy); // Fill 2D GL1 goodmap histogram for each trigger bit
            h2_good_jet_gl1[b]->Fill(goodmap_jet[jet_ebin][jet_pbin], max_energy_jet); // Fill 2D GL1 goodmap jet histogram for each trigger bit
            h_cluster_energy_gl1[b]->Fill(max_energy_clus); // Fill GL1 cluster energy histogram for each trigger bit
            h_emcal_energy_gl1[b]->Fill(max_energy); // Fill GL1 EMCal energy histogram for each trigger bit -- I think this is 8x8
            h_hcal_energy_gl1[b]->Fill(max_energy_hcal); // Fill GL1 HCal energy histogram for each trigger bit
            h_jet_energy_gl1[b]->Fill(max_energy_jet); // Fill GL1 jet energy histogram for each trigger bit
            h_jet_emcal_energy_gl1[b]->Fill(energymap_jet_emcal[jet_ebin][jet_pbin]); // Fill GL1 jet EMCal energy histogram for each trigger bit
            h_jet_hcalin_energy_gl1[b]->Fill(energymap_jet_hcalin[jet_ebin][jet_pbin]); // Fill GL1 jet HCal inner energy histogram for each trigger bit
            h_jet_hcalout_energy_gl1[b]->Fill(energymap_jet_hcalout[jet_ebin][jet_pbin]); // Fill GL1 jet HCal outer energy histogram for each trigger bit

            // Fill jet EMCal fraction histogram
            h_jet_emcal_fraction_gl1[b]->Fill(max_energy_jet, energymap_jet_emcal[jet_ebin][jet_pbin] / max_energy_jet);
            // Fill jet HCal inner fraction histogram
            h_jet_hcalin_fraction_gl1[b]->Fill(max_energy_jet, energymap_jet_hcalin[jet_ebin][jet_pbin] / max_energy_jet);
            // Fill jet HCal outer fraction histogram
            h_jet_hcalout_fraction_gl1[b]->Fill(max_energy_jet, energymap_jet_hcalout[jet_ebin][jet_pbin] / max_energy_jet);

            // Fill imbalance histograms if energy conditions are met
            if (energymap_jet_emcal[jet_ebin][jet_pbin] < 1) {
                h_hcalout_imbalance[b]->Fill(energymap_jet_hcalout[jet_ebin][jet_pbin]);
            }
            if (energymap_jet_hcalout[jet_ebin][jet_pbin] < 1) {
                h_emcal_imbalance[b]->Fill(energymap_jet_emcal[jet_ebin][jet_pbin]);
            }
        }
        std::cout << "Filled scaled GL1 histograms for entry " << i << std::endl;
    }

    TFile *outfile = new TFile(output_filename, "RECREATE");
    TTree *outputTree = new TTree("SummaryTree", "Summary data including accumulated counts");
    outputTree->Branch("accumulated_live", accumulated_live.data(), "accumulated_live[64]/F");
    outputTree->Branch("accumulated_scaled", accumulated_scaled.data(), "accumulated_scaled[64]/F");
    
    h_mbd_triggers->Write();
    h_triggers->Write();
    std::string name;
    h_emcal_spread->Write();
    
    outputTree->Fill();

    for (int i = 0; i < 64; ++i) {
        histograms_ecore[i]->Write();
        histograms_energy[i]->Write();
        histograms_ecore_liveCounter[i]->Write();
        histograms_energy_liveCounter[i]->Write();
        
        
        h2_good_gl1[i]->Write();
        h2_good_jet_gl1[i]->Write();
        h_emcal_energy_gl1[i]->Write();
        h_emcal_energy_raw_gl1[i]->Write();
        h_emcal_energy_scaled_gl1[i] = (TH1D*) h_emcal_energy_gl1[i]->Clone();
        name = "h_emcal_energy_scaled_gl1_" + to_string(i);
        h_emcal_energy_scaled_gl1[i]->SetName(name.c_str());
        h_emcal_energy_scaled_gl1[i]->Write();

        h_cluster_energy_gl1[i]->Write();
        h_cluster_energy_raw_gl1[i]->Write();
        h_cluster_energy_scaled_gl1[i] = (TH1D*) h_cluster_energy_gl1[i]->Clone();
        name = "h_cluster_energy_scaled_gl1_" + to_string(i);
        h_cluster_energy_scaled_gl1[i]->SetName(name.c_str());
        h_cluster_energy_scaled_gl1[i]->Write();

        h_hcal_energy_gl1[i]->Write();
        h_jet_energy_gl1[i]->Write();
        h_jet_energy_raw_gl1[i]->Write();
        h_jet_energy_scaled_gl1[i] = (TH1D*) h_jet_energy_gl1[i]->Clone();
        name = "h_jet_energy_scaled_gl1_" + to_string(i);
        h_jet_energy_scaled_gl1[i]->SetName(name.c_str());
        h_jet_energy_scaled_gl1[i]->Write();

        h_jet_emcal_energy_gl1[i]->Write();
        h_jet_hcalin_energy_gl1[i]->Write();
        h_jet_hcalout_energy_gl1[i]->Write();
        h_jet_emcal_fraction_gl1[i]->Write();
        h_jet_hcalin_fraction_gl1[i]->Write();
        h_jet_hcalout_fraction_gl1[i]->Write();
        h_hcalout_imbalance[i]->Write();
        h_emcal_imbalance[i]->Write();
    }
    
    outputTree->Write();
    h2_good->Write();
    h2_good_jet->Write();
    h_emcal_energy_ref->Write();
    h_jet_energy_ref->Write();
    h_emcal_energy->Write();
    h_cluster_energy->Write();
    h_jet_energy->Write();
    h_hcal_energy->Write();
    
    outfile->Close();
    file->Close();
    std::cout << "Histograms and accumulated counts written to file and files closed." << std::endl;
    
    // Print the tabulated results
    std::cout << "Tabulated Results:" << std::endl;
    std::cout << "Entry\tmaxEcoreMyCalc\tmaxEcoreDansCode\t8x8MyCalc\t8x8DansCalc" << std::endl;
    for (size_t i = 0; i < maxEcoreMyCalc.size(); ++i) {
        std::cout << i << "\t\t" << maxEcoreMyCalc[i] << "\t\t" << maxEcoreDansCode[i] << "\t\t" << max8x8MyCalc[i] << "\t\t" << max8x8DansCalc[i] << std::endl;
    }
}
