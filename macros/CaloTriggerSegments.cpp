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

    ULong64_t gl1_scaled[64], gl1_live[64], gl1_scaledvec, gl1_livevec;
    std::vector<float> *emcal_energy = nullptr, *cluster_ecore = nullptr;
    std::vector<int> *emcal_phibin = nullptr, *emcal_etabin = nullptr;
    std::vector<bool> *emcal_good = nullptr;
    float mbd_vertex_z;
    UInt_t trigger_sum_emcal[6144], trigger_sumkey_emcal[6144];


    tree->SetBranchAddress("gl1_scaled", gl1_scaled);
    tree->SetBranchAddress("gl1_live", gl1_live);
    tree->SetBranchAddress("gl1_scaledvec", &gl1_scaledvec);
    tree->SetBranchAddress("gl1_livevec", &gl1_livevec); //equivalent to TriggerVector
    
    /*
     Vectors of length 24576 for 96 * 216 towers
     */
    tree->SetBranchAddress("emcal_energy", &emcal_energy);
    tree->SetBranchAddress("emcal_good",&emcal_good);
    tree->SetBranchAddress("cluster_ecore", &cluster_ecore);
    tree->SetBranchAddress("emcal_phibin", &emcal_phibin);
    tree->SetBranchAddress("emcal_etabin", &emcal_etabin);
    tree->SetBranchAddress("trigger_sum_emcal", &trigger_sum_emcal);
    tree->SetBranchAddress("trigger_sumkey_emcal", &trigger_sumkey_emcal);
    
    /*
     Use for z vertex cut?
     */
    tree->SetBranchAddress("mbd_vertex_z", &mbd_vertex_z);

    TH1F *histograms_ecore[64];
    TH1F *histograms_energy[64];
    TH1F *histograms_ecore_liveCounter[64];
    TH1F *histograms_energy_liveCounter[64];
    std::vector<float> accumulated_live(64, 0);
    std::vector<float> accumulated_scaled(64, 0);

    for (int i = 0; i < 64; ++i) {
        histograms_ecore[i] = new TH1F(Form("trigger_%d_ecore", i), Form("Trigger %d Maximum Cluster ECore", i), 100, 0, 20);
        histograms_energy[i] = new TH1F(Form("trigger_%d_energy", i), Form("Trigger %d Maximum EMCAL Energy", i), 100, 0, 20);
        histograms_ecore_liveCounter[i] = new TH1F(Form("trigger_%d_ecore_liveCounter", i), Form("Trigger %d Maximum Cluster ECore, Live Counter", i), 100, 0, 20);
        histograms_energy_liveCounter[i] = new TH1F(Form("trigger_%d_energy_liveCounter", i), Form("Trigger %d Maximum EMCAL Energy, Live Counter", i), 100, 0, 20);
    }

    Long64_t nentries = std::min(tree->GetEntries(), (nEvents > 0 ? nEvents : tree->GetEntries()));
    std::cout << "Number of entries to process: " << nentries << std::endl;

    for (Long64_t i = 0; i < nentries; ++i) {
        tree->GetEntry(i);
        
        // Apply the z vertex cut
        if (std::abs(mbd_vertex_z) >= 10) {
            continue;
        }

        std::bitset<64> bits(gl1_scaledvec);
        std::bitset<64> bits_live(gl1_livevec);
        std::cout << "Processing entry " << i << ", gl1_scaledvec (bits): " << bits.to_string() << std::endl;
        std::cout << "Processing entry " << i << ", trigger_vector (bits): " << bits_live.to_string() << std::endl;
        
        std::unordered_map<EtaPhiKey, float> energy_map;
        std::map<int, std::map<int, float>> emcal_energies; // Declare emcal_energies here

        // Preprocess the emcal_energy vector to store energy values in a map for quick access
        for (size_t k = 0; k < emcal_energy->size(); ++k) {
            if ((*emcal_good)[k]) {
                EtaPhiKey key = {(*emcal_etabin)[k], (*emcal_phibin)[k]};
                energy_map[key] += (*emcal_energy)[k];
            }
        }
        
        // Loop over all trigger sums
        for (int idx = 0; idx < 6144; ++idx) {
            unsigned int sumk = trigger_sumkey_emcal[idx];
            uint16_t sum_phi = TriggerDefs::getSumPhiId(sumk) + 2 * TriggerDefs::getPrimitivePhiId_from_TriggerSumKey(sumk);
            uint16_t sum_eta = TriggerDefs::getSumEtaId(sumk);

            float energy_sum = 0.0;

            // Calculate the 8x8 energy sum for each tower in the 8x8 grid
            for (int itower = 0; itower < 64; ++itower) {
                int ieta = sum_eta * 8 + itower % 8;
                int iphi = sum_phi * 8 + itower / 8;

                EtaPhiKey key = {ieta, iphi};
                if (energy_map.find(key) != energy_map.end()) {
                    energy_sum += energy_map[key];
                }
            }

            // Store the calculated energy sum in the map
            emcal_energies[sum_eta][sum_phi] = energy_sum;
        }


        for (int j = 0; j < 64; ++j) {
            if (bits.test(j)) {
                accumulated_live[j] += gl1_live[j];
                accumulated_scaled[j] += gl1_scaled[j];
                std::cout << "Trigger index " << j << " fired. Accumulated live: " << accumulated_live[j] << ", Accumulated scaled: " << accumulated_scaled[j] << std::endl;

                if (!cluster_ecore->empty()) {
                    float maxEcore = *std::max_element(cluster_ecore->begin(), cluster_ecore->end());
                    histograms_ecore[j]->Fill(maxEcore);
                    std::cout << "Max ECore: " << maxEcore << " added to histogram." << std::endl;
                }

                float max8x8Energy = 0.0;
                for (const auto& entry : emcal_energies) {
                    for (const auto& sub_entry : entry.second) {
                        if (sub_entry.second > max8x8Energy) {
                            max8x8Energy = sub_entry.second;
                        }
                    }
                }
                histograms_energy[j]->Fill(max8x8Energy);
                std::cout << "Max 8x8 Energy: " << max8x8Energy << " added to histogram." << std::endl;
            }

            if (bits_live.test(j)) {
                if (!cluster_ecore->empty()) {
                    float maxEcore = *std::max_element(cluster_ecore->begin(), cluster_ecore->end());
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
    }

    TFile *outfile = new TFile(output_filename, "RECREATE");
    TTree *outputTree = new TTree("SummaryTree", "Summary data including accumulated counts");
    outputTree->Branch("accumulated_live", accumulated_live.data(), "accumulated_live[64]/F");
    outputTree->Branch("accumulated_scaled", accumulated_scaled.data(), "accumulated_scaled[64]/F");
    outputTree->Fill();

    for (int i = 0; i < 64; ++i) {
        histograms_ecore[i]->Write();
        histograms_energy[i]->Write();
        histograms_ecore_liveCounter[i]->Write();
        histograms_energy_liveCounter[i]->Write();
    }
    outputTree->Write();
    outfile->Close();
    file->Close();
    std::cout << "Histograms and accumulated counts written to file and files closed." << std::endl;
}
