#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TKey.h>
#include <TSystem.h>
#include <TList.h>
#include <TDirectory.h>
#include <iostream>
#include <map>
#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <sys/stat.h>
#include <sys/types.h>

void CombineHistograms_oneRun(int RunNumber) {
    TString directoryPath = TString::Format("/sphenix/u/patsfan753/scratch/analysis/calotriggeremulator/outputHistRootFiles/%d", RunNumber);
    TString combinedDirectoryPath = directoryPath + "/combined";
    TString outputPath = combinedDirectoryPath + "/Combined_TriggerHists_" + TString::Format("%d", RunNumber) + ".root";
    std::string scaleDownFilePath = std::string(combinedDirectoryPath.Data()) + "/Scaledowns_" + std::to_string(RunNumber) + ".txt";

    // Create the combined directory if it does not exist
    struct stat info;
    if (stat(combinedDirectoryPath.Data(), &info) != 0 || !(info.st_mode & S_IFDIR)) {
        if (gSystem->MakeDirectory(combinedDirectoryPath.Data()) != 0) {
            std::cerr << "Failed to create combined directory!" << std::endl;
            return;
        }
    }

    TFile *outputFile = new TFile(outputPath, "RECREATE");
    if (!outputFile->IsOpen()) {
        std::cerr << "Failed to create output file!" << std::endl;
        return;
    }
    std::cout << "Output file created: " << outputPath << std::endl;

    void *dirp = gSystem->OpenDirectory(directoryPath);
    const char *filename;
    std::map<std::string, TH1*> histogramMap;
    std::vector<float> accumulated_live(64, 0.0);
    std::vector<float> accumulated_scaled(64, 0.0);

    while ((filename = gSystem->GetDirEntry(dirp))) {
        TString filepath = TString(directoryPath) + "/" + TString(filename);
        if (!filepath.EndsWith(".root")) continue;

        TFile *file = TFile::Open(filepath);
        if (!file || file->IsZombie()) {
            std::cerr << "Error opening file: " << filepath << std::endl;
            continue;
        }
        std::cout << "Processing file: " << filepath << std::endl;

        TTree *summaryTree;
        file->GetObject("SummaryTree", summaryTree);
        if (!summaryTree) {
            std::cerr << "SummaryTree not found in " << filepath << std::endl;
            file->Close();
            continue;
        }

        float live[64], scaled[64];
        summaryTree->SetBranchAddress("accumulated_live", live);
        summaryTree->SetBranchAddress("accumulated_scaled", scaled);
        summaryTree->GetEntry(0);

        for (int i = 0; i < 64; i++) {
            accumulated_live[i] += live[i];
            accumulated_scaled[i] += scaled[i];
            std::cout << "Accumulated live[" << i << "]: " << accumulated_live[i] << ", Accumulated scaled[" << i << "]: " << accumulated_scaled[i] << std::endl;
        }

        TList* list = file->GetListOfKeys();
        TIter next(list);
        TKey* key;
        while ((key = (TKey*)next())) {
            TH1* hist = (TH1*)file->Get(key->GetName());
            if (!hist) continue;
            if (!hist->InheritsFrom("TH1")) continue;

            if (histogramMap.find(key->GetName()) == histogramMap.end()) {
                histogramMap[key->GetName()] = (TH1*)hist->Clone();
                histogramMap[key->GetName()]->SetDirectory(0);
                std::cout << "New histogram added: " << key->GetName() << std::endl;
            } else {
                histogramMap[key->GetName()]->Add(hist);
                std::cout << "Existing histogram updated: " << key->GetName() << std::endl;
            }
        }

        file->Close();
    }

    std::ofstream scaleDownFile(scaleDownFilePath);
    if (!scaleDownFile.is_open()) {
        std::cerr << "Failed to create scale-down factor file!" << std::endl;
        return;
    }

    outputFile->cd(); // Ensure the output file is the current directory
    for (int i = 0; i < 64; i++) {
        if (accumulated_scaled[i] > 0) { // Avoid division by zero
            float scaleFactor = accumulated_live[i] / accumulated_scaled[i];
            scaleDownFile << i << " " << scaleFactor << std::endl;
            std::string histName = "trigger_" + std::to_string(i) + "_ecore";
            if (histogramMap.find(histName) != histogramMap.end()) {
                TH1* unscaledHist = (TH1*)histogramMap[histName]->Clone((histName + "_unscaled").c_str());
                unscaledHist->SetDirectory(outputFile);
                unscaledHist->Write();
                histogramMap[histName]->Scale(scaleFactor);
                std::cout << "Scaled " << histName << " by factor " << scaleFactor << std::endl;
            }
            histName = "trigger_" + std::to_string(i) + "_energy";
            if (histogramMap.find(histName) != histogramMap.end()) {
                TH1* unscaledHist = (TH1*)histogramMap[histName]->Clone((histName + "_unscaled").c_str());
                unscaledHist->SetDirectory(outputFile);
                unscaledHist->Write();
                histogramMap[histName]->Scale(scaleFactor);
                std::cout << "Scaled " << histName << " by factor " << scaleFactor << std::endl;
            }
            // Add live counter histograms
            histName = "trigger_" + std::to_string(i) + "_ecore_liveCounter";
            if (histogramMap.find(histName) != histogramMap.end()) {
                TH1* unscaledHist = (TH1*)histogramMap[histName]->Clone((histName + "_unscaled").c_str());
                unscaledHist->SetDirectory(outputFile);
                unscaledHist->Write();
                histogramMap[histName]->Scale(scaleFactor);
                std::cout << "Scaled " << histName << " by factor " << scaleFactor << std::endl;
            }
            histName = "trigger_" + std::to_string(i) + "_energy_liveCounter";
            if (histogramMap.find(histName) != histogramMap.end()) {
                TH1* unscaledHist = (TH1*)histogramMap[histName]->Clone((histName + "_unscaled").c_str());
                unscaledHist->SetDirectory(outputFile);
                unscaledHist->Write();
                histogramMap[histName]->Scale(scaleFactor);
                std::cout << "Scaled " << histName << " by factor " << scaleFactor << std::endl;
            }
        }
    }

    for (auto& entry : histogramMap) {
        outputFile->cd();
        entry.second->Write();
        std::cout << "Histogram written to file: " << entry.first << std::endl;
    }

    scaleDownFile.close();
    outputFile->Close();
    std::cout << "Histograms combined and saved to " << outputPath << std::endl;
}

std::vector<int> ParseRunNumbers(const std::string& runNumbersStr) {
    std::vector<int> runNumbers;
    std::stringstream ss(runNumbersStr);
    std::string runNumber;
    while (std::getline(ss, runNumber, ',')) {
        runNumbers.push_back(std::stoi(runNumber));
    }
    return runNumbers;
}

void CombineHistograms(const std::vector<int>& runNumbers) {
    for (const int& runNumber : runNumbers) {
        CombineHistograms_oneRun(runNumber);
    }
}

void CombineHistogramsWrapper(const std::string& runNumbersStr) {
    std::vector<int> runNumbers = ParseRunNumbers(runNumbersStr);
    CombineHistograms(runNumbers);
}

// Ensure main is declared after the above functions
void CombineHistograms(const char* runNumbersStr) {
    CombineHistogramsWrapper(runNumbersStr);
}
