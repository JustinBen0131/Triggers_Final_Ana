void CombineHistograms() {
    // List of ROOT files to combine
    std::vector<std::string> inputFiles = {
        "44686_combined_histograms.root",
        "44687_combined_histograms.root",
        "44689_combined_histograms.root",
        "44690_combined_histograms.root",
        "44691_combined_histograms.root",
        "44702_combined_histograms.root",
        "44703_combined_histograms.root",
        "44707_combined_histograms.root"
    };

    // Map to store the final combined histograms
    std::map<std::string, TH2F*> combinedHistograms;

    // Loop over all files
    for (const auto& file : inputFiles) {
        TFile* inputFile = TFile::Open(file.c_str(), "READ");
        if (!inputFile || inputFile->IsZombie()) {
            std::cerr << "Failed to open file: " << file << std::endl;
            continue;
        }
        std::cout << "Processing file: " << file << std::endl;

        // Iterate over all keys (objects) in the file
        TIter next(inputFile->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)next())) {
            TH2F* hist = dynamic_cast<TH2F*>(key->ReadObj());
            if (!hist) continue;  // Skip if not a TH2F

            std::string histName = hist->GetName();
            std::cout << "Found histogram: " << histName << " with entries: " << hist->GetEntries() << std::endl;

            if (combinedHistograms.find(histName) == combinedHistograms.end()) {
                // Clone histogram if it does not exist in the map
                TH2F* clonedHist = (TH2F*)hist->Clone();
                clonedHist->SetDirectory(0); // Detach from file
                combinedHistograms[histName] = clonedHist;
                std::cout << "Cloned new histogram: " << histName << std::endl;
            } else {
                // Add to the existing histogram if it already exists
                std::cout << "Adding to existing histogram: " << histName << std::endl;
                combinedHistograms[histName]->Add(hist);
            }
        }
        inputFile->Close(); // Close the input file
    }

    // Create a new ROOT file to save the combined histograms
    TFile* outputFile = new TFile("combined_histograms_DansRunList.root", "RECREATE");
    for (auto& histPair : combinedHistograms) {
        histPair.second->Write(); // Write each histogram to the output file
        std::cout << "Writing combined histogram: " << histPair.first << std::endl;
    }

    outputFile->Close(); // Close the output file
    std::cout << "All histograms have been combined and saved to combined_histograms_PhotonTrigAna.root" << std::endl;
}
