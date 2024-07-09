void CombineOnlineHistograms_CEMC() {
    // Array of input ROOT file paths
    const char* inputFiles[] = {
    "/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/CEMCMON/Run_47722-CEMCMON_0.root",
    "/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/CEMCMON/Run_47722-CEMCMON_1.root",
    "/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/CEMCMON/Run_47722-CEMCMON_10.root",
    "/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/CEMCMON/Run_47722-CEMCMON_11.root",
    "/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/CEMCMON/Run_47722-CEMCMON_12.root",
    "/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/CEMCMON/Run_47722-CEMCMON_13.root",
    "/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/CEMCMON/Run_47722-CEMCMON_14.root",
    "/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/CEMCMON/Run_47722-CEMCMON_15.root",
    "/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/CEMCMON/Run_47722-CEMCMON_2.root",
    "/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/CEMCMON/Run_47722-CEMCMON_3.root",
    "/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/CEMCMON/Run_47722-CEMCMON_4.root",
    "/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/CEMCMON/Run_47722-CEMCMON_5.root",
    "/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/CEMCMON/Run_47722-CEMCMON_6.root",
    "/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/CEMCMON/Run_47722-CEMCMON_7.root",
    "/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/CEMCMON/Run_47722-CEMCMON_8.root",
    "/sphenix/lustre01/sphnxpro/commissioning/online_monitoring/histograms/CEMCMON/Run_47722-CEMCMON_9.root"
    };
    int numFiles = sizeof(inputFiles) / sizeof(char*);

    std::map<std::string, TH2*> histogramMap;
    std::regex pattern("^h2_cemc_hits_trig_bit_([0-9]|[1-5][0-9]|6[0-3])$");

    for (int i = 0; i < numFiles; i++) {
        TFile* inputFile = TFile::Open(inputFiles[i], "READ");
        if (!inputFile || inputFile->IsZombie()) {
            std::cerr << "Failed to open file: " << inputFiles[i] << std::endl;
            continue;
        }
        std::cout << "Processing file: " << inputFiles[i] << std::endl;

        TIter next(inputFile->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)next())) {
            TH2* hist = dynamic_cast<TH2*>(key->ReadObj());
            if (!hist) {
                std::cerr << "Not a TH2 object or failed to cast: " << key->GetName() << std::endl;
                continue;
            }
            std::string histName = hist->GetName();
            std::cout << "Checking histogram: " << histName << " with entries: " << hist->GetEntries() << std::endl;
            if (std::regex_match(histName, pattern)) {
                if (hist->GetEntries() == 0) {
                    std::cout << "Skipping empty histogram: " << histName << std::endl;
                    continue; // Skip empty histograms
                }
                if (histogramMap.find(histName) == histogramMap.end()) {
                    histogramMap[histName] = (TH2*)hist->Clone();
                    histogramMap[histName]->SetDirectory(0); // Detach from file
                    std::cout << "Cloned new histogram: " << histName << std::endl;
                } else {
                    histogramMap[histName]->Add(hist);
                    std::cout << "Added to existing histogram: " << histName << std::endl;
                }
            }
            delete hist; // Prevent memory leaks by deleting the histogram
        }
        inputFile->Close();
    }
    TFile* outputFile = new TFile("47722_combined_histograms.root", "RECREATE");
    for (auto& it : histogramMap) {
        it.second->Write();
        std::cout << "Writing combined histogram: " << it.first << std::endl;
    }
    outputFile->Close();
    std::cout << "All histograms have been combined and written to combined_histograms.root" << std::endl;
}
