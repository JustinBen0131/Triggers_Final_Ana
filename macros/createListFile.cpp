#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>
#include <vector>
#include <sstream>

namespace fs = std::filesystem;

void createListFile(const std::vector<std::string>& runNumbers) {
    std::string baseOutputPath = "/sphenix/user/patsfan753/analysis/calotriggeremulator/";

    for (const auto& runNumber : runNumbers) {
        std::string runDirectory = baseOutputPath + "output/" + runNumber;
        std::string listFilePath = baseOutputPath + runNumber + "_segmentList.list";

        std::ofstream listFile(listFilePath);
        if (!listFile.is_open()) {
            std::cerr << "Error: Unable to open file " << listFilePath << std::endl;
            continue;
        }

        for (const auto& entry : fs::directory_iterator(runDirectory)) {
            if (entry.path().extension() == ".root") {
                listFile << entry.path().string() << std::endl;
            }
        }

        listFile.close();
        std::cout << "List file created at " << listFilePath << std::endl;
    }
}

std::vector<std::string> parseRunNumbers(const std::string& runNumbersStr) {
    std::vector<std::string> runNumbers;
    std::stringstream ss(runNumbersStr);
    std::string runNumber;
    while (std::getline(ss, runNumber, ',')) {
        runNumbers.push_back(runNumber);
    }
    return runNumbers;
}

void createListFileWrapper(const std::string& runNumbersStr) {
    std::vector<std::string> runNumbers = parseRunNumbers(runNumbersStr);
    createListFile(runNumbers);
}

// ROOT-compatible main function
void createListFile(const char* runNumbersStr) {
    createListFileWrapper(runNumbersStr);
}
