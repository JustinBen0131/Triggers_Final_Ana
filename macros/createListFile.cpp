#include <iostream>
#include <fstream>
#include <filesystem>
#include <string>

namespace fs = std::filesystem;

void createListFile(const std::string& directory, const std::string& runNumber) {
    std::string outputPath = "output/" + runNumber + "/";
    std::string listFilePath = outputPath + runNumber + "_segmentList.list";

    // Create the output directory if it does not exist
    fs::create_directories(outputPath);

    std::ofstream listFile(listFilePath);
    if (!listFile.is_open()) {
        std::cerr << "Error: Unable to open file " << listFilePath << std::endl;
        return;
    }

    for (const auto& entry : fs::directory_iterator(directory)) {
        if (entry.path().extension() == ".root") {
            listFile << entry.path().string() << std::endl;
        }
    }

    listFile.close();
    std::cout << "List file created at " << listFilePath << std::endl;
}

int main(int argc, char* argv[]) {
    if (argc != 3) {
        std::cerr << "Usage: " << argv[0] << " <directory> <runNumber>" << std::endl;
        return 1;
    }

    std::string directory = argv[1];
    std::string runNumber = argv[2];

    createListFile(directory, runNumber);

    return 0;
}
