#include <TFile.h>
#include <TH2.h>
#include <TH1D.h>
#include <TKey.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <regex>
#include <map>

int custom_sector_mapping(unsigned int eta, unsigned int phi) {
    int sector = -1;
    if (eta >= 48 && eta < 96) {
        sector = (phi / 8);
    } else if (eta < 48) {
        sector = 32 + (phi / 8);
    }
    return sector;
}
int custom_ib_board(int eta, int phi) {
    int ib_board;
    if (eta >= 48 && eta < 56) {
        ib_board = 0;
    } else if (eta >= 56 && eta < 64) {
        ib_board = 1;
    } else if (eta >= 64 && eta < 72) {
        ib_board = 2;
    } else if (eta >= 72 && eta < 80) {
        ib_board = 3;
    } else if (eta >= 80 && eta < 88) {
        ib_board = 4;
    } else if (eta >= 88 && eta < 96) {
        ib_board = 5;
    } else if (eta >= 40 && eta < 48) {
        ib_board = 0;
    } else if (eta >= 32 && eta < 40) {
        ib_board = 1;
    } else if (eta >= 24 && eta < 32) {
        ib_board = 2;
    } else if (eta >= 16 && eta < 24) {
        ib_board = 3;
    } else if (eta >= 8 && eta < 16) {
        ib_board = 4;
    } else if (eta >= 0 && eta < 8) {
        ib_board = 5;
    } else {
        ib_board = -1;  // For out-of-range values, if any
    }
    return ib_board;
}
void plotEmptyGrid() {
    // Define the output directory
    std::string outputDirectory = "/Users/patsfan753/Desktop/Desktop/TriggerAna/OnlineMonitoringAna/";

    // Create an empty 2D histogram
    TH2D* hist = new TH2D("emptyHist", "Template EMCal;#eta;#phi", 96, 0, 96, 256, 0, 256);

    // Configure the histogram appearance
    hist->SetStats(kFALSE); // Remove the statistics box
    gStyle->SetPalette(kBird); // Set the color map to kBird
    hist->GetXaxis()->SetTitle("#eta");
    hist->GetYaxis()->SetTitle("#phi");
    hist->GetZaxis()->SetTitle("Number of events");
    hist->GetZaxis()->SetTitleOffset(1.3);

    // Create a canvas to draw the histogram
    TCanvas* canvas = new TCanvas("canvas", "Canvas", 1200, 1600);
    hist->Draw("COLZ");
    gPad->SetRightMargin(0.145);

    // Add custom visual elements
    for (int i = 0; i <= hist->GetNbinsX(); i += 8) {
        TLine *line = new TLine(i, 0, i, hist->GetNbinsY());
        line->SetLineStyle(1);
        line->SetLineColor(kBlack);
        line->SetLineWidth(2);
        line->Draw();
    }
    for (int j = 0; j <= hist->GetNbinsY(); j += 8) {
        TLine *line = new TLine(0, j, hist->GetNbinsX(), j);
        line->SetLineStyle(1);
        line->SetLineColor(kBlack);
        line->SetLineWidth(2);
        line->Draw();
    }

    TLine *middleLine = new TLine(hist->GetNbinsX() / 2, 0, hist->GetNbinsX() / 2, hist->GetNbinsY());
    middleLine->SetLineWidth(4);
    middleLine->Draw();

    for (int i = 0; i < hist->GetNbinsX(); i += 8) {
        for (int j = 0; j < hist->GetNbinsY(); j += 8) {
            int ib_board = custom_ib_board(i, j);
            TLatex *text = new TLatex(i + 4, j + 4, Form("%d", ib_board));
            text->SetTextAlign(22);
            text->SetTextSize(0.022);
            text->SetTextFont(62);
            text->Draw();
        }
    }
    for (int i = 0; i < hist->GetNbinsX(); i += 48) {
        for (int j = 0; j < hist->GetNbinsY(); j += 8) {
            int sectorNumber = custom_sector_mapping(i, j);
            TText *text = new TText(i + 24, j + 4, Form("%d", sectorNumber));
            text->SetTextAlign(22);
            text->SetTextSize(0.03);
            text->SetTextFont(62);
            text->Draw();
        }
    }
    hist->GetXaxis()->SetNdivisions(-512);
    hist->GetYaxis()->SetNdivisions(-32);

    // Save the canvas as an image file
    std::string outputFileName = outputDirectory + "EmptyGrid.png";
    canvas->SaveAs(outputFileName.c_str());

    delete canvas; // Clean up
    delete hist; // Clean up
}

void plotCombinedOnlineTriggerHists() {
    // Define the output directory
    std::string outputDirectory = "/Users/patsfan753/Desktop/Desktop/TriggerAna/OnlineMonitoringAna/";

    // Mapping of trigger indices to descriptions
    std::map<int, std::string> triggerMap = {
        {1, "ZDC South"}, {2, "ZDC North"}, {3, "ZDC Coincidence"}, {4, "HCAL Singles"},
        {5, "HCAL Coincidence"}, {8, "MBD S >= 1"}, {9, "MBD N >= 1"}, {10, "MBD N&S >= 1"},
        {11, "MBD N&S >= 2"}, {12, "MBD N&S >= 1, vtx < T1"}, {13, "MBD N&S >= 1, vtx < T2"},
        {14, "MBD N&S >= 1, vtx < T3"}, {15, "HCAL Singles + MBD NS >= 1"}, {16, "Jet 1 + MBD NS >= 1 -4"},
        {17, "Jet 2 + MBD NS >= 1 -6"}, {18, "Jet 3 + MBD NS >= 1 -8"}, {19, "Jet 4 + MBD NS >= 1 -10 GeV"},
        {20, "Jet 1"}, {21, "Jet 2"}, {22, "Jet 3"}, {23, "Jet 4"}, {24, "Photon 1 + MBD NS >= 1"},
        {25, "Photon 2 + MBD NS >= 1"}, {26, "Photon 3 + MBD NS >= 1"}, {27, "Photon 4 + MBD NS >= 1"},
        {28, "Photon 1"}, {29, "Photon 2"}, {30, "Photon 3"}, {31, "Photon 4"}
    };

    // Open the ROOT file
    TFile* file = TFile::Open("combined_histograms_DansRunList.root", "READ");
    if (!file || file->IsZombie()) {
        std::cerr << "Error opening file" << std::endl;
        return;
    }

    // Create a 1D histogram for trigger indices
    TH1D* triggerHist = new TH1D("triggerHist", "Total Number of Events by Trigger Index;Trigger Index;Total Number of Events", 64, 0, 64);

    // Iterate over all keys in the file
    TIter next(file->GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)next())) {
        // Attempt to read the object
        TObject* obj = key->ReadObj();
        if (TH2* hist = dynamic_cast<TH2*>(obj)) {
            // Extract the trigger index from the histogram name using regex
            std::regex rgx(".*_(\\d+)$");
            std::smatch match;
            std::string histName = hist->GetName();
            if (std::regex_search(histName, match, rgx) && match.size() > 1) {
                int triggerIndex = std::stoi(match[1].str());
                if (triggerIndex >= 0 && triggerIndex < 64) {
                    // Sum all entries and assign to the corresponding bin in triggerHist
                    triggerHist->SetBinContent(triggerIndex + 1, hist->GetEntries());

                    // Set the histogram title to include the trigger description and run number
                    if (triggerMap.find(triggerIndex) != triggerMap.end()) {
                        std::string title = triggerMap[triggerIndex] + ", Run 45161, 45162, 45164, 45167, 45174, 45175, 45177, 45178";
                        hist->SetTitle(title.c_str());
                    }
                }
            }

            // Configure the histogram appearance
            hist->SetStats(kFALSE); // Remove the statistics box
            gStyle->SetPalette(kBird); // Set the color map to kBird
            hist->GetXaxis()->SetTitle("#eta");
            hist->GetYaxis()->SetTitle("#phi");
            hist->GetZaxis()->SetTitle("Number of events");
            hist->GetZaxis()->SetTitleOffset(1.3);



            // Calculate percentiles for color scaling
            std::vector<double> contents;
            for (int i = 1; i <= hist->GetNbinsX(); ++i) {
                for (int j = 1; j <= hist->GetNbinsY(); ++j) {
                    double content = hist->GetBinContent(i, j);
                    if (content > 0) {
                        contents.push_back(content);
                    }
                }
            }
            std::sort(contents.begin(), contents.end());
            double lowPercentile = contents.size() > 100 ? contents[contents.size() / 100] : contents.front();
            double highPercentile = contents.size() > 100 ? contents[contents.size() * 99 / 100] : contents.back();
            hist->GetZaxis()->SetRangeUser(lowPercentile, highPercentile);

            // Create a canvas to draw the histogram
            TCanvas* canvas = new TCanvas("canvas", "Canvas", 1200, 1600);
            hist->Draw("COLZ");
            gPad->SetRightMargin(0.145);

            // Add custom visual elements
              for (int i = 0; i <= hist->GetNbinsX(); i += 8) {
                  TLine *line = new TLine(i, 0, i, hist->GetNbinsY());
                  line->SetLineStyle(1);
                  line->SetLineColor(kBlack);
                  line->SetLineWidth(2);
                  line->Draw();
              }
              for (int j = 0; j <= hist->GetNbinsY(); j += 8) {
                  TLine *line = new TLine(0, j, hist->GetNbinsX(), j);
                  line->SetLineStyle(1);
                  line->SetLineColor(kBlack);
                  line->SetLineWidth(2);
                  line->Draw();
              }

              TLine *middleLine = new TLine(hist->GetNbinsX() / 2, 0, hist->GetNbinsX() / 2, hist->GetNbinsY());
              middleLine->SetLineWidth(4);
              middleLine->Draw();

              for (int i = 0; i < hist->GetNbinsX(); i += 8) {
                  for (int j = 0; j < hist->GetNbinsY(); j += 8) {
                      int ib_board = custom_ib_board(i, j);
                      TLatex *text = new TLatex(i + 4, j + 4, Form("%d", ib_board));
                      text->SetTextAlign(22);
                      text->SetTextSize(0.022);
                      text->SetTextFont(62);
                      text->Draw();
                  }
              }
              for (int i = 0; i < hist->GetNbinsX(); i += 48) {
                  for (int j = 0; j < hist->GetNbinsY(); j += 8) {
                      int sectorNumber = custom_sector_mapping(i, j);
                      TText *text = new TText(i + 24, j + 4, Form("%d", sectorNumber));
                      text->SetTextAlign(22);
                      text->SetTextSize(0.03);
                      text->SetTextFont(62);
                      text->Draw();
                  }
              }
              hist->GetXaxis()->SetNdivisions(-512);
              hist->GetYaxis()->SetNdivisions(-32);

            
            
            // Construct the output file name
            std::string outputFileName = outputDirectory + histName + ".png";
            canvas->SaveAs(outputFileName.c_str());
            
            TH1D* projX = hist->ProjectionX((histName + "_px").c_str());

            // Customize the projection histogram if needed:
            projX->SetTitle((hist->GetTitle() + std::string(" (X Projection)")).c_str());
            projX->SetStats(kFALSE); // Remove statistics box
            projX->GetXaxis()->SetTitle("#eta");
            projX->GetYaxis()->SetTitle("Number of events");
            projX->GetYaxis()->SetTitleOffset(1.6);

            // Save the projection histogram
            TCanvas* canvasProjX = new TCanvas("canvasProjX", "ProjectionX Canvas", 800, 600);
            projX->Draw();
            std::string projXFileName = outputDirectory + histName + "_projX.png";
            canvasProjX->SaveAs(projXFileName.c_str());
            delete canvasProjX;
            delete projX; // Clean up projection histogram

            delete canvas; // Clean up
        }
        delete obj; // Clean up the object
    }

    // Plot the 1D histogram as a BAR CHART
    TCanvas* canvasTrigger = new TCanvas("canvasTrigger", "Trigger Canvas", 800, 600);
    triggerHist->SetStats(kFALSE);

    // Key change: Set bar chart style
    triggerHist->SetFillColor(kAzure-9); // Light blue fill color for bars
    triggerHist->SetBarWidth(0.8);      // Adjust bar width if needed
    triggerHist->SetBarOffset(0.1);     // Space between bars

    triggerHist->Draw("BAR");           // Draw as bar chart
    
    triggerHist->SetStats(kFALSE); // Remove the statistics box
    
    // Print out the trigger index and total events
    std::cout << "\nFinal Total Number of Events by Trigger Index:" << std::endl;
    for (int i = 1; i <= 64; ++i) {
        std::cout << "Trigger Index " << i-1 << ": " << triggerHist->GetBinContent(i) << " events" << std::endl;
    }
    
//    TLatex latex;
//    latex.SetTextSize(0.02);
//    latex.SetTextAlign(13);  // Align at top right
//    double xPos = 0.6, yPos = 0.8;
//    for (int i = 0; i <= triggerHist->GetNbinsX(); ++i) {
//         int counts = static_cast<int>(triggerHist->GetBinContent(i));
//         if (counts > 0 && triggerMap.count(i - 1) > 0) {  // Adjust index for map access
//             std::string triggerName = triggerMap[i - 1];  // Accessing map with i-1
//             std::string label = Form("%s (%d): %d", triggerName.c_str(), i - 1, counts);
//             latex.DrawLatexNDC(xPos, yPos, label.c_str());
//             yPos -= 0.03;  // Move down for the next entry
//             if (yPos < 0.2) {  // Prevent drawing below the canvas
//                 xPos = 0.55;  // Move to the second column if needed
//                 yPos = 0.95;
//             }
//         }
//     }
    std::string triggerHistFilename = outputDirectory + "TriggerIndexHistogram.png";
    canvasTrigger->SaveAs(triggerHistFilename.c_str());

    // Close the file
    file->Close();
    delete file;
    delete triggerHist; // Clean up the 1D histogram
}
