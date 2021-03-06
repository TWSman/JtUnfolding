from defs import createResponseInverse
from os.path import splitext, split
from rootpy.io import root_open
from rootpy.plotting import Hist
import JtUnfolder
import os
import sys


def get_hists(in_file, base_folder, histname, iFinder, Njets=0):
    if Njets == 0:
        return in_file.Get(
            "{}/{}/{}NFin{:02d}".format(base_folder, histname, histname, iFinder)
        )
    else:
        return [
            in_file.Get(
                "{}/{}/{}NFin{:02d}JetPt{:02d}".format(
                    base_folder, histname, histname, iFinder, iJet
                )
            )
            for iJet in range(Njets)
        ]


def get_numbers(in_file, base_folder, histname, iFinder, Njets):
    return [
        x.Integral() for x in get_hists(in_file, base_folder, histname, iFinder, Njets)
    ]


def main():
    print("Number of arguments: ", len(sys.argv), "arguments.")
    filename = sys.argv[1]
    print("Input file: {}".format(filename))

    jetBinBorders = [5, 10, 20, 30, 40, 60, 80, 100, 150, 500]
    jetPtBins = [(a, b) for a, b in zip(jetBinBorders, jetBinBorders[1:])]
    JetPtCenter = [7.5, 15, 25, 35, 50, 70, 90, 125, 325]
    JetPtError = [2.5, 5, 5, 5, 10, 10, 10, 25, 175]
    Njets = len(jetBinBorders) - 1
    Njets = 9
    f = root_open(filename, "read")
    if len(sys.argv) > 2:
        filename_test = sys.argv[2]
        f_test = root_open(filename_test)
        print("Open {} for test data".format(filename_test))
    else:
        f_test = root_open(filename)
    if len(sys.argv) < 4:
        f_data = None
    else:
        filename_data = sys.argv[3]
        f_data = root_open(filename_data)
        if len(sys.argv) > 4:
            if sys.argv[4] == "MB":
                do_mb = True
                do_triggered = False
            elif sys.argv[4] == "Triggered":
                do_mb = False
                do_triggered = True
            else:
                do_mb = False
                do_triggered = False

    nR = 3
    iFinder = 0
    iMCFinder = iFinder + nR * 2

    base_folder = "AliJJetJtTask/AliJJetJtHistManager"
    base_folder_trig = "AliJJetJtTask_kEMCEJE/AliJJetJtHistManager"
    base_folder_MC = "AliJJetJtTask/AliJJetJtMCHistManager"

    if f_data:
        numberJetsDataMB = get_numbers(f_data, base_folder, "JetPtBin", iFinder, Njets)

        numberJetsDataTriggered = get_numbers(
            f_data, base_folder_trig, "JetPtBin", iFinder, Njets
        )

        hTrackJtDataMB = get_hists(
            f_data, base_folder, "JetConeJtWeightBin", iFinder, Njets
        )
        hTrackJtDataTriggered = get_hists(
            f_data, base_folder_trig, "JetConeJtWeightBin", iFinder, Njets
        )

        hTrackJt2DDataMB = get_hists(f_data, base_folder, "JtWeight2D", iFinder)
        hTrackJt2DDataTriggered = get_hists(
            f_data, base_folder_trig, "JtWeight2D", iFinder
        )

        hJetPtDataMBCoarse = Hist(jetBinBorders)
        for n, i in zip(numberJetsDataMB, range(1, Njets + 1)):
            hJetPtDataMBCoarse.SetBinContent(i, n)

        hJetPtDataTriggeredCoarse = Hist(jetBinBorders)
        for n, i in zip(numberJetsDataTriggered, range(1, Njets + 1)):
            hJetPtDataTriggeredCoarse.SetBinContent(i, n)

        hJetPtDataMB = get_hists(f_data, base_folder, "JetPt", iFinder)
        hJetPtDataTriggered = get_hists(f_data, base_folder_trig, "JetPt", iFinder)

        # Get Background distributions
        hBgJtDataMB = get_hists(f_data, base_folder, "BgJtWeightBin", iFinder, Njets)
        hBgJtDataTriggered = get_hists(
            f_data, base_folder_trig, "BgJtWeightBin", iFinder, Njets
        )

        # Get number of background jets
        hBgNumbersDataMB = get_numbers(
            f_data, base_folder, "BgTrkNumberBin", iFinder, Njets
        )

        hBgNumbersDataTriggered = get_numbers(
            f_data, base_folder_trig, "BgTrkNumberBin", iFinder, Njets
        )

    numberJetsMeas = get_numbers(f_test, base_folder, "JetPtBin", iFinder, Njets)

    numberJetsTrue = get_numbers(f_test, base_folder, "JetPtBin", iMCFinder, Njets)

    hTrackJtMeas = get_hists(f_test, base_folder, "JetConeJtWeightBin", iFinder, Njets)
    hTrackJtTrue = get_hists(
        f_test, base_folder, "JetConeJtWeightBin", iMCFinder, Njets
    )

    numberJetsMeasTrain = get_numbers(f, base_folder, "JetPtBin", iFinder, Njets)

    hTrackJtCorr2D = get_hists(f, base_folder_MC, "TrackJtCorr2D", iFinder, Njets)
    hTrackJtMeas2D = get_hists(f_test, base_folder, "JtWeight2D", iFinder)
    hTrackJtTrue2D = get_hists(f_test, base_folder, "JtWeight2D", iMCFinder)
    hTrackJtMisses2D = get_hists(f, base_folder_MC, "TrackJtMisses2D", iFinder)
    hTrackJtFakes2D = get_hists(f, base_folder, "JetConeJtUnfBg2D", iFinder)
    hTrackJtBg = get_hists(f, base_folder, "BgJtWeightBin", iFinder, Njets)

    hTrackJtBgNumber = get_numbers(f, base_folder, "BgTrkNumberBin", iFinder, Njets)
    # Get number of background jets

    hTrackMatchSuccess = get_hists(
        f, base_folder_MC, "TrackMatchSuccess", iFinder, Njets
    )

    LogBinsJt = [
        hTrackJtMeas2D.GetXaxis().GetBinLowEdge(iBin)
        for iBin in range(1, hTrackJtMeas2D.GetNbinsX() + 2)
    ]

    for h in hTrackJtCorr2D:
        print("{}".format(h.GetTitle()))
    hJetPtMeas = get_hists(f_test, base_folder, "JetPt", iFinder)
    hJetPtMeasCoarse = Hist(jetBinBorders)
    for n, i in zip(numberJetsMeas, range(1, Njets + 1)):
        hJetPtMeasCoarse.SetBinContent(i, n)

    hJetPtTrue = get_hists(f_test, base_folder, "JetPt", iMCFinder)
    hJetPtTrueCoarse = Hist(jetBinBorders)
    for n, i in zip(numberJetsTrue, range(1, Njets + 1)):
        hJetPtTrueCoarse.SetBinContent(i, n)
    LogBinsPt = [
        hJetPtTrue.GetXaxis().GetBinLowEdge(iBin)
        for iBin in range(1, hJetPtTrue.GetNbinsX() + 2)
    ]
    hJetPtResponse = get_hists(f, base_folder_MC, "JetPtCorr", iFinder)
    hJetPtResponseCoarse = get_hists(f, base_folder_MC, "JetPtCorrCoarse", iFinder)

    if False:
        TrackJtUnfolder = JtUnfolder.JtUnfolder(
            "TrackJtUnfolder",
            jetBinBorders=jetBinBorders,
            Njets=Njets,
            Iterations=5,
            Data=True,
        )
        TrackJtUnfolder.setTrackMatch(hTrackMatchSuccess)
        # TrackJtUnfolder.drawTrackMatch("TrackMatch",'single')
        TrackJtUnfolder.setPtBins(LogBinsPt)
        TrackJtUnfolder.setJtBins(LogBinsJt)
        TrackJtUnfolder.setJtMeas2D(hTrackJtMeas2D)
        TrackJtUnfolder.setJtTestMeas2D(hTrackJtMeas2D)
        TrackJtUnfolder.setJtTrue2D(hTrackJtTrue2D)
        TrackJtUnfolder.setJtTestTrue2D(hTrackJtTrue2D)
        TrackJtUnfolder.setMisses2D(hTrackJtMisses2D)
        TrackJtUnfolder.setFakes2D(hTrackJtFakes2D)
        TrackJtUnfolder.setJetPtMeas(hJetPtMeas)
        TrackJtUnfolder.setJetPtTrue(hJetPtTrue)
        TrackJtUnfolder.setJetPtMeasCoarse(hJetPtMeasCoarse)
        TrackJtUnfolder.setJetPtTrueCoarse(hJetPtTrueCoarse)
        TrackJtUnfolder.setJetPtResponse(
            createResponseInverse(hJetPtMeas, hJetPtResponse)
        )
        TrackJtUnfolder.setJetPtResponseCoarse(
            createResponseInverse(hJetPtMeasCoarse, hJetPtResponseCoarse)
        )
        TrackJtUnfolder.setNumberJetsMeas(numberJetsMeas)
        TrackJtUnfolder.setNumberJetsTrue(numberJetsTrue)
        TrackJtUnfolder.setNumberJetsTestMeas(numberJetsMeas)
        TrackJtUnfolder.setNumberJetsTestTrue(numberJetsTrue)
        TrackJtUnfolder.setNumberJetsMeasTrain(sum(numberJetsMeasTrain))
        TrackJtUnfolder.set2Dresponse(hTrackJtCorr2D)
        TrackJtUnfolder.setJtBackground(hTrackJtBg)
        TrackJtUnfolder.setJtBackgroundNumbers(hTrackJtBgNumber)
        TrackJtUnfolder.unfold()
        TrackJtUnfolder.write_files("dataUnfolded.root")
        # TrackJtUnfolder.plotResponse()
        # TrackJtUnfolder.plotJetPt()
        # TrackJtUnfolder.plotJt("PythiaTest", Rebin=4)
        return

    if f_data:
        split_file = splitext(split(filename_data)[1])
        suffix = "_unfolded"
        if do_mb:
            suffix += "_MB"
        if do_triggered:
            suffix += "_triggered"
        output_file = "output/" + split_file[0] + suffix + split_file[1]

        if not os.path.exists('output'):
            os.makedirs('output')

        if do_mb:
            MBDataJtUnfolder = JtUnfolder.JtUnfolder(
                "MBDataUnfolder",
                jetBinBorders=jetBinBorders,
                Njets=Njets,
                Data=True,
                Iterations=5,
            )
            MBDataJtUnfolder.setPtBins(LogBinsPt)
            MBDataJtUnfolder.setJtBins(LogBinsJt)
            MBDataJtUnfolder.setJtMeas2D(hTrackJt2DDataMB)
            MBDataJtUnfolder.setJetPtMeas(hJetPtDataMB)
            MBDataJtUnfolder.setJetPtMeasCoarse(hJetPtDataMBCoarse)
            MBDataJtUnfolder.setNumberJetsMeas(numberJetsDataMB)
            MBDataJtUnfolder.setJtBackground(hBgJtDataMB)
            MBDataJtUnfolder.setJtBackgroundNumbers(hBgNumbersDataMB)

            MBDataJtUnfolder.setJtTrue2D(hTrackJtTrue2D)
            MBDataJtUnfolder.setJtTestTrue2D(hTrackJtTrue2D)
            MBDataJtUnfolder.setFakes2D(hTrackJtFakes2D)
            MBDataJtUnfolder.setMisses2D(hTrackJtMisses2D)
            MBDataJtUnfolder.setJetPtResponse(createResponseInverse(hJetPtMeas,
                                                                    hJetPtResponse))
            MBDataJtUnfolder.setJetPtResponseCoarse(
                createResponseInverse(hJetPtMeasCoarse, hJetPtResponseCoarse)
            )
            MBDataJtUnfolder.setNumberJetsMeasTrain(sum(numberJetsMeasTrain))
            MBDataJtUnfolder.set2Dresponse(hTrackJtCorr2D)
            MBDataJtUnfolder.unfold()
            MBDataJtUnfolder.write_files(output_file)
            # MBDataJtUnfolder.plotJt("MBDataUnfolded", Rebin=4)

        if do_triggered:
            TriggeredDataJtUnfolder = JtUnfolder.JtUnfolder(
                "TriggeredDataUnfolder",
                jetBinBorders=jetBinBorders,
                Njets=Njets,
                Data=True,
                Iterations=5,
            )
            TriggeredDataJtUnfolder.setPtBins(LogBinsPt)
            TriggeredDataJtUnfolder.setJtBins(LogBinsJt)
            TriggeredDataJtUnfolder.setJtMeas2D(hTrackJt2DDataTriggered)
            TriggeredDataJtUnfolder.setJetPtMeas(hJetPtDataTriggered)
            TriggeredDataJtUnfolder.setJetPtMeasCoarse(hJetPtDataTriggeredCoarse)
            TriggeredDataJtUnfolder.setNumberJetsMeas(numberJetsDataTriggered)
            TriggeredDataJtUnfolder.setJtBackground(hBgJtDataTriggered)
            TriggeredDataJtUnfolder.setJtBackgroundNumbers(hBgNumbersDataTriggered)

            TriggeredDataJtUnfolder.setJtTrue2D(hTrackJtTrue2D)
            TriggeredDataJtUnfolder.setJtTestTrue2D(hTrackJtTrue2D)
            TriggeredDataJtUnfolder.setFakes2D(hTrackJtFakes2D)
            TriggeredDataJtUnfolder.setMisses2D(hTrackJtMisses2D)
            TriggeredDataJtUnfolder.setJetPtResponse(createResponseInverse(hJetPtMeas,
                                                                           hJetPtResponse))
            TriggeredDataJtUnfolder.setJetPtResponseCoarse(
                createResponseInverse(hJetPtMeasCoarse, hJetPtResponseCoarse)
            )
            TriggeredDataJtUnfolder.setNumberJetsMeasTrain(sum(numberJetsMeasTrain))
            TriggeredDataJtUnfolder.set2Dresponse(hTrackJtCorr2D)
            TriggeredDataJtUnfolder.unfold()
            TriggeredDataJtUnfolder.write_files(output_file)
            # TriggeredDataJtUnfolder.plotJt("TriggeredDataUnfolded", Rebin=4)


if __name__ == "__main__":
    main()
