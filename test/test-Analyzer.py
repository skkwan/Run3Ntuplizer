import os
import FWCore.ParameterSet.Config as cms
from Configuration.StandardSequences.Eras import eras
process = cms.Process("L1TTauSummary", eras.Run2_2018)

#import EventFilter.L1TXRawToDigi.util as util

from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing()
options.register('runNumber', 0, VarParsing.multiplicity.singleton, VarParsing.varType.int, 'Run to analyze')
options.register('lumis', '1-max', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Lumis')
options.register('dataStream', '/ExpressPhysics/Run2015D-Express-v4/FEVT', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Dataset to look for run in')
options.register('inputFiles', [], VarParsing.multiplicity.list, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.register('inputFileList', '', VarParsing.multiplicity.singleton, VarParsing.varType.string, 'Manual file list input, will query DAS if empty')
options.register('useORCON', False, VarParsing.multiplicity.singleton, VarParsing.varType.bool, 'Use ORCON for conditions.  This is necessary for very recent runs where conditions have not propogated to Frontier')
options.parseArguments()

def formatLumis(lumistring, run) :
    lumis = (lrange.split('-') for lrange in lumistring.split(','))
    runlumis = (['%d:%s' % (run,lumi) for lumi in lrange] for lrange in lumis)
    return ['-'.join(l) for l in runlumis]

print 'Getting files for run %d...' % options.runNumber
#if len(options.inputFiles) is 0 and options.inputFileList is '' :
#    inputFiles = util.getFilesForRun(options.runNumber, options.dataStream)
#elif len(options.inputFileList) > 0 :
#    with open(options.inputFileList) as f :
#        inputFiles = list((line.strip() for line in f))
#else :
#    inputFiles = cms.untracked.vstring(options.inputFiles)
#if len(inputFiles) is 0 :
#    raise Exception('No files found for dataset %s run %d' % (options.dataStream, options.runNumber))
#print 'Ok, time to analyze'


# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.Geometry.GeometryExtended2016Reco_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration.StandardSequences.RawToDigi_Data_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_data', '')

# To get L1 CaloParams
#process.load('L1Trigger.L1TCalorimeter.caloStage2Params_cfi')
# To get CaloTPGTranscoder
#process.load('SimCalorimetry.HcalTrigPrimProducers.hcaltpdigi_cff')
#process.HcalTPGCoderULUT.LUTGenerationMode = cms.bool(False)

process.load('L1Trigger.Configuration.SimL1Emulator_cff')
process.load('L1Trigger.Configuration.CaloTriggerPrimitives_cff')

process.load('EventFilter.L1TXRawToDigi.caloLayer1Stage2Digis_cfi')

process.load('L1Trigger.L1TCaloSummary.uct2016EmulatorDigis_cfi')

process.load("L1Trigger.Run3Ntuplizer.l1TRun3Ntuplizer_cfi")
process.l1NtupleProducer.isData = cms.bool(False)

process.uct2016EmulatorDigis.useECALLUT = cms.bool(False)
process.uct2016EmulatorDigis.useHCALLUT = cms.bool(False)
process.uct2016EmulatorDigis.useHFLUT = cms.bool(False)
process.uct2016EmulatorDigis.useLSB = cms.bool(True)
process.uct2016EmulatorDigis.verbose = cms.bool(False)
process.uct2016EmulatorDigis.ecalToken = cms.InputTag("l1tCaloLayer1Digis")
process.uct2016EmulatorDigis.hcalToken = cms.InputTag("l1tCaloLayer1Digis")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

#process.source = cms.Source("PoolSource")

#process.source.fileNames = cms.untracked.vstring( 
# "root://cmsxrootd.fnal.gov//store/mc/PhaseIITDRSpring19MiniAOD/GluGluHToTauTau_M125_14TeV_powheg_pythia8_TuneCP5/MINIAODSIM/NoPU_106X_upgrade2023_realistic_v3-v2/130000/26D2C453-C4D7-214D-902B-0BE730BA38D2.root")


#process.source.secondaryFileNames = cms.untracked.vstring( 
# "root://cmsxrootd.fnal.gov//store/mc/PhaseIITDRSpring19DR/GluGluHToTauTau_M125_14TeV_powheg_pythia8_TuneCP5/GEN-SIM-DIGI-RAW/NoPU_106X_upgrade2023_realistic_v3-v2/130000/534564BB-36D1-EE4C-9C20-6275568C1C09.root","root://cmsxrootd.fnal.gov//store/mc/PhaseIITDRSpring19DR/GluGluHToTauTau_M125_14TeV_powheg_pythia8_TuneCP5/GEN-SIM-DIGI-RAW/NoPU_106X_upgrade2023_realistic_v3-v2/130000/A9B3A43E-20FB-F54B-BC79-1E1EB489BA02.root")

#process.source.inputCommands = cms.untracked.vstring("keep *", 
#                                                     "drop patHcalDepthEnergyFractionsedmValueMap_packedPFCandidates_hcalDepthEnergyFractions_RECO")

#process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("1:49","1:42")

process.source = cms.Source("PoolSource",
                            #fileNames = cms.untracked.vstring(inputFiles)#,
                            #secondaryFileNames = cms.untracked.vstring(secondaryMap[options.inputFiles[0]])
                            fileNames = cms.untracked.vstring('root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18MiniAOD/DYToLL_M-50_14TeV_pythia8_pilot1/MINIAODSIM/PU200_pilot_103X_upgrade2023_realistic_v2-v1/80000/A40F504F-88AA-BA44-822B-2FF02ADFACF3.root'),
                            secondaryFileNames = cms.untracked.vstring(
                                'root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_pythia8_pilot1/FEVT/PU200_pilot_103X_upgrade2023_realistic_v2-v1/80000/257D1965-FB09-1148-9DA4-B7DAD07BB6D3.root',
                                'root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_pythia8_pilot1/FEVT/PU200_pilot_103X_upgrade2023_realistic_v2-v1/80000/1F1DD211-18B8-EF4C-A9DF-A79ED5E48983.root',
                                'root://cmsxrootd.fnal.gov//store/mc/PhaseIIMTDTDRAutumn18DR/DYToLL_M-50_14TeV_pythia8_pilot1/FEVT/PU200_pilot_103X_upgrade2023_realistic_v2-v1/80000/1BA14D04-C386-5D41-859F-818864B04297.root'
                                                                   ),
                                inputCommands = cms.untracked.vstring("keep *", 
                                                                      "drop patHcalDepthEnergyFractionsedmValueMap_packedPFCandidates_hcalDepthEnergyFractions_RECO")
)


process.source.lumisToProcess = cms.untracked.VLuminosityBlockRange("1:46", "1:5", "1:77")

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    annotation = cms.untracked.string('step2 nevts:1'),
    name = cms.untracked.string('Applications'),
    version = cms.untracked.string('$Revision: 1.19 $')
)



process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string("l1TFullEvent.root"),
    outputCommands = cms.untracked.vstring('keep *') #'keep *_*_*_L1TCaloSummaryTest')
    #outputCommands = cms.untracked.vstring('drop *', 'keep *_l1tCaloLayer1Digis_*_*, keep *_*_*_L1TCaloSummaryTest' )
)

#Output
process.TFileService = cms.Service(
    "TFileService",
    fileName = cms.string("analyzer.root")
)


#process.e = cms.EndPath(process.out)

process.p = cms.Path(process.l1tCaloLayer1Digis*process.uct2016EmulatorDigis*process.l1NtupleProducer)
process.schedule = cms.Schedule(process.p)

from PhysicsTools.PatAlgos.tools.helpers import associatePatAlgosToolsTask
associatePatAlgosToolsTask(process)

# Add early deletion of temporary data products to reduce peak memory need
from Configuration.StandardSequences.earlyDeleteSettings_cff import customiseEarlyDelete
process = customiseEarlyDelete(process)
# End adding early deletion

dump_file = open('dump.py','w')
dump_file.write(process.dumpPython())
