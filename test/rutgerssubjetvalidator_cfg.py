###############################
####### Parameters ############
###############################
from FWCore.ParameterSet.VarParsing import VarParsing

options = VarParsing ('python')

options.register('reportEvery',
    10,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.int,
    "Report every N events (default is N=10)"
)
options.register('wantSummary',
    False,
    VarParsing.multiplicity.singleton,
    VarParsing.varType.bool,
    "Print out trigger and timing summary"
)
## 'maxEvents' is already registered by the Framework, changing default value
options.setDefault('maxEvents', 10)

options.parseArguments()

import FWCore.ParameterSet.Config as cms

process = cms.Process("USER")

process.load("FWCore.MessageService.MessageLogger_cfi")
############## IMPORTANT ########################################
# If you run over many samples and you save the log, remember to reduce
# the size of the output by prescaling the report of the event number
process.MessageLogger.cerr.FwkReport.reportEvery = options.reportEvery
process.MessageLogger.cerr.default.limit = 10
#################################################################

## Options and output report
process.options   = cms.untracked.PSet(
    wantSummary = cms.untracked.bool(options.wantSummary)
)

## Events to process
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(options.maxEvents) )

## Input files
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/mc/Summer12/WWtoAnything_ptmin500_TuneZ2Star_8TeV-pythia6-tauola/AODSIM/PU_S7_START52_V9-v1/0000/FC07A5F9-8495-E111-9E8E-00304867903E.root'
    )
)

## Load module that produces GenParticlesForJets
process.load('RecoJets.Configuration.GenJetParticles_cff')

## Load RutgersSubJetProducer
process.load("RutgersSandbox.RutgersSubJetAlgorithm.ak5GenJetsRU_cfi")
process.ak5GenJetsRU.doAreaFastjet = True
#process.ak5GenJetsRU.useExplicitGhosts = False
process.ak5GenJetsRU.jetAlgorithm = 'Kt'

## Define RutgersSubJetValidator
process.subJetValidator = cms.EDAnalyzer('RutgersSubJetValidator')

## Output Module Configuration
process.out = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string(options.outputFile),
    # save only events passing the full path
    SelectEvents   = cms.untracked.PSet( SelectEvents = cms.vstring('p') ),
    outputCommands = cms.untracked.vstring('keep *' )
)

# Path definition
process.p = cms.Path(process.genParticlesForJets*process.ak5GenJetsRU*process.subJetValidator)

## Endpath
process.outpath = cms.EndPath(process.out)
