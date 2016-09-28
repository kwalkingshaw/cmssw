import FWCore.ParameterSet.VarParsing as VarParsing

process = cms.Process("ProcessOne")

options = VarParsing.VarParsing()
options.register('destinationDatabase',
                '',
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "the destination database connection string")
options.register('destinationTag',
                '',
                VarParsing.VarParsing.multiplicity.singleton,
                VarParsing.VarParsing.varType.string,
                "the destination tag name")
options.parseArguments()
