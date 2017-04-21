from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_("General")
config.General.requestName = 'Bfinder_2012B'
config.General.workArea = 'crab_projects_Bfinder'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Bfinder_MEX_cfg.py'

config.section_("Data")
config.Data.inputDataset = '/MuOnia/Run2012B-22Jan2013-v1/AOD'
config.Data.inputDBS =	'global'
config.Data.splitting = 'LumiBased'
config.Data.unitsPerJob = 4
#config.Data.totalUnits = 90000
config.Data.publication = False
config.Data.publishDBS = 'phys03'
config.Data.outputDatasetTag = 'CRAB3_Bfinder'
config.Data.outLFNDirBase = '/store/user/'+getUsernameFromSiteDB()+'/Bfinder/'

config.section_("Site")
config.Site.storageSite = 'T2_RU_JINR'

if __name__ == '__main__':
    print 'multisubmit.\nunitsPerJob ~ 4 for maximum splitting\nHave you done scram b -j8?!'
    import sys
    from CRABAPI.RawCommand import crabCommand
    from httplib import HTTPException
    #from copy import deepcopy
    global_lumi_mask_2012 = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.txt'
    dset = ''
    task = 'Bc_v2'
    units_per_job = 90
    run_range=''
    # print 'aaa', sys.argv, len(sys.argv)
    
    if len(sys.argv) >= 2:
        n = int(sys.argv[1])
        year = '2012'
        #
        if n == 1:
            dset = '/DoubleMuParked/Run2012A'
            config.Site.ignoreGlobalBlacklist = True
            # config.Site.whitelist += ['T2_IN_TIFR', 'T2_US_Wisconsin', 'T3_TW_NTU_HEP']
        if n == 2:
            dset = '/MuOniaParked/Run2012B'
        if n == 3:
            dset = '/MuOniaParked/Run2012C'
        if n == 4:
            dset = '/MuOniaParked/Run2012D'
    
    def submit(cfg):
		try:
			print crabCommand('submit', config = cfg)
		except HTTPException, hte:
			print hte.headers
    
    # print 'Bfinder_' + task + '_parked_' + dset[-1]
    #
    config.General.requestName = 'Bfinder_' + task + '_parked_' + dset[-1]
    config.General.workArea = 'crab_projects_Bfinder_' + task
    config.Data.inputDataset = dset + '-22Jan2013-v1/AOD'
    #
    lumi_mask = global_lumi_mask_2012 #if year == '2012' else global_lumi_mask_2011
    #
    config.Data.unitsPerJob = units_per_job
    config.Data.lumiMask = lumi_mask
    if run_range: config.Data.runRange = run_range
    submit(config)


