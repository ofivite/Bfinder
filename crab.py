from CRABClient.UserUtilities import config, getUsernameFromSiteDB
config = config()

config.section_("General")
config.General.requestName = 'Bfinder_2012B'
config.General.workArea = 'crab_projects_Bfinder'
config.General.transferLogs = True

config.section_("JobType")
config.JobType.pluginName = 'Analysis'
config.JobType.psetName = 'Bfinder_MEX_cfg.py'
#config.JobType.psetName = 'Bfinder_cfg.py'

config.section_("Data")
config.Data.inputDataset = '/MuOnia/Run2012B-22Jan2013-v1/AOD'
##config.Data.inputDataset = '/BcToJpsiPi_EtaPtFilter_8TeV-bcvegpy2/Summer12_DR53X-PU_RD2_START53_V19F-v1/AODSIM'
config.Data.inputDBS =	'global'
config.Data.splitting = 'LumiBased'
##config.Data.splitting = 'FileBased'
config.Data.unitsPerJob = 30
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
	#import argparse
	from CRABAPI.RawCommand import crabCommand
	from httplib import HTTPException
	#from copy import deepcopy

	#print sys.argv

	global_lumi_mask_2012 = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions12/8TeV/Reprocessing/Cert_190456-208686_8TeV_22Jan2013ReReco_Collisions12_JSON_MuonPhys.txt'
	global_lumi_mask_2011 = 'https://cms-service-dqm.web.cern.ch/cms-service-dqm/CAF/certification/Collisions11/7TeV/Reprocessing/Cert_160404-180252_7TeV_ReRecoNov08_Collisions11_JSON_MuonPhys_v2.txt' 

#	parser = argparse.ArgumentParser(description='Multisubmit tool for crab with useful default config')
#	parser.add_argument('dset', choices='A B C D'.split(), help='the letter for the 2012 dataset to take')
	dset = ''
#	parser.add_argument('task', help='The task name. For example PiPi_recovery')
	task = 'Xinc_v0'
#	parser.add_argument('-u', '--units_per_job', default=7, type=int)
	units_per_job = 150
#	parser.add_argument('-r', '--run-range')
	run_range=''

	if len(sys.argv) == 2:
		n = int(sys.argv[1])
		year = ''
		if n <= 7: year = '2012'
		else: year = '2011'
		if n == 1:
			dset = 'A'
			run_range = ''
		if n == 2:
			dset = 'B'
			run_range = '193834-195182'
		if n == 3:
			dset = 'B'
			run_range = '195183-196531'
		if n == 4:
			dset = 'C'
			run_range = '198002-200882'
		if n == 5:
			dset = 'C'
			run_range = '200883-203742'
		if n == 6:
			dset = 'D'
			run_range = '203777-206231'
		if n == 7:
			dset = 'D'
			run_range = '206232-208686'
		if n == 8:
			dset = 'A'
		if n == 9:
			dset = 'B'


#
#	args = parser.parse_args()

	def submit(cfg):
		try:
			print crabCommand('submit', config = cfg)
		except HTTPException, hte:
			print hte.headers

	config.General.requestName = 'Bfinder_' + task + '_'+ year + dset + (('_' + run_range) if run_range else '') 
	config.General.workArea = 'crab_projects_Bfinder_' + task
	if year == '2012': config.Data.inputDataset = '/MuOnia/Run2012' + dset + '-22Jan2013-v1/AOD'
	if year == '2011': config.Data.inputDataset = '/MuOnia/Run2011' + dset + '-12Oct2013-v1/AOD'
	
	#	parser.add_argument('-l','--lumi_mask', default=global_lumi_mask, help='local or global file to use as a lumi mask. (default: %(default)s)')
	lumi_mask = global_lumi_mask_2012 if year == '2012' else global_lumi_mask_2011

	config.Data.unitsPerJob = units_per_job
	config.Data.lumiMask = lumi_mask
	if run_range: config.Data.runRange = run_range
	submit(config)

