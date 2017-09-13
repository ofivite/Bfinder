import glob 
##      folder names array 
jnames  = [
'crab_Bfinder_B0_parked_parked_A', 
'crab_Bfinder_B0_parked_parked_B', 
'crab_Bfinder_B0_parked_parked_C', 
'crab_Bfinder_B0_parked_parked_D', 
]

##      corresponding short names
jnas    = ['A', 'B', 'C', 'D']

##      number of subjobs for each job (respectively)
jNs         = [191, 610, 902, 855]   

jSplit      = 1 ## if all jobs are finished, set to f.e. 200 for maximum speed,
##              else jSplit = 1 is recommended to try each file separetely.
##              if you use RemoveExistingFiles, and almost all jobs are finished 
##              you can firstly run with = 10-40 and then recover remaining 
##              files with = 1 (several times, unless "0 to download")
##
parallel    = True  ## parallelize, using parallN processes. 
parallN     = 6     ## 5-6-7-8 is ok, > 10 leads to instability and crashes

RemoveExistingFiles     = True
PrintOutSkippedExisting = False
RemoveLogsAutomatically = True
AutoUpdateScripts       = True

####################

doall = ''

if len(jnames) != len(jNs):
    print 'LENGTH (job names) != LENGTH(subjob count)'
    exit(0)

for i in range(len(jNs)):
    #
    jname   = jnames[i]
    jna     = jnas[i]
    jN      = jNs[i]
    ss = ''
    nskipped = 0

    f=open('my_%s.run'%jna, 'w')

    for jj in range(1, jN+1):
        ## firstly check if time to remove logs came 
        if jj % 800 == 0:
            if RemoveLogsAutomatically:
                f.write('\nrm ' + jname + '/crab.log\n\n')
        #
        ## then check if the file is already in 
        if  glob.glob(jname+'/results/Bfinder_2012X_%i.*'%(jj)):
            if not(RemoveExistingFiles): 
                if PrintOutSkippedExisting: print 'adding existing file %s %i'%(jname, jj)
            else: 
                if PrintOutSkippedExisting: print 'skipping existing file %s %i'%(jname, jj)
                nskipped += 1
                continue
        #
        ## if we need to download file number jj 
        ss += str(jj) + ','
        if (jj-nskipped) % jSplit == 0:
            #
            if parallel:
                f.write('crab getoutput -d ' + jname + ' --jobids=' + ss[:-1] + ' & \n')
                if int((jj - nskipped)/jSplit) % parallN == 0:
                    f.write('\nwait\n\n')
            #
            else:
                f.write('crab getoutput -d ' + jname + ' --jobids=' + ss[:-1] + '\n\n')
            #
            ss = ''
    #
    if ss:
        f.write('crab getoutput -d ' + jname + ' --jobids=' + ss[:-1] + '\n\n')
        ss = ''
    #
    f.write('\necho "\ndataset %s DONE \n"\n'%jna)
    #
    if parallel: 
        f.write('wait\n')
    #
    if RemoveLogsAutomatically:
        f.write('rm ' + jname + '/crab.log\n')
    #
    if AutoUpdateScripts:
        f.write('python MakeGetoutputScript.py\n')
    #
    if (RemoveExistingFiles and not(AutoUpdateScripts)): 
        f.write('\necho "\n\nexecution finished; you may want to update my_X file(s) by rerunning\npython MakeGetoutputScript.py\n"\n')
    #
    f.close()
    print jna, ': skipped %i files, will download remaining %i files'%(nskipped, jN-nskipped)
    # print 'execute by (without quotes) ". my_%s.run"'%jna
    doall += '. my_%s.run\n'%jna


print '\nNOTE: do not forget to do cmsenv and setup crab and grid proxy !'
if not(RemoveLogsAutomatically): 
    print 'NOTE: crab logs will be quite large in the end; you may want to delete them !'
    print 'NOTE: you may do it automatically if you want, set RemoveLogsAutomatically = True'
else:
    print 'NOTE: logs will be automatically deleted'
#
print '\n****  Commands to RUN downloading:  ****\n%s'%(doall)
