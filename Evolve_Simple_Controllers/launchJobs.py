import os;
import time;
import math;
import sys;

def ComputeRunIndex(turningAmount,targetDistance,totalObjects,expRegime,currRun):

        runIndex = 0        + turningAmount*100000;
        runIndex = runIndex + int(targetDistance/10)*10000;
        runIndex = runIndex + totalObjects*1000;
        runIndex = runIndex + expRegime*100;
        runIndex = runIndex + currRun;
        runIndex = int(runIndex);

        return( runIndex );

def AdvanceState(expRegime,currRun,turningAmount):

        expRegime = expRegime + 1;
        if ( expRegime>5 ):
                expRegime = 0;
                currRun = currRun + 1;

        if ( currRun==100 ):
                currRun = 0;
                turningAmount = turningAmount + 1;

        if ( turningAmount > 5 ):
                turningAmount = -1;

        return expRegime, currRun, turningAmount;

def NumJobsOnVACC():

	cmd = "showq | grep jbongard > allMyJobs.listing";
	os.system(cmd);

	f = open('allMyJobs.listing','r');

	numJobs = 0;
	for line in f.readlines():
		numJobs = numJobs + 1;

	f.close();

	cmd = 'rm allMyJobs.listing';
	os.system(cmd);

	time.sleep(1);

	return numJobs;

def SomeJobsAreIdle():

        cmd = "showq | grep jbongard | grep 'Deferred' > allMyJobs.listing";
        os.system(cmd);

        f = open('allMyJobs.listing','r');

        numJobs = 0;
        for line in f.readlines():
                numJobs = numJobs + 1;
	
        f.close();

        cmd = 'rm allMyJobs.listing';
        os.system(cmd);

	time.sleep(1);

        return ( numJobs>0 );

def SomeJobsAreDeferred():

        cmd = "showq | grep jbongard | grep 'Idle' > allMyJobs.listing";
        os.system(cmd);

        f = open('allMyJobs.listing','r');

        numJobs = 0;
        for line in f.readlines():
                numJobs = numJobs + 1;

        f.close();

        cmd = 'rm allMyJobs.listing';
        os.system(cmd);

        time.sleep(1);

        return ( numJobs>0 );

def Write_Preamble(f):

	numNodes = '#PBS -l nodes=1:ppn=1 \n';
	f.write(numNodes);
	wallTime = "#PBS -l walltime=10:00:00\n";
	f.write(wallTime);
	f.write("#PBS -N myjob\n");
	f.write("#PBS -j oe\n");

def Write_Modularity(f,regimeName,runIndex):

	f.write("cd $HOME/"+regimeName+"/Target \n");
	writeCmd = "./Modularity -r "+str(runIndex)+" \n";
	f.write(writeCmd);

def Write_GA(f,regimeName,numSlots,expRegime,runIndex):

		f.write("cd $HOME/"+regimeName+"/EE\n");
		commandStr =   "./GA_Camera -r "+str(runIndex);
		commandStr = commandStr + " -s "+str(numSlots);
                commandStr = commandStr + " -m 0.05";
                commandStr = commandStr + " -h 0";
                commandStr = commandStr + " -core 1";
                commandStr = commandStr + " -left 1";
                commandStr = commandStr + " -right 1";
		commandStr = commandStr + " -ua 20 -es -co ";
		commandStr = commandStr + " -to 4";
		commandStr = commandStr + " -el 1000";
		commandStr = commandStr + " -td 20";
                if ( expRegime==0 ):
                        commandStr = commandStr + " ";
                elif ( expRegime==1 ):
                        commandStr = commandStr + " -ws3 ";
                elif ( expRegime==2 ):
                        commandStr = commandStr + " -esc ";
                elif ( expRegime==3 ):
                        commandStr = commandStr + " -ws3 -esc ";

                commandStr = commandStr + " \n";
                f.write(commandStr);

def NotYetDeployed(regimeName,runIndex):

	randIndex = runIndex - int(float(runIndex)/100.0)*100;

	runFile = "/users/j/b/jbongard/"+regimeName+"/Data/"+str(randIndex)+"/runData_"+str(runIndex)+".dat";

	return( os.path.exists(runFile)==False );

def DeployJob(fileName,regimeName):

	cmd = 'chmod 777 ' + fileName + '\n';
	print cmd;
	os.system(cmd);

        cmd = 'cd /users/j/b/jbongard/'+regimeName+'/Batch_Files \n';
	cmd = cmd + "qsub myjob.script\n";
        print cmd;
        os.system(cmd);

def RemoveMyJobScripts():

        cmd = 'rm Scaff_GECCO/Batch_Files/myjob.o* \n';
        print cmd;
        os.system(cmd);

def RemoveIdleJobs():

	cmd = 'python removeIdle.py';
	print cmd;
	os.system(cmd);

def SaveOutState(expRegime,currRun,turningAmount):

	f = open('envState.dat','w');

	line = str(expRegime)+'\n';
	f.write(line);
	line = str(currRun)+'\n';
	f.write(line);
        line = str(turningAmount)+'\n';
        f.write(line);

	f.close();

def SpawnNewJob(njov):

	regimeName = 'Quad';
	numSlots = 3;

	jobsDeployed = 0;

	for r in range(0,99+1):
	
		for regime in range(0,3+1):

			fileName = "/users/j/b/jbongard/"+regimeName+"/Batch_Files/myjob.script";
		        f = open(fileName,"w");
		        Write_Preamble(f,numSlots);
			Write_MorphEngine(f,regimeName,numSlots,regime*100+r);
		        Write_GA(f,regimeName,numSlots,regime,regime*100+r);
		      	f.close();

			if ( NotYetDeployed(regimeName,regime*100+r) & (jobsDeployed<(200-njov)) ):

				DeployJob(fileName,regimeName);
				jobsDeployed = jobsDeployed + 1;

# -------------------------------- Main function

regimeName = 'Modularity';

for r in range(0,9+1):

	fileName = "/users/j/b/jbongard/"+regimeName+"/Batch_Files/myjob.script";
	f = open(fileName,"w");
	Write_Preamble(f);
	Write_Modularity(f,regimeName,r);
	f.close();

	DeployJob(fileName,regimeName);

