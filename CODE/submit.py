import csv, subprocess
parameter_file_full_path = "/storage/home/l/lfi5013/ROBUSTnew/CODE/file.csv"

with open(parameter_file_full_path,"rb") as csvfile:
        reader = csv.reader(csvfile)
        for job in reader:
		print job
            	qsub_command = """qsub -v N={0},d={1},true_s={2},k_s={3},R={4},k_n={5},sim_type={6},rep_tot={7},sseed={8}, juliaMIP_rob.pbs""".format(*job)
            	exit_status=subprocess.call(qsub_command,shell=True)
            	if exit_status is 1:
                    print "Job {0} failed to submit.".format(qsub_command)
print "Done submitting jobs!"



