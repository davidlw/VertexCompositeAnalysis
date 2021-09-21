import subprocess
name = "submitMC_v0.jdl"
f = open(name, "w")

command_lines = '''Universe   = vanilla
GetEnv     = True
Executable = submitMC_v0.sh
Arguments  = 000
Log        = log/submitMC_v0.$(Process).log
Output     = out/submitMC_v0.$(Process).out
Error      = err/submitMC_v0.$(Process).err

+JobFlavour           = "microcentury"
#+JobFlavour           = "workday" 
#+JobFlavour           = "longlunch" 

Queue
'''

for i in range(1, 58):
   temp = '''
Arguments  = %03d
Queue
   ''' % i
   command_lines += temp

f.write(command_lines)
f.close()
subprocess.call(["condor_submit", name]);
