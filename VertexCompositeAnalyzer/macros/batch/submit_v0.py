import subprocess
name = "submit_v0.jdl"
f = open(name, "w")

command_lines = '''Universe   = vanilla
GetEnv     = True
Executable = submit_v0.sh
Arguments  = 000
Log        = log/submit_v0.$(Process).log
Output     = out/submit_v0.$(Process).out
Error      = err/submit_v0.$(Process).err

+JobFlavour           = "microcentury"
#+JobFlavour           = "workday" 
#+JobFlavour           = "longlunch" 

Queue
'''

for i in range(1, 304):
   temp = '''
Arguments  = %03d
Queue
   ''' % i
   command_lines += temp

f.write(command_lines)
f.close()
subprocess.call(["condor_submit", name]);
