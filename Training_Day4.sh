### Training Day 4 ### 
## By Cole Nawrocki ##

## Topic: What happens if you run out of RAM? -- Our workstation and Rivanna ##


## Rivanna
# This is UVA's research high performance computing cluster (HPC). I am pretty sure the CS department has something similar, 
# so maybe you are familiar with this stuff through class. The cluster runs linux. It can be connected to in two ways: from 
# terminal via ssh or through a web browser. 

# Terminal via ssh

# Requesting login
ssh -Y ccn7wn@rivanna.hpc.virginia.edu

# I will be prompted for my password. This is my UVA NetBadge password.

# Now I am in a "login node" and can navigate through the filesystem. I will automatically be logged into the following location: 
# /home/ccn7wn 
pwd

# This is my "home" directory and is supposedly where I can store things long-term. The location where you are supposed to place 
# files that you are actively working on is "scratch" and these files will be deleted every 90 days. So, you work in scratch, then
# when you have a completed file, you move it over to home. But, what about sharing files? No one else can go into your home or
# scratch. We have a shared folder at the following location: 
# /project/Grainger-Lab

# Some basic commands

# "change directory" allows you to change where you are. Use ".."" as the path to go back a level.
cd <relative-path>

# "list" allows you to see the location that you are in.
ls 

# "print working directory" gives you the path to your location.
pwd 

# "remove" allows you to delete things.
rm <file> 

# "move" allows you to put files in new places. Note that this doubles as the renaming tool. 
mv <current-relative-path> <new-absolute-path>

# "make directory" allows you to make a new folder.
mkdir <relative-path>

# This function allows you to view a text file, like a README, in the terminal. 
less <relative-path>

# Example of me navigating to my scratch foler: 
cd ../.. 
cd /scratch/ccn7wn
ls 

# Note: in bash, which is the language we are using, you can add arguments after the function with a "--" or a "-". The big dash means
# you have to write the whole arguments name afterwards. The little dash means you can write the shortened one. For example, in order 
# to remove a directory, you have to add the "-r" arguemnt. 
rm -r test-dir

# Via web browser

# Just go to the following page and login: https://rivanna-portal.hpc.virginia.edu/pun/sys/dashboard/
# Now you can navigate the way you do on your laptop. 

# IMPORTANT: To use the terminal login off-grounds, you need to be connected to the UVA Anywhere VPN. You do not need to be connected 
# to use the web browser login though. 
# UVA Anywhere VPN: https://virginia.service-now.com/its?id=itsweb_kb_article&sys_id=f24e5cdfdb3acb804f32fb671d9619d0

# COMPARING THE TWO METHODS: 
# - The ssh method is better for running jobs
# - Once you get quick with the commands, the ssh method is faster
# - Uploading files and downloading with the ssh method is MUCH faster
# - The web browser method allows you to use GUIs like RStudio Server and Jupyter Lab, which makes interactive sessions MUCH easier

# Running Jobs
# The whole point of this is to be able to utilize huge amounts of RAM when necessary. We do this by running jobs on the HPC. When we 
# log in, we are just on a "login node" not a "compute node." You cannot run stuff on the login nodes. So, you have to request a compute 
# node and specify the allocations that you want by writing a script that the Slurm workload manager can recognize and then submit that 
# script. I have made a tutorial with pictures that is in the "Learning-Resources" folder on OneDrive. We will walk through that now. 
# You can also run interactive jobs, which do not require you to to submit a script. Instead, you can type line by line into the console. 
# Often, this is not as useful as it seems. 


# Using scp to upload files. 
# This is a really useful tool. Here is an example:
scp <start-location> ccn7wn@rivanna.hpc.virginia.edu:/<end-location>

# Note that to move directories, you will have to again provide the "-r" argument. 

# Sublime Text
# If you are going to be writing some scripts to submit to Rivanna, then I suggest downloading a better text editor. The macOS default is 
# TextEdit. I use Sublime Text. It can recognice code based on the file estension, which I like. 


## Our Workstation
# I built this workstation after receiving a small grant of behalf of the lab. It is for the entire lab's use, not just mine. So, if you 
# want to use it, then you are welcome to. Just make your own account. All of the Bioinformatics data is currently backed up to its 
# hard drive. When I leave, someone will have to take over the responsibility of backing up the data to the workstation. 

# Here are the specifications, so you can get an idea of its capability if you are interested in using it. 
# - OS: Ubuntu
# - Processor: Intel Core i9 13th Generation CPU
# - RAM: 128 GB
# - Storage: 1 TB

# To use the workstation remotely, you will need to download Microsoft Remote Desktop (MRD) on you Mac from the app store. Next, you go to 
# Settings -> Sharing -> Remote Desktop and turn everything on. Back on MRD, add a PC and enter the IP address: 172.27.175.55 
# Note: to do this, your Mac must be connected to the UVA Anywhere VPN if you are off-grounds. 


