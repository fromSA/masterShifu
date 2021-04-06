# Change the configuration here.
# Include your useid/name as part of IMAGENAME to avoid conflicts
IMAGENAME = docker-fhe009-masters
COMMAND   = bash
JUP	  = ./jupscript
DISKS     = -v $(PWD)/code:/project
USERID    = $(shell id -u)
GROUPID   = $(shell id -g)
USERNAME  = $(shell whoami)
PORT      = -p 9998:9998 -p 9997:9997
NETWORK   = --network host
CONTNAME  = --name fhe009-docker-container
RUNTIME   =
# --runtime=nvidia
# No need to change anything below this line

# Allows you to use sshfs to mount disks
SSHFSOPTIONS = --cap-add SYS_ADMIN --device /dev/fuse

USERCONFIG   = --build-arg user=$(USERNAME) --build-arg uid=$(USERID) --build-arg gid=$(GROUPID)

# Check for detached mode
ifeq ($(DETACHED), 1)
	DETACH = -d
else
	DETACH =
endif

.docker: Dockerfile
	docker build $(USERCONFIG) -t $(USERNAME)-$(IMAGENAME) $(NETWORK) -f Dockerfile .

# Using -it for interactive use
RUNCMD=docker run $(RUNTIME) $(NETWORK) --rm --user $(USERID):$(GROUPID) $(DETACH) $(PORT) $(SSHFSOPTIONS) $(DISKS) $(CONTNAME) -it $(USERNAME)-$(IMAGENAME)

# Replace 'bash' with the command you want to do
default: .docker
	$(RUNCMD) $(COMMAND) $(JUP)

# requires CONFIG=jupyter
jupyter:
	$(RUNCMD) jupyter notebook --ip $(hostname -I) --port 9998
