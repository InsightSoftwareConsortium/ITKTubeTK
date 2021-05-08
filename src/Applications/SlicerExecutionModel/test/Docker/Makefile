#
# Parameters
#

# Name of the docker executable
DOCKER = docker

# DockerHub organization and repository to pull/push the images from/to
ORG = slicer
REPO = slicerexecutionmodel

DIRECTORIES = \
  Docker-ITK-master_USE_SYSTEM_LIBRARIES-OFF \
  Docker-ITK-v4.13.0_USE_SYSTEM_LIBRARIES-OFF \
  Docker-ITK-v4.10.1_USE_SYSTEM_LIBRARIES-OFF \
  Docker-ITK-v4.8.0_USE_SYSTEM_LIBRARIES-OFF

# On CircleCI, do not attempt to delete container
# See https://circleci.com/docs/docker-btrfs-error/
RM = --rm
ifeq ("$(CIRCLECI)", "true")
	RM =
endif

#
# images: This target builds all IMAGES (because it is the first one, it is built by default)
#
images: $(DIRECTORIES)

#
# display
#
display_images:
	for directory in $(DIRECTORIES); do echo $$directory | rev | cut -d'/' -f 1 | rev | cut -d'-' -f 2- | tr '[:upper:]' '[:lower:]'; done

$(VERBOSE).SILENT: display_images

#
# build implicit rule
#
$(DIRECTORIES): %: ../%/Dockerfile
	$(eval DIR := $@)
	$(eval TAG := $(shell echo $(DIR) | rev | cut -d'/' -f 1 | rev | cut -d'-' -f 2- | tr '[:upper:]' '[:lower:]'))
	$(eval IMAGEID := $(shell $(DOCKER) images -q $(ORG)/$(REPO):$(TAG)))
	$(DOCKER) build -t $(ORG)/$(REPO):$(TAG) \
		--build-arg IMAGE=$(ORG)/$(REPO):$(TAG) \
		--build-arg VCS_REF=`git rev-parse --short HEAD` \
	  --build-arg VCS_URL=`git config --get remote.origin.url` \
		--build-arg BUILD_DATE=`date -u +"%Y-%m-%dT%H:%M:%SZ"` \
		../$@
	CURRENT_IMAGEID=$$($(DOCKER) images -q $(ORG)/$(REPO)) && \
	if [ -n "$(IMAGEID)" ] && [ "$(IMAGEID)" != "$$CURRENT_IMAGEID" ]; then $(DOCKER) rmi "$(IMAGEID)"; fi

#
# run implicit rule
#
.SECONDEXPANSION:
$(addsuffix .run,$(DIRECTORIES)):
	$(eval DIR := $(basename $@))
	$(eval TAG := $(shell echo $(DIR) | rev | cut -d'/' -f 1 | rev | cut -d'-' -f 2- | tr '[:upper:]' '[:lower:]'))
	$(DOCKER) run -v $$(pwd)/../..:/usr/src/SlicerExecutionModel -ti $(RM) $(ORG)/$(REPO):$(TAG) bash

#
# push implicit rule
#
.SECONDEXPANSION:
$(addsuffix .push,$(DIRECTORIES)): $$(basename $$@)
	$(eval DIR := $(basename $@))
	$(eval TAG := $(shell echo $(DIR) | rev | cut -d'/' -f 1 | rev | cut -d'-' -f 2- | tr '[:upper:]' '[:lower:]'))
	$(DOCKER) push $(ORG)/$(REPO):$(TAG)

push: $(addsuffix .push,$(DIRECTORIES))

.PHONY: images display_images $(DIRECTORIES) $(addsuffix .run,$(DIRECTORIES)) $(addsuffix .push,$(DIRECTORIES))
