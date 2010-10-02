# Define the task name
TASK = pyger

# Set Flight environment to be SKA.  The other choice is TST.  Include the
# Makefile.FLIGHT make file that does most of the hard work
FLIGHT = SKA
include /proj/sot/ska/include/Makefile.FLIGHT

WWW  = $(INSTALL)/www
# Define the installed executables for the task.  This directory is reserved
# for documented tools and not for dedicated scripts, one-off codes etc
# BIN = fetch fetch_server
DOC = doc/_build/html/*

# Installed data.  These are template RDB files, mostly relevant for testing
# DATA = archfiles_def.sql filetypes.dat filetypes_all.dat task_schedule.cfg msid_bad_times.dat

# telem_archive uses a number of dedicated perl and IDL scripts
SHARE = examples 

.PHONY: doc

doc:
	cd doc; make html

install: $(TEST_DEPS)
	mkdir -p $(INSTALL_DATA)
	mkdir -p $(INSTALL_SHARE)
	mkdir -p $(INSTALL_DOC)
#
#	rsync --times --cvs-exclude $(BIN) $(INSTALL_BIN)/
	cp pyger.sh $(INSTALL_BIN)/pyger
	rsync -a --times --cvs-exclude $(SHARE) $(INSTALL_SHARE)/
	rsync --archive --times $(DOC)   $(INSTALL_DOC)/
