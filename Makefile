all:
	cd src && $(MAKE)

clean: clean-code clean-work

clean-code:
	cd src && $(MAKE) clean
	cd doc && $(MAKE) clean
	cd work/latex && $(MAKE) clean

clean-work:
	rm -f work/lw/config_longwave_*.nam
	rm -f work/sw/config_shortwave_*.nam
	rm -f work/cloud/config_rt_*.nam

clean-autosaves:
	rm -f *~ */*~ */*/*~

# Not compatible with git:
# Convert working directory name from /my/dir/ckdmip-1.0 to '#define
# CKDMIP_VERSION "1.0"' and put it in src/ckdmip_version.h
#ckdmip_version.h:
#	basename `pwd` | awk -F- '{print "#define CKDMIP_VERSION \""$$NF"\""}' > src/ckdmip_version.h

.PHONY: clean-code clean-work clean clean-autosaves
