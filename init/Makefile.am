AUTOMAKE_OPTIONS = foreign
bin_SCRIPTS = Multiwfn
CLEANFILES  = $(bin_SCRIPTS)
EXTRA_DIST  = Multiwfn.01 Multiwfn.02 settings.ini 

do_subst = ( sed \
  -e 's,[@]pkglibdir[@],$(pkglibdir),g' \
  )

Multiwfn: Multiwfn.01 Multiwfn.02 settings.ini Makefile
	$(do_subst) < Multiwfn.01 > Multiwfn
	cat settings.ini >> Multiwfn
	cat Multiwfn.02 >> Multiwfn
	chmod +x Multiwfn

SUBDIRS = src init
