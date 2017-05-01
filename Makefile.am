AUTOMAKE_OPTIONS = foreign
SUBDIRS = src
#VPATH = src
bin_SCRIPTS = Multiwfn
CLEANFILES  = $(bin_SCRIPTS)
EXTRA_DIST  = Multiwfn.01 Multiwfn.02

do_subst = ( sed \
  -e 's,[@]pkglibdir[@],$(pkglibdir),g' \
  )

Multiwfn: Multiwfn.in Makefile
	$(do_subst) < Multiwfn.01 > Multiwfn
	cat settings.ini >> Multiwfn
	cat Multiwfn.02 >> Multiwfn
	chmod +x Multiwfn
