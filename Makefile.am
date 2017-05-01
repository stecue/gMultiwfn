AUTOMAKE_OPTIONS = foreign
SUBDIRS = src

bin_SCRIPTS = Multiwfn
CLEANFILES  = $(bin_SCRIPTS) Multiwfn.01 Multiwfn.02 settings.ini
EXTRA_DIST  = init

do_subst = ( sed \
  -e 's,[@]pkglibdir[@],$(pkglibdir),g' \
  )


Multiwfn:
	cp $(top_srcdir)/init/settings.ini ./
	cp $(top_srcdir)/init/Multiwfn.01 ./
	cp $(top_srcdir)/init/Multiwfn.02 ./
	$(do_subst) < Multiwfn.01 > Multiwfn
	cat settings.ini >> Multiwfn
	cat Multiwfn.02 >> Multiwfn
	chmod +x Multiwfn
