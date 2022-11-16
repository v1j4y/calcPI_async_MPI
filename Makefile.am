AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = -I m4 --install

VERSION_MAJOR   = @VERSION_MAJOR@
VERSION_MINOR   = @VERSION_MINOR@
VERSION_PATCH   = @VERSION_PATCH@

pkgconfigdir = $(libdir)/pkgconfig
pkgconfig_DATA = pkgconfig/mpiexamplelib.pc

src_main_c = src/main.c src/calculate.c
src_main_h = src/calculate.h

src/calculate.tangled: $(abs_top_builddir)/org/calculate.org
	top_builddir=$(abs_top_builddir) srcdir=$(abs_srcdir) $(top_builddir)/tools/missing bash $(top_builddir)/tools/tangle.sh $(top_srcdir)/org/calculate.org

src/main.tangled: $(abs_top_builddir)/org/main.org 
	top_builddir=$(abs_top_builddir) srcdir=$(abs_srcdir) $(top_builddir)/tools/missing bash $(top_builddir)/tools/tangle.sh $(top_srcdir)/org/main.org

src/main.c: src/main.tangled src/calculate.c src/calculate.h
src/calculate.c: src/calculate.tangled
src/calculate.h: src/calculate.tangled

bin_PROGRAMS = mpiexample
mpiexample_SOURCES = $(src_main_c) $(src_main_h)
mpiexample_CPPFLAGS = $(AM_CPPFLAGS)
mpiexample_CFLAGS = $(AM_CFLAGS)
mpiexample_CFLAGS += "-I${top_srcdir}/include"  "-I${top_srcdir}/src" "-lm"
mpiexample_LDFLAGS = $(AM_LDFLAGS)

lib_LTLIBRARIES = lib/libmpiexample.la
include_HEADERS = $(src_main_h)
lib_libmpiexample_la_SOURCES = $(src_main_c) $(src_main_h)
lib_libmpiexample_la_CFLAGS = $(AM_CFLAGS)
lib_libmpiexample_la_CPPFLAGS = "-I${top_srcdir}/src"
lib_libmpiexample_la_CPPFLAGS += "-I${top_srcdir}/include"
lib_libmpiexample_la_LDFLAGS = -version-info 0:0:0
