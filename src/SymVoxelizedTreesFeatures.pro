TEMPLATE = app
CONFIG += console c++20
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        createfeatures.cpp \
        main.cpp \
        r_createfeatures.cpp

## comment this out if you need a different version of R,
## and set set R_HOME accordingly as an environment variable
R_HOME = $$system(R RHOME)

## include headers and libraries for R
RCPPFLAGS =         $$system($$R_HOME/bin/R CMD config --cppflags)
RLDFLAGS =      $$system($$R_HOME/bin/R CMD config --ldflags)
RBLAS =         $$system($$R_HOME/bin/R CMD config BLAS_LIBS)
RLAPACK =       $$system($$R_HOME/bin/R CMD config LAPACK_LIBS)

## if you need to set an rpath to R itself, also uncomment
RRPATH =        -Wl,-rpath,$$R_HOME/lib

## include headers and libraries for Rcpp interface classes
## note that RCPPLIBS will be empty with Rcpp (>= 0.11.0) and can be omitted
RCPPINCL =      $$system($$R_HOME/bin/Rscript -e \"Rcpp:::CxxFlags\(\)\")
#RCPPLIBS =      $$system($$R_HOME/bin/Rscript -e \"Rcpp:::LdFlags\(\)\")

## include headers and libraries for RInside embedding classes
RINSIDEINCL =       $$system($$R_HOME/bin/Rscript -e \"RInside:::CxxFlags\(\)\")
RINSIDELIBS =       $$system($$R_HOME/bin/Rscript -e \"RInside:::LdFlags\(\)\")

RINC = $$system($$R_HOME/bin/Rscript -e \"cat(paste0(Rcpp:::Rcpp.system.file(),\'/include/\'))\")
INCLUDEPATH += $$RINC

## compiler etc settings used in default make rules
QMAKE_CXXFLAGS +=   $$RCPPWARNING $$RCPPFLAGS $$RCPPINCL $$RINSIDEINCL
#QMAKE_LIBS +=           $$RLDFLAGS $$RBLAS $$RLAPACK $$RINSIDELIBS $$RCPPLIBS
QMAKE_LIBS +=           $$RLDFLAGS $$RBLAS $$RLAPACK $$RINSIDELIBS

HEADERS += \
    createfeatures.h
