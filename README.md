# dune-PFASST

a time-parallel solver for partial differential equations using the finite element method for spacial descretisation

# installation

for dune-PFASST you need the following dune-modules which you can download via git:

git clone https://gitlab.dune-project.org/core/dune-common.git

git clone https://gitlab.dune-project.org/core/dune-geometry.git

git clone https://gitlab.dune-project.org/core/dune-grid.git

git clone https://gitlab.dune-project.org/core/dune-localfunctions.git

git clone https://gitlab.dune-project.org/core/dune-istl.git

git clone https://gitlab.dune-project.org/staging/dune-functions.git

git clone https://gitlab.dune-project.org/staging/dune-typetree.git

git clone https://git.imp.fu-berlin.de/agnumpde/dune-fufem.git

git clone https://git.imp.fu-berlin.de/agnumpde/dune-matrix-vector.git

to install them type 

./dune-common/bin/dunecontrol all

sudo ./dune-common/bin/dunecontrol make install

to build dune-PFASST just run:

./dune/dune-common/bin/dunecontrol --opts=dune.opts --only=dune-PFASST all

in dune.opts you need to specify the location of Eigen3 and PFASST++, you can use the PFASST folder includet in dune-PFASST:

eg: 

CMAKE_FLAGS="-Dpfasst_include=something/dune-PFASST/PFASST -Deigen_include=/something/eigen"

