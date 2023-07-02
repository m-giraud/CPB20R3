// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#include "config.h"

#include <array>
#include <vector>

#define DEBUG_GRIDGLUE_PARALLELMERGE 1

#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid-glue/common/ringcomm.hh>

using namespace Dune;
// using namespace Dune::GridGlue;

#if HAVE_MPI
void eh( MPI_Comm *comm, int *err, ... )
{
  int len = 1024;
  char error_txt[len];

  MPI_Error_string(*err, error_txt, &len);
  assert(len <= 1024);
  DUNE_THROW(Dune::Exception, "MPI ERROR -- " << error_txt);
}
#endif // HAVE_MPI

int main(int argc, char *argv[])
{
  auto & mpihelper = Dune::MPIHelper::instance(argc, argv);
  Dune::dinfo.attach(std::cout);
  int rank = 0;

#if HAVE_MPI
  MPI_Errhandler errhandler;
  MPI_Comm_create_errhandler(eh, &errhandler);
  MPI_Comm_set_errhandler(MPI_COMM_WORLD, errhandler);
  rank = mpihelper.rank();
#endif

  std::vector<bool> seen(mpihelper.size(), false);
  std::vector<int> data1({mpihelper.rank(), 1, 2, 3});
  using Vec = Dune::FieldVector<double,2>;
  std::vector<Vec> data2(1000*(mpihelper.rank()+1), Vec(0.0));
  auto op =
    [&](int remote, const std::vector<int>&v1,
      const std::vector<Vec>&v2){
      std::cout << rank << " received " << v1[0] << "/" << v2.size()
                << " from " << remote
                << std::endl;
      assert(v1[0] == remote);
      assert(v2.size() == 1000*(remote+1));
      seen[remote] = true;
  };

  int sz; MPI_Comm_size(MPI_COMM_WORLD, &sz);
  std::cout << "SIZE: " << sz << std::endl;
  std::cout << rank << ": COMM " << MPI_COMM_WORLD << std::endl;

  Dune::Parallel::MPI_AllApply(MPI_COMM_WORLD, op, data1, data2);

  for (int i=0; i<mpihelper.size(); i++)
    assert(seen[i] == true);

  return 0;
}
