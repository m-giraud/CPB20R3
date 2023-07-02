// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
/*   IMPLEMENTATION OF CLASS   G R I D  G L U E   */

#include "intersection.hh"
#include <vector>
#include <iterator>
#include "../gridglue.hh"
#if HAVE_MPI
#include "../common/ringcomm.hh"
#endif

#include <dune/common/unused.hh>

namespace Dune {
namespace GridGlue {

template<typename P0, typename P1>
GridGlue<P0, P1>::GridGlue(const std::shared_ptr< const GridPatch<0> >& gp0, const std::shared_ptr< const GridPatch<1> >& gp1, const std::shared_ptr<Merger>& merger)
  : patches_{gp0, gp1}, merger_(merger)
{
#if HAVE_MPI
  // if we have only seq. meshes don't use parallel glueing
  if (gp0->gridView().comm().size() == 1
      && gp1->gridView().comm().size() == 1)
    mpicomm_ = MPI_COMM_SELF;
  else
    mpicomm_ = MPI_COMM_WORLD;
#endif // HAVE_MPI
  std::cout << "GridGlue: Constructor succeeded!" << std::endl;
}

template<typename P0, typename P1>
void GridGlue<P0, P1>::build()
{
  int myrank = 0;
#if HAVE_MPI
  int commsize = 1;
  MPI_Comm_rank(mpicomm_, &myrank);
  MPI_Comm_size(mpicomm_, &commsize);
#endif // HAVE_MPI

  // clear the contents from the current intersections array
  {
    std::vector<IntersectionData> dummy(1); // we need size 1, as we always store data for the end-intersection
    intersections_.swap(dummy);
  }

  std::vector<Dune::FieldVector<ctype, dimworld> > patch0coords;
  std::vector<unsigned int> patch0entities;
  std::vector<Dune::GeometryType> patch0types;
  std::vector<Dune::FieldVector<ctype,dimworld> > patch1coords;
  std::vector<unsigned int> patch1entities;
  std::vector<Dune::GeometryType> patch1types;

  /*
   * extract global surface patchs
   */

  // retrieve the coordinate and topology information from the extractors
  // and apply transformations if necessary
  extractGrid(patch<0>(), patch0coords, patch0entities, patch0types);
  extractGrid(patch<1>(), patch1coords, patch1entities, patch1types);

  std::cout << ">>>> rank " << myrank << " coords: "
            << patch0coords.size() << " and " << patch1coords.size() << std::endl;
  std::cout << ">>>> rank " << myrank << " entities: "
            << patch0entities.size() << " and " << patch1entities.size() << std::endl;
  std::cout << ">>>> rank " << myrank << " types: "
            << patch0types.size() << " and " << patch1types.size() << std::endl;

#ifdef WRITE_TO_VTK
  const char prefix[] = "GridGlue::Builder::build() : ";
  char patch0surf[256];
  sprintf(patch0surf, "/tmp/vtk-patch0-test-%i", myrank);
  char patch1surf[256];
  sprintf(patch1surf, "/tmp/vtk-patch1-test-%i", myrank);

  // std::cout << prefix << "Writing patch0 surface to '" << patch0surf << ".vtk'...\n";
  // VtkSurfaceWriter vtksw(patch0surf);
  // vtksw.writeSurface(patch0coords, patch0entities, grid0dim, dimworld);
  // std::cout << prefix << "Done writing patch0 surface!\n";

  // std::cout << prefix << "Writing patch1 surface to '" << patch1surf << ".vtk'...\n";
  // vtksw.setFilename(patch1surf);
  // vtksw.writeSurface(patch1coords, patch1entities, grid1dim, dimworld);
  // std::cout << prefix << "Done writing patch1 surface!\n";
#endif // WRITE_TO_VTK

  // we start with an empty set
  index__sz = 0;

#if HAVE_MPI
  if (commsize > 1)
  {
    // setup parallel indexset
    patch0_is_.beginResize();
    patch1_is_.beginResize();
  }

  auto op =
    [&](
      const int mergingrank,
      const std::vector<Dune::FieldVector<ctype,dimworld> >& remotePatch0coords,
      const std::vector<unsigned int>& remotePatch0entities,
      const std::vector<Dune::GeometryType>& remotePatch0types,
      const std::vector<Dune::FieldVector<ctype,dimworld> >& remotePatch1coords,
      const std::vector<unsigned int>& remotePatch1entities,
      const std::vector<Dune::GeometryType>& remotePatch1types
      )
    {
      if (remotePatch1entities.size() > 0 && patch0entities.size() > 0)
        mergePatches(patch0coords, patch0entities, patch0types, myrank,
          remotePatch1coords, remotePatch1entities, remotePatch1types, mergingrank);
      if (mergingrank != myrank &&
        remotePatch0entities.size() > 0 && patch1entities.size() > 0)
        mergePatches(remotePatch0coords, remotePatch0entities, remotePatch0types, mergingrank,
          patch1coords, patch1entities, patch1types, myrank);
    };
  Parallel::MPI_AllApply(mpicomm_, op,
    patch0coords, patch0entities, patch0types,
    patch1coords, patch1entities, patch1types
    );

  if (commsize > 1)
  {
    // finalize ParallelIndexSet & RemoteIndices
    patch0_is_.endResize();
    patch1_is_.endResize();

    // setup remote index information
    remoteIndices_.setIncludeSelf(true);
    // #warning add list of neighbors ...
    remoteIndices_.setIndexSets(patch0_is_, patch1_is_, mpicomm_) ;
    remoteIndices_.rebuild<true/* all indices are public */>();

    // DEBUG Print all remote indices
#ifdef DEBUG_GRIDGLUE_PARALLELMERGE
    for (auto it = remoteIndices_.begin(); it != remoteIndices_.end(); it++)
    {
      std::cout << myrank << "\tri-list\t" << it->first << std::endl;
      for (auto xit = it->second.first->begin(); xit != it->second.first->end(); ++xit)
        std::cout << myrank << "\tri-list 1 \t" << it->first << "\t" << *xit << std::endl;
      for (auto xit = it->second.second->begin(); xit != it->second.second->end(); ++xit)
        std::cout << myrank << "\tri-list 2 \t" << it->first << "\t" << *xit << std::endl;
    }
#endif
  }
#else // HAVE_MPI

  if (patch1entities.size() > 0 && patch0entities.size() > 0)
  {
    mergePatches(patch0coords, patch0entities, patch0types, myrank,
      patch1coords, patch1entities, patch1types, myrank);
#ifdef CALL_MERGER_TWICE
    mergePatches(patch0coords, patch0entities, patch0types, myrank,
      patch1coords, patch1entities, patch1types, myrank);
#endif
  }

#endif // HAVE_MPI

}

template<typename T>
void printVector(const std::vector<T> & v, std::string name, int rank)
{
  std::cout << rank << ": " << name << std::endl;
  for (size_t i=0; i<v.size(); i++)
  {
    std::cout << v[i] << "   ";
  }
  std::cout << std::endl;
}

template<typename P0, typename P1>
void GridGlue<P0, P1>::mergePatches(
  const std::vector<Dune::FieldVector<ctype,dimworld> >& patch0coords,
  const std::vector<unsigned int>& patch0entities,
  const std::vector<Dune::GeometryType>& patch0types,
  const int patch0rank,
  const std::vector<Dune::FieldVector<ctype,dimworld> >& patch1coords,
  const std::vector<unsigned int>& patch1entities,
  const std::vector<Dune::GeometryType>& patch1types,
  const int patch1rank)
{

  // howto handle overlap etc?

  int myrank = 0;
#if HAVE_MPI
  int commsize = 1;
  MPI_Comm_rank(mpicomm_, &myrank);
  MPI_Comm_size(mpicomm_, &commsize);
#endif // HAVE_MPI

  // which patches are local?
  const bool patch0local = (myrank == patch0rank);
  const bool patch1local = (myrank == patch1rank);

  // remember the number of previous remote intersections
  const unsigned int offset = intersections_.size()-1;

  std::cout << myrank
            << " GridGlue::mergePatches : rank " << patch0rank << " / " << patch1rank << std::endl;

  // start the actual build process
  merger_->build(patch0coords, patch0entities, patch0types,
                 patch1coords, patch1entities, patch1types);

  // append to intersections list
  intersections_.resize(merger_->nSimplices() + offset + 1);
  for (unsigned int i = 0; i < merger_->nSimplices(); ++i)
    intersections_[offset + i] = IntersectionData(*this, i, offset, patch0local, patch1local);

  index__sz = intersections_.size() - 1;

  std::cout << myrank
            << " GridGlue::mergePatches : "
            << "The number of remote intersections is " << intersections_.size()-1 << std::endl;

#if 0
  printVector(patch0coords,"patch0coords",myrank);
  printVector(patch0entities,"patch0entities",myrank);
  printVector(patch0types,"patch0types",myrank);
  printVector(patch1coords,"patch1coords",myrank);
  printVector(patch1entities,"patch1entities",myrank);
  printVector(patch1types,"patch1types",myrank);
#endif

#if HAVE_MPI
  if (commsize > 1)
  {
    // update remote index sets
    assert(Dune::RESIZE == patch0_is_.state());
    assert(Dune::RESIZE == patch1_is_.state());

    for (unsigned int i = 0; i < merger_->nSimplices(); i++)
    {
      // #warning only handle the newest intersections / merger info
      const IntersectionData & it = intersections_[i];
      GlobalId gid(patch0rank, patch1rank, i);
      if (it.template local<0>())
      {
        Dune::PartitionType ptype = patch<0>().element(it.template index<0>()).partitionType();
        patch0_is_.add (gid, LocalIndex(offset+i, ptype) );
      }
      if (it.template local<1>())
      {
        Dune::PartitionType ptype = patch<1>().element(it.template index<1>()).partitionType();
        patch1_is_.add (gid, LocalIndex(offset+i, ptype) );
      }
    }
  }
#endif // HAVE_MPI

  // cleanup the merger
  merger_->clear();
}

template<typename P0, typename P1>
template<typename Extractor>
void GridGlue<P0, P1>::extractGrid (const Extractor & extractor,
                                    std::vector<Dune::FieldVector<ctype, dimworld> > & coords,
                                    std::vector<unsigned int> & entities,
                                    std::vector<Dune::GeometryType>& geometryTypes) const
{
  std::vector<typename Extractor::Coords> tempcoords;
  std::vector<typename Extractor::VertexVector> tempentities;

  extractor.getCoords(tempcoords);
  coords.clear();
  coords.reserve(tempcoords.size());

  for (unsigned int i = 0; i < tempcoords.size(); ++i)
  {
    assert(int(dimworld) == int(Extractor::dimworld));
    coords.push_back(Dune::FieldVector<ctype, dimworld>());
    for (size_t j = 0; j <dimworld; ++j)
      coords.back()[j] = tempcoords[i][j];
  }

  extractor.getFaces(tempentities);
  entities.clear();

  for (unsigned int i = 0; i < tempentities.size(); ++i) {
    for (unsigned int j = 0; j < tempentities[i].size(); ++j)
      entities.push_back(tempentities[i][j]);
  }

  // get the list of geometry types from the extractor
  extractor.getGeometryTypes(geometryTypes);

}

template<typename P0, typename P1>
template<class DataHandleImp, class DataTypeImp>
void GridGlue<P0, P1>::communicate(
  Dune::GridGlue::CommDataHandle<DataHandleImp,DataTypeImp> & data,
  Dune::InterfaceType iftype, Dune::CommunicationDirection dir) const
{
  typedef Dune::GridGlue::CommDataHandle<DataHandleImp,DataTypeImp> DataHandle;
  typedef typename DataHandle::DataType DataType;

#if HAVE_MPI

  if (mpicomm_ != MPI_COMM_SELF)
  {
    /*
     * P A R A L L E L   V E R S I O N
     */
    // setup communication interfaces
    Dune::dinfo << "GridGlue: parallel communication" << std::endl;
    typedef Dune::EnumItem <Dune::PartitionType, Dune::InteriorEntity> InteriorFlags;
    typedef Dune::EnumItem <Dune::PartitionType, Dune::OverlapEntity>  OverlapFlags;
    typedef Dune::EnumRange <Dune::PartitionType, Dune::InteriorEntity, Dune::GhostEntity>  AllFlags;
    Dune::Interface interface;
    assert(remoteIndices_.isSynced());
    switch (iftype)
    {
    case Dune::InteriorBorder_InteriorBorder_Interface :
      interface.build (remoteIndices_, InteriorFlags(), InteriorFlags() );
      break;
    case Dune::InteriorBorder_All_Interface :
      if (dir == Dune::ForwardCommunication)
        interface.build (remoteIndices_, InteriorFlags(), AllFlags() );
      else
        interface.build (remoteIndices_, AllFlags(), InteriorFlags() );
      break;
    case Dune::Overlap_OverlapFront_Interface :
      interface.build (remoteIndices_, OverlapFlags(), OverlapFlags() );
      break;
    case Dune::Overlap_All_Interface :
      if (dir == Dune::ForwardCommunication)
        interface.build (remoteIndices_, OverlapFlags(), AllFlags() );
      else
        interface.build (remoteIndices_, AllFlags(), OverlapFlags() );
      break;
    case Dune::All_All_Interface :
      interface.build (remoteIndices_, AllFlags(), AllFlags() );
      break;
    default :
      DUNE_THROW(Dune::NotImplemented, "GridGlue::communicate for interface " << iftype << " not implemented");
    }

    // setup communication info (class needed to tunnel all info to the operator)
    typedef Dune::GridGlue::CommInfo<GridGlue,DataHandleImp,DataTypeImp> CommInfo;
    CommInfo commInfo;
    commInfo.dir = dir;
    commInfo.gridglue = this;
    commInfo.data = &data;

    // create communicator
    Dune::BufferedCommunicator bComm ;
    bComm.template build< CommInfo >(commInfo, commInfo, interface);

    // do communication
    // choose communication direction.
    if (dir == Dune::ForwardCommunication)
      bComm.forward< Dune::GridGlue::ForwardOperator >(commInfo, commInfo);
    else
      bComm.backward< Dune::GridGlue::BackwardOperator >(commInfo, commInfo);
  }
  else
#endif // HAVE_MPI
  {
    /*
     * S E Q U E N T I A L   V E R S I O N
     */
    Dune::dinfo << "GridGlue: sequential fallback communication" << std::endl;

    // get comm buffer size
    int ssz = size() * 10;       // times data per intersection
    int rsz = size() * 10;

    // allocate send/receive buffer
    auto sendbuffer = std::make_unique<DataType[]>(ssz);
    auto receivebuffer = std::make_unique<DataType[]>(rsz);

    // gather
    Dune::GridGlue::StreamingMessageBuffer<DataType> gatherbuffer(sendbuffer.get());
    for (const auto& in : intersections(*this))
    {
      /*
         we need to have to variants depending on the communication direction.
       */
      if (dir == Dune::ForwardCommunication)
      {
        /*
           dir : Forward (grid0 -> grid1)
         */
        if (in.self())
        {
          data.gather(gatherbuffer, in.inside(), in);
        }
      }
      else         // (dir == Dune::BackwardCommunication)
      {
        /*
           dir : Backward (grid1 -> grid0)
         */
        if (in.neighbor())
        {
          data.gather(gatherbuffer, in.outside(), in.flip());
        }
      }
    }

    assert(ssz == rsz);
    for (int i=0; i<ssz; i++)
      receivebuffer[i] = sendbuffer[i];

    // scatter
    Dune::GridGlue::StreamingMessageBuffer<DataType> scatterbuffer(receivebuffer.get());
    for (const auto& in : intersections(*this))
    {
      /*
         we need to have to variants depending on the communication direction.
       */
      if (dir == Dune::ForwardCommunication)
      {
        /*
           dir : Forward (grid0 -> grid1)
         */
        if (in.neighbor())
          data.scatter(scatterbuffer, in.outside(), in.flip(),
                       data.size(in));
      }
      else         // (dir == Dune::BackwardCommunication)
      {
        /*
           dir : Backward (grid1 -> grid0)
         */
        if (in.self())
          data.scatter(scatterbuffer, in.inside(), in,
                       data.size(in));
      }
    }
  }
}

} // end namespace GridGlue
} // end namespace Dune
