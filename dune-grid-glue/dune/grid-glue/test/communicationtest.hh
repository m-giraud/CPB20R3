// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
#ifndef DUNE_GRIDGLUE_TEST_COMMUNICATIONTEST_HH
#define DUNE_GRIDGLUE_TEST_COMMUNICATIONTEST_HH

template <typename ctype, int dimw>
class CheckGlobalCoordDataHandle :
  public Dune::GridGlue::CommDataHandle< CheckGlobalCoordDataHandle<ctype, dimw>, Dune::FieldVector<ctype,dimw> >
{
public:
  template<class RISType>
  size_t size (RISType& i) const
  {
    if (i.self())
      return i.geometry().corners();
    else
      return i.geometryOutside().corners();
  }

  template<class MessageBuffer, class EntityType, class RISType>
  void gather (MessageBuffer& buff, const EntityType& e, const RISType & i) const
  {
    assert(i.self());
    for (size_t n=0; n<size(i); n++)
      buff.write(i.geometry().corner(n));
  }

  template<class MessageBuffer, class EntityType, class RISType>
  void scatter (MessageBuffer& buff, const EntityType& e, const RISType & i, size_t n)
  {
    assert(i.self());
    assert(n == size(i));
    for (size_t n=0; n<size(i); n++)
    {
      Dune::FieldVector<ctype,dimw> x;
      buff.read(x);
      assert( (x - i.geometry().corner(n)).two_norm() < 1e-6 );
    }
  }
};

template <class GlueType>
void testCommunication (const GlueType& glue)
{
  typedef typename GlueType::ctype ctype;
  static constexpr int dimw = GlueType::dimworld;
  CheckGlobalCoordDataHandle<ctype, dimw> dh;
  glue.communicate(dh, Dune::All_All_Interface, Dune::ForwardCommunication);
  glue.communicate(dh, Dune::All_All_Interface, Dune::BackwardCommunication);
}

#endif // DUNE_GRIDGLUE_TEST_COMMUNICATIONTEST_HH
