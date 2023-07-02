#ifndef DUNE_ALU3DGRIDINDEXSETS_HH
#define DUNE_ALU3DGRIDINDEXSETS_HH

#include <vector>

#include <dune/common/stdstreams.hh>
#include <dune/common/bigunsignedint.hh>
#include <dune/common/hash.hh>

#include <dune/grid/common/grid.hh>
#include <dune/grid/common/indexidset.hh>

#include "alu3dinclude.hh"
#include "topology.hh"
#include "alu3diterators.hh"

namespace Dune
{

  // External Forward Declarations
  // -----------------------------

  template<int dim, int dimworld, ALU3dGridElementType, class >
  class ALU3dGrid;

  template<int cd, int dim, class GridImp>
  class ALU3dGridEntity;



  // ALU3dGridHierarchicIndexSet
  // ---------------------------

  //! hierarchic index set of ALU3dGrid
  template<int dim, int dimworld, ALU3dGridElementType elType, class Comm >
  class ALU3dGridHierarchicIndexSet
  : public IndexSet< ALU3dGrid< dim, dimworld, elType, Comm >, ALU3dGridHierarchicIndexSet< dim, dimworld, elType, Comm > >
  {
    typedef ALU3dGridHierarchicIndexSet< dim, dimworld, elType, Comm > This;

    typedef ALU3dGrid< dim, dimworld, elType, Comm > GridType;

    friend class ALU3dGrid<dim, dimworld, elType, Comm >;

    // constructor
    ALU3dGridHierarchicIndexSet( const GridType &grid )
    : grid_( grid )
    {}

  public:
    typedef typename GridType::Traits::template Codim<0>::Entity EntityCodim0Type;

    //! return hierarchic index of given entity
    template <class EntityType>
    int index (const EntityType & ep) const
    {
      enum { cd = EntityType :: codimension };
      return index<cd>(ep);
    }

    //! return hierarchic index of given entity
    template< int codim >
    int index ( const typename GridType::Traits::template Codim< codim >::Entity &entity ) const
    {
      return GridType::getRealImplementation( entity ).getIndex();
    }

    template< class Entity >
    int subIndex ( const Entity &entity, int i, unsigned int codim ) const
    {
      return subIndex< Entity::codimension >( entity, i, codim );
    }

    //! return subIndex i of given entity for subEntity with codim
    template< int cd >
    int subIndex ( const typename GridType::Traits::template Codim< cd >::Entity &e, int i, unsigned int codim ) const
    {
      // call method subIndex on real implementation
      return GridType::getRealImplementation( e ).subIndex( i, codim );
    }

    //! return size of indexset, i.e. maxindex+1
    //! for given type, if type is not exisiting within grid 0 is returned
    int size ( GeometryType type ) const
    {
      if( elType == tetra && !type.isSimplex() ) return 0;
      if( elType == hexa  && !type.isCube() ) return 0;
      // return size of hierarchic index set
      return this->size(GridType::dimension-type.dim());
    }

    //! return size of indexset, i.e. maxindex+1
    int size ( int codim ) const
    {
      // return size of hierarchic index set
      return grid_.hierSetSize(codim);
    }

    //! deliver all geometry types used in this grid
    const std::vector<GeometryType>& geomTypes (int codim) const
    {
      return grid_.geomTypes(codim);
    }

    //! return true because all entities are contained in this set
    template <class EntityType>
    bool contains (const EntityType &) const { return true; }

  private:
    // our Grid
    const GridType & grid_;
  };

  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////
  //////////////////////////////////////////////////////////////////////

  class ALUMacroKey : public ALU3DSPACE Key4<int>
  {
    typedef int A;
    typedef ALUMacroKey ThisType;
    typedef ALU3DSPACE Key4<A> BaseType;

    public:
    ALUMacroKey() : BaseType(-1,-1,-1,-1) {}
    ALUMacroKey(const A&a,const A&b,const A&c,const A&d) : BaseType(a,b,c,d) {}
    ALUMacroKey(const ALUMacroKey & org ) : BaseType(org) {}
    ALUMacroKey & operator = (const ALUMacroKey & org )
    {
      BaseType::operator = (org);
      return *this;
    }

    bool operator == (const ALUMacroKey & org) const
    {
      return ( (this->_a == org._a) &&
               (this->_b == org._b) &&
               (this->_c == org._c) &&
               (this->_d == org._d) );
    }

    // operator < is already implemented in BaseType
    bool operator > (const ALUMacroKey & org) const
    {
      return ( (!this->operator == (org)) && (!this->operator <(org)) );
    }


    void extractKey(std::vector<int> &key) const
    {
      alugrid_assert ( key.size() == 4 );
      key[0] = this->_a;
      key[1] = this->_b;
      key[2] = this->_c;
      key[3] = this->_d;
    }

    void print(std::ostream & out) const
    {
      out << "[" << this->_a << "," << this->_b << "," << this->_c << "," << this->_d << "]";
    }

    inline friend std::size_t hash_value(const ALUMacroKey& arg)
    {
      std::size_t seed = 0;
      hash_combine(seed,arg._a);
      hash_combine(seed,arg._b);
      hash_combine(seed,arg._c);
      hash_combine(seed,arg._d);
      return seed;
    }
  };

  template <class MacroKeyImp>
  class ALUGridId
  {
    MacroKeyImp key_;
    int nChild_;
    int codimLevel_;

    // this means that only up to 300000000 entities are allowed
    static constexpr int codimOffset = 300000000 ;

  public:
    ALUGridId() : key_()
                , nChild_(-1)
                , codimLevel_(-1)
    {}

    explicit ALUGridId(const MacroKeyImp & key, const int nChild , const int codim, const int level)
      : key_(key) , nChild_(nChild)
      , codimLevel_( codim * codimOffset + level )
    {}

    ALUGridId(const ALUGridId & org )
      : key_(org.key_)
      , nChild_(org.nChild_)
      , codimLevel_(org.codimLevel_)
    {}

    ALUGridId & operator = (const ALUGridId & org )
    {
      key_         = org.key_;
      nChild_      = org.nChild_;
      codimLevel_  = org.codimLevel_;
      return *this;
    }

    bool operator == (const ALUGridId & org) const
    {
      return equals(org);
    }

    bool operator != (const ALUGridId & org) const
    {
      return ! equals(org);
    }

    bool operator <= (const ALUGridId & org) const
    {
      if(equals(org)) return true;
      else return lesser(org);
    }

    bool operator >= (const ALUGridId & org) const
    {
      if(equals(org)) return true;
      else return ! lesser(org);
    }

    bool operator < (const ALUGridId & org) const
    {
      return lesser(org);
    }

    bool operator > (const ALUGridId & org) const
    {
      return (!equals(org) && ! lesser(org));
    }

    const MacroKeyImp & getKey() const { return key_; }
    int nChild() const { return nChild_; }
    int codim() const  { return codimLevel_ / codimOffset ; }
    int level() const  { return codimLevel_ % codimOffset ; }

    bool isValid () const
    {
      return ( (nChild_ >= 0) && (codimLevel_  >= 0) );
    }

    void reset()
    {
      nChild_ = -1;
      codimLevel_  = -1;
    }

    void print(std::ostream & out) const
    {
      out << "(" << getKey() << "," << nChild_ << "," << codimLevel_ << ")";
    }

    inline friend std::size_t hash_value(const ALUGridId& arg)
    {
      std::size_t seed = hash<MacroKeyImp>()(arg.getKey());
      hash_combine(seed,arg.nChild_);
      hash_combine(seed,arg.codimLevel_);
      return seed;
    }

  protected:
    // returns true is the id is lesser then org
    bool lesser(const ALUGridId & org) const
    {
      if(getKey() < org.getKey() ) return true;
      if(getKey() > org.getKey() ) return false;
      if(getKey() == org.getKey() )
      {
        if(nChild_ == org.nChild_)
        {
          return codimLevel_ < org.codimLevel_;
        }
        else
          return nChild_ < org.nChild_;
      }
      alugrid_assert ( equals(org) );
      return false;
    }

    // returns true if this id equals org
    bool equals(const ALUGridId & org) const
    {
      return ( (getKey() == org.getKey() ) && (nChild_ == org.nChild_)
            && (codimLevel_ == org.codimLevel_) );
    }
  };

} // drop out of namespace Dune, as hash definitions have to be done in global namespace

DUNE_DEFINE_HASH(DUNE_HASH_TEMPLATE_ARGS(),DUNE_HASH_TYPE(Dune::ALUMacroKey))
DUNE_DEFINE_HASH(DUNE_HASH_TEMPLATE_ARGS(typename MacroKeyImp),DUNE_HASH_TYPE(Dune::ALUGridId<MacroKeyImp>))

namespace Dune {

  inline std::ostream& operator<< (std::ostream& s, const ALUMacroKey & key)
  {
    key.print(s);
    return s;
  }

  template <class KeyImp>
  inline std::ostream& operator<< (std::ostream& s, const ALUGridId<KeyImp> & id)
  {
    id.print(s);
    return s;
  }



  // ALU3dGlobalIdSet
  // ----------------

  template<int dim, int dimworld, ALU3dGridElementType elType, class Comm >
  class ALU3dGridGlobalIdSet
  : public IdSet< ALU3dGrid< dim, dimworld, elType, Comm >, ALU3dGridGlobalIdSet< dim, dimworld, elType, Comm >,
                  typename ALU3dGrid< dim, dimworld, elType, Comm >::Traits::GlobalIdType >,
    public ALU3DSPACE AdaptRestrictProlongType
  {
    typedef ALU3dGrid< dim, dimworld, elType, Comm > GridType;
    typedef typename GridType::HierarchicIndexSet  HierarchicIndexSetType;

    typedef ALU3dImplTraits< elType, Comm > ImplTraitsType;
    typedef typename ImplTraitsType::IMPLElementType IMPLElementType;
    typedef typename ImplTraitsType::GEOElementType GEOElementType;
    typedef typename ImplTraitsType::GEOFaceType GEOFaceType;
    typedef typename ImplTraitsType::GEOEdgeType GEOEdgeType;

    typedef typename ImplTraitsType::GitterImplType GitterImplType;

    typedef typename ImplTraitsType::HElementType HElementType;
    typedef typename ImplTraitsType::HFaceType HFaceType;
    typedef typename ImplTraitsType::HEdgeType HEdgeType;
    typedef typename ImplTraitsType::VertexType VertexType;
    typedef typename ImplTraitsType::HBndSegType HBndSegType;

    typedef EntityCount< elType > EntityCountType;

    using ALU3DSPACE AdaptRestrictProlongType::postRefinement;
    using ALU3DSPACE AdaptRestrictProlongType::preCoarsening;

  public:
    typedef typename GridType::Traits::GlobalIdType IdType;

  private:
    typedef ALUMacroKey MacroKeyType;

    typedef ALUGridId <  MacroKeyType > MacroIdType;
    enum { numCodim = 4 }; // we are always using the 3d grid here

    typedef typename GridType::Traits::template Codim<0>::Entity EntityCodim0Type;

  private:
    mutable std::map< int , IdType > ids_[ numCodim ];

    // our Grid
    const GridType & grid_;

    // the hierarchicIndexSet
    const HierarchicIndexSetType & hset_;

    int vertexKey_[4];

    enum { startOffSet_ = 0 };

  public:

    //! import default implementation of subId<cc>
    //! \todo remove after next release
    using IdSet < GridType , ALU3dGridGlobalIdSet, IdType > :: subId;

    //! create id set, only allowed for ALU3dGrid
    ALU3dGridGlobalIdSet(const GridType & grid)
      : grid_(grid), hset_(grid.hierarchicIndexSet())
    {
      if(elType == hexa)
      {
        // see dune/alugrid/impl/serial/gitter_mgb.cc
        // InsertUniqueHexa
        const int vxKey[4] = {0,1,3,4};
        for(int i=0; i<4; i++) vertexKey_[i] = vxKey[i];
      }
      else
      {
        alugrid_assert ( elType == tetra );
        // see dune/alugrid/impl/serial/gitter_mgb.cc
        // InsertUniqueTetra
        const int vxKey[4] = {0,1,2,3};
        for(int i=0; i<4; i++) vertexKey_[i] = vxKey[i];
      }

      // setup the id set
      buildIdSet();
    }

    virtual ~ALU3dGridGlobalIdSet() {}

    // update id set after adaptation
    void updateIdSet()
    {
      // to be revised
      buildIdSet();
    }

    // print all ids
    void print () const
    {
      for(int i=0 ;i<numCodim; ++i)
      {
        std::cout << "*****************************************************\n";
        std::cout << "Ids for codim " << i << "\n";
        std::cout << "*****************************************************\n";
        for(unsigned int k=0; k<ids_[i].size(); ++k)
        {
          std::cout << "Item[" << i << "," << k <<"] has id " << ids_[i][k] << "\n";
        }
        std::cout << "\n\n\n";
      }
    }

    template <class IterType>
    void checkId(const IdType & macroId, const IterType & idIter) const //int codim , unsigned int num ) const
    {

      IdType id = getId(macroId);
      for(int i=0 ;i<numCodim; ++i)
      {
        typedef typename std::map<int,IdType>::iterator IteratorType;
        IteratorType end = ids_[i].end();
        for(IteratorType it = ids_[i].begin(); it != end; ++it)
        {
          if(idIter == it) continue;
          const IdType & checkMId = (*it).second;
          IdType checkId = getId(checkMId);
          if( id == checkId )
          {
            //std::cout << "Check(codim,num = " << codim<< "," << num <<") failed for k="<<k << " codim = " << i << "\n";
            std::cout << id << " equals " << checkId << std::endl;
            std::cout << idIter->first << " != " << it->first << std::endl;
            alugrid_assert ( id != checkId );
            DUNE_THROW(GridError," " << id << " equals " << checkId << "\n");
          }
          else
          {
            bool lesser  = (id < checkId);
            bool greater = (id > checkId);
            alugrid_assert ( lesser != greater );
            if( lesser == greater )
            {
              alugrid_assert ( lesser != greater );
              DUNE_THROW(GridError," lesser equals greater of one id ");
            }
          }
        }
      }
    }

    // check id set for uniqueness
    void uniquenessCheck() const
    {
      for(int i=0 ;i<numCodim; i++)
      {
        typedef typename std::map<int,IdType>::iterator IteratorType;
        IteratorType end = ids_[i].end();
        for(IteratorType it = ids_[i].begin(); it != end; ++it)
        {
          const IdType & id = (*it).second;
          if( id.isValid() )
            checkId(id,it);
        }
      }
    }

    void setChunkSize( int chunkSize )
    {
    }

    // creates the id set
    void buildIdSet ()
    {
      for(int i=0; i<numCodim; ++i)
      {
        ids_[i].clear();
      }

      GitterImplType &gitter = grid_.myGrid();

      // all interior and border vertices
      {
        typename ALU3DSPACE AccessIterator< VertexType >::Handle fw( gitter.container() );
        for( fw.first (); !fw.done(); fw.next() )
        {
          int idx = fw.item().getIndex();
          ids_[3][idx] = buildMacroVertexId( fw.item() );
        }
      }

      // all ghost vertices
      {
        typedef typename ALU3DSPACE ALU3dGridLevelIteratorWrapper< 3, Ghost_Partition, Comm > IteratorType;
        IteratorType fw (grid_ , 0 , grid_.nlinks() );
        typedef typename IteratorType :: val_t val_t;
        for (fw.first () ; ! fw.done () ; fw.next ())
        {
          val_t & item = fw.item();
          alugrid_assert ( item.first );
          VertexType & vx = * (item.first);
          int idx = vx.getIndex();
          ids_[3][idx] = buildMacroVertexId( vx );
        }
      }

      {
        // create ids for all macro edges
        {
          typename ALU3DSPACE AccessIterator< HEdgeType >::Handle w( gitter.container() );
          for (w.first(); !w.done(); w.next())
          {
            int idx = w.item().getIndex();
            ids_[2][idx] = buildMacroEdgeId( w.item() );
            buildEdgeIds( w.item() , ids_[2][idx] , startOffSet_ );
          }
        }

        // all ghost edges
        {
          typedef typename ALU3DSPACE ALU3dGridLevelIteratorWrapper< 2, Ghost_Partition, Comm > IteratorType;
          IteratorType fw( grid_, 0, grid_.nlinks() );
          typedef typename IteratorType :: val_t val_t;
          for (fw.first () ; ! fw.done () ; fw.next ())
          {
            val_t & item = fw.item();
            alugrid_assert ( item.first );
            HEdgeType & edge = * (item.first);
            int idx = edge.getIndex();

            ids_[2][idx] = buildMacroEdgeId( edge );
            buildEdgeIds( edge , ids_[2][idx] , startOffSet_ );
          }
        }
      }


      // for all macro faces and all children
      {
        typename ALU3DSPACE AccessIterator< HFaceType >::Handle w( gitter.container() );
        for (w.first () ; ! w.done () ; w.next ())
        {
          int idx = w.item().getIndex();
          ids_[1][idx] = buildMacroFaceId( w.item() );
          buildFaceIds( w.item() , ids_[1][idx] , startOffSet_ );
        }
      }

      // all ghost faces
      if( grid_.comm().size() > 1 )
      {
        typedef typename ALU3DSPACE ALU3dGridLevelIteratorWrapper< 1, Ghost_Partition, Comm > IteratorType;
        IteratorType fw (grid_ , 0 , grid_.nlinks() );
        typedef typename IteratorType :: val_t val_t;
        for (fw.first () ; ! fw.done () ; fw.next ())
        {
          val_t & item = fw.item();
          alugrid_assert ( item.first );
          HFaceType & face = * (item.first);
          int idx = face.getIndex();
          ids_[1][idx] = buildMacroFaceId( face );
          buildFaceIds( face , ids_[1][idx] , startOffSet_ );
        }
      }

      // for all macro elements and all internal entities
      {
        typename ALU3DSPACE AccessIterator< HElementType >::Handle w( gitter.container() );
        for (w.first () ; ! w.done () ; w.next ())
        {
          int idx = w.item().getIndex();
          ids_[0][idx] = buildMacroElementId( w.item() );
          buildElementIds( w.item() , ids_[0][idx] , startOffSet_ );
        }
      }

      // all ghost elements
      if( grid_.comm().size() > 1 )
      {
        typedef typename ALU3DSPACE ALU3dGridLevelIteratorWrapper< 0, Ghost_Partition, Comm > IteratorType;
        IteratorType fw (grid_ , 0 , grid_.nlinks() );
        typedef typename IteratorType :: val_t val_t;
        for (fw.first () ; ! fw.done () ; fw.next ())
        {
          val_t & item = fw.item();
          alugrid_assert ( item.second );
          HElementType & elem = * ( item.second->getGhost().first );
          int idx = elem.getIndex();
          ids_[0][idx] = buildMacroElementId( elem );
          buildElementIds( elem , ids_[0][idx] , startOffSet_ );
        }
      }

      // check uniqueness of id only in serial, because
      // in parallel some faces and edges of ghost exists more than once
      // but have the same id, but not the same index, there for the check
      // will fail for ghost elements
      // be carefull with this check, it's complexity is O(N^2)
#ifdef ALUGRID_CHECK_GLOBALIDSET_UNIQUENESS
#warning "GlobalIdSet uniqueness check enabled!"
      uniquenessCheck();
#endif
    }

    IdType buildMacroVertexId(const VertexType & item )
    {
      int vx[4] = { item.ident(), -1, -1, -1};
      enum { codim = 3 };
      MacroKeyType key(vx[0],vx[1],vx[2],vx[3]);
      return MacroIdType(key, 1, codim, startOffSet_ );
    }

    IdType buildMacroEdgeId(const HEdgeType & item )
    {
      const GEOEdgeType & edge = static_cast<const GEOEdgeType &> (item);
      int vx[4] = {-1,-1,-1,-1};
      for(int i=0; i<2; ++i)
      {
        vx[i] = edge.myvertex(i)->ident();
      }

      enum { codim = 2 };
      MacroKeyType key(vx[0],vx[1],vx[2],vx[3]);
      return MacroIdType(key, 1, codim, startOffSet_ );
    }

    IdType buildMacroFaceId(const HFaceType & item )
    {
      const GEOFaceType & face = static_cast<const GEOFaceType &> (item);
      int vx[4] = {-1,-1,-1,-1};
      for(int i=0; i<3; ++i)
      {
        vx[i] = face.myvertex(i)->ident();
      }

      enum { codim = 1 };
      MacroKeyType key(vx[0],vx[1],vx[2],vx[3]);
      return MacroIdType(key,1, codim, startOffSet_ );
    }

    IdType buildMacroElementId(const HElementType & item )
    {
      const GEOElementType & elem = static_cast<const GEOElementType &> (item);
      int vx[4] = {-1,-1,-1,-1};
      for(int i=0; i<4; ++i)
      {
        vx[i] = elem.myvertex(vertexKey_[i])->ident();
      }
      enum { codim = 0 };
      MacroKeyType key(vx[0],vx[1],vx[2],vx[3]);
      return MacroIdType(key,1, codim, startOffSet_ );
    }

    template <int cd, class Item>
    IdType createId(const Item& item , const IdType& creatorId , int nChild )
    {
      alugrid_assert ( creatorId.isValid() );

      // we have up to 12 internal hexa faces, therefore need 100 offset
      enum { childOffSet = (dim == 2) ? 4 : ((cd == 1) && (elType == hexa)) ? 16 : 8 };
      alugrid_assert ( nChild < childOffSet );

      const int newChild   = ((creatorId.codim()+1) * creatorId.nChild() * childOffSet ) + nChild;

      IdType newId( creatorId.getKey() , newChild , cd, creatorId.level() + 1  );
      alugrid_assert( newId != creatorId );
      return newId;
    }

    // build ids for all children of this element
    void buildElementIds(const HElementType & item , const IdType & macroId , int nChild)
    {
      enum { codim = 0 };
      ids_[codim][item.getIndex()] = createId<codim>(item,macroId,nChild);

      const IdType & itemId = ids_[codim][item.getIndex()];

      buildInteriorElementIds(item,itemId);
    }

    // build ids for all children of this element
    void buildInteriorElementIds(const HElementType & item , const IdType & fatherId)
    {
      alugrid_assert ( fatherId.isValid() );

      // build id for inner vertex
      {
        const VertexType * v = item.innerVertex() ;
        // for tetras there is no inner vertex, therefore check
        if(v) buildVertexIds(*v,fatherId );
      }

      // build edge ids for all inner edges
      {
        int inneredge = startOffSet_;
        for(const HEdgeType * e = item.innerHedge () ; e ; e = e->next ())
        {
          buildEdgeIds(*e,fatherId,inneredge);
          ++inneredge;
        }
      }

      // build face ids for all inner faces
      {
        int innerface = startOffSet_;
        for(const HFaceType * f = item.innerHface () ; f ; f = f->next ())
        {
          buildFaceIds(*f,fatherId,innerface);
          ++innerface;
        }
      }

      // build ids of all children
      {
        int numChild = startOffSet_;
        for(const HElementType * child = item.down(); child; child =child->next() )
        {
          alugrid_assert ( numChild == child->nChild() );
          buildElementIds(*child, fatherId, numChild);
          ++numChild;
        }
      }
    }

    // build ids for all children of this face
    void buildFaceIds(const HFaceType & face, const IdType & fatherId , int innerFace )
    {
      enum { codim = 1 };
      ids_[codim][face.getIndex()] = createId<codim>(face,fatherId,innerFace);
      const IdType & faceId = ids_[codim][face.getIndex()];

      buildInteriorFaceIds(face,faceId);
    }

    // build ids for all children of this face
    void buildInteriorFaceIds(const HFaceType & face, const IdType & faceId)
    {
      alugrid_assert ( faceId.isValid () );

      // build id for inner vertex
      {
        const VertexType * v = face.innerVertex() ;
        if(v) buildVertexIds(*v,faceId );
      }

      // build ids for all inner edges
      {
        int inneredge = startOffSet_;
        for (const HEdgeType * e = face.innerHedge () ; e ; e = e->next ())
        {
          buildEdgeIds(*e,faceId ,inneredge );
          ++inneredge;
        }
      }

      // build ids for all child faces
      {
        int child = startOffSet_;
        for(const HFaceType * f = face.down () ; f ; f = f->next ())
        {
          alugrid_assert ( child == f->nChild()+startOffSet_);
          buildFaceIds(*f,faceId,child);
          ++child;
        }
      }
    }

    // build ids for all children of this edge
    void buildEdgeIds(const HEdgeType & edge, const IdType & fatherId , int inneredge)
    {
      enum { codim = 2 };
      ids_[codim][edge.getIndex()] = createId<codim>(edge,fatherId,inneredge);
      const IdType & edgeId = ids_[codim][edge.getIndex()];
      buildInteriorEdgeIds(edge,edgeId);
    }

    void buildInteriorEdgeIds(const HEdgeType & edge, const IdType & edgeId)
    {
      alugrid_assert ( edgeId.isValid() );

      // build id for inner vertex
      {
        const VertexType * v = edge.innerVertex() ;
        if(v) buildVertexIds(*v,edgeId );
      }

      // build ids for all inner edges
      {
        int child = startOffSet_;
        for (const HEdgeType * e = edge.down () ; e ; e = e->next ())
        {
          alugrid_assert ( child == e->nChild()+ startOffSet_ );
          buildEdgeIds(*e,edgeId , child );
          ++child;
        }
      }
    }

    // build id for this vertex
    void buildVertexIds(const VertexType & vertex, const IdType & fatherId )
    {
      enum { codim = 3 };
      // inner vertex number is 1
      ids_[codim][vertex.getIndex()] = createId<codim>(vertex,fatherId,1);
      alugrid_assert ( ids_[codim][vertex.getIndex()].isValid() );
    }

    friend class ALU3dGrid< dim, dimworld, elType, Comm >;

    const IdType & getId(const IdType & macroId) const
    {
      return macroId;
    }

  public:
    //! return global id of given entity
    template <class EntityType>
    IdType id (const EntityType & ep) const
    {
      enum { codim = ( dim == EntityType :: codimension ) ? 3 : EntityType :: codimension };
      alugrid_assert ( ids_[codim].find( hset_.index(ep) ) != ids_[codim].end() );
      const IdType & macroId = ids_[codim][hset_.index(ep)];
      alugrid_assert ( macroId.isValid() );
      return getId(macroId);
    }

    //! return global id of given entity
    template <int cd>
    IdType id (const typename GridType:: template Codim<cd> :: Entity & ep) const
    {
      const unsigned int codim = ( dim == cd ) ? 3 : cd ;
      alugrid_assert ( ids_[codim].find( hset_.index(ep) ) != ids_[codim].end() );
      const IdType & macroId = ids_[codim][hset_.index(ep)];
      alugrid_assert ( macroId.isValid() );
      return getId(macroId);
    }

    //! return subId of given entity
    IdType subId ( const EntityCodim0Type &e, int i, unsigned int codim ) const
    {
      const int hIndex = hset_.subIndex( e, i, codim );
      // idCodim is the codim used in the id storage which relates to the 3d grid
      const unsigned int idCodim = ( dim == codim ) ? 3 : codim ;
      alugrid_assert ( ids_[ idCodim ].find( hIndex ) != ids_[ idCodim ].end() );
      const IdType &macroId = ids_[ idCodim ][ hIndex ];
      alugrid_assert ( macroId.isValid() );
      return getId( macroId );
    }

    // create ids for refined elements
    int postRefinement( HElementType & item )
    {
      {
        enum { elCodim = 0 };
        const IdType & fatherId = ids_[elCodim][item.getIndex()];
        alugrid_assert ( fatherId.isValid() );
        buildInteriorElementIds(item, fatherId );
      }

      const IMPLElementType & elem = static_cast<const IMPLElementType &> (item);
      for(int i=0; i<EntityCountType::numFaces; ++i)
      {
        enum { faceCodim = 1 };
        const HFaceType & face  = *( elem.myhface( i ) );
        const IdType & id = ids_[faceCodim][face.getIndex()];
        alugrid_assert ( id.isValid() );
        buildInteriorFaceIds( face, id);
      }

      {
        for(int i=0; i<EntityCountType::numEdges; ++i)
        {
          enum { edgeCodim = 2 };
          const HEdgeType & edge  = *( elem.myhedge(i));
          const IdType & id = ids_[edgeCodim][edge.getIndex()];
          alugrid_assert ( id.isValid() );
          buildInteriorEdgeIds(edge,id);
        }
      }
      return 0;
    }

    // dummy functions
    int preCoarsening( HElementType & elem )
    {
      /*
      const IdType & fatherId = ids_[0][item.getIndex()];

      removeElementIds(item,fatherId,item.nChild());

      for(int i=0; i<EntityCountType::numFaces; ++i)
        BuildIds<dim,elType>::buildFace(*this,item,i,ids_[1]);

      for(int i=0; i<EntityCountType::numEdges; ++i)
      {
        const IMPLElementType & elem = static_cast<const IMPLElementType &> (item);
        const HEdgeType & edge  = *( elem.myhedge(i));
        const HEdgeType * child = edge.down();
        alugrid_assert ( child );
        if( ids_[2][child->getIndex() ] > zero_ ) continue;
        buildEdgeIds(edge,ids_[2][edge.getIndex()],0);
      }
#ifdef ALUGRIDDEBUG
      //uniquenessCheck();
#endif
      */
      return 0;
    }

    // dummy functions
    int preCoarsening ( HBndSegType & el ) { return 0; }

    //! prolong data, elem is the father
    int postRefinement ( HBndSegType & el ) { return 0; }

  };

  //***********************************************************
  //
  //  --LocalIdSet
  //
  //***********************************************************

  //! hierarchic index set of ALU3dGrid
    template< int dim, int dimworld, ALU3dGridElementType elType, class Comm >
  class ALU3dGridLocalIdSet
  : public IdSet< ALU3dGrid< dim, dimworld, elType, Comm >, ALU3dGridLocalIdSet< dim, dimworld, elType, Comm >, int >,
    public ALU3DSPACE AdaptRestrictProlongType
  {
    typedef ALU3dGridLocalIdSet< dim, dimworld, elType, Comm > This;

    typedef ALU3dImplTraits< elType, Comm > ImplTraitsType;
    typedef typename ImplTraitsType::HElementType HElementType;
    typedef typename ImplTraitsType::HBndSegType HBndSegType;

    typedef ALU3dGrid< dim, dimworld, elType, Comm > GridType;
    typedef typename GridType::HierarchicIndexSet HierarchicIndexSetType;

    // this means that only up to 300000000 entities are allowed
    enum { codimOffSet = 300000000 };
    typedef typename GridType::Traits::template Codim<0>::Entity EntityCodim0Type;

    // create local id set , only for the grid allowed
    ALU3dGridLocalIdSet(const GridType & grid) : hset_(grid.hierarchicIndexSet())
    {
      for( int codim = 0; codim <= GridType::dimension; ++codim )
        codimStart_[ codim ] = codim * codimOffSet;
    }

    friend class ALU3dGrid< dim, dimworld, elType, Comm >;

    // fake method to have the same method like GlobalIdSet
    void updateIdSet() {}

    using ALU3DSPACE AdaptRestrictProlongType :: postRefinement ;
    using ALU3DSPACE AdaptRestrictProlongType :: preCoarsening ;

  public:
    //! export type of id
    typedef int IdType;

    //! import default implementation of subId<cc>
    //! \todo remove after next release
    using IdSet < GridType , ALU3dGridLocalIdSet, IdType > :: subId;

    //! return global id of given entity
    template <class EntityType>
    int id (const EntityType & ep) const
    {
      enum { cd = EntityType :: codimension };
      alugrid_assert ( hset_.size(cd) < codimOffSet );
      return codimStart_[cd] + hset_.index(ep);
    }

    //! return global id of given entity
    template <int codim>
    int id (const typename GridType:: template Codim<codim> :: Entity & ep) const
    {
      //enum { cd = EntityType :: codimension };
      alugrid_assert ( hset_.size(codim) < codimOffSet );
      return codimStart_[codim] + hset_.index(ep);
    }

    //! return subId of given entity
    IdType subId ( const EntityCodim0Type &e, int i, unsigned int codim ) const
    {
      alugrid_assert ( hset_.size( codim ) < codimOffSet );
      return codimStart_[ codim ] + hset_.subIndex( e, i, codim );
    }

    // dummy functions
    int preCoarsening( HElementType & elem )  { return 0; }
    // create ids for refined elements
    int postRefinement( HElementType & item )  { return 0; }

    // dummy functions
    int preCoarsening ( HBndSegType & el ) { return 0; }

    //! prolong data, elem is the father
    int postRefinement ( HBndSegType & el ) { return 0; }

    void setChunkSize( int chunkSize ) {}

  private:
    // our HierarchicIndexSet
    const HierarchicIndexSetType & hset_;

    // store start of each codim numbers
    int codimStart_[ GridType::dimension+1 ];
  };

} // end namespace Dune

#endif // #ifndef DUNE_ALU3DGRIDINDEXSETS_HH
