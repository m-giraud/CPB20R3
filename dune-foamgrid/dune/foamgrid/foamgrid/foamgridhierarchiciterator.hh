#ifndef DUNE_FOAMGRID_HIERARCHIC_ITERATOR_HH
#define DUNE_FOAMGRID_HIERARCHIC_ITERATOR_HH

/** \file
* \brief The FoamGridHierarchicIterator class
*/

#include <stack>

#include <dune/foamgrid/foamgrid/foamgridentitypointer.hh>

namespace Dune {


//**********************************************************************
//
/** \brief Iterator over the descendants of an entity.
* \ingroup FoamGrid
Mesh entities of codimension 0 ("elements") allow to visit all entities of
codimension 0 obtained through nested, hierarchic refinement of the entity.
Iteration over this set of entities is provided by the HierarchicIterator,
starting from a given entity.
*/
template<class GridImp>
class FoamGridHierarchicIterator :
        public Dune::FoamGridEntityPointer <0,GridImp>
{
    enum {dimworld = GridImp::dimensionworld};
    enum {dimgrid  = GridImp::dimension};

    friend class FoamGridEntity<0, dimgrid,GridImp>;

    public:

        typedef typename GridImp::template Codim<0>::Entity Entity;
        typedef const FoamGridEntityImp<dimgrid, dimgrid, dimworld>* StackEntry;

    //! Constructor
    FoamGridHierarchicIterator(int maxlevel)
        : FoamGridEntityPointer<0, GridImp>(nullptr),
          maxlevel_(maxlevel), elemStack()
    {}

        //! \todo Please doc me !
        void increment()
        {
            if (elemStack.empty())
                return;

            StackEntry old_target = elemStack.top();
            elemStack.pop();

            // Traverse the tree no deeper than maxlevel
            if (old_target->level_ < maxlevel_) {

                // Load sons of old target onto the iterator stack
                if (!old_target->isLeaf()) {

                    for (size_t i=0; i<old_target->nSons(); i++)
                        elemStack.push(old_target->sons_[i]);

                }

            }

            GridImp::getRealImplementation(this->virtualEntity_).setToTarget((elemStack.empty())
                                             ? nullptr : elemStack.top());
        }


private:

    //! max level to go down
    int maxlevel_;

    /** \brief For depth-first search */
    std::stack<StackEntry> elemStack;
};


}  // end namespace Dune

#endif
