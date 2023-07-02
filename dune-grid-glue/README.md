The `dune-grid-glue` module
===========================

The `dune-grid-glue` module provides infrastructure for the coupling of two unrelated Dune grids.
The coupling may be overlapping or nonoverlapping, conforming or nonconforming.
The two grids are not requested to be of the same type, and they may even be of different dimensions.

Couplings are described as sets of remote intersections.
Conceptually, these remote intersections are very close to what the regular intersections in the Dune grid interface are, with the difference that the inside and outside entities are taken from different grids.

Installation
------------

`dune-grid-glue` requires the DUNE core modules, version 2.3 or later.

Please see the [general instructions for building DUNE modules](https://www.dune-project.org/doc/installation-notes.html) for detailed instructions on how to build the module.

Development
-----------

The [development version of `dune-grid-glue`](https://gitlab.dune-project.org/extensions/dune-grid-glue) can be obtained from the DUNE project's Gitlab installation.
At the same place an [issue tracker](https://gitlab.dune-project.org/extensions/dune-grid-glue/issues) can be found.

Publications
------------

* [P. Bastian, G. Buse, O. Sander: Infrastructure for the Coupling of Dune Grids, In 'Proceedings of ENUMATH 2009', Springer, 2010, pp. 107-114](https://dx.doi.org/10.1007/978-3-642-11795-4_10)
* [C. Engwer, S. Müthing: Concepts for flexible parallel multi-domain simulations, In 'Domain Decomposition Methods in Science and Engineering XXII', Springer](https://dx.doi.org/10.1007/978-3-319-18827-0_17)

License
-------

The `dune-grid-glue` module is licensed under the GNU Lesser General Public License, version 3 or later, or the GNU General Public License, version 2, with a special runtime exception.

Please see the COPYING file for details.
