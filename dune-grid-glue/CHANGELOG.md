Changes in dune-grid-glue v2.7.0
================================

Major changes in dune-grid-glue v2.7.0
--------------------------------------

* Constants and type aliases including the "side" in their names like
    `Grid0View` have been deprecated in favor of template constexpr
    functions or template aliases.  Old code like

    ```c++
    using Glue = Dune::GridGlue::GridGlue<...>;
    ... Glue::grid0dim ...
    ... Glue::Grid1View ...
    ```
    should be replaced by
    ```c++
    using Glue = Dune::GridGlue::GridGlue<...>;
    ... Glue::griddim<0>() ...
    ... Glue::GridView<1> ...
    ```

Changes in dune-grid-glue v2.5.0
================================

Incompatible changes in dune-grid-glue v2.5.0
---------------------------------------------

* The autotools-based build system has been removed.

* Non-default projection directions in ContactMerge are now set and stored via
  std::function instead of Dune::VirtualFunction.

Changes in dune-grid-glue v2.4.0
================================

Incompatible changes in dune-grid-glue v2.4.0
---------------------------------------------

* All interfaces have been moved into the `Dune::GridGlue` namespace.

* The `contains()` method of the `Codim0Extractor` and `Codim1Extractor`
  classes now take a codim-0-`Entity` instead of a codim-0-`EntityPointer`.
  New code can also use a `std::function` instead of the extractor classes.

* Methods that returned an `EntityPointer` now return an `Entity` instead
  when `dune-grid-glue` is built against version 2.4 or later of the DUNE
  core modules.

Major changes in dune-grid-glue v2.4.0
--------------------------------------

* This is the first release that supports version 2.4 of the DUNE core modules.
  The older 2.3 release is also still supported.

* `dune-grid-glue` now requires C++11 support.

* Both `ContactMerge` and `OverlappingMerge` had large changes and were
  partially rewritten.  The newer version should both be faster and hopefully
  have less bugs as well.

* Support for the range-based for statement has been added.  It is
  possible to iterate over all intersections of a `GridGlue` object by
  code like
  ```
  GridGlue<...> glue;
  for (const auto& in : intersections(glue)) ...;
  ```
  See the documentation for further details.

* `PSurfaceMerge` has been replaced by a wrapper around `ContactMerge` and
  `OverlappingMerge`.  It is now also deprecated and will be removed in the
  next release.
