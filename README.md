# BdfCGAL
Nastran (BDF) and Medit (MESH) mesh generator based on CGAL library.

Requires [CGAL](https://www.cgal.org/).

Build with CMake:

```
cmake .
make
```

Type `help` to see possible options.

```
BdfCGAL help
```

Basic usage:

```
BdfCGAL generate example/ahu.off
```

The command above should generate `ahu.bdf` and `ahu.mesh` files.

Boundary conditions can be squashed with `squash_from` and `squash_to` options.

