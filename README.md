# BdfCGAL
Nastran (BDF) and Medit (MESH) mesh generator based on CGAL library.

Requires [CGAL](https://www.cgal.org/).

Type `help` to see possible options.

```
BdfCGAL help
```

Basic usage:

```
BdfCGAL model.off
```

The command above should generate `model.bdf` and `model.mesh` files.

Boundary conditions can be squashed with `squash_from` and `squash_to` options.

