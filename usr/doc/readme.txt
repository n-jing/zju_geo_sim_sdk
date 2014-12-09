# -*- mode: org; -*-
#+AUTHOR:    Jin Huang
#+DATE:      2010-09-08

Geometry Simulation Group SDK

* Goal

For cross compiling between Linux and Windows, Mingw and Visual C++,
32 bit and 64 bit.

* Implementation

Like Linux package management, a typical package includes four parts:
include, bin, lib and doc.  The file names and path should consistent
with Linux package layouts.  We use $HOME/usr as the root of the sdk.

To distinguish incompatible tool chain, such as 64 bit gnu tool chain,
32 bit Visual C++ tool chain, and incompatible runtime libraries
(especially for Visual C++), we describe a build with keywords of (OS,
[BIT], [compiler], [version], [configuration]), and some keywords are
optional if it will not cause problem.  For example, (Windows, 64, vc,
9sp1, double_bounding_box), (Windows, 32).  Two builds are compatible
if and if only their keywords are not conflict.

To make it easy, we organize the builds of a package into a
hierarchical structure. We require the builds in the upper level are
compatible with all the lower ones.  For example: The header in
/usr/include can be used for both Windows and Linux, and the dll put
in /usr/Windows/bin can be used for both Mingw and Visual C++.

The order of hierarchical is OS/BITS/COMPILER+VERSION:
|------------------+-------------------------|
| OS               | Windows, Linux, Dawin   |
| BITS             | 32, 64                  |
| COMPILER+VERSION | gcc3, gcc4, vc9sp1      |
|------------------+-------------------------|

Debug is distinguished from release by a suffix "d".

We do not use COMPILER/VERSION for simple, since in most of the cases,
only one version of a compiler is used.

* Shortage

The hierarchical structure does not reflect the above rule.  For
example, if a header file (e.g. inttypes.h) is compatible and only
compatible for both 32/vc9sp1 and 64/vc9sp1.  In current
configuration, we have to copy these files in two places:
32/vc9sp1/include and 64/vc9sp1/include.  Another choice is to use
another structure OS/COMPILER+VERSION/BITS.  However, if so, it's a
problem to put BLAS.dll, which is compatible for all the compiler with
the same BITS. If in Linux, such problem can be easily solved by
symbolic link, but Windows has no such feature.

* References
The name of OS is follow cmake, debug suffix is follow osg and many
other packages.
