---
project: fortran_cubatures
summary: Gaussian cubature rules on the reference elements of the finite element method: the derived type `cubature`, selected by degree of exactness.
author: Will Klausler
github: https://github.com/willklausler
project_github: https://github.com/willklausler/fortran_cubatures
license: mit
src_dir: ./src
output_dir: ./build/doc
page_dir: ./build/pages
preprocess: false
display: public
source: true
search: true
graph: false
---

`use cubatures` gives the type `cubature` and the element constants `CUB_LIN`,
`CUB_TRI`, `CUB_QUA`, `CUB_TET`, `CUB_HEX` and `CUB_WED`.  Start with the
[user manual](page/index.html), which is also available as plain text in
[llms.txt](llms.txt).

The specific procedures behind the constructor and `set` are private; call
them through the generic names.

```fortran
use cubatures
type(cubature) :: q
q = cubature(CUB_TET, 4)
print *, q%npoints, sum(q%weights)   ! 27  0.16666666666666669
```
