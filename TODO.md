## TODO

- [x] doublecheck and fix plot aspect ratio
  - [ ] fix padding issues, ideally dont need to call `ax.set_aspect('equal')` and can just handle it through calculation
- [x] ~~offset of initial wakepoint for vortex blob mechanism~~ Limit vortex blobs to wake votices
- [ ] vectorize rollup loops
- [ ] fix use of trapezoid in lift calculations to account for coincident x-coords from bodies with upper and lower surfaces
- [ ] vectorization self influence terms
- [ ] vortex core merging
- [ ] moveable wake separation point
- [ ] rotation and offset point for body geometry definition
- [ ] refactor code into modules
- [ ] update motion section as function and reqork how positioning is handled
- [ ] (possibly) handle upper and lower surfaces of closed bodies as separate bodies to iterate over
- [ ] (possibly) change format of how wake info is indexed. currently newest[0] to oldest[-1]. Slightly unintuivitve in some sections and is not consistent with other arrays, but changing could be more trouble than its worth
- [x] Make sure kutta condition is being enforced correctly for closed bodies
