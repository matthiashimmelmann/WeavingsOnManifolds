# WeavingsOnManifolds
 
This Julia package contains test classes and the nonlinear optimization problem corresponding to a tensegrity model for weavings on surfaces. The following functionality is supplied:

+ The class `WeavingsOnManifolds` contains information about the `manifold` $\mathcal{M}$, the over-under information of the weave (`offsetList`), the `offset` of the $\varepsilon$-offset surface $\bigcup_{p\in \mathcal{M}}\{p\pm \varepsilon \cdot n_\mathcal{M}(p)\}$ for the unit normal field $n_\mathcal{M}(p)$, the framework's `bars` and `cables` and the variable and pinned vertex positions (`variablePositions` and `pinnedVertices`).
+ When the tensegrity model for weavings is set up, we can equilibrate it using `computeOptimalWeaving`. 
+ The weave can be visualized using the routine `plotWeaving` and there is a functionality to turn it into a `.poly` file as a smooth weave on the surface (e.g. `toSphere`).

The following test cases are implemented:

+ Sphere: `test_sphere_truncated_dodecahedron`, `test_sphere_truncated_icosihedron`, `test_sphere_truncated_octahedron`, `test_sphere_antiprism_penta`, `test_sphere_antiprism_square`, `test_sphere_trunctetra`, `test_sphere_truncatedcube`, `test_sphere_dodecahedron`, `test_sphere_octahedron`, `test_sphere_tetrahedron` and `test_sphere_rhombicuboctahedron`
+ Spheroid: `test_spheroid_tetrahedron`, `test_spheroid_octahedron`, `test_spheroid_rhombioctahedron`, `test_spheroid_dodecahedron`, `test_spheroid_trunctetra`, `test_spheroid_truncatedcube`, `test_spheroid_truncatedoctahedron` and `test_spheroid_truncatedicosihedron`
+ Cylinder: `test_cylinder_kagome` and `test_cylinder_kagome2`
+ Hyperboloid: `test_hyperboloid`
+ Torus of Revolution: `test_torus_trefoil`
+ Flat Torus: `test_flattorus_plain`

Moreover, we can compute deformations paths of the weavings from the sphere into the spheroid and record the corresponding hysteresis loops using `plot_hystereces`. 
