# ActiveMatterCyclicStretch

Self-propelled particle model with uniaxial cyclic stretch proposed in T. Tagaki et al. [1]

## Description

The source code for simulating the model.

<!--
## include Folder

The `include` folder contains the header files used in the simulation. Each header file is responsible for specific functions and calculations.

- `LJ.hpp`: Calculates the Lennard-Jones potential.
- `nematic.hpp`: Calculates nematic interactions.
- `observe.hpp`: Handles data observation and analysis.
- `pair.hpp`: Manages pair lists and mesh generation.
- `systemparam.hpp`: Defines system parameters and constants.
- `vicsek.hpp`: Calculates the Vicsek model.

## src Folder

The `src` folder contains the implementations of the functions declared in the `include` folder.

- `LJ.cpp`: Implements the calculations for the Lennard-Jones potential.
- `nematic.cpp`: Implements the calculations for $m=2$.
- `pair.cpp`: Implements the pair list and mesh management.
- `vicsek.cpp`: Implements the calculations for $m=1$.
-->

## Requirement

* [Required Library or Tool]

## Usage

1. Set the parameters in `./lib/include/systemparam.hpp`.
2. Compile and run the main program `./cpp/stretch.cpp`.
3. The results of the simulation will be saved in the directory specified by `DATA_PATH`.

## Compilation

Use the following command to compile the program:

```sh
g++ -std=c++17 ./cpp/stretch.cpp ./lib/src/pair.cpp ./lib/src/nematic.cpp ./lib/src/LJ.cpp 
```

## References
 [1] Tagaki T., Nishikawa S., Ishihara S., Active Matter under Cyclic Stretch: Modeling Microtubule Alignment and Bundling. XXXX. 2025;XX. doi:......
