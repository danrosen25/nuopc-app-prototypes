# ESMX_MetadataExchange

[ESMX](https://github.com/esmf-org/esmf/tree/develop/src/addon/ESMX) is used to implement a coupled system with two components (CompA and CompO). Both components are provided as NUOPC-compliant models. 

## Primary Artifacts

Files and sub-directories that implement the fundamental concept demonstrated by the prototype. These are the primary artifacts to look at and to pattern actual user code after.

- `esmxBuild.yaml` - Standard ESMX YAML file describing the build dependencies of the `example.exe` (the executable) on components CompO and CompA.
- `esmxRun.yaml`   - Standard ESMX YAML file describing the run configuration: CompO is C2, CompA is C1, and defining the run sequence.
- `CompO`          - Example C2 component. It uses a simple CMake based build system.
- `CompA`          - Example C1 component. It uses a GNU Make based build system.

### Usage

1. Build the ESMX executable by using the command line tool:
   ```
   ESMX_Builder
   ```
   This assumes that the `bin` directory of the desired ESMF installation is present in the user's `PATH` environemnt variable. The `ESMX_Builder` tool first compiles all the required sub-components, and then links it into the final executable.
2. Run the `./install/bin/example.exe` executable on 4 PETs using the appropriate MPI launch procedure. E.g.:
   ```
   mpirun -np 4 ./install/bin/example.exe
   ```

## Secondary Artifacts

Files that are needed for the integration into ESMF's automated testing infrastructure for regression testing. These artifacts might be interesting to look at, but generally should *not* be used as patterns to follow in actual projects.

- `Makefile`        - GNU Makefile that defines targets that are used by the automated ESMF regression testing script.

### Usage

1. The default target of the `Makefile` calls the `ESMX_Builder` command line tool to build the ESMX executable:
   ```
   make
   ```
2. The `run` target of the `Makefile` uses the MPI launch procedure - identified by ESMF - to run `./install/bin/example.exe` on 4 PETs.
   ```
   make run
   ```
3. The `distclean` target of the `Makefile` removes all of the generated files.
   ```
   make distclean
   ```
