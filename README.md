DA-MDS
=========
Deterministic Annealing Multidimensional Scaling (DA-MDS) is a high performance
implementation of the [WDA-SMACOF](http://grids.ucs.indiana.edu/ptliupages/publications/WDA-SMACOF_v1.02.pdf)
algorithm.

Success Stories
-----
* Million Sequence Clustering at http://salsahpc.indiana.edu/millionseq/
* The Fungi Phylogenetic Project at http://salsafungiphy.blogspot.com/

# Prerequisites
-----
1. Operating System
  * SPDIAL is extensively tested and known to work on,
    *  Red Hat Enterprise Linux Server release 6.7 (Santiago)
    *  Red Hat Enterprise Linux Server release 5.10 (Tikanga)
    *  Ubuntu 12.04.3 LTS
    *  Ubuntu 12.10
  * This may work in Windows systems depending on the ability to setup OpenMPI properly, however, this has not been tested and we recommend choosing a Linux based operating system instead.

2. Java
  * Download Oracle JDK 8 from http://www.oracle.com/technetwork/java/javase/downloads/index.html
  * Extract the archive to a folder named `jdk1.8.0`
  * Set the following environment variables.
  ```
    JAVA_HOME=<path-to-jdk1.8.0-directory>
    PATH=$JAVA_HOME/bin:$PATH
    export JAVA_HOME PATH
  ```
3. Apache Maven
  * Download latest Maven release from http://maven.apache.org/download.cgi
  * Extract it to some folder and set the following environment variables.
  ```
    MVN_HOME=<path-to-Maven-folder>
    $PATH=$MVN_HOME/bin:$PATH
    export MVN_HOME PATH
  ```

5. OpenMPI
  * We recommend using `OpenMPI 1.10.1` although it work with the previous 1.8 versions. Note, if using a version other than 1.10.1 please remember to set Maven dependency appropriately in the `pom.xml`.

  * Download OpenMPI 1.10.1 from http://www.open-mpi.org/software/ompi/v1.10/downloads/openmpi-1.10.1.tar.gz
  * Extract the archive to a folder named `openmpi-1.10.1`
  * Also create a directory named `build` in some location. We will use this to install OpenMPI
  * Set the following environment variables
  ```
    BUILD=<path-to-build-directory>
    OMPI_1101=<path-to-openmpi-1.10.1-directory>
    PATH=$BUILD/bin:$PATH
    LD_LIBRARY_PATH=$BUILD/lib:$LD_LIBRARY_PATH
    export BUILD OMPI_1101 PATH LD_LIBRARY_PATH
  ```
  * The instructions to build OpenMPI depend on the platform. Therefore, we highly recommend looking into the `$OMPI_1101/INSTALL` file. Platform specific build files are available in `$OMPI_1101/contrib/platform` directory.
  * In general, please specify `--prefix=$BUILD` and `--enable-mpi-java` as arguments to `configure` script. If Infiniband is available (highly recommended) specify `--with-verbs=<path-to-verbs-installation>`. Usually, the path to verbs installation is `/usr`. In summary, the following commands will build OpenMPI for a Linux system.
  ```
    cd $OMPI_1101
    ./configure --prefix=$BUILD --enable-mpi-java
    make;make install
  ```
  * If everything goes well `mpirun --version` will show `mpirun (Open MPI) 1.10.1`. Execute the following command to instal `$OMPI_1101/ompi/mpi/java/java/mpi.jar` as a Maven artifact.
  ```
    mvn install:install-file -DcreateChecksum=true -Dpackaging=jar -Dfile=$OMPI_1101/ompi/mpi/java/java/mpi.jar -DgroupId=ompi -DartifactId=ompijavabinding -Dversion=1.10.1
  ```
  * Few examples are available in `$OMPI_1101/examples`. Please use `mpijavac` with other parameters similar to `javac` command to compile OpenMPI Java programs. Once compiled `mpirun [options] java -cp <classpath> class-name arguments` command with proper values set as arguments will run the MPI Java program.

Building DA-MDS
-----
* Check all prerequisites are satisfied before building DA-MDS
* Clone this git repository from `git@github.com:DSC-SPIDAL/damds.git` Let's call this directory `damdshome`
* Once above two steps are completed, building DA-MDS requires only one command, `mvn install`, issued within `damdshome`.

Running DA-MDS
-----
The following shell script may be used with necessary modifications to run the program.
```sh
#!/bin/bash

# Java classpath. This should include paths to damds dependent jar files and the damds-1.0.jar
# The dependent jar files may be obtained by running mvn dependency:build-classpath command within damdshome
cp=<classpath>

# Obtain working directory
wd=`pwd`
# Character x as a variable
x='x'

# A text file listing available nodes
hosts=<path-to-hostfile>
# Number of nodes
nodes=<num-nodes>
# Number of cores per node
cpn=8

# Options for Java runtime
jopts="-Xms64M -Xmx64M"

# Number of threads to use within one dapwc process
tpn=<threads-per-process>
# Number of processes per node
ppn=<processes-per-node>
# Total parallelism expressed as a pattern TxPxN
# where T is number of threads per process, P is processes per node, and N is total nodes
pat=$tpn$x$ppn$x$nodes

# Number of computing units assigned per process, assuming $cpn is divisible by $ppn
bw=$(($cpn/$ppn))

mmaps=1
# Directory to memory map. Ideally, set this to where tmpfs is in Linux
mmapdir=/dev/shm/$USER

echo "Running $pat on `date`" >> status.txt
# Invoke MPI to run dapwc
mpirun --report-bindings --mca btl ^tcp --hostfile $hostfile --map-by ppr:$ppn:node:PE=$bw --rank-by core -np $(($nodes*$ppn)) java $jopts -cp $cp edu.indiana.soic.spidal.damds.Program -c config$pat.properties -n $nodes -t $tpn -mmaps $mmaps -mmapdir $mmapdir | tee $pat/mds-out.txt
echo "Finished $pat on `date`" >> status.txt
```
The arguments listed in the `mpirun` command fall into three categories.
* OpenMPI Runtime Parameters
  * `--report-bindings` requests OpenMPI runtime to output how processes are mapped to processing elements (cores) in the allocated nodes.
  * `--mca btl ^tcp` instructs to disable tcp, which is useful when running on Infiniband.
  * `--hostfile` indicates the file listing available nodes. Each node has to be a in a separate line.
  * `--map-by ppr:$ppn:node:PE=$bw` controls process mapping and binding. This is a topic on its own right, but the specific values in this example requests processes to be mapped by node while binding each to `bw` number of processing elements. A good set of slides on this topic is available at http://www.slideshare.net/jsquyres/open-mpi-explorations-in-process-affinity-eurompi13-presentation
  * `-np $(($nodes*$ppn))` determines the total number of processes to run and in this case it is equal to `nodes*ppn`
* Java Runtime Parameters
  * `$jopts` in this case lists initial and maximum heap sizes for a JVM instance.
  * `-cp` indicates paths to find required classes where each entry is separated by a `:` (in Linux)
* Program (dapwc) Parameters
  * `-c` points to the configuration file. This is a Java properties files listing values for each parameter that damds requires. Details on these parameters will follow in a later section.
  * `-n` indicates the total number of nodes
  * `-t` denotes the number of threads to use within one instance of damds
  * `-mmaps` is the number of memory maps to use. Set this to 1 for the best performance
  * `-mmapdir` points to the directory where memory maps should be created. Ideally, it should point to `tmpfs` directory in Linux

Configuring damds
-----
The following table summarizes the parameters used in dapwc.

| Parameter | Description | Default Value | Type |
|----------------------------------|-----------------------------------------------------------------|---------------|-----------------------------------------|
| DistanceMatrixFile | Path of the pairwise distance file. | n/a | String |
| WeightMatrixFile | Path of the pairwise weight matrix file or the simple linear weights text file. | n/a | String |
| LabelDataFile | Path of the points' labels file. | n/a | String |
| InitialPoints File | Path of the initial points file. | n/a | String |
| PointsFile | Path of the output points file. | n/a | String |
| TimingFile | Path of the output timing information. | n/a | String |
| SummaryFile | Path of the output summary file. | n/a | String |
| NumberDataPoints | Total number of data points. | n/a | Integer |
| TargetDimension  | Target dimension of the output points. | 3 | Integer |
| DistanceTransform | Raise distances to the power of DistanceTransform. | 1.0 | Double |
| Threshold  | Stress threshold. | 0.000001 | Double |
| Alpha | Cooling factor. | 0.95 | Double |
| TminFactor | The minimum temperature factor. | 0.5 | Double |
| StressIterations | The maximum number of stress loops to run. | 10000 | Integer |
| CGIterations | The maximum number of conjugate gradient loops | 20 | Integer |
| IsSammon | The flag to determine if sammon distances should be used | false | Boolean |
| BlockSize | The block size to use in block matrix multiplication | 64 | Integer |
| IsBigEndian | The flat to indicate the endianness of the binary distance file.  | false | Boolean |
| IsMemoryMapped | The flag to indicate if memory mapped files should be used to load data. | true | Boolean |
| TransformationFunction | The path of the jar file containing additional distance transformations. | null | String |
| WeightTransformationFunction | The path of the jar file containing additional weight transformations. | null | String |
| Repetitions | The number of repetitions (see below). | 1 | Integer |
| MaxTempLoops | The maximum number of temperature loops (see below). | 0 | Integer |
| IsSimpleWeights | The flag to indicate if weights are read from a simple linear file. | false | Boolean |

`Repetitions` is a quick way to test large data sizes using a smaller original
distances and weights files. For example with a NxN data set and a `Repetitions=2`,
DA-MDS will do a 2Nx2N run. It does so by tiling the NxN matrix 4 times
(2 horizontally and 2 vertically).

`MaxTempLoops` allows to test the program for performance without running for
full convergence by allowing it to run only the specified number of temperature loops.
Setting this to `0` will disable it and will do the full run.

Acknowledgement
-----
We like to express our sincere gratitude to Prof. Vivek Sarkar
and his team at Rice University for giving us access and
continuous support for the HJ library. We are also
thankful to FutureSystems project and its support team for their support with
HPC systems. Also, we thank Intel  for their support of the Juliet cluster system
that we used to test DA-MDS. Last but not least OpenMPI community deserves equal
recognition for their valuable support.

We also like to thank the following companies for providing us Open Source licences for their profiler software.

![alt text][jprofiler logo]
 - ej-technologies the creator of JProfiler  (http://www.ej-technologies.com/products/jprofiler/overview.html)

![alt text][yourkit logo]
 - YourKit supports open source projects with its full-featured Java Profiler. YourKit, LLC is the creator of [YourKit Java Profiler](https://www.yourkit.com/java/profiler/index.jsp) and [YourKit .NET Profiler](https://www.yourkit.com/.net/profiler/index.jsp), innovative and intelligent tools for profiling Java and .NET applications


[jprofiler logo]: https://www.ej-technologies.com/images/product_banners/jprofiler_medium.png
[yourkit logo]: https://www.yourkit.com/images/yklogo.png
