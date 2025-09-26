
# Running with docker: two minimal working examples 

In these two examples, we download the image
[scottperkins/gwat](https://hub.docker.com/r/scottperkins/gwat),
and run an example usage of `gwat` that performs parameter
estimation on a gravitational wave injection.
The first example shows how to do this on your local computer
using docker,
the second shows how to do this with `singularity`, which the
container application that runs on clusters (instead of docker).

## Running on your local computer

1. Assuming you already have docker desktop installed,
start docker desktop, e.g. in the terminal type
```
systemctl --user start docker-desktop
``` 

2. Make a directory, e.g. `GW-Analysis`, and cd to that directory. 

3. Clone [gw_analysis_tools](https://github.com/scottperkins/gw_analysis_tools)
to the directory. cd to `gw_analysis_tools/examples/injection_IMRPhenomD`.

4. There should be bash script called `docker_launch.sh`.

Running this file downloads and runs the docker image 
[scottperkins/gwat](https://hub.docker.com/r/scottperkins/gwat)
from Docker Hub. 

5. Run the script 
```
bash docker_launch.sh
```
You should now be in the docker image, and you should be in a directory called `/injection/`.

6. To make and run the binary, execute the following commands
```
make 
./bin/injection.exe
```
The output of the code will be saved
under `/data` (which also contains the injection parameters).
It may take some time to run on your computer. 

7. To examine the output, there are routines already packaged 
in `gw_analysis_tools` and `BayesShip` to unpack and process the 
raw data coming from the sampler. Included in the `python/` directory 
of this example is an example notebook that loads the data and plots a 
corner plot, trace plots, and some useful convergence diagnostics. 
To launch it, you can run 
```
./jupyter_docker_launch.sh
```
from your computer or laptop to launch another docker container hosting a jupyter server.
The server can now be accessed by typing the following URL into the browser on your computer: "localhost:8990".
The token to login to the session is printed in the output of the script, but can also be seen at any time by running 
```
jupyter server list
```
if you login to the container through another session.

If you are having trouble opening jupyter, remember to copy and past
`localhost:8990` (or the other value) in your browser. The token to access
the notebook can be found in the terminal output, e.g. 
```
[I 2023-02-22 20:37:09.988 ServerApp] http://e84cb84deac2:8989/lab?token=[some long string of letters and numbers: this is the token]
```
To shutdown the server, press "ctrl-c" twice in the terminal running the docker container. 
The server will shutdown and the docker container will stop and remove itself.

Alternatively, if you do not want to run docker to run the jupyter notebook, under
`examples/data_analysis`, you'll find the necessary source code to run an example
jupyter notebook analysis. 
Essentially, that directory contains the source code to read from the hdf5 output the
code produces, with the right formatting.

## Running on a cluster

1. Make sure singularity is available, e.g. through `module load singularity`

2. As above, make and cd to directory, e.g. `GW-Analysis`.

3. Clone [gw_analysis_tools](https://github.com/scottperkins/gw_analysis_tools)
to the directory. cd to `gw_analysis_tools/examples/injection_IMRPhenomD`.

4. There should be bash script called `singularity_launch.slurm`.
You may need to change some of the SBATCH commands for your local cluster.

Running this file downloads and runs the docker image 
[scottperkins/gwat](https://hub.docker.com/r/scottperkins/gwat)
from Docker Hub, and then runs the `runscript_for_slurm.sh`
script, which makes and runs the executable within the image.

5. Launch with 
```
slurm singularity_launch.slurm
```

Output should be saved in the relevant `data` directory
in your scratch folder (or elsewhere if you changed the .slurm file). 
