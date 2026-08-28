0 Install docker
1 docker pull acculix/solvers:latest. 
  OR
  docker pull docker.1ms.run/acculix/solvers:latest 
2 Check the image:
  docker run -it --rm -v $(pwd):/work -w /work acculix/solvers:latest /bin/bash
  OR
  docker run -it --rm -v $(pwd):/work -w /work docker.1ms.run/acculix/solvers:latest /bin/bash】.
  Then 'exit'.
3 Run
  [a] mkdir work; cd work
  [b] Create script file:run.sh
      #!/bin/sh 
     sudo docker run -it --rm -v $(pwd):/work -w /work docker.1ms.run/acculix/solvers:latest /bin/bash "$@"      
  [c] submit job：./run.sh mpirun --allow-run-as-root -n 2 ccFEM_mpi ./4pointbending.txt
  [d] Run paraview for postprocessing.
