# Docker

Dockerfile and ancillary files to build a Docker version of the NGS pipeline contained in the repository mri-ihg/ngs_pipeline.

- The github mri-ihg/ngs_pipeline repository content must be copied into build/pipeline
- In the configuration_test folder i's necessary to:

  - Provide MYSQL parameters into the file "fields"
  - Then execute the entryconfig.sh script in that directory to create the config files from the templates provided.

The Docker container is then built with build.sh and run with ./run.sh 
The two scripts can be customized to define custom mount points.

In the current status individual pipeline steps can be run from inside the Docker.



