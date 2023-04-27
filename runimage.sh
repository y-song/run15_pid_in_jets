#!/bin/bash

#singularity run -B /gpfs01 -B /star -B /gpfs01/star/pwg/youqi/run15 star_test.simg
#singularity run -B /gpfs01 -B /star -B /gpfs01/star/pwg/youqi/run15 /gpfs01/star/pwg/elayavalli/star_test.simg
singularity run -B /gpfs01 -B /star -B /gpfs01/star/pwg/youqi/run15 star_star.simg
