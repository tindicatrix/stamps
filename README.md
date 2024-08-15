# stamps
Takes SExtractor catalog files and makes stamps from Roman sims for Real/Bogus Classification

Workflow
-
1. Run SExtractor on the image w/ (sex [input.fits] -c config.sex -PARAMETER_NAME config.sex.param -FILTER_NAME config.sex.conv -CATALOG_NAME [output_name.cat])
    - sample configuration parameters are provided, but can be modified
    - ensure that the output parameters include ALPHA_J2000 and DELTA_J2000 values
2. Run either just fullstamps.py OR (stamps.py THEN makeCSV.py)
    - there's no difference in either of the two methods, fullstamps.py is just the other two combined into one python file.
  
Notes
-
Make sure you edit every file you run to your specifications

A sample SLURM job script is provided (job.sh)

Minor edits have been made without testing before upload, there may be bugs (hopefully not)
