# Streptococcus-Project (Summer, 2016)
Scripts Repository

Pipeline: 
Sequencing (to obtian short fragmented DNA sequences which are stored in .fastq files) --> Cortex Var Process (to compare sample genomes with reference genome so that the variances will be summarized in a .vcf file) --> Upload Data to Database (parse the .vcf file into tables based on the framework of database) --> Biostatistical Analysis and Data Representation --> Webcoding (retrieve data from database and convert data into .json files)

Database Framework: genome (TABLE, Primary Key: gid): stores metadata of reference and sample genomes , varlist (TABLE, Primary Key: gid, var_id): the primary column to store information of variances in comparison between reference genome and sample genomes, reforf (TABLE, Primary Key: var_id, gid): to store reference genome information, var (TABLE Primary Key: var_id): to store the statistical information of variances.
