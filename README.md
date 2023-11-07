
# PANDA2: Protein function prediction using graph neural network

This repository includes code and a pre-trained model of PANDA2 for protein function prediction.

## Stand-alone PANDA2
1. Download the trained model
```
wget http://dna.cs.miami.edu/PANDA2/download_files/gcn37onlyesm_cafa.cnn ./
wget http://dna.cs.miami.edu/PANDA2/download_files/gcn37onlyesm_cafa.gcn ./
```
3. Run PANDA2 with a fasta format file
```
# perl parse_seq_cut.pl $input_fasta $output_dir
>> perl parse_seq_cut.pl example/test_5_samples.fasta example/test_5_samples.
# python panda2_features.py $output_dir
# a $output_dir.2 will automaticly created.
>> python panda2_features.py example/test_5_samples
# make prediction with the feature file
>> python panda2_prediction.py example/test_5_samples.2/panda2_features.pkl
```
## Output explaination:
The output is saved in [example/test_5_samples.2/panda2_prediction.txt](example/test_5_samples.2/panda2_prediction.txt). The file format is as follows:  
- The first two lines contain model information.  
- Subsequent lines follow this format:  
  - The first column is the name of the PDB structure file.  
  - The second column is the GO term ID.  
  - The third column is the confidence score predicted by PANDA2.  
- The last line contains only "END."
```txt
AUTHOR PANDA2
MODEL 1
KEYWORDS graph network, sequence alignment.
T100900000046   GO:0097458      0.49
T100900000046   GO:0003674      0.91
...
T100900000026   GO:0006887      0.16
T100900000026   GO:0046578      0.09
END
```

>**Dependencies**   
>The conda environment is shared via "[environment.yml](environment.yml)".  
>In order to locally run PANDA2, you also need to install blast-2.2.23 and then update $blast_path in panda2_psiblast.py.

## Web-server PANDA2
Submit jobs at
http://dna.cs.miami.edu/PANDA2/

## Citation
```
@article{10.1093/nargab/lqac004,
    author = {Zhao, Chenguang and Liu, Tong and Wang, Zheng},
    title = "{PANDA2: protein function prediction using graph neural networks}",
    journal = {NAR Genomics and Bioinformatics},
    volume = {4},
    number = {1},
    pages = {lqac004},
    year = {2022},
    month = {02},
    issn = {2631-9268},
    doi = {10.1093/nargab/lqac004},
}
```