# Example run 

1) Generate input file (input_example.txt) varying ```gamma``` between 1-3: 

```
WHO_H5N1_2009.newick
4,0.1,1-3(1) 
```

2) Run PhyCLIP. Note the following flags:  
  i. ```--tree_outgroup a/goose/guangdong/1/1996_{0}```: Reroot tree using 'a/goose/guangdong/1/1996_{0}' as outgroup.  
  ii. ```--collapse_zero_branch_length 1```: Collapse internal nodes with zero branch length.  
  iii. ```--subsume_subclusters 1```: Subsume any sub-clusters into their parent cluster to avoid over-delineation.  
  iv. ```--pdf_tree```: Generate pdf tree files with cluster annotations.  
  v. ```--optimise intermediate```: Programatically determine the optimal parameter set for intermediate-resolution clustering.  

```
phyclip.py -i input_example.txt --tree_outgroup a/goose/guangdong/1/1996_{0} --collapse_zero_branch_length 1 --subsume_subclusters 1  --pdf_tree --optimise intermediate 
```
