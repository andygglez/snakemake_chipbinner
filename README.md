
Start by running:

``` bash
python genconfig.py --mode chipbinner
```

This line will create a `config/config.yaml` with the instructions to run the pipeline, and `config/sample_list.csv` file. Some configurations should be modified according to the analysis. For example:
`
``` yaml
organism: hg38
MS_normalize: True
merge_replicates: True

MS_yaml: references/MS_coefficients.yaml

comparisons:
- treatment1_control1
```

Specify `MS_normalize = True` if your analysis implies MS normalizing `bigWig` files, and set `merge_replicates = False` if you want to use a single replicate as `bigWig`.

If you want to use 2 replicates, specify those replicates in the `config/sample_list.csv` file, or just 1 if no merging is needed. Each replicate/pair-of-replicates should have a condition and a mark associated, and specified in the `config/sample_list.csv`. After each condition is defined in the `config/sample_list.csv` file, the comparisons should be set in the `config/config.yaml` in the `comparisons` level, and separated by an `_`. For example:

![[Pasted image 20241125145716.png]]

![[Pasted image 20241125145758.png]]

The file `references/MS_coefficients.yaml` should specify the MS coefficients used to normalize the `bigWig` files, and each entry should refer to a replicate file, with the `".bam"` suffix removed. For example:

![[Pasted image 20241125150259.png]]


##### Limitations:

- ==Don't specify conditions that only contain one replicate and other with two replicates. If you use this pipelines, either all conditions should have one replicate or all have both replicates.==
- ==Due to the common use of the folders to deposit the final bw files, after each run, you should erase these output files before rerunning the script!==