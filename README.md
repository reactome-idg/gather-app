#### Hierarchy to allow for dataOps parallel docker services 

```
└── gather-app
      │
      ├── python
      │     ├── switch
      │     ├── docker                 
      │     └── scripts
      │
      └── data 
            ├── curated
            │     |── gk_central 
            │     |── FI_reactome
            │     └── TCRD_Pharos 
            │
            ├── source
            │     |── GTEx
            │     │    |── ui-fetch 
            │     │    |── raw-data
            │     │    └── meta-data 
            │     │ 
            |     └── ...
            │ 
            └──  processed
                  |── GTEx
                  │     |── raw-count
                  │     |── cpm
                  │     │    |── cpm 
                  │     │    |── wo-outliers
                  │     │    └── eda
                  │     |── fi
                  │     |    |── gene-ensg-fi.csv
                  │     |    |── cpm-spearman
                  │     |    |── wel-spearman
                  │     |    |── mcl-spearman
                  │     |    └── ea-spearman
                  │     |    
                  │     └──  adjacency 
                  │           └── spearman
                  │
                  └── ...

```

