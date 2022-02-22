# KEGG over representative analysis

this pipeline is still under development. 

1. run emapper
```bash
./my_enrichment.py annotate -o example/EP1 example/EP1.faa
```

2. gene gene2ko.xls file
``` bash
./my_enrichment.py gene2x -o gene2x example/EP1.emapper.annotations
```

3. build database
``` bash
./my_enrichment.py kegg_db -o kegg_db -g Ralstonia
```

4. run kegg_enrich
``` bash
my_enrichment.py kegg_enrich \
    -l example/interested.list \
    -a gene2x/gene2ko.xls \
    -d kegg_db \
    -s kegg_db/slim/Ralstonia_pathway.list \
    -o example/test_out.xls
```
