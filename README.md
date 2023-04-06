# Build Detection VCF/BCF

Python script to detect genomic build from a VCF/BCF file. 
Use the coordinates of common variants to identify the build / assembly of a genotype file that is being loaded.

## Dependencies 
```python
import os 
from cyvcf2 import VCF
import pandas as pd
import argparse
```

## Usage

```python
tabix <file_name.vcf.gz>
python3 build_detection_VCF.py <file_name.vcf.gz>
```

## Output

The script returns:

|Message|Genome Build|
|:-:|:-:|
|hg19/GRCh37|build 37|
|hg38/GRCh38|build 38|
|No target variants found to detect the build|no detection|
|Error in build detection|Error|
